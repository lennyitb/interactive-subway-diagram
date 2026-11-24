"""Collapse the 4D diagram's parallel service lanes into a line graph with one centreline per trunk.

Reads  data/reference/nyc-subway-4d.json   (source-faithful fit: one lane per service, as drawn)
Writes data/nyc-subway-diagram.json        (nodes, edges carrying ordered lanes, routes as edge walks)

Method
  1. Sample every route pattern path every STEP units.
  2. At each sample, look across the line: other routes running parallel at the source lane pitch
     (10 units, chained outward) belong to the same trunk, unless trunks_def.NEVER_MERGE says otherwise.
     Shift the sample to the middle of that trunk.
  3. Smooth the shift along each route so branches peel off the trunk centreline gradually.
  4. Snap the shifted samples of all routes onto shared graph nodes, so services in one trunk share
     one centreline, then contract chains into edges between junctions, platforms and ends.
Needs numpy only.
"""
import json, math, os, sys, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from trunks_def import NEVER_MERGE, ALWAYS_MERGE, LANE_OF
from layout_def import SHARE, PATHS, ORDER, POINTS, MOVE, WALK, GROUP, PLATFORM, APART, ANCHOR

REF = os.path.join(HERE, '../../data/reference/nyc-subway-4d.json')
OUT = os.path.join(HERE, '../../data/nyc-subway-diagram.json')

STEP = 2.0            # sample spacing along route paths
PITCH = 10.0          # source lane pitch
PITCH_TOL = 1.6       # accepted deviation from the pitch when chaining lanes
SAME_LANE = 1.5       # lateral distance under which two routes share a drawn lane
PAR_COS = 0.96        # parallel test (about 16 degrees)
MAX_HALF = 60.0       # widest trunk looked at, either side
MEDIAN_WIN = 9        # samples: removes single-sample flicker at sharp corners
MIN_RUN = 12          # samples: shorter trunk memberships are treated as noise
GAP_FILL = 20         # samples: bridge lanes that stop looking parallel for this long at jogs
RAMP = 10             # samples: a branch peels off its trunk over this many samples either side
SNAP = 3.0            # a route joins a chain laid down by another route within this distance
SNAP_OUT = 5.0        # ... and leaves it beyond this one
TIGHT = 1.0           # distance counted as sitting on the chain (joins and departures happen there)
PLATFORM_JOIN = 7.0   # stop nodes of one station closer than this form one platform
PLATFORM_ALONG = 30.0 # ... as do stop nodes up to this far apart along one chain
PLATFORM_NEAR = 11.0  # ... and stops of services that share track close by
SIMPLIFY = 0.35       # Douglas-Peucker tolerance for edge polylines
BEND_SHORT = 9.0      # segments shorter than this between straight pieces are part of one bend
JOG_MAX = 8.0         # sideways nudges smaller than this between parallel straights are dropped
STRAIGHT_TOL = 5.0    # track within this of a straight line ...
STRAIGHT_MIN = 60.0   # ... for at least this long is made exactly straight

ref = json.load(open(REF))
routes = ref['routes']
stations = {s['id']: s for s in ref['stations']}
name_xy = {s['name']: (s['x'], s['y']) for s in ref['stations']}
color = {r['id']: r['color'] for r in routes}
lane_of = lambda rid: LANE_OF.get(rid, rid)

def in_scope(scope, p):
    if scope is None: return True
    if len(scope) == 2:
        (x, y), r = name_xy[scope[0]], scope[1]
        return math.hypot(p[0] - x, p[1] - y) <= r
    x0, y0, x1, y1 = scope
    return x0 <= p[0] <= x1 and y0 <= p[1] <= y1
rules = [(set(a), set(b), scope) for a, b, scope in NEVER_MERGE]
forced = [(set(a), set(b), gap, scope) for a, b, gap, scope in ALWAYS_MERGE]

def conflict(ra, rb, p):
    for a, b, scope in rules:
        if in_scope(scope, p) and ((ra & a and rb & b) or (ra & b and rb & a)):
            return True
    return False
def joins(ra, rb, gap, p):
    # are two neighbouring lanes, gap apart, part of one trunk?
    if conflict(ra, rb, p): return False
    if abs(gap - PITCH) <= PITCH_TOL: return True
    for a, b, wide, scope in forced:
        if gap <= wide and in_scope(scope, p) and ((ra & a and rb & b) or (ra & b and rb & a)):
            return True
    return False

# ---------------------------------------------------------------- 0. clean input paths
# The source fit occasionally runs a route past a vertex and back again. Drop those reversal spikes.
def despike(path, stop_idx):
    pts = [tuple(p) for p in path]; keep_stop = set(stop_idx); idx = list(range(len(pts))); removed = 0
    changed = True
    while changed:
        changed = False
        for k in range(1, len(idx) - 1):
            a, b, c = np.array(pts[idx[k - 1]]), np.array(pts[idx[k]]), np.array(pts[idx[k + 1]])
            u, v = b - a, c - b
            nu, nv = np.hypot(*u), np.hypot(*v)
            if idx[k] in keep_stop or nu < 1e-9 or nv < 1e-9: continue
            if (u @ v) / (nu * nv) < -0.5:
                del idx[k]; removed += 1; changed = True; break
    remap = {old: new for new, old in enumerate(idx)}
    return [list(pts[i]) for i in idx], [remap[i] for i in stop_idx], removed
spikes = 0
for r in routes:
    for p in r['patterns']:
        p['path'], p['stopPathIndex'], n = despike(p['path'], p['stopPathIndex']); spikes += n
print('reversal spikes removed from source paths', spikes)

# ---------------------------------------------------------------- 1. samples
pats = []          # (route id, pattern dict)
P, U, OWN, VTX = [], [], [], []
pat_slice = []
for r in routes:
    for p in r['patterns']:
        pi = len(pats); pats.append((r['id'], p))
        path = np.array(p['path'], float)
        start = len(P); vtx = []
        for i in range(len(path) - 1):
            a, b = path[i], path[i + 1]
            L = float(np.hypot(*(b - a)))
            n = max(1, int(math.ceil(L / STEP)))
            vtx.append(len(P))
            for k in range(n):
                P.append(a + (b - a) * (k / n)); OWN.append(pi)
        vtx.append(len(P)); P.append(path[-1]); OWN.append(pi)
        VTX.append([v for v in vtx])
        pat_slice.append((start, len(P)))
P = np.array(P); OWN = np.array(OWN)
U = np.zeros_like(P)
for pi, (s, e) in enumerate(pat_slice):
    seg = P[s:e]
    d = np.zeros_like(seg)
    d[1:-1] = seg[2:] - seg[:-2]; d[0] = seg[1] - seg[0]; d[-1] = seg[-1] - seg[-2]
    U[s:e] = d / np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-9)[:, None]
NRM = np.stack([-U[:, 1], U[:, 0]], 1)
ROUTE = np.array([pats[o][0] for o in OWN])
print('samples', len(P), 'patterns', len(pats))

class Grid:
    def __init__(self, cell):
        self.cell = cell; self.cells = collections.defaultdict(list)
    def key(self, p):
        return (int(math.floor(p[0] / self.cell)), int(math.floor(p[1] / self.cell)))
    def add(self, i, p):
        self.cells[self.key(p)].append(i)
    def near(self, p, r):
        cx, cy = self.key(p); k = int(math.ceil(r / self.cell)); out = []
        for x in range(cx - k, cx + k + 1):
            for y in range(cy - k, cy + k + 1):
                out.extend(self.cells.get((x, y), ()))
        return out

# ---------------------------------------------------------------- 2. cross-sections
# LATQ[i, q]: lateral offset (own normal frame) from sample i to pattern q where q runs parallel alongside.
g = Grid(20.0)
for i, p in enumerate(P): g.add(i, p)
NPAT = len(pats)
LATQ = np.full((len(P), NPAT), np.nan)
for i in range(len(P)):
    J = np.array(g.near(P[i], MAX_HALF), int)
    J = J[OWN[J] != OWN[i]]
    d = P[J] - P[i]
    lon = d @ U[i]; lat = d @ NRM[i]; par = np.abs(U[J] @ U[i])
    m = (np.abs(lon) < STEP * 0.6) & (np.abs(lat) < MAX_HALF) & (par > PAR_COS)
    best = {}
    for j, lo_, la in zip(J[m].tolist(), lon[m].tolist(), lat[m].tolist()):
        q = OWN[j]
        if q not in best or abs(lo_) < best[q][0]: best[q] = (abs(lo_), la)
    for q, (_, la) in best.items(): LATQ[i, q] = la
# bridge short gaps where the parallel test fails at staggered corners and jogs
for s, e in pat_slice:
    for q in range(NPAT):
        col = LATQ[s:e, q]; k = 0; n = len(col)
        while k < n:
            if not np.isnan(col[k]): k += 1; continue
            j = k
            while j < n and np.isnan(col[j]): j += 1
            if k > 0 and j < n and j - k <= GAP_FILL and abs(col[k - 1] - col[j]) < 2.5:
                col[k:j] = np.linspace(col[k - 1], col[j], j - k + 2)[1:-1]
            k = j
MID = np.zeros(len(P))          # lateral offset from the sample to its trunk centre (own normal frame)
SIG = [None] * len(P)           # the routes in the sample's trunk
for i in range(len(P)):
    # variants (<6>, <7>, <F>) count as their parent service when deciding what forms a trunk
    items = [(float(LATQ[i, q]), lane_of(pats[q][0])) for q in np.nonzero(~np.isnan(LATQ[i]))[0]]
    items.append((0.0, lane_of(ROUTE[i]))); items.sort()
    lanes = []                    # [first offset, {routes}, last offset]
    for o, rid in items:
        if lanes and o - lanes[-1][2] < SAME_LANE:
            lanes[-1][1].add(rid); lanes[-1][2] = o
        else:
            lanes.append([o, {rid}, o])
    lanes = [((l[0] + l[2]) / 2, l[1]) for l in lanes]
    own = min(range(len(lanes)), key=lambda k: abs(lanes[k][0]))
    lo = hi = own
    while hi + 1 < len(lanes) and joins(lanes[hi][1], lanes[hi + 1][1], lanes[hi + 1][0] - lanes[hi][0], P[i]):
        hi += 1
    while lo - 1 >= 0 and joins(lanes[lo][1], lanes[lo - 1][1], lanes[lo][0] - lanes[lo - 1][0], P[i]):
        lo -= 1
    MID[i] = (lanes[lo][0] + lanes[hi][0]) / 2
    SIG[i] = frozenset().union(*[lanes[k][1] for k in range(lo, hi + 1)])

# ---------------------------------------------------------------- 3. smoothing
def median(a, w):
    h = w // 2; out = a.copy()
    for k in range(len(a)):
        out[k] = np.median(a[max(0, k - h):k + h + 1])
    return out
LAT = np.zeros(len(P))          # own lane position relative to trunk centre, own normal frame
C = np.zeros_like(P)            # sample moved onto the trunk centreline
def centre_path(s, e):
    """Centreline for one route: its samples moved onto their trunk centre.
    Cross-sections are unreliable where lanes turn corners at staggered points, swap sides or peel off, so
    runs of a trunk membership or offset lasting under MIN_RUN samples are rebuilt from the stable runs
    either side: with the same offset (a trunk turning a corner) the offset is held, keeping the corner
    concentric; with a different one (lanes crossing, a branch leaving) the centre follows a straight chord,
    ramped over up to RAMP samples of straight path either side when it steps sideways."""
    n = e - s; own = P[s:e]
    mid = median(MID[s:e], MEDIAN_WIN)
    state = [(SIG[s + k], int(round(mid[k] / 2.0))) for k in range(n)]
    runs = []; k = 0
    while k < n:
        j = k
        while j < n and state[j] == state[k]: j += 1
        runs.append((k, j)); k = j
    ok = np.zeros(n, bool)
    for a_, b_ in runs:
        if b_ - a_ >= MIN_RUN or a_ == 0 or b_ == n: ok[a_:b_] = True
    blocks = []; k = 1
    while k < n:
        if not ok[k]:
            j = k
            while j < n and not ok[j]: j += 1
            blocks.append((k, j)); k = j
        elif not ok[k - 1] or state[k] == state[k - 1]: k += 1
        else: blocks.append((k, k)); k += 1      # two stable runs meet
    chords = []
    for k, j in blocks:
        if k == 0 or j >= n: continue
        lo, hi = k - 1, j
        if abs(mid[hi] - mid[lo]) < 1.0:
            mid[k:j] = np.linspace(mid[lo], mid[hi], j - k + 2)[1:-1]
        else:
            chords.append((lo, hi))
    sharp = own + NRM[s:e] * mid[:, None]
    c = sharp.copy()
    arc = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(own, axis=0).T))])
    for lo, hi in chords:
        if float(np.hypot(*(sharp[hi] - sharp[lo]))) - (arc[hi] - arc[lo]) > 0.5:
            u0, u1 = U[s + lo], U[s + hi]           # sideways step: ramp, along straight path only
            for _ in range(RAMP):
                if lo > 0 and U[s + lo - 1] @ u0 > 0.966: lo -= 1
                if hi < n - 1 and U[s + hi + 1] @ u1 > 0.966: hi += 1
        t = (arc[lo:hi + 1] - arc[lo]) / max(1e-9, arc[hi] - arc[lo])
        c[lo:hi + 1] = sharp[lo] + (sharp[hi] - sharp[lo]) * t[:, None]
    return mid, c
for s, e in pat_slice:
    mid, c = centre_path(s, e)
    LAT[s:e] = -mid
    C[s:e] = c

# ---------------------------------------------------------------- 4. build the graph
# Routes are added longest first. A route that runs within SNAP of the graph built so far follows it
# (projected onto its segments); elsewhere it lays down new nodes. Nodes are only inserted on existing
# chains where a route joins, leaves or stops, so routes sharing a trunk share one chain.
NP = []                                  # node positions
adj = collections.defaultdict(dict)      # node -> {node: set(route ids)}
sg = Grid(6.0)                           # segment grid: cell -> set((a, b))
seg_cells = {}
def new_node(p):
    NP.append(np.array(p, float)); return len(NP) - 1
def seg_add(a, b):
    k = (a, b) if a < b else (b, a)
    pa, pb = NP[a], NP[b]; cells = set()
    for t in np.linspace(0, 1, max(2, int(np.hypot(*(pb - pa)) / 3) + 2)):
        cells.add(sg.key(pa + (pb - pa) * t))
    seg_cells[k] = cells
    for c in cells: sg.cells[c].append(k)
def seg_del(a, b):
    k = (a, b) if a < b else (b, a)
    for c in seg_cells.pop(k, ()):
        sg.cells[c].remove(k)
def link(a, b, rid):
    if b not in adj[a]: seg_add(a, b)
    adj[a].setdefault(b, set()).add(rid); adj[b].setdefault(a, set()).add(rid)
MINE = set()                             # nodes laid down by the route being added
def project(c, u, r=None):
    r = r or SNAP; best = None
    for (a, b) in set(sg.near(c, r)):
        if a in MINE or b in MINE: continue
        pa, pb = NP[a], NP[b]; ab = pb - pa; L2 = float(ab @ ab)
        if L2 < 1e-12: continue
        if abs(ab @ u) / math.sqrt(L2) < 0.94: continue     # only follow lines running alongside, never ones crossing
        t = min(1.0, max(0.0, float((c - pa) @ ab) / L2)); q = pa + ab * t
        dd = float(np.hypot(*(q - c)))
        if dd < r and (best is None or dd < best[0]): best = (dd, a, b, t)
    return best
SPLIT = {}                               # (a, b) -> node inserted between them
def split(a, b, t):
    if t < 0.02 or np.hypot(*((NP[b] - NP[a]) * t)) < 0.3: return a
    if t > 0.98 or np.hypot(*((NP[b] - NP[a]) * (1 - t))) < 0.3: return b
    q = new_node(NP[a] + (NP[b] - NP[a]) * t); rs = adj[a].pop(b); adj[b].pop(a); seg_del(a, b)
    SPLIT[(a, b)] = q; SPLIT[(b, a)] = q
    for x, y in ((a, q), (q, b)):
        adj[x][y] = set(rs); adj[y][x] = set(rs); seg_add(x, y)
    return q
def walk(a, b, budget):
    # shortest path a -> b along the graph, no longer than budget
    if a == b: return [a]
    dist = {a: 0.0}; prev = {}; pq = [(0.0, a)]
    import heapq
    while pq:
        dcur, x = heapq.heappop(pq)
        if x == b: break
        if dcur > dist.get(x, 1e18) or dcur > budget: continue
        for y in adj[x]:
            if y in MINE: continue
            nd = dcur + float(np.hypot(*(NP[y] - NP[x])))
            if nd < dist.get(y, 1e18) and nd <= budget:
                dist[y] = nd; prev[y] = x; heapq.heappush(pq, (nd, y))
    if b not in dist: return None
    out = [b]
    while out[-1] != a: out.append(prev[out[-1]])
    return out[::-1]
stop_samples = collections.defaultdict(set)
for pi in range(len(pats)):
    for vi in pats[pi][1]['stopPathIndex']: stop_samples[pi].add(VTX[pi][vi])
SNODE = {}                               # stop sample -> node
pat_nodes = {}
warn_walk = 0
order = sorted(range(len(pats)), key=lambda pi: -(pat_slice[pi][1] - pat_slice[pi][0]))
for pi in order:
    s, e = pat_slice[pi]; rid = pats[pi][0]; MINE.clear()
    seq = []; on = False; tight = None; arc = 0.0; prev_d = None
    def go_to(n, follow):
        # follow: n is further along the chain we are on, so walk there over the existing graph
        global warn_walk
        if not seq: seq.append(n); return
        if n == seq[-1]: return
        w = walk(seq[-1], n, arc + 6.0) if follow else None
        if w is None or len(w) < 2:
            if follow: warn_walk += 1
            w = [seq[-1], n]
        for x, y in zip(w[:-1], w[1:]): link(x, y, rid)
        seq.extend(w[1:])
    def off_node(i):
        n = new_node(C[i]); MINE.add(n); go_to(n, False)
        if i in stop_samples[pi]: SNODE[i] = n
    for i in range(s, e):
        if i > s: arc += float(np.hypot(*(C[i] - C[i - 1])))
        stop = i in stop_samples[pi]; last = i == e - 1
        if on:
            pr = project(C[i], U[i], SNAP_OUT)
            if pr is not None:
                if stop or last:
                    n = split(pr[1], pr[2], pr[3]); go_to(n, True); arc = 0.0; tight = i
                    if stop: SNODE[i] = n
                elif pr[0] < TIGHT:
                    tight = i
                continue
            # leaving: depart from the last sample that sat tight on the chain
            pr = project(C[tight], U[tight], SNAP_OUT)
            if pr is not None:
                n = split(pr[1], pr[2], pr[3]); go_to(n, True)
            on = False; prev_d = None
            for j in range(tight + 1, i + 1): off_node(j)
            continue
        pr = project(C[i], U[i], SNAP)
        settled = pr is not None and (pr[0] < TIGHT or stop or last or not seq or (prev_d is not None and pr[0] >= prev_d - 0.05))
        if settled:
            n = split(pr[1], pr[2], pr[3]); go_to(n, False); arc = 0.0; on = True; tight = i
            if stop: SNODE[i] = n
            continue
        prev_d = pr[0] if pr is not None else None
        off_node(i)
    pat_nodes[pi] = seq
def expand(a, b):
    if b in adj[a]: return [b]
    q = SPLIT[(a, b)]
    return expand(a, q) + expand(q, b)
for pi, seq in pat_nodes.items():
    full = seq[:1]
    for a, b in zip(seq[:-1], seq[1:]): full.extend(expand(a, b))
    pat_nodes[pi] = full
pos = np.array(NP)
print('raw nodes', len(pos), 'walk fallbacks', warn_walk)

# ---------------------------------------------------------------- 5. platforms
def nearest_marker(st, xy):
    return min(st['markers'], key=lambda m: math.hypot(m['x'] - xy[0], m['y'] - xy[1]))
stop_raw = collections.defaultdict(list)    # station id -> [(pi, k, node)]
for pi, (rid, p) in enumerate(pats):
    for k, (sid, vi) in enumerate(zip(p['stops'], p['stopPathIndex'])):
        stop_raw[sid].append((pi, k, int(SNODE[VTX[pi][vi]])))
platform_of = {}                             # (pi, k) -> node
platforms = collections.defaultdict(list)    # station id -> [node]
pat_set = {pi: set(seq) for pi, seq in pat_nodes.items()}
MINE.clear()
def shared_near(a, ra, rb, radius):
    # do routes ra and rb run together on some stretch of graph within `radius` of node a?
    dist = {a: 0.0}; q = [(0.0, a)]
    import heapq
    while q:
        dcur, x = heapq.heappop(q)
        if dcur > dist.get(x, 1e18): continue
        for y, rs in adj[x].items():
            if ra in rs and rb in rs: return True
            nd = dcur + float(np.hypot(*(pos[y] - pos[x])))
            if nd <= radius and nd < dist.get(y, 1e18): dist[y] = nd; heapq.heappush(q, (nd, y))
    return False
def same_platform(ia, ib):
    a, b = ia[2], ib[2]
    d = float(np.hypot(*(pos[a] - pos[b])))
    if d < PLATFORM_JOIN: return True
    # dots the source staggered along one line: close, and joined by a straight run of the chain
    if d < PLATFORM_ALONG and walk(a, b, 1.25 * d + 2.0) is not None: return True
    # services of one trunk whose stops landed either side of where they part (local and express
    # sides of a junction station)
    ra, rb = pats[ia[0]][0], pats[ib[0]][0]
    return d < PLATFORM_NEAR and ra != rb and (shared_near(a, ra, rb, 25.0) or shared_near(b, ra, rb, 25.0))
for sid, lst in stop_raw.items():
    groups = []
    for item in lst:
        hits = [gr for gr in groups if any(same_platform(item, o) for o in gr)]
        if not hits: groups.append([item]); continue
        hits[0].append(item)
        for gr in hits[1:]: hits[0].extend(gr); groups.remove(gr)
    for gr in groups:
        todo = list(gr)
        while todo:
            cands = {o[2] for o in todo}
            # also consider every node on the members' walks near the group's stops
            for o in todo:
                cands.update(n for n in pat_nodes[o[0]] if np.hypot(*(pos[n] - pos[o[2]])) < PLATFORM_JOIN)
            ctr = np.mean([pos[o[2]] for o in todo], 0)
            best = max(cands, key=lambda n: (sum(n in pat_set[o[0]] for o in todo), -np.hypot(*(pos[n] - ctr))))
            # routes that end or branch off just short of the platform: route them through it
            for o in todo:
                if best in pat_set[o[0]] or np.hypot(*(pos[best] - pos[o[2]])) >= PLATFORM_NEAR: continue
                seq = pat_nodes[o[0]]; j = seq.index(o[2])
                if j == 0: seq = [best] + seq[1:]                 # starts here: start on the shared platform
                elif j == len(seq) - 1: seq = seq[:-1] + [best]   # ends here: end on it
                else:
                    w1 = walk(seq[j - 1], best, 4 * PLATFORM_JOIN) or [seq[j - 1], best]
                    w2 = walk(best, seq[j + 1], 4 * PLATFORM_JOIN) or [best, seq[j + 1]]
                    seq = seq[:j - 1] + w1 + w2[1:] + seq[j + 2:]
                pat_nodes[o[0]] = seq; pat_set[o[0]] = set(seq)
            took = [o for o in todo if best in pat_set[o[0]]]
            for o in took: platform_of[(o[0], o[1])] = best
            platforms[sid].append(best)
            todo = [o for o in todo if best not in pat_set[o[0]]]
platform_station = {}
for sid, ns in platforms.items():
    for n in ns:
        if n in platform_station and platform_station[n] != sid:
            print('WARN node shared by stations', sid, platform_station[n])
        platform_station[n] = sid

# ---------------------------------------------------------------- 6. contract chains into edges
def rebuild_adj():
    adj.clear()
    for pi, seq in pat_nodes.items():
        rid = pats[pi][0]
        for x, y in zip(seq[:-1], seq[1:]):
            adj[x].setdefault(y, set()).add(rid); adj[y].setdefault(x, set()).add(rid)
# named points that the stretch rules (ORDER, GROUP, APART, ANCHOR) start or end at stay nodes of the graph
RULE_POINTS = {s[1:] for rule in ORDER + GROUP + APART + ANCHOR for s in rule[1:3] if s.startswith('@') and ',' not in s}
def is_key(n):
    if n in platform_station: return True
    if any(named.get(p) == n for p in RULE_POINTS): return True
    nb = list(adj[n])
    if len(nb) != 2: return True
    return adj[n][nb[0]] != adj[n][nb[1]]
def contract():
    key = {n for n in adj if is_key(n)}
    edges = []; edge_of = {}; seen = set()
    for a in sorted(key):
        for nb in sorted(adj[a]):
            if (a, nb) in seen: continue
            chain = [a, nb]
            while chain[-1] not in key:
                x = chain[-1]; nxt = [y for y in adj[x] if y != chain[-2]]
                chain.append(nxt[0])
            for u_, v_ in zip(chain[:-1], chain[1:]): seen.add((u_, v_)); seen.add((v_, u_))
            ei = len(edges)
            edges.append(dict(a=chain[0], b=chain[-1], chain=chain, routes=set(adj[chain[0]][chain[1]]),
                              length=float(sum(np.hypot(*(pos[y] - pos[x])) for x, y in zip(chain[:-1], chain[1:])))))
            edge_of[(chain[0], chain[1])] = (ei, 1)
            edge_of[(chain[-1], chain[-2])] = (ei, -1)
    return key, edges, edge_of
def replace_run(seq, old, new):
    # replace every occurrence of node run `old` (either direction) in seq with `new` (same direction)
    out = []; k = 0; n = len(old); ro = old[::-1]
    while k < len(seq):
        if seq[k:k + n] == old: out.extend(new); k += n
        elif seq[k:k + n] == ro: out.extend(new[::-1]); k += n
        else: out.append(seq[k]); k += 1
    dedup = [out[0]]
    for x in out[1:]:
        if x != dedup[-1]: dedup.append(x)
    return dedup
def despike_seq(seq):
    """Drop out-and-back steps and sharp reversals from a route's node walk.
    - (a, b, a) at a node that is not a platform: the step to b goes.
    - a route that runs up to a platform and straight back the same way (it joined its trunk just past the
      station): it now goes from its branch straight to the platform, so the junction moves to the station.
    - a sharp reversal at a node: a pass-through node goes; a platform is moved back onto the line."""
    while True:
        changed = False
        for k in range(1, len(seq) - 1):
            if seq[k - 1] != seq[k + 1] or seq[k] not in platform_station: continue
            m = 1
            while k - m - 1 >= 0 and k + m + 1 < len(seq) and seq[k - m - 1] == seq[k + m + 1]: m += 1
            seq = seq[:k - m] + seq[k:]; changed = True; break
        out = [seq[0]]
        for n in seq[1:]:
            if len(out) >= 2 and n == out[-2] and out[-1] not in platform_station: out.pop(); continue
            if n != out[-1]: out.append(n)
        k = 1; changed = changed or len(out) != len(seq)
        while k < len(out) - 1:
            a, b, c = pos[out[k - 1]], pos[out[k]], pos[out[k + 1]]
            u, v = b - a, c - b; nu, nv = np.hypot(*u), np.hypot(*v)
            if nu > 1e-9 and nv > 1e-9 and (u @ v) / (nu * nv) < -0.5:
                if out[k] not in platform_station:
                    del out[k]; changed = True; k = max(1, k - 1); continue
                ac = c - a; L2 = float(ac @ ac)          # platform sitting past a hairpin: put it back on the line
                t = 0.5 if L2 < 1e-9 else min(1.0, max(0.0, float((b - a) @ ac) / L2))
                if out[k] not in FIXED: pos[out[k]] = a + ac * t; changed = True
            k += 1
        seq = out
        if not changed: return seq
# ---------------------------------------------------------------- 5b. hand placement: shared track
def stop_node(pi, station):
    # platform node where pattern pi stops at a station given by display name or id
    for k, sid in enumerate(pats[pi][1]['stops']):
        if sid == station or stations[sid]['name'] == station: return platform_of[(pi, k)]
    return None
def between(seq, a, b):
    i, j = seq.index(a), seq.index(b)
    return (i, j, False) if i <= j else (j, i, True)
for rid, follow, sa, sb in SHARE:
    src_seq = None
    for pi, (r, _) in enumerate(pats):
        if r != follow: continue
        a, b = stop_node(pi, sa), stop_node(pi, sb)
        if a in pat_nodes[pi] and b in pat_nodes[pi]:
            i, j, _ = between(pat_nodes[pi], a, b); src_seq = (a, b, pat_nodes[pi][i:j + 1] if pat_nodes[pi].index(a) <= pat_nodes[pi].index(b) else pat_nodes[pi][i:j + 1][::-1]); break
    if not src_seq: print('SHARE: no track for', follow, sa, sb); continue
    a, b, run = src_seq; n_done = 0
    for pi, (r, _) in enumerate(pats):
        if r != rid: continue
        seq = pat_nodes[pi]
        if a not in seq or b not in seq: print('SHARE:', rid, 'does not stop at both', sa, sb); continue
        i, j, rev = between(seq, a, b)
        pat_nodes[pi] = seq[:i] + (run[::-1] if rev else run) + seq[j + 1:]; n_done += 1
    print('SHARE', rid, 'on', follow, sa, '->', sb, 'patterns', n_done)

# ---------------------------------------------------------------- 5c. hand placement: pins and walks
FIXED = set()                            # nodes placed by hand: later passes leave them where they are
def walk_node(pi, station):
    # node of a station on pattern pi's walk: where it stops, or a platform it passes through;
    # '@name' is a named point (POINTS), '@x,y' the walk's node nearest those coordinates
    if station.startswith('@'):
        if ',' in station:
            q = np.array([float(v) for v in station[1:].split(',')])
            return min(pat_nodes[pi], key=lambda n: float(np.hypot(*(pos[n] - q))))
        return named.get(station[1:])
    n = stop_node(pi, station)
    if n is not None: return n
    for sid, ns in platforms.items():
        if sid == station or stations[sid]['name'] == station:
            hit = [x for x in ns if x in pat_set_now(pi)]
            if hit: return hit[0]
    return None
def pat_set_now(pi): return set(pat_nodes[pi])
def patterns_of(rid): return [pi for pi, (r, _) in enumerate(pats) if lane_of(r) == lane_of(rid) and (r == rid or LANE_OF.get(r) == rid)]
for st, keep_r, other_r in PLATFORM:     # one platform node for services stopping side by side
    keep = next((stop_node(pi, st) for pi in patterns_of(keep_r) if stop_node(pi, st) is not None), None)
    for pi in patterns_of(other_r):
        for k, sid in enumerate(pats[pi][1]['stops']):
            if sid != st and stations[sid]['name'] != st: continue
            old_n = platform_of[(pi, k)]
            if old_n == keep: continue
            platform_of[(pi, k)] = keep
            pat_nodes[pi] = [keep if n == old_n else n for n in pat_nodes[pi]]
            if old_n in platforms[sid]: platforms[sid].remove(old_n)
            platform_station.pop(old_n, None)
    print('PLATFORM', st, other_r, 'stops at', keep_r + "'s platform")
for rid, st, x, y in MOVE:
    done_ = False
    for pi in patterns_of(rid):
        n = stop_node(pi, st)
        if n is None: continue
        pos[n] = [pos[n][0] if x is None else x, pos[n][1] if y is None else y]; FIXED.add(n); done_ = True
    if not done_: print('MOVE: no stop for', rid, st)
named = {}
def point_node(name):
    global pos
    if name not in named:
        pos = np.vstack([pos, np.array(POINTS[name], float)]); named[name] = len(pos) - 1; FIXED.add(named[name])
    return named[name]
for rid, sa, sb, mids in WALK:
    n_done = 0
    for pi in patterns_of(rid):
        seq = pat_nodes[pi]; a, b = walk_node(pi, sa), walk_node(pi, sb)
        if a is None or b is None or a not in seq or b not in seq: continue
        mid = []
        for m in mids:
            node = point_node(m[1:]) if m.startswith('@') else walk_node(pi, m)
            if node is None:
                cand = [x for sid, ns in platforms.items() if sid == m or stations[sid]['name'] == m for x in ns]
                node = cand[0] if cand else None
            if node is None: print('WALK: unknown stop', m); break
            mid.append(node)
        else:
            i, j = seq.index(a), seq.index(b)
            pat_nodes[pi] = seq[:i] + [a] + mid + [b] + seq[j + 1:] if i < j else seq[:j] + [b] + mid[::-1] + [a] + seq[i + 1:]
            n_done += 1
    print('WALK', rid, sa, '->', sb, 'patterns', n_done)

JMERGE = 16.0       # junctions closer than this along an edge become one junction
JPLAT = 6.0         # a junction this close to a platform moves onto it
for it in range(20):
    for pi in pat_nodes: pat_nodes[pi] = despike_seq(pat_nodes[pi])
    pat_set = {pi: set(seq) for pi, seq in pat_nodes.items()}
    rebuild_adj(); key, edges, edge_of = contract()
    changed = 0
    groups = collections.defaultdict(list)
    for ei, e in enumerate(edges):
        if e['a'] != e['b']: groups[tuple(sorted((e['a'], e['b'])))].append(ei)
    for (x, y), eis in groups.items():
        if len(eis) < 2: continue
        keep = max(eis, key=lambda ei: (len(edges[ei]['routes']), -edges[ei]['length']))
        kc = edges[keep]['chain'] if edges[keep]['a'] == x else edges[keep]['chain'][::-1]
        for ei in eis:
            if ei == keep: continue
            oc = edges[ei]['chain'] if edges[ei]['a'] == x else edges[ei]['chain'][::-1]
            for pi in pat_nodes: pat_nodes[pi] = replace_run(pat_nodes[pi], oc, kc)
            changed += 1
    merge = {}; collapse_runs = []
    for e in edges:
        a, b = e['a'], e['b']
        if a == b: continue
        plats = (a in platform_station) + (b in platform_station)
        if e['length'] >= (JMERGE if plats == 0 else JPLAT if plats == 1 else 0): continue
        if a in merge or b in merge or a in merge.values() or b in merge.values(): continue
        # never if some route reaches both ends without using this edge: it would loop back on itself
        ch = e['chain']
        def via_edge(seq):
            i, j = seq.index(a), seq.index(b)
            run = seq[min(i, j):max(i, j) + 1]
            return run == ch or run == ch[::-1]
        if any(a in pat_set[pi] and b in pat_set[pi] and not via_edge(seq) for pi, seq in pat_nodes.items()): continue
        if b in platform_station: a, b = b, a; ch = ch[::-1]
        elif plats == 0: pos[a] = (pos[a] + pos[b]) / 2      # two junctions become one, halfway
        merge[b] = a; collapse_runs.append(ch if ch[0] == a else ch[::-1])
    if merge:
        for pi in pat_nodes:
            for ch in collapse_runs: pat_nodes[pi] = replace_run(pat_nodes[pi], ch, [ch[0]])
        for pi, seq in pat_nodes.items():
            out = []
            for n in seq:
                n = merge.get(n, n)
                if not out or out[-1] != n: out.append(n)
            pat_nodes[pi] = out
        changed += len(merge)
    if not changed: break
    print('cleanup pass', it, 'bubbles + stubs folded', changed)
rebuild_adj(); key, edges, edge_of = contract()
for e in edges: e.pop('length')

def dp(pts, tol):
    if len(pts) < 3: return pts
    a, b = np.array(pts[0]), np.array(pts[-1]); ab = b - a; L = np.hypot(*ab)
    best, bi = -1, 0
    for i in range(1, len(pts) - 1):
        q = np.array(pts[i]) - a
        d = abs(ab[0] * q[1] - ab[1] * q[0]) / L if L > 1e-9 else np.hypot(*q)
        if d > best: best, bi = d, i
    if best <= tol: return [pts[0], pts[-1]]
    return dp(pts[:bi + 1], tol)[:-1] + dp(pts[bi:], tol)

# lateral order of lanes on each edge, from the source lane positions: project every sample onto the
# nearest edge segment that carries its route and read its lane offset in that edge's frame
lat_acc = [collections.defaultdict(list) for _ in edges]
eg = Grid(10.0); eseg = []
for ei, e in enumerate(edges):
    for k in range(len(e['chain']) - 1):
        pa, pb = pos[e['chain'][k]], pos[e['chain'][k + 1]]
        si = len(eseg); eseg.append((ei, pa, pb))
        for t in np.linspace(0, 1, max(2, int(np.hypot(*(pb - pa)) / 5) + 2)): eg.add(si, pa + (pb - pa) * t)
for i in range(len(P)):
    best = None
    for si in set(eg.near(C[i], 4.0)):
        ei, pa, pb = eseg[si]
        if ROUTE[i] not in edges[ei]['routes']: continue
        ab = pb - pa; L2 = float(ab @ ab)
        if L2 < 1e-12: continue
        t = min(1.0, max(0.0, float((C[i] - pa) @ ab) / L2)); dd = float(np.hypot(*(pa + ab * t - C[i])))
        if dd < 4.0 and (best is None or dd < best[0]): best = (dd, ei, ab)
    if best is None: continue
    sgn = 1.0 if best[2] @ U[i] >= 0 else -1.0
    lat_acc[best[1]][lane_of(ROUTE[i])].append(sgn * LAT[i])
rank = {r['id']: k for k, r in enumerate(routes)}
TIE = 1.0           # lanes this close in the source share one drawn lane, i.e. one track
for ei, e in enumerate(edges):
    lanes = sorted({lane_of(r) for r in e['routes']})
    med = {l: (float(np.median(lat_acc[ei][l])) if lat_acc[ei].get(l) else 0.0) for l in lanes}
    lanes.sort(key=lambda l: (med[l], rank.get(l, 0)))
    ties = []
    for l in lanes:
        if ties and med[l] - med[ties[-1][-1]] < TIE: ties[-1].append(l)
        else: ties.append([l])
    e['lanes'] = lanes; e['ties'] = ties

# Lanes sharing a source lane have no order of their own. Pick one and carry it across nodes so it
# stays the same whichever way each edge is stored.
def travel(ei, x, into):
    # e's lanes left to right when travelling into node x (into=True) or out of it
    e = edges[ei]; fwd = (e['b'] == x) if into else (e['a'] == x)
    return e['lanes'] if fwd else e['lanes'][::-1]
node_edges = collections.defaultdict(list)
for ei, e in enumerate(edges):
    node_edges[e['a']].append(ei); node_edges[e['b']].append(ei)
# 1. A plain stretch (edges joined by pass-through nodes, same lanes throughout) keeps one lane order:
#    the order of its longest edge. The source sometimes swaps two lanes mid-stretch; that is dropped.
elen = [float(sum(np.hypot(*(pos[y] - pos[x])) for x, y in zip(e['chain'][:-1], e['chain'][1:]))) for e in edges]
def plain(x):
    es = node_edges[x]
    return len(es) == 2 and es[0] != es[1] and set(edges[es[0]]['lanes']) == set(edges[es[1]]['lanes'])
sid_of = {}; stretches = []
for ei in sorted(range(len(edges)), key=lambda k: -elen[k]):
    if ei in sid_of: continue
    sid = len(stretches); members = [ei]; sid_of[ei] = sid; stack = [ei]
    while stack:
        e1 = stack.pop()
        for x in (edges[e1]['a'], edges[e1]['b']):
            if not plain(x): continue
            e2 = [y for y in node_edges[x] if y != e1][0]
            if e2 in sid_of: continue
            inbound = list(travel(e1, x, True))
            edges[e2]['lanes'] = inbound if edges[e2]['a'] == x else inbound[::-1]
            sid_of[e2] = sid; members.append(e2); stack.append(e2)
    stretches.append(members)
stretch_ties = []
for members in stretches:
    gid = {l: l for l in edges[members[0]]['lanes']}
    def root(l):
        while gid[l] != l: l = gid[l]
        return l
    for m in members:
        for t in edges[m]['ties']:
            for l in t[1:]: gid[root(l)] = root(t[0])
    groups = collections.defaultdict(set)
    for l in gid: groups[root(l)].add(l)
    stretch_ties.append(list(groups.values()))
# 2. Lanes sharing a source lane have no order of their own. Carry one across junctions from stretch to
#    stretch by swapping names within each tied group, which keeps every stretch consistent.
done = set()
for s0 in range(len(stretches)):
    if s0 in done: continue
    done.add(s0); q = collections.deque([s0])
    while q:
        s1 = q.popleft()
        for e1 in stretches[s1]:
            for x in (edges[e1]['a'], edges[e1]['b']):
                for e2 in node_edges[x]:
                    s2 = sid_of[e2]
                    if s2 in done: continue
                    done.add(s2); q.append(s2)
                    pos_in = {l: k for k, l in enumerate(travel(e1, x, True))}
                    out = travel(e2, x, False); swap = {}
                    for t in stretch_ties[s2]:
                        if len(t) < 2: continue
                        cur = [l for l in out if l in t]
                        swap.update(zip(cur, sorted(cur, key=lambda l: (pos_in.get(l, 1e9), rank.get(l, 0)))))
                    if swap:
                        for m in stretches[s2]: edges[m]['lanes'] = [swap.get(l, l) for l in edges[m]['lanes']]
# 3. Services that carry on from a trunk into a branch keep the trunk's left-to-right order, so lanes only
#    cross where the routes really have to. Trunks are handled largest first; a branch stretch is reordered
#    as a whole (by swapping lane names, which keeps it consistent along its length) at most once.
through = collections.defaultdict(set)        # (node, edge, edge) -> lanes that run from one edge onto the other
for pi, seq in pat_nodes.items():
    walk_e = []; k = 0
    while k < len(seq) - 1:
        if seq[k] in key and (seq[k], seq[k + 1]) in edge_of:
            ei, dr = edge_of[(seq[k], seq[k + 1])]; walk_e.append((ei, seq[k])); k += len(edges[ei]['chain']) - 1
        else: k += 1
    for (e1, _), (e2, x) in zip(walk_e, walk_e[1:]):
        through[(x, e1, e2)].add(lane_of(pats[pi][0])); through[(x, e2, e1)].add(lane_of(pats[pi][0]))
reordered = set()
def walk_edges(seq):
    out = []; k = 0
    while k < len(seq) - 1:
        if seq[k] in key and (seq[k], seq[k + 1]) in edge_of:
            ei, dr = edge_of[(seq[k], seq[k + 1])]; out.append((ei, dr)); k += len(edges[ei]['chain']) - 1
        else: k += 1
    return out
for rid, sa, sb, want in ORDER:     # hand-set lane orders (layout_def.ORDER) go first and are not changed
    for pi, (r, _) in enumerate(pats):
        if r != rid: continue
        a, b = walk_node(pi, sa), walk_node(pi, sb)
        walk = walk_edges(pat_nodes[pi]); nodes = [edges[walk[0][0]]['a'] if walk[0][1] > 0 else edges[walk[0][0]]['b']]
        for ei, dr in walk: nodes.append(edges[ei]['b'] if dr > 0 else edges[ei]['a'])
        if a not in nodes or b not in nodes: continue
        i, j = nodes.index(a), nodes.index(b)
        sub = walk[i:j] if i < j else [(ei, -dr) for ei, dr in walk[j:i][::-1]]
        for ei, dr in sub:
            e = edges[ei]; cur = e['lanes'] if dr > 0 else e['lanes'][::-1]   # current order in the direction of travel
            travel_order = [l for l in want if l in cur] + [l for l in cur if l not in want]
            lanes = travel_order if dr > 0 else travel_order[::-1]
            swap = dict(zip(e['lanes'], lanes))
            for m in stretches[sid_of[ei]]: edges[m]['lanes'] = [swap.get(l, l) for l in edges[m]['lanes']]
            reordered.add(sid_of[ei])
        print('ORDER', rid, sa, '->', sb, 'edges', len(sub))
        break
for x in sorted(node_edges, key=lambda n: -max(len(edges[ei]['lanes']) for ei in node_edges[n])):
    es = sorted(set(node_edges[x]), key=lambda ei: -len(edges[ei]['lanes']))
    for i1, e1 in enumerate(es):
        for e2 in es[i1 + 1:]:
            s2 = sid_of[e2]
            if s2 in reordered or s2 == sid_of[e1]: continue
            common = through.get((x, e1, e2), set())
            if len(common) < 2: continue
            want = [l for l in travel(e1, x, True) if l in common]
            cur = [l for l in travel(e2, x, False) if l in common]
            if want == cur: continue
            swap = dict(zip(cur, want))
            for m in stretches[s2]: edges[m]['lanes'] = [swap.get(l, l) for l in edges[m]['lanes']]
            reordered.add(s2)
for ei, e in enumerate(edges):
    e['ties'] = [[l for l in e['lanes'] if l in t] for t in stretch_ties[sid_of[ei]]]

# ---------------------------------------------------------------- 7. tracks (express / local groups)
# Start from shared source lanes, then join neighbouring groups that make the same stops.
stops_at = collections.defaultdict(set)   # node -> lanes that stop there
for (pi, k), n in platform_of.items():
    stops_at[n].add(lane_of(pats[pi][0]))
def stretch(ei):
    # edges either side that carry the same lanes, through simple pass-through nodes
    ls = set(edges[ei]['lanes']); out = [ei]
    for start in (edges[ei]['a'], edges[ei]['b']):
        cur_e, cur_n = ei, start
        for _ in range(400):
            nxt = [x for x in node_edges[cur_n] if x != cur_e and set(edges[x]['lanes']) == ls]
            if len(node_edges[cur_n]) != 2 or len(nxt) != 1: break
            cur_e = nxt[0]; out.append(cur_e)
            cur_n = edges[cur_e]['b'] if edges[cur_e]['a'] == cur_n else edges[cur_e]['a']
    return out
fam = {r['id']: r['color'] for r in routes}
done = set(); blind = []
for ei, e in enumerate(edges):
    if ei in done: continue
    es = stretch(ei); done.update(es)
    inner = {n for x in es for n in (edges[x]['a'], edges[x]['b']) if len(node_edges[n]) == 2}
    plats = sorted(n for n in inner if n in platform_station and any(l in stops_at[n] for l in e['lanes']))
    if not plats:   # no stations inside the stretch: judge express and local by the stations at its ends
        ends = {n for x in es for n in (edges[x]['a'], edges[x]['b'])} - inner
        plats = sorted(n for n in ends if n in platform_station and any(l in stops_at[n] for l in e['lanes']))
    if not plats: blind.append(es)
    sig = {l: tuple(l in stops_at[n] for n in plats) for l in e['lanes']}
    # one grouping for the whole stretch: lanes tied anywhere along it share a track, then
    # neighbouring groups with the same stops join
    gid = {l: l for l in e['lanes']}
    def root(l):
        while gid[l] != l: l = gid[l]
        return l
    for x in es:
        for t in edges[x]['ties']:
            for l in t[1:]: gid[root(l)] = root(t[0])
    order = e['lanes']
    for l1, l2 in zip(order, order[1:]):
        if root(l1) != root(l2) and sig[l1] == sig[l2] and (any(sig[l1]) or fam[l1] == fam[l2]):
            gid[root(l2)] = root(l1)
    for x in es:
        tracks = []
        for l in edges[x]['lanes']:
            if tracks and root(tracks[-1][-1]) == root(l): tracks[-1].append(l)
            else: tracks.append([l])
        edges[x]['tracks'] = tracks
# A stretch with no stations at all says nothing about express and local: take the grouping of the biggest
# trunk it connects to that carries all its lanes, so lanes keep their spacing through the junctions.
for es in blind:
    ends = {n for x in es for n in (edges[x]['a'], edges[x]['b'])}
    lanes0 = set(edges[es[0]]['lanes'])
    nbrs = [y for n in ends for y in node_edges[n] if y not in es and lanes0 <= set(edges[y]['lanes'])]
    if not nbrs: continue
    big = max(nbrs, key=lambda y: len(edges[y]['lanes']))
    gid = {l: gi for gi, t in enumerate(edges[big]['tracks']) for l in t}
    for x in es:
        tracks = []
        for l in edges[x]['lanes']:
            if tracks and gid.get(tracks[-1][-1]) == gid.get(l): tracks[-1].append(l)
            else: tracks.append([l])
        edges[x]['tracks'] = tracks

# hand-set separate tracks (layout_def.APART): listed lanes never share a track group
for rid, sa, sb, apart in APART:
    # lanes: plain names (each its own track) or lists (lanes that may share a track with each other)
    side = {}
    for gi, g in enumerate(apart):
        for l in ([g] if isinstance(g, str) else g): side[l] = gi
    n_done = 0
    for pi in patterns_of(rid):
        a, b = walk_node(pi, sa), walk_node(pi, sb)
        walk = walk_edges(pat_nodes[pi]); nodes = [edges[walk[0][0]]['a'] if walk[0][1] > 0 else edges[walk[0][0]]['b']]
        for ei, dr in walk: nodes.append(edges[ei]['b'] if dr > 0 else edges[ei]['a'])
        if a not in nodes or b not in nodes: continue
        i, j = sorted((nodes.index(a), nodes.index(b)))
        for ei, dr in walk[i:j]:
            tracks = []
            for t in edges[ei]['tracks']:
                cur = []
                for l in t:
                    if cur and l in side and any(x in side and side[x] != side[l] for x in cur): tracks.append(cur); cur = []
                    cur.append(l)
                tracks.append(cur)
            if len(tracks) != len(edges[ei]['tracks']): n_done += 1
            edges[ei]['tracks'] = tracks
        break
    print('APART', rid, sa, '->', sb, apart, 'edges', n_done)

# hand-set track groups (layout_def.GROUP): listed lanes running side by side become one group
for rid, sa, sb, together in GROUP:
    n_done = 0
    for pi in patterns_of(rid):
        a, b = walk_node(pi, sa), walk_node(pi, sb)
        walk = walk_edges(pat_nodes[pi]); nodes = [edges[walk[0][0]]['a'] if walk[0][1] > 0 else edges[walk[0][0]]['b']]
        for ei, dr in walk: nodes.append(edges[ei]['b'] if dr > 0 else edges[ei]['a'])
        if a not in nodes or b not in nodes: continue
        i, j = sorted((nodes.index(a), nodes.index(b)))
        for ei, dr in walk[i:j]:
            e = edges[ei]; present = [l for l in e['lanes'] if l in together]
            if len(present) < 2: continue
            idx = sorted(e['lanes'].index(l) for l in present)
            if idx[-1] - idx[0] != len(idx) - 1: continue           # not side by side
            tracks, merged = [], None
            for t in e['tracks']:
                if any(l in together for l in t):
                    if merged is None: merged = list(t); tracks.append(merged)
                    else: merged.extend(t)
                else: tracks.append(t)
            e['tracks'] = tracks; e['regroup'] = True; n_done += 1
        break
    print('GROUP', rid, sa, '->', sb, together, 'edges', n_done)

# hand-set reference lanes (layout_def.ANCHOR): the lane runs on the centreline, others are offset from it
for rid, sa, sb, lane in ANCHOR:
    n_done = 0
    for pi in patterns_of(rid):
        a, b = walk_node(pi, sa), walk_node(pi, sb)
        walk = walk_edges(pat_nodes[pi]); nodes = [edges[walk[0][0]]['a'] if walk[0][1] > 0 else edges[walk[0][0]]['b']]
        for ei, dr in walk: nodes.append(edges[ei]['b'] if dr > 0 else edges[ei]['a'])
        if a not in nodes or b not in nodes: continue
        i, j = sorted((nodes.index(a), nodes.index(b)))
        for ei, dr in walk[i:j]: edges[ei]['anchor'] = lane; edges[ei]['anchor_from'] = a; n_done += 1
        break
    print('ANCHOR', rid, sa, '->', sb, lane, 'edges', n_done)

# ---------------------------------------------------------------- 7b. one corner per bend
def consolidate(pts):
    """Each bend becomes a single corner. A run of short segments between two longer straight pieces
    (a curve sampled from the source, or leftover wobble) is replaced by the point where those straight
    pieces meet; a sideways jog between parallel pieces becomes one straight connector. The renderer
    then draws one sweeping arc at each corner instead of rounding every little vertex."""
    pts = [list(p) for p in pts]
    while True:
        n = len(pts)
        if n < 4: return pts
        L = [math.dist(pts[i], pts[i + 1]) for i in range(n - 1)]
        done = True
        i = 1
        while i < n - 2:
            if L[i] >= BEND_SHORT: i += 1; continue
            j = i
            while j < n - 1 and L[j] < BEND_SHORT: j += 1
            # short segments i..j-1 join vertex i to vertex j; segments i-1 and j are the straight pieces
            if j >= n - 1: break
            A, B, Cc, Dd = (np.array(pts[k], float) for k in (i - 1, i, j, j + 1))
            u, v = B - A, Dd - Cc
            nu, nv = np.hypot(*u), np.hypot(*v)
            run = sum(L[i:j])
            cr = (u[0] * v[1] - u[1] * v[0]) / (nu * nv)
            # a small sideways nudge between two nearly parallel straights: drop it, the straights meet directly
            if (u @ v) / (nu * nv) > 0.9 and abs(((Cc - B)[0] * u[1] - (Cc - B)[1] * u[0]) / nu) < JOG_MAX:
                pts = pts[:i] + pts[j + 1:]; done = False; break
            if j - i >= 1 and abs(cr) > 0.08:
                # meet the two straight pieces at their intersection
                t = ((Cc - B)[0] * v[1] - (Cc - B)[1] * v[0]) / (u[0] * v[1] - u[1] * v[0])
                X = B + u * t
                if np.hypot(*(X - B)) < 2 * run + 2 and np.hypot(*(X - Cc)) < 2 * run + 2 and (X - A) @ u > 0 and (Dd - X) @ v > 0:
                    pts = pts[:i] + [list(X)] + pts[j + 1:]; done = False; break
            elif j - i >= 2 and (u @ v) > 0:
                pts = pts[:i + 1] + pts[j:]; done = False; break   # a jog: one straight connector
            i = j + 1
        if done: return pts

# ---------------------------------------------------------------- 7c. straight runs
# Where a route's track stays within STRAIGHT_TOL of a straight line for at least STRAIGHT_MIN, it is made
# exactly straight: nodes, junctions and edge vertices along it are moved onto that line. This removes the
# small wanders the source's lane shuffles leave in a centreline. Longest runs first; nothing moves twice.
for e in edges:
    e['pts'] = dp([[float(pos[n][0]), float(pos[n][1])] for n in e['chain']], SIMPLIFY)
def walk_edges(seq):
    out = []; k = 0
    while k < len(seq) - 1:
        if seq[k] in key and (seq[k], seq[k + 1]) in edge_of:
            ei, dr = edge_of[(seq[k], seq[k + 1])]; out.append((ei, dr)); k += len(edges[ei]['chain']) - 1
        else: k += 1
    return out
def point_refs(walk):
    refs = []
    for ei, dr in walk:
        e = edges[ei]; n = len(e['pts'])
        idx = list(range(n)) if dr > 0 else list(range(n - 1, -1, -1))
        for m, i in enumerate(idx):
            ref = ('n', e['a']) if i == 0 else ('n', e['b']) if i == n - 1 else ('v', ei, i)
            if refs and m == 0: continue
            refs.append(ref)
    return refs
def ref_xy(r):
    return np.array(pos[r[1]], float) if r[0] == 'n' else np.array(edges[r[1]]['pts'][r[2]], float)
def set_ref(r, q):
    if r[0] == 'v': edges[r[1]]['pts'][r[2]] = [float(q[0]), float(q[1])]; return
    pos[r[1]] = q
    for ei in node_edges[r[1]]:
        e = edges[ei]
        if e['a'] == r[1]: e['pts'][0] = [float(q[0]), float(q[1])]
        if e['b'] == r[1]: e['pts'][-1] = [float(q[0]), float(q[1])]
runs = []
for pi, seq in pat_nodes.items():
    refs = point_refs(walk_edges(seq)); xy = [ref_xy(r) for r in refs]; n = len(refs); i = 0
    while i < n - 2:
        best = None; j = i + 2
        while j < n:
            a, b = xy[i], xy[j]; ab = b - a; L = float(np.hypot(*ab))
            if L < 1e-6: break
            dev = max(abs(ab[0] * (xy[k] - a)[1] - ab[1] * (xy[k] - a)[0]) / L for k in range(i + 1, j))
            # stay a line: every point within tolerance, and progressing along it
            if dev > STRAIGHT_TOL or any((xy[k + 1] - xy[k]) @ ab < -1e-6 for k in range(i, j)): break
            best = (j, dev, L); j += 1
        if best and best[2] >= STRAIGHT_MIN and best[1] > 0.3:
            runs.append((best[2], refs[i:best[0] + 1], xy[i], xy[best[0]]))
        i = best[0] if best else i + 1
moved = {('n', n) for n in FIXED}; nmoved = 0
for L, refs, a, b in sorted(runs, key=lambda r: -r[0]):
    ab = b - a; L2 = float(ab @ ab)
    for r in refs[1:-1]:
        if r in moved: continue
        q = ref_xy(r); t = float((q - a) @ ab) / L2
        set_ref(r, a + ab * t); moved.add(r); nmoved += 1
    moved.add(refs[0]); moved.add(refs[-1])
print('straight runs', len(runs), 'points moved', nmoved)
for e in edges:     # a short stub at either end that turns into its node: run straight into the node instead
    for _ in range(3):
        P_ = e['pts']
        for end in (0, -1):
            if len(P_) < 3: break
            n0, p1, p2 = (P_[0], P_[1], P_[2]) if end == 0 else (P_[-1], P_[-2], P_[-3])
            u = np.subtract(p1, n0); v = np.subtract(p2, p1); nu, nv = np.hypot(*u), np.hypot(*v)
            if nu < BEND_SHORT and nv > 1e-6 and nu > 1e-6 and (u @ v) / (nu * nv) < math.cos(math.radians(8)):
                del P_[1 if end == 0 else -2]
for e in edges:     # moving shared points can leave a vertex behind its neighbour: drop such reversals
    P_ = e['pts']; k = 1
    while k < len(P_) - 1:
        u = np.subtract(P_[k], P_[k - 1]); v = np.subtract(P_[k + 1], P_[k]); nu, nv = np.hypot(*u), np.hypot(*v)
        if nu < 0.05 or nv < 0.05 or (u @ v) / (nu * nv) < -0.2: del P_[k]; k = max(1, k - 1); continue
        k += 1

# ---------------------------------------------------------------- 7d. hand placement: paths
def place_path(rid, sa, sb, waypoints):
    for pi, (r, _) in enumerate(pats):
        if r != rid: continue
        a, b = stop_node(pi, sa), stop_node(pi, sb)
        if a is None or b is None: continue
        walk = walk_edges(pat_nodes[pi]); nodes = [edges[walk[0][0]]['a'] if walk[0][1] > 0 else edges[walk[0][0]]['b']]
        for ei, dr in walk: nodes.append(edges[ei]['b'] if dr > 0 else edges[ei]['a'])
        if a not in nodes or b not in nodes: continue
        i, j = nodes.index(a), nodes.index(b)
        sub = walk[i:j] if i < j else [(ei, -dr) for ei, dr in walk[j:i][::-1]]
        seq = nodes[i:j + 1] if i < j else nodes[j:i + 1][::-1]
        break
    else:
        print('PATH: no walk for', rid, sa, sb); return
    T = [np.array(pos[a], float)]
    snaps = []
    for w in waypoints:
        if w[0] == 'J':
            q = np.array(w[1:], float); T.append(q)
            jn = min((n for n in seq[1:-1] if n not in platform_station), key=lambda n: np.hypot(*(pos[n] - q)), default=None)
            if jn is not None: snaps.append((jn, q))
        else: T.append(np.array(w, float))
    T.append(np.array(pos[b], float))
    cum = [0.0]
    for p0, p1 in zip(T[:-1], T[1:]): cum.append(cum[-1] + float(np.hypot(*(p1 - p0))))
    def proj(q, k0):        # nearest point on T at or after segment k0: (arc length, segment, point)
        best = None
        for k in range(k0, len(T) - 1):
            p0, p1 = T[k], T[k + 1]; d_ = p1 - p0; L2 = float(d_ @ d_)
            t = 0.0 if L2 < 1e-12 else min(1.0, max(0.0, float((q - p0) @ d_) / L2))
            x = p0 + d_ * t; dd = float(np.hypot(*(x - q)))
            if best is None or dd < best[0] - 1e-9: best = (dd, cum[k] + t * float(np.sqrt(L2)), k, x)
        return best[1:]
    snapped = dict(snaps); s_of = {a: 0.0, b: cum[-1]}; k0 = 0
    for n in seq[1:-1]:
        if n in snapped:
            s, k0, x = proj(snapped[n], 0)
        else:
            s, k0, x = proj(np.array(pos[n], float), k0)
        set_ref(('n', n), x); s_of[n] = s
    for (ei, dr), u, v in zip(sub, seq[:-1], seq[1:]):
        s0, s1 = s_of[u], s_of[v]
        inner = [list(map(float, T[k])) for k in range(1, len(T) - 1) if s0 + 1e-6 < cum[k] < s1 - 1e-6]
        pts = [list(map(float, pos[u]))] + inner + [list(map(float, pos[v]))]
        edges[ei]['pts'] = pts if dr > 0 else pts[::-1]
        edges[ei]['locked'] = True
    print('PATH', rid, sa, '->', sb, 'nodes', len(seq), 'edges', len(sub))
for rid, sa, sb, wps in PATHS: place_path(rid, sa, sb, wps)

# ---------------------------------------------------------------- 8. output
key_list = sorted(key, key=lambda n: (round(pos[n][1], 1), round(pos[n][0], 1)))
nid = {n: f'n{k}' for k, n in enumerate(key_list)}
eid = {ei: f'e{ei}' for ei in range(len(edges))}
r2 = lambda v: round(float(v), 2)
node_out = []
for n in key_list:
    deg = len(adj[n])
    kind = 'platform' if n in platform_station else ('end' if deg == 1 else 'junction')
    o = dict(id=nid[n], x=r2(pos[n][0]), y=r2(pos[n][1]), kind=kind)
    if n in platform_station: o['station'] = platform_station[n]
    node_out.append(o)
edge_out = []
for ei, e in enumerate(edges):
    pinned = all(n in FIXED for n in e['chain'][1:-1]) and len(e['chain']) > 2
    pts = e['pts'] if e.get('locked') or pinned else consolidate(e['pts'])
    pts = [[r2(x), r2(y)] for x, y in pts]
    routes_ord = sorted(e['routes'], key=lambda r: (e['lanes'].index(lane_of(r)), r))
    edge_out.append(dict(id=eid[ei], **{'from': nid[e['a']], 'to': nid[e['b']]}, path=pts,
                         lanes=e['lanes'], tracks=e['tracks'], routes=routes_ord))
    if e.get('anchor'): edge_out[-1]['anchor'] = e['anchor']; edge_out[-1]['anchorFrom'] = nid[e['anchor_from']]
    if e.get('regroup'): edge_out[-1]['regroup'] = True
REG = {'always', 'interchange'}
route_out = []
for r in routes:
    pout = []
    for pi, (rid, p) in enumerate(pats):
        if rid != r['id']: continue
        seq = pat_nodes[pi]
        walk = []; k = 0
        while k < len(seq) - 1:
            if seq[k] not in key: k += 1; continue
            ei, dr = edge_of[(seq[k], seq[k + 1])]
            walk.append([eid[ei], dr]); k += len(edges[ei]['chain']) - 1
        styles = []
        for sid, vi in zip(p['stops'], p['stopPathIndex']):
            m = nearest_marker(stations[sid], p['path'][vi])
            styles.append('regular' if m['style'] in REG else m['style'].replace('interchange-part-time', 'part-time'))
        stop_nodes = [nid[platform_of[(pi, k)]] for k in range(len(p['stops']))]
        pout.append(dict(id=p['id'], stops=p['stops'], stopNodes=stop_nodes, stopStyles=styles, edges=walk))
    rr = {k: v for k, v in r.items() if k != 'patterns'}
    rr['lane'] = lane_of(r['id']); rr['patterns'] = pout
    route_out.append(rr)
st_out = []
for s in ref['stations']:
    o = {k: s[k] for k in ('id', 'name', 'system', 'x', 'y', 'routes')}
    o['platforms'] = [nid[n] for n in platforms.get(s['id'], [])]
    if s.get('note'): o['note'] = s['note']
    st_out.append(o)
transfers = [dict(stations=c['stations'], type=c['type']) for c in ref['connectors'] if c['type'] != 'in-station' and len(c['stations']) == 2]
for t in ref['transfers']:
    if not any(set(t['stations']) == set(x['stations']) for x in transfers): transfers.append(t)
others = [dict(id=b['id'], system=b['system'], color=b['color'], lines=b['lines'])
          for b in ref['bands'] if b['system'] in ('path', 'airtrain')]
meta = dict(ref['meta'])
meta['title'] = 'New York City Subway diagram (4D) - line graph'
meta['notes'] = [
    'One centreline per trunk: services that share a line are lanes on one edge, not separate geometry.',
    '`nodes` are platforms (one per station per trunk), junctions and line ends. `edges` run between them.',
    'Edge `lanes` are ordered across the edge, left to right when walking `from` to `to` with y down. '
    '`tracks` groups adjacent lanes that make the same stops (express vs local). Variants such as <6> share their parent lane.',
    'Route `patterns[].edges` is the walk as [edge id, +1 forward / -1 reverse]; `stopNodes[i]` is the platform for `stops[i]`.',
    '`stopStyles` is regular or a part-time style from the source legend. Stop sequences reflect regular weekday daytime service.',
    'The source-faithful lanes and markers are in data/reference/nyc-subway-4d.json.']
data = dict(meta=meta, routes=route_out, stations=st_out, nodes=node_out, edges=edge_out,
            transfers=transfers, ferries=ref['ferries'], others=others)
json.dump(data, open(OUT, 'w'), indent=1)

# ---------------------------------------------------------------- diagnostics
kinds = collections.Counter(o['kind'] for o in node_out)
print('nodes', dict(kinds), 'edges', len(edge_out))
short = [(e['id'], round(sum(math.dist(a, b) for a, b in zip(e['path'][:-1], e['path'][1:])), 1)) for e in edge_out]
print('edges under 4 units:', [x for x in short if x[1] < 4])
pairs = collections.Counter(tuple(sorted((e['from'], e['to']))) for e in edge_out)
print('parallel edges (bubbles):', [p for p, c in pairs.items() if c > 1])
multi = {s['id']: len(s['platforms']) for s in st_out if len(s['platforms']) > 1}
print('stations with several platforms:', len(multi))
junc_near_plat = []
for o in node_out:
    if o['kind'] != 'junction': continue
    for q in node_out:
        if q['kind'] == 'platform' and math.hypot(o['x'] - q['x'], o['y'] - q['y']) < 4:
            junc_near_plat.append((o['id'], q['station'])); break
print('junctions within 4 of a platform:', junc_near_plat)
