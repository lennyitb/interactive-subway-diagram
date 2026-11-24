"""Sanity checks for data/nyc-subway-diagram.json (the collapsed line graph)."""
import json, os, collections
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../data/nyc-subway-diagram.json')))
E = {e['id']: e for e in d['edges']}; N = {n['id']: n for n in d['nodes']}; S = {s['id']: s for s in d['stations']}
bad = 0; used = collections.defaultdict(set)
for r in d['routes']:
    for p in r['patterns']:
        nodes = []
        for eid, dr in p['edges']:
            e = E[eid]; a, b = (e['from'], e['to']) if dr == 1 else (e['to'], e['from'])
            if nodes and nodes[-1] != a: print('BROKEN WALK', r['id'], p['id'], eid); bad += 1
            if not nodes: nodes.append(a)
            nodes.append(b); used[eid].add(r['id'])
            if r['id'] not in e['routes']: print('ROUTE NOT ON EDGE', r['id'], eid); bad += 1
        k = 0
        for sid, sn in zip(p['stops'], p['stopNodes']):
            if N[sn].get('station') != sid: print('STOP NODE STATION MISMATCH', r['id'], sid, sn); bad += 1
            while k < len(nodes) and nodes[k] != sn: k += 1
            if k == len(nodes): print('STOP NOT ON WALK IN ORDER', r['id'], p['id'], sid); bad += 1; k = 0
for e in d['edges']:
    if set(e['routes']) != used[e['id']]: print('EDGE ROUTES MISMATCH', e['id'], e['routes'], sorted(used[e['id']])); bad += 1
    if sorted(sum(e['tracks'], [])) != sorted(e['lanes']): print('TRACKS != LANES', e['id']); bad += 1
import math
for e in d['edges']:
    P = e['path']
    for a, b, c in zip(P, P[1:], P[2:]):
        u = (b[0] - a[0], b[1] - a[1]); v = (c[0] - b[0], c[1] - b[1])
        nu, nv = math.hypot(*u), math.hypot(*v)
        if nu > 1e-6 and nv > 1e-6 and (u[0] * v[0] + u[1] * v[1]) / (nu * nv) < -0.5:
            print('EDGE DOUBLES BACK', e['id'], e['lanes'], b); bad += 1; break
for n in d['nodes']:
    if n['kind'] == 'end':
        es = [e for e in d['edges'] if n['id'] in (e['from'], e['to'])]
        other = [x for e in es for x in (e['from'], e['to']) if x != n['id']]
        print('tail end', n['id'], n['x'], n['y'], 'next to', [S[N[o]['station']]['name'] if N[o].get('station') else N[o]['kind'] for o in other])
# report (not a failure): wobbles, two direction changes within 15 units along a route that aren't one bend
wob = set()
for r in d['routes']:
    for pt in r['patterns']:
        P = []
        for eid, dr in pt['edges']:
            seg = E[eid]['path'] if dr > 0 else E[eid]['path'][::-1]
            P.extend(seg if not P else seg[1:])
        turns = []
        for i in range(1, len(P) - 1):
            u = (P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1]); v = (P[i + 1][0] - P[i][0], P[i + 1][1] - P[i][1])
            if math.hypot(*u) < 1e-6 or math.hypot(*v) < 1e-6: continue
            ang = math.degrees(math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]))
            if abs(ang) > 4: turns.append((i, ang))
        for (i, a1), (j, a2) in zip(turns, turns[1:]):
            if sum(math.dist(P[k], P[k + 1]) for k in range(i, j)) < 15 and (a1 * a2 < 0 or abs(a1) + abs(a2) < 60):
                wob.add((round(P[i][0] / 20), round(P[i][1] / 20)))
print('wobble spots (report only):', len(wob))
print('problems:', bad)
