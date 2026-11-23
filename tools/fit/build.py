import json, re, heapq, math, collections, sys
import numpy as np
sys.path.insert(0,'.')
from routes_def import R
C=json.load(open('clusters_frozen.json')); B=json.load(open('bands.json'))
dots={d['k']:d for d in C['dots']}
# ---------- stations ----------
names={}
for ln in open('names.txt'):
    ln=ln.rstrip('\n')
    if ln.strip(): k,v=ln.split(' ',1); names[int(k)]=v
groups={}  # key -> list of dot ks
for c in C['clusters']:
    v=names[c['c']]
    if v=='X': continue
    ks=[k for k in c['dots'] if not (dots[k]['d']>=9.5 and dots[k]['cls']=='st1')]
    if v.startswith('='):
        tgt=names[int(v[1:])]; groups.setdefault(tgt,[]).extend(ks); continue
    if v.startswith('#'):  # '#x,y+x,y=Name;#rest=Name'
        used=set()
        for part in v.split(';'):
            sel,nm=part[1:].split('=')
            if sel=='rest': pick=[k for k in ks if k not in used]
            else:
                pts=[tuple(map(float,p.split(','))) for p in sel.split('+')]
                pick=[k for k in ks if any(abs(dots[k]['x']-x)<1.5 and abs(dots[k]['y']-y)<1.5 for x,y in pts)]
            used|=set(pick); groups.setdefault(nm,[]).extend(pick)
        continue
    if v.startswith('@'):
        for part in v.split(';'):
            col,nm=part[1:].split('=')
            groups.setdefault(nm,[]).extend(k for k in ks if dots[k]['bcol']==col)
        continue
    groups.setdefault(v,[]).extend(ks)
def slug(s): return re.sub(r'[^a-z0-9]+','-',s.lower().replace("'",'')).strip('-')
AIR={'AirTrain Jamaica','Lefferts Blvd (AirTrain)','Federal Circle','Terminal 1','Terminal 2','Terminal 4','Terminal 5','Terminal 7','Terminal 8'}
stations={}
for key,ks in groups.items():
    if not ks: print('EMPTY station',key); continue
    disp=key.split('|')[0]
    kind='path' if disp.startswith('PATH ') else 'airtrain' if key in AIR else 'subway'
    if kind=='path': disp=disp[5:]
    if kind=='airtrain' and key=='Lefferts Blvd (AirTrain)': disp='Lefferts Blvd'
    if kind=='airtrain' and key=='AirTrain Jamaica': disp='Jamaica'
    xs=[dots[k]['x'] for k in ks]; ys=[dots[k]['y'] for k in ks]
    sid=slug(key.replace('|','-'))
    if kind!='subway': sid=kind+'-'+slug(disp)
    stations[key]=dict(id=sid,name=disp,kind=kind,x=round(float(np.mean(xs)),2),y=round(float(np.mean(ys)),2),
        markers=[dict(x=dots[k]['x'],y=dots[k]['y'],band=B[dots[k]['band']]['id'],style={'st1':'regular','st52':'interchange','st49':'no-regular-service','st61':'path','st62':'airtrain'}.get(dots[k]['classes'][0] if len(dots[k]['classes'])==1 else [c for c in dots[k]['classes'] if c!='st1'][0],'part-time')) for k in ks],_ks=ks)
ids=collections.Counter(s['id'] for s in stations.values()); assert max(ids.values())==1,[i for i,c in ids.items() if c>1]
dot2st={k:key for key,s in stations.items() for k in s['_ks']}
for s in stations.values():
    if s['kind']=='subway' and all(B[dots[k]['band']]['col']=='SIR' for k in s['_ks']): s['kind']='sir'; 
# ---------- graph ----------
FAM={'red':{'red123'},'green':{'green456'},'purple':{'purple7'},'blue':{'blueACE'},'orange':{'orangeBDFM'},'lime':{'limeG'},'brown':{'brownJZ'},'greyL':{'greyL'},'yellow':{'yellowNQRW','yellow2'},'shuttle':{'shuttleS'},'sir':{'SIR'},'path':{'PATH'},'airtrain':{'airtrain'}}
col2fam={c:f for f,cs in FAM.items() for c in cs}
LEG=lambda b: 2300<=b['i']<2700
bands=[(bi,b) for bi,b in enumerate(B) if b['col'] in col2fam and not LEG(b) and b['subs']]
# polylines list: (bi, subidx, pts)
polys=[(bi,si,np.array(s,float)) for bi,b in bands for si,s in enumerate(b['subs'])]
def project(pts,p):
    a=pts[:-1]; ab=pts[1:]-a; L2=np.maximum((ab*ab).sum(1),1e-12)
    t=np.clip(((p-a)*ab).sum(1)/L2,0,1); q=a+ab*t[:,None]; d=np.hypot(*(q-p).T); j=int(d.argmin())
    return j,float(t[j]),q[j],float(d[j])
inserts=collections.defaultdict(list)   # poly idx -> [(segj,t,tag)]
extra_edges=[]
# dots
for k,d in dots.items():
    if k not in dot2st: continue
    bi=d['band']; best=None
    for pi,(b2,si,pts) in enumerate(polys):
        if b2!=bi: continue
        j,t,q,dist=project(pts,np.array([d['x'],d['y']]))
        if best is None or dist<best[3]: best=(pi,j,t,dist)
    if best and best[3]<12: inserts[best[0]].append((best[1],best[2],('dot',k)))
    fam0=col2fam.get(B[bi]['col'])
    for pi,(b2,si,pts) in enumerate(polys):
        if b2==bi or col2fam[B[b2]['col']]!=fam0: continue
        j,t,q,dist=project(pts,np.array([d['x'],d['y']]))
        if dist<22: inserts[pi].append((j,t,('vdot',k,pi)))
# junctions
for pi,(bi,si,pts) in enumerate(polys):
    fam=col2fam[B[bi]['col']]
    for end,ep in ((0,pts[0]),(1,pts[-1])):
        if len(pts)>2 and np.hypot(*(pts[0]-pts[-1]))<0.5: continue
        for pj,(bj,sj,pts2) in enumerate(polys):
            if pj==pi or col2fam[B[bj]['col']]!=fam: continue
            j,t,q,dist=project(pts2,ep)
            if dist<6:
                inserts[pj].append((j,t,('jn',pi,end)))
                extra_edges.append((('end',pi,end),('jn',pi,end,pj),dist))
# build nodes
nodes=[]; nid={}; adj=collections.defaultdict(list)
def node(key,xy):
    if key not in nid: nid[key]=len(nodes); nodes.append(dict(xy=(float(xy[0]),float(xy[1])),key=key))
    return nid[key]
def edge(a,b,w,fam):
    adj[a].append((b,w)); adj[b].append((a,w))
for pi,(bi,si,pts) in enumerate(polys):
    fam=col2fam[B[bi]['col']]
    seq=[]
    for j in range(len(pts)):
        seq.append((j,0.0,('v',pi,j),pts[j]))
    for (j,t,tag) in inserts[pi]:
        q=pts[j]+(pts[j+1]-pts[j])*t
        key=tag if tag[0] in('dot','vdot') else ('jn',tag[1],tag[2],pi)
        seq.append((j,t if t>0 else 1e-9,key,q))
    seq.sort(key=lambda s:(s[0],s[1]))
    ids_=[]
    for s in seq:
        k=s[2]
        if k==('v',pi,0): k=('end',pi,0)
        if k==('v',pi,len(pts)-1): k=('end',pi,1)
        n=node(k,s[3]); nodes[n]['fam']=fam; nodes[n]['band']=bi; ids_.append(n)
    for a,b in zip(ids_[:-1],ids_[1:]):
        if a!=b: edge(a,b,math.dist(nodes[a]['xy'],nodes[b]['xy']),fam)
for a,b,w in extra_edges:
    if a in nid and b in nid: edge(nid[a],nid[b],w,None)
print('graph nodes',len(nodes))
dotnode={k[1]:n for k,n in nid.items() if k[0]=='dot'}
vdotnodes=collections.defaultdict(list)
for k,n in nid.items():
    if k[0]=='vdot': vdotnodes[k[1]].append(n)
_f=json.load(open('flat.json'))
_c=[e for e in _f if (e['bbox'][2]-e['bbox'][0])<=7.5 and (e['bbox'][3]-e['bbox'][1])<=7.5 and ((e['fill']=='#000000' and e['cls'] is None) or e['cls']=='st80')]
_cx=np.array([[(e['bbox'][0]+e['bbox'][2])/2,(e['bbox'][1]+e['bbox'][3])/2] for e in _c])
SOFT=set()
for k,d in dots.items():
    if (np.hypot(*(_cx-[d['x'],d['y']]).T)<3.2).any() or set(d['classes'])-{'st1','st52'}: SOFT.add(k)
def pen(n,allowed):
    k=nodes[n]['key']
    if k[0]!='dot' or dot2st[k[1]] in allowed: return 0
    st=stations[dot2st[k[1]]]; d=dots[k[1]]
    return 40 if k[1] in SOFT else 150
def dijkstra(srcs,tgts,fam,allowed):
    dist={};prev={};pq=[]
    for n,c in srcs.items(): dist[n]=c; heapq.heappush(pq,(c,n))
    while pq:
        c,n=heapq.heappop(pq)
        if n==-1:
            path=[prev[-1]]
            while path[-1] in prev: path.append(prev[path[-1]])
            return c,path[::-1]
        if c>dist.get(n,1e18): continue
        if n in tgts:
            nc=c+tgts[n]
            if nc<dist.get(-1,1e18): dist[-1]=nc; prev[-1]=n; heapq.heappush(pq,(nc,-1))
        for m,w in adj[n]:
            if nodes[m]['fam'] not in fam: continue
            nc=c+w*(0.95 if nodes[m].get('band') in CURPREF else 1.0)+pen(m,allowed)
            if nc<dist.get(m,1e18): dist[m]=nc; prev[m]=n; heapq.heappush(pq,(nc,m))
    return None,None
def stnodes(key,fams):
    real=[dotnode[k] for k in stations[key]['_ks'] if k in dotnode and nodes[dotnode[k]]['fam'] in fams]
    virt=[n for k in stations[key]['_ks'] for n in vdotnodes.get(k,[]) if nodes[n]['fam'] in fams]
    return real,virt
PIN={}; CURPREF=set(); PREF={}
def run(rid,stops,fams):
    global CURPREF
    CURPREF={bi for bi,b in enumerate(B) if b['id'] in PREF.get(rid,())}
    path=[];stopidx=[];warn=[]
    cur=None
    for i,(a,b) in enumerate(zip(stops[:-1],stops[1:])):
        Ar,Av=stnodes(a,fams); Br,Bv=stnodes(b,fams)
        if (rid,b) in PIN:
            px,py=PIN[(rid,b)]; Br=[n for n in Br if math.dist(nodes[n]['xy'],(px,py))<3]; Bv=[]
        if not (Ar or Av) or not (Br or Bv): warn.append(f'no dots on family for {a if not (Ar or Av) else b}'); continue
        srcs={n:(0 if n==cur else 30 if cur is not None else 0) for n in Ar}
        for n in Av: srcs.setdefault(n,(0 if n==cur else 50))
        T={n:40 for n in Bv}; T.update({n:0 for n in Br})
        c,p=dijkstra({cur:0},T,fams,{a,b}) if cur is not None else (None,None)
        if p is None: c,p=dijkstra(srcs,T,fams,{a,b})
        elif cur is not None:
            c2,p2=dijkstra(srcs,T,fams,{a,b})
            if p2 is not None and c2+60<c: c,p=c2,p2
        if p is None: warn.append(f'NO PATH {a} -> {b}'); continue
        passed=[dot2st[nodes[n]['key'][1]] for n in p[1:-1] if nodes[n]['key'][0]=='dot' and dot2st[nodes[n]['key'][1]] not in (a,b)]
        passed=[s for s in passed]
        if passed: warn.append(f'{a} -> {b} passes {sorted(set(passed))} (cost {c:.0f})')
        seg=[nodes[n]['xy'] for n in p]
        if not path: stopidx.append(0); path.extend(seg)
        else:
            if math.dist(path[-1],seg[0])>0.01:
                warn.append(f'JUMP {math.dist(path[-1],seg[0]):.1f} at {a}'); path.append(seg[0])
            path.extend(seg[1:])
        stopidx.append(len(path)-1); cur=p[-1]
    return path,stopidx,warn
PIN[('E','Chambers St/WTC/Park Pl/Cortlandt St')]=(453,1733)
PREF['N']=['polyline2630']; PREF['R']=['polyline2636']; PREF['W']=['polyline2636']; PREF['Q']=['path43216']
PREF.update({'B':['polyline2672','polyline2674'],'D':['polyline2664'],'6':['polyline2646'],'6X':['polyline2696'],'7':['polyline2758'],'7X':['polyline2760'],'J':['polyline2698'],'Z':['polyline2742']})
