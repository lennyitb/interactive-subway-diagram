from build import *
import numpy as np
f=json.load(open('flat.json'))
# ---- marker styles ----
cand=[e for e in f if (e['bbox'][2]-e['bbox'][0])<=7.5 and (e['bbox'][3]-e['bbox'][1])<=7.5]
cx=np.array([[(e['bbox'][0]+e['bbox'][2])/2,(e['bbox'][1]+e['bbox'][3])/2] for e in cand])
RING={'st64','st63','st59','st48','st60','st45','st69'}
def style(k):
    d=dots[k]; cl=set(d['classes'])
    dd=np.hypot(*(cx-[d['x'],d['y']]).T); near=[cand[j] for j in np.where(dd<3.2)[0]]
    if any(e['fill']=='#000000' and e['cls'] is None for e in near): return 'weekdays-only'
    if any(e['cls']=='st80' for e in near): return 'nights-weekends-only'
    if cl&{'st58','st65','st67','st68'}: return 'interchange-part-time'
    if 'st52' in cl: return 'interchange'
    if 'st49' in cl: return 'except-rush-peak'
    if cl&RING: return 'rush-peak-only'
    if 'st61' in cl: return 'path'
    if 'st62' in cl: return 'airtrain'
    return 'always'
# ---- Z pattern from diagram ----
zb=[bi for bi,b in enumerate(B) if b['id']=='polyline2742'][0]
Z=[]
for key in R['J']['stops']:
    ks=stations[key]['_ks']
    onz=any(dots[k]['band']==zb for k in ks)
    if onz or R['J']['stops'].index(key)>=R['J']['stops'].index('Marcy Av'): Z.append(key)
R['Z']['stops']=Z
# ---- routes ----
LABEL={'6X':'6','7X':'7','FX':'F','GS':'S','FS':'S','H':'S'}
routes=[]; bandroutes=collections.defaultdict(set); stroutes=collections.defaultdict(set)
for rid,r in R.items():
    pats=r.get('branches') or {rid:r['stops']}
    out=[]
    for pn,stops in pats.items():
        path,idx,warn=run(rid,stops,{r['fam']})
        assert len(idx)==len(stops),(rid,pn,warn)
        # record bands used: nearest band to midpoints
        keep=set(idx)
        simp=[]; sidx=[]
        for i,p in enumerate(path):
            if simp and math.dist(simp[-1],p)<0.05 and i not in keep:
                continue
            if len(simp)>=2 and i not in keep and (len(simp)-1) not in sidx_set if False else False: pass
            simp.append(p)
            if i in keep: sidx.append(len(simp)-1)
        # remove collinear interior points that aren't stops
        pts=simp; stopset=set(sidx); res=[];m=[]
        for i,p in enumerate(pts):
            if 0<i<len(pts)-1 and i not in stopset:
                a=np.array(res[-1]); b=np.array(p); c=np.array(pts[i+1])
                u=b-a;v=c-b
                if abs(u[0]*v[1]-u[1]*v[0])<1e-3*max(1,np.linalg.norm(c-a)) and np.dot(b-a,c-b)>=0: continue
            res.append(p); m.append(i)
        remap={old:new for new,old in enumerate(m)}
        sidx=[remap[i] for i in sidx]
        ids=[stations[s]['id'] for s in stops]
        for s in stops: stroutes[s].add(rid)
        out.append(dict(id=pn,stops=ids,stopPathIndex=sidx,path=[[round(x,2),round(y,2)] for x,y in res]))
    routes.append(dict(id=rid,label=LABEL.get(rid,rid),name=r['name'],color=r['color'],textColor=r.get('text','#FFFFFF'),
        bullet=r.get('bullet','circle'),system='sir' if rid=='SIR' else 'subway',patterns=out))
# which drawn bands each route uses: sample route paths and match to band polylines
polyarr=[(bi,np.array(s,float)) for bi,b in bands for s in b['subs']]
def nearest_band(p):
    best=(1e9,None)
    for bi,pts in polyarr:
        j,t,q,dd=project(pts,np.array(p))
        if dd<best[0]: best=(dd,bi)
    return best
for r in routes:
    fam=R[r['id']]['fam']
    for pat in r['patterns']:
        P=pat['path']
        for a,b in zip(P[:-1],P[1:]):
            L=math.dist(a,b); n=max(1,int(L//15))
            for t in np.linspace(0,1,n+2)[1:-1]:
                p=(a[0]+(b[0]-a[0])*t,a[1]+(b[1]-a[1])*t)
                dd,bi=nearest_band(p)
                if dd<1.0 and col2fam[B[bi]['col']]==fam: bandroutes[bi].add(r['id'])
# ---- stations ----
order=list(R)
st_out=[]
for key,s in sorted(stations.items(),key=lambda kv:(kv[1]['y'],kv[1]['x'])):
    rs=sorted(stroutes.get(key,[]),key=order.index)
    sys_='airtrain' if s['kind']=='airtrain' else 'path' if s['kind']=='path' else 'sir' if s['kind']=='sir' else 'subway'
    st_out.append(dict(id=s['id'],name=s['name'],system=sys_,x=s['x'],y=s['y'],routes=rs,
        markers=[dict(x=round(dots[k]['x'],2),y=round(dots[k]['y'],2),style=style(k),routes=sorted(bandroutes.get(dots[k]['band'],[]),key=order.index)) for k in s['_ks']]))
# notes
notes={'aqueduct-racetrack':'Served by the A northbound only.','75-av':'E stops here late nights only; F always.'}
for s in st_out:
    if s['id'] in notes: s['note']=notes[s['id']]
# ---- drawn bands ----
DASH={'st14','st17','st19','st25','st27','st28'}
band_out=[]
for bi,b in bands:
    if col2fam[B[bi]['col']] in ('airtrain','path'): rr=[]
    else: rr=sorted(bandroutes.get(bi,[]),key=order.index)
    band_out.append(dict(id=b['id'],system={'path':'path','airtrain':'airtrain','sir':'sir'}.get(col2fam[B[bi]['col']],'subway'),
        color=next(e['stroke'] for e in f if e['i']==b['i']).upper(),width=b['sw'],dashed=b['cls'] in DASH,routes=rr,
        lines=[[[round(x,2),round(y,2)] for x,y in s] for s in b['subs']]))
# ---- connectors / transfers ----
kpos={k:(dots[k]['x'],dots[k]['y']) for k in dot2st}
def st_near(p):
    best=(1e9,None)
    for k,(x,y) in kpos.items():
        dd=math.hypot(x-p[0],y-p[1])
        if dd<best[0]: best=(dd,dot2st[k])
    return best
conn=[]; transfers=[]
for bi,b in enumerate(B):
    if b['col'] not in ('connector','grey999') or LEG(b) or not b['subs']: continue
    if b['id'] in ('polyline2883',): continue
    for s in b['subs']:
        hits=set()
        for p in s:
            dd,key=st_near(p)
            if dd<7: hits.add(key)
        typ='part-time' if b['col']=='grey999' else 'out-of-station' if b['cls']=='st32' else 'in-station'
        ids_=sorted(stations[h]['id'] for h in hits)
        conn.append(dict(type=typ,stations=ids_,line=[[round(x,2),round(y,2)] for x,y in s]))
# merge part-time dashed pieces (3 pieces) by station pairs
tset={}
for c in conn:
    if len(c['stations'])==2: tset[tuple(c['stations'])]=c['type']
pt=[c for c in conn if c['type']=='part-time']
if pt:
    hits=set(); [hits.update(c['stations']) for c in pt]
    for c in pt:
        for p in c['line']:
            dd,key=st_near(p)
            if dd<12: hits.add(stations[key]['id'])
    if len(hits)==2: tset[tuple(sorted(hits))]='part-time'
transfers=[dict(stations=list(k),type=v) for k,v in sorted(tset.items())]
# Livonia/Junius etc. are covered by out-of-station connectors.
ferry=[b for b in B if b['id']=='polyline2883'][0]
data=dict(
  meta=dict(title='New York City Subway diagram (4D) - fitted network data',
    source='NYC_subway-4D.svg by CountZ (Jake Berman, 53studio.com) at English Wikipedia, via Wikimedia Commons (https://commons.wikimedia.org/wiki/File:NYC_subway-4D.svg), licensed CC BY-SA 3.0. Station positions, lines and stop sequences were fitted to it; this derived data is shared under CC BY-SA 3.0.',license='CC-BY-SA-3.0',
    coordinateSpace=dict(units='SVG user units of the source diagram',viewBox=[0,0,2500,2700],yAxis='down'),
    generated='2026-11-23',
    notes=['Stations are station complexes; `markers` are the individual dots drawn on each service band.',
      'Route `patterns[].path` follows the drawn bands; `stopPathIndex[i]` is the path vertex of `stops[i]`.',
      'Stop sequences reflect regular weekday daytime service. Marker `style` reflects the diagram legend (rush-peak-only, weekdays-only, etc).',
      '`bands` reproduce the source strokes (one per service group); parallel same-colour services are separate bands, as in the source.']),
  routes=routes, stations=st_out, transfers=transfers, connectors=conn, bands=band_out,
  ferries=[dict(id='staten-island-ferry',name='Staten Island Ferry',color='#5696C1',stations=['south-ferry-whitehall-st','st-george'],line=[[round(x,2),round(y,2)] for x,y in ferry['subs'][0]])])
import os
os.makedirs('../../data',exist_ok=True)
json.dump(data,open('../../data/nyc-subway-diagram.json','w'),indent=1)
print('routes',len(routes),'stations',len(st_out),'bands',len(band_out),'connectors',len(conn),'transfers',transfers)
print('Z:',Z)
print('unserved:',[s['id'] for s in st_out if not s['routes']])
print('bands with no routes:',[ (b['id'],b['color']) for b in band_out if not b['routes'] and b['system'] in('subway','sir')])
