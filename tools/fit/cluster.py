import json, numpy as np, collections
D=json.load(open('dots.json')); bands=json.load(open('bands.json'))
# dedupe
u={}
for d in D:
    key=(round(d['x']),round(d['y']))
    if key in u: u[key]['classes'].add(d['cls']); continue
    d['classes']={d['cls']}; u[key]=d
dots=list(u.values())
for k,d in enumerate(dots): d['k']=k
P=np.array([[d['x'],d['y']] for d in dots])
par=list(range(len(dots)))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def union(a,b): par[find(a)]=find(b)
TH=11.5
for a in range(len(dots)):
    dd=np.hypot(*(P-P[a]).T)
    for b in np.where((dd<TH)&(dd>0))[0]: union(a,int(b))
def pdist(poly,pt):
    best=1e9
    for p,q in zip(poly[:-1],poly[1:]):
        p=np.array(p);q=np.array(q);ab=q-p;t=np.clip(np.dot(pt-p,ab)/max(np.dot(ab,ab),1e-9),0,1)
        best=min(best,np.hypot(*(p+ab*t-pt)))
    return best
conn=[b for b in bands if b['col']=='connector' and not (2300<=b['i']<2700)]
links=[]
for b in conn:
    for s in b['subs']:
        hit=[k for k in range(len(dots)) if pdist(s,P[k])<6]
        for h in hit[1:]: union(hit[0],h)
        if len(hit)<2: print('connector with <2 dots',b['id'],[round(v) for v in s[0]],[round(v) for v in s[-1]],hit)
dashed=[b for b in bands if b['col']=='grey999' and not (2300<=b['i']<2700)]
cl=collections.defaultdict(list)
for k in range(len(dots)): cl[find(k)].append(k)
clusters=sorted(cl.values(),key=lambda ks:(min(P[k][1] for k in ks)))
out=[]
for ci,ks in enumerate(clusters):
    for k in ks: dots[k]['c']=ci
    out.append(dict(c=ci,dots=ks,x=float(P[ks,0].mean()),y=float(P[ks,1].mean())))
for d in dots: d['classes']=sorted(d['classes'])
json.dump(dict(dots=dots,clusters=out),open('clusters.json','w'))
print('dots',len(dots),'clusters',len(out))
print(collections.Counter(len(c['dots']) for c in out))
