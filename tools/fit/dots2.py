import json, numpy as np, collections, re
f=json.load(open('flat.json')); bands=json.load(open('bands.json'))
segs=[]
for bi,b in enumerate(bands):
    if b['col'] in ('connector','grey999','e00034','5696c1'): continue
    if 2300<=b['i']<2700: continue
    for s in b['subs']:
        for p,q in zip(s[:-1],s[1:]): segs.append((bi,*p,*q))
S=np.array(segs)
def near(x,y):
    a=S[:,1:3]; b=S[:,3:5]; ab=b-a; t=np.clip(((np.array([x,y])-a)*ab).sum(1)/np.maximum((ab*ab).sum(1),1e-9),0,1)
    p=a+ab*t[:,None]; d=np.hypot(*(p-[x,y]).T); o=np.argsort(d)
    res=[]
    for k in o:
        bi=int(S[k,0])
        if bi not in [r[0] for r in res]: res.append((bi,float(d[k])))
        if len(res)==3: break
    return res
OK={'st1','st52','st49','st64','st63','st59','st62','st48','st60','st61','st45','st46','st58','st65','st66','st67','st68','st69'}
dots=[]
for e in f:
    pid=re.sub(r'[\d_-].*','',e['id'] or 'None')
    b=e['bbox'];w=b[2]-b[0];h=b[3]-b[1]
    if e['cls'] in OK and 5<=w<=7.5 and 5<=h<=7.5 and (pid in('circle','path')):
        x=(b[0]+b[2])/2;y=(b[1]+b[3])/2
        if 1950<x and 200<y<620: continue
        if pid=='path' and not (abs(w-h)<0.1 and min(abs(w-6),abs(w-7))<0.06 and min(abs(h-6),abs(h-7))<0.06):
            # make sure it's circle-like: has curve commands
            continue
            if not re.search('[CcAa]',e['d']): continue
        nb=near(x,y)
        dots.append(dict(k=len(dots),i=e['i'],x=round(x,2),y=round(y,2),cls=e['cls'],band=nb[0][0],bcol=bands[nb[0][0]]['col'],d=round(nb[0][1],2),alt=[(bands[a]['col'],bands[a]['id'],round(d,1)) for a,d in nb[1:]]))
print(len(dots))
print(collections.Counter(d['d']<2 for d in dots))
for d in dots:
    if d['d']>=2: print('FAR',d['k'],d['x'],d['y'],d['cls'],d['bcol'],d['d'])
json.dump(dots,open('dots.json','w'))
