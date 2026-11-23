import json, numpy as np
from svgelements import Path, Move, Line, Close
f=json.load(open('flat.json'))
COL={'#ee352e':'red123','#00933c':'green456','#b933ad':'purple7','#0039a6':'blueACE','#ff6319':'orangeBDFM','#6cbe45':'limeG','#996633':'brownJZ','#a7a9ac':'greyL','#fccc0a':'yellowNQRW','#fecb00':'yellow2','#808183':'shuttleS','#006bb6':'SIR','#e00034':'e00034','#c3a3cd':'PATH','#03705f':'airtrain','#4d4d4d':'connector','#999999':'grey999','#5696c1':'5696c1'}
def poly(d):
    p=Path(d); subs=[]; cur=[]
    for seg in p:
        if isinstance(seg,Move):
            if len(cur)>1: subs.append(cur)
            cur=[(seg.end.x,seg.end.y)]
        elif isinstance(seg,Close):
            if cur: cur.append(cur[0])
        else:
            n=1 if isinstance(seg,Line) else 8
            for k in range(1,n+1):
                q=seg.point(k/n); cur.append((q.x,q.y))
    if len(cur)>1: subs.append(cur)
    return subs
bands=[]
for e in f:
    if e['fill'] is None and e['stroke'] in COL and (e['sw'] or 0)>=4:
        subs=poly(e['d'])
        L=sum(np.hypot(*np.diff(np.array(s),axis=0).T).sum() for s in subs)
        bands.append(dict(i=e['i'],id=e['id'],cls=e['cls'],col=COL[e['stroke']],sw=e['sw'],subs=subs,len=L))
json.dump(bands,open('bands.json','w'))
for b in bands:
    s=b["subs"]
    if not s: print(b["i"],b["id"],"EMPTY",b["col"]); continue
    print(b['i'],b['id'],b['cls'],b['col'],b['sw'],'nsub',len(s),'len',round(b['len']),'start',[round(v) for v in s[0][0]],'end',[round(v) for v in s[-1][-1]])
