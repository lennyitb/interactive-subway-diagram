from build import *
allkeys=set(stations)
for rid,r in R.items():
    pats=r.get('branches') or {'main':r['stops']}
    for pn,stops in pats.items():
        if stops=='AUTO': continue
        bad=[s for s in stops if s not in allkeys]
        if bad: print(rid,pn,'UNKNOWN',bad); continue
        path,idx,warn=run(rid,stops,{r['fam']})
        L=sum(math.dist(a,b) for a,b in zip(path[:-1],path[1:]))
        print(f'== {rid}/{pn}: {len(stops)} stops, len {L:.0f}')
        for w in warn: print('   ',w)
