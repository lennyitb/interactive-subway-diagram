import json, sys
from svgelements import SVG, Path, Shape, Group, Circle, Ellipse, Text
svg=SVG.parse('../../NYC_subway-4D.svg', reify=True)
out=[]
i=0
for e in svg.elements():
    if isinstance(e,(Group,SVG)) or not isinstance(e,Shape): continue
    try: bb=e.bbox()
    except Exception: bb=None
    if bb is None: continue
    out.append(dict(i=i,id=e.values.get('id'),cls=e.values.get('class'),tag=type(e).__name__,
        fill=str(e.fill) if e.fill is not None and e.fill.value is not None else None,
        stroke=str(e.stroke) if e.stroke is not None and e.stroke.value is not None else None,
        sw=e.stroke_width, bbox=[round(v,2) for v in bb], d=Path(e).d() if not isinstance(e,Path) else e.d()))
    i+=1
json.dump(out,open('flat.json','w'))
print(len(out))
