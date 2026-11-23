#!/bin/sh
# Re-fit the network from ../../NYC_subway-4D.svg. Needs: pip install numpy lxml svgelements
# names.txt maps cluster ids (from clusters_frozen.json) to station names; routes_def.py holds stop sequences.
set -e; cd "$(dirname "$0")"
[ -f flat.json ] || python3 flat.py      # flatten SVG shapes (slow, ~1 min)
python3 bands.py >/dev/null               # coloured line strokes
python3 dots2.py >/dev/null               # station dots
python3 check.py                          # route path-finding report
python3 export.py                         # writes ../../data/nyc-subway-diagram.json
