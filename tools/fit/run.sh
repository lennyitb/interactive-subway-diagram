#!/bin/sh
# Re-fit the network from ../../NYC_subway-4D.svg. Needs: pip install numpy lxml svgelements
# names.txt maps cluster ids (from clusters_frozen.json) to station names; routes_def.py holds stop sequences;
# trunks_def.py lists lines that are drawn side by side but must not be merged into one trunk.
set -e; cd "$(dirname "$0")"
[ -f flat.json ] || python3 flat.py      # flatten SVG shapes (slow, ~1 min)
python3 bands.py >/dev/null               # coloured line strokes
python3 dots2.py >/dev/null               # station dots
python3 check.py                          # route path-finding report
python3 export.py                         # writes ../../data/reference/nyc-subway-4d.json (one lane per service, as drawn)
python3 collapse.py                       # writes ../../data/nyc-subway-diagram.json (one centreline per trunk)
python3 check_graph.py                    # consistency checks on the line graph
