# Interactive subway diagram data

NYC subway network data fitted to a schematic (non-geographic) diagram, for building a pan/zoom diagram explorer where lines and markers stay a constant on-screen size and labels are placed at runtime with collision avoidance.

## Contents

| Path | What it is |
| --- | --- |
| `data/nyc-subway-diagram.json` | The network as a line graph: one centreline per trunk, with the services on it as ordered lanes. Diagram coordinate space (viewBox `0 0 2500 2700`, y down) |
| `data/reference/nyc-subway-4d.json` | The source-faithful fit: one drawn lane per service and one marker dot per lane, as in the source diagram. Input to the line graph |
| `data/basemap.json` | Water, land, park and borough-boundary shapes as SVG path data in the same space |
| `preview.html` | Canvas map renderer: lanes offset from the centrelines at a constant on-screen pitch, platform pills, transfer links, greedy labels; pan, zoom and rotate |
| `tools/fit/` | Scripts that extracted and fitted the data from the source SVG, and collapse it to the line graph |
| `tools/fit/preview-reference.html` | Viewer for the reference fit (open `/tools/fit/preview-reference.html` from the same server) |
| `NYC_subway-4D.svg` | Source diagram (see attribution below) |

### Data model (`nyc-subway-diagram.json`)

Parallel services are not separate geometry. A trunk such as Queens Blvd is one centreline; the E, F, M and R are lanes on it, and the renderer decides how far apart to draw them at each zoom level.

- **nodes**: `platform` (one per station per trunk, with its `station`), `junction` (where lanes join or leave a trunk) and `end`. Each has `x`/`y` on the centreline.
- **edges**: a centreline `path` between two nodes (`from`, `to`). `lanes` lists the services on it across the edge, left to right when walking `from` to `to` (y down). `tracks` groups neighbouring lanes that share a track or make the same stops (express vs local), e.g. `[["E","F"],["M","R"]]`. `routes` also lists variants (`6X`, `7X`, `FX`), which ride their parent's lane.
- **routes**: 29 services with `color`, `label`, `bullet` and `lane` (the lane they are drawn in). Each pattern has `stops` (station ids), `stopNodes` (the platform node for each stop), `stopStyles` (`regular` or a part-time style from the source legend, such as `weekdays-only` or `rush-peak-only`) and `edges`, the walk as `[edge id, 1 forward / -1 reverse]`. Stop sequences reflect regular weekday daytime service.
- **stations**: complexes with `id`, `name`, `system`, an anchor `x`/`y`, the `routes` that stop there and their `platforms` (node ids). A complex with several platforms is an in-station transfer.
- **transfers**: out-of-station and part-time transfers between complexes.
- **ferries**, **others**: the Staten Island Ferry, and PATH and AirTrain strokes as drawn.

## Running the preview

```sh
python3 -m http.server 8765
# open http://127.0.0.1:8765/preview.html
```

Drag to pan, scroll to zoom, shift- or alt-scroll (or two-finger twist) to rotate. The camera is kept in the URL hash. The checkboxes overlay the land, the source lanes and the graph's nodes and centrelines.

## Re-fitting

Station names are in `tools/fit/names.txt` (cluster id to name) and stop sequences in `tools/fit/routes_def.py`. After editing either, run `tools/fit/run.sh`. It needs `numpy`, `lxml` and `svgelements`.

To rebuild only the line graph from the reference fit, run `python3 tools/fit/collapse.py` and then `python3 tools/fit/check_graph.py` (numpy only). `collapse.py` treats services drawn side by side at the source's lane pitch (10 units) as one trunk and moves them onto its centreline. Separate lines that happen to be drawn at that spacing, such as 6 Av and 8 Av at W 4 St, are kept apart by `NEVER_MERGE` in `tools/fit/trunks_def.py`; one line drawn wider apart than that, such as Eastern Pkwy, is joined by `ALWAYS_MERGE`.

## Attribution and license

The source diagram `NYC_subway-4D.svg` is by [CountZ (Jake Berman)](https://en.wikipedia.org/wiki/User:CountZ) at English Wikipedia ([53studio.com](https://53studio.com)), via [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:NYC_subway-4D.svg), licensed [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).

The data in `data/` is derived from it (positions, geometry and shapes extracted and fitted; names and stop sequences added) and is shared under the same **CC BY-SA 3.0** license.

The code (`preview.html`, `tools/`) is released under the MIT License.

This project is not affiliated with or endorsed by the Metropolitan Transportation Authority.
