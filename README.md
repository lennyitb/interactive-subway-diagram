# Interactive subway diagram data

NYC subway network data fitted to a schematic (non-geographic) diagram, for building a pan/zoom diagram explorer where lines and markers stay a constant on-screen size and labels are placed at runtime with collision avoidance.

## Contents

| Path | What it is |
| --- | --- |
| `data/nyc-subway-diagram.json` | Stations, routes, drawn line bands, transfers, connectors and ferry, in the diagram's coordinate space (viewBox `0 0 2500 2700`, y down) |
| `data/basemap.json` | Water, land, park and borough-boundary shapes as SVG path data in the same space |
| `preview.html` | Canvas viewer used to check the data (pan, zoom, route highlight, greedy label placement) |
| `tools/fit/` | Scripts that extracted and fitted the data from the source SVG |
| `NYC_subway-4D.svg` | Source diagram (see attribution below) |

### Data model

- **stations**: `id`, `name`, `system` (`subway`, `sir`, `airtrain`, `path`), anchor `x`/`y`, the `routes` that stop there, and `markers`, which are the individual dots drawn on each service line. Each marker has a `style` taken from the diagram legend (`always`, `interchange`, `weekdays-only`, `rush-peak-only`, `except-rush-peak`, `nights-weekends-only`, ...). Stations linked by solid connectors are merged into one complex.
- **routes**: 29 service patterns (1-7, A-Z, the <6>/<7>/<F> express variants, three shuttles and SIR). Each pattern has ordered `stops`, a `path` polyline traced along the drawn lines, and `stopPathIndex`, the path vertex for each stop. Stop sequences reflect regular weekday daytime service.
- **bands**: the source strokes as drawn (one per service group), with colour, width, dash flag and the routes using them. Draw these to reproduce the original look, or draw `routes[].patterns[].path` per service.
- **connectors / transfers**: in-station, out-of-station and part-time transfer links.

## Running the preview

```sh
python3 -m http.server 8765
# open http://127.0.0.1:8765/preview.html
```

## Re-fitting

Station names are in `tools/fit/names.txt` (cluster id to name) and stop sequences in `tools/fit/routes_def.py`. After editing either, run `tools/fit/run.sh`. It needs `numpy`, `lxml` and `svgelements`.

## Attribution and license

The source diagram `NYC_subway-4D.svg` is by [CountZ (Jake Berman)](https://en.wikipedia.org/wiki/User:CountZ) at English Wikipedia ([53studio.com](https://53studio.com)), via [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:NYC_subway-4D.svg), licensed [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).

The data in `data/` is derived from it (positions, geometry and shapes extracted and fitted; names and stop sequences added) and is shared under the same **CC BY-SA 3.0** license.

The code (`preview.html`, `tools/`) is released under the MIT License.

This project is not affiliated with or endorsed by the Metropolitan Transportation Authority.
