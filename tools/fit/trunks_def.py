# Hand rules for collapsing the 4D diagram's parallel service lanes into one centreline per trunk.
#
# The source draws every service on its own lane, 10 units apart. collapse.py treats lanes that run
# parallel at that spacing as one trunk. Where two separate lines (different tunnels or levels) happen
# to be drawn at lane spacing, list them here so they stay separate.
#
# Each rule is (routes on one side, routes on the other side, scope). Scope is None for everywhere, or
# (station display name, radius in diagram units) to apply only near that station.

NEVER_MERGE = [
    ({'2', '3'}, {'N', 'Q', 'R', 'W'}, None),            # 7 Av vs Broadway at 49 St / Times Sq
    ({'2', '3', '4', '5'}, {'B', 'Q'}, None),             # Eastern Pkwy vs Brighton, Atlantic Av to Prospect Park
    ({'7', '7X'}, {'GS'}, None),                          # Flushing vs 42 St Shuttle
    ({'7', '7X'}, {'N', 'W'}, None),                      # Flushing vs Astoria at Queensboro Plaza
    ({'E'}, {'B', 'D'}, None),                            # 53 St vs 6 Av at 7 Av, 8 Av vs 6 Av at W 4 St
    ({'E'}, {'J', 'Z'}, None),                            # Archer Av upper vs lower level
    ({'A', 'C'}, {'B', 'D'}, ('W 4 St-Washington Sq', 40)),  # 8 Av vs 6 Av levels (they do share Central Park West)
    ({'G'}, {'A', 'C'}, None),                            # Crosstown vs Fulton St
    ({'A', 'C'}, {'J', 'Z'}, None),                       # Fulton St vs Jamaica at Broadway Jct
    ({'F', 'FX'}, {'Q', 'N'}, None),                      # 63 St / Culver vs Brighton / Sea Beach (Lex/63, W 8 St, Coney Island)
]

# Services that are drawn as their own lane in the data but share a lane with their parent when rendered.
LANE_OF = {'6X': '6', '7X': '7', 'FX': 'F'}

# Lines drawn further apart than the lane pitch that are nonetheless one trunk. Each rule is
# (routes on one side, routes on the other side, widest gap between their lanes, scope); scope is a
# station (name, radius) or a box (x0, y0, x1, y1) in diagram units.
ALWAYS_MERGE = [
    ({'2', '3'}, {'4', '5'}, 34.0, (810, 1860, 1300, 2080)),   # Eastern Pkwy, Hoyt St to Franklin Av
]
