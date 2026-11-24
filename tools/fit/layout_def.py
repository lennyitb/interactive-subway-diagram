# Hand placement, applied after collapse.py has built the line graph. Use it for the handful of places the
# automatic fit gets awkward. Stations are named as in the data (display name, or station id when a name is
# shared, such as the several Canal Sts or 7 Avs). A rule for a route also applies to its variants (<F>, <6>).
import math

def arc(cx, cy, r, a0, a1, n=6):
    """Points along a circular arc (degrees, y down: 0 = east, 90 = south, 180 = west, 270 = north)."""
    return [(cx + r * math.cos(math.radians(a0 + (a1 - a0) * k / n)),
             cy + r * math.sin(math.radians(a0 + (a1 - a0) * k / n))) for k in range(1, n + 1)]

# POINTS: named positions used by WALK. Routes that list the same point share it: it becomes a junction
# exactly there.
POINTS = {
    # Queens Blvd: locals leave at qJ1 and go straight up, express E is joined by the F at qJ2, all meet at qJ3
    'qJ1': (895, 1055), 'qJ2': (929.2, 1017.5), 'qC': (895, 897), 'qJ3': (1039, 897),
    # top of 6 Av: B/D join the F at sBD, the M leaves the E at sE and joins at sM
    'sBD': (462, 1128.5), 'sE': (474, 1137), 'sMc': (462, 1137), 'sM': (462, 1149),
    # bottom of Central Park West: B/D join the 8 Av trunk at one corner
    'cpw': (258, 1128.5),
    # Broadway at 57 St: Q carries on north, N/R/W turn east along 60 St
    'b57': (366, 1058), 'b63': (366, 1017.5),
    # 63 St: F and Q share the tunnel from where the F turns in (j63) to where the Q turns north (jq)
    'j63': (462, 1017.5), 'jq': (678, 1017.5),
    # Downtown Brooklyn: A/C join the F at jA above Jay St and leave it at j1 with one 90-degree turn; the G
    # meets the F at jFG; A/C and G share Hoyt-Schermerhorn and run together east to jG, where the G turns
    # NE; A/C turn NE once at k1 onto the straight line to Broadway Jct
    'jA': (845, 1813), 'j1': (845, 1954), 'jFG': (845, 1981.5), 'jG': (1033, 1954),
    # A/C part company (A express, C local) at dAC, just past Hoyt-Schermerhorn
    'dAC': (887, 1954),
    # 8 Av: the locals' pairing changes from C/E to C/A at cAC, halfway between Canal St and Chambers St
    'cAC': (443.25, 1651.0),
    # R: middle of the 45-degree climb from Court St to Jay St, where it moves across to its 4 Av lane
    'rCJ': (827.5, 1858.5),
    # DeKalb Av: N/D leave the bridge at dJ1 and go round a jug handle (east to dC, then south down 4 Av) that
    # crosses the Q/B diagonal at dX; the R joins the diagonal at dJR, stops at DeKalb with Q/B, and turns
    # south onto 4 Av at dX
    # (all on one straight line from the bridge at (740.4, 1628) to Atlantic Av-Barclays Ctr)
    'dJ1': (936.2, 1825.3), 'dC': (984.3, 1825.3), 'dX': (984.3, 1873.6), 'dJR': (950.9, 1840),
    # Queens Plaza: the R leaves N/W at rJ, drops at 45 degrees and joins E/M at qJR just before the platform;
    # E/M turn once at ecq from 53 St onto the Queens Blvd diagonal
    'rJ': (771, 1058), 'qJR': (834.6, 1121.6), 'ecq': (820.6, 1137),
}

# MOVE: (route, station, x, y) puts the platform where the route stops at that station at (x, y); None keeps
# that coordinate.
MOVE = [
    # Queens Blvd locals: 36 St and Steinway on one vertical, 46 St and Northern Blvd on one horizontal
    ('M', '36 St', 895, 990), ('M', 'Steinway St', 895, 935), ('M', '46 St', 970, 897), ('M', 'Northern Blvd', None, 897),
    # 6 Av: the F comes straight down from 57 St
    ('A', 'Hoyt-Schermerhorn Sts', 872.5, 1954), ('C', 'Lafayette Av', 1119, 1954),
    ('F', 'York St', 845, None), ('F', 'Jay St-MetroTech', 845, None), ('F', 'Bergen St', 845, None), ('F', 'Carroll St', 845, None), ('F', 'Smith-9 Sts', 880, 2086.9),
    ('Q', 'DeKalb Av', 967.6, 1856.8),
    ('N', 'Atlantic Av-Barclays Ctr', 984.3, None), ('R', '59 St', 984.3, None),   # 4 Av: one vertical from the jug handle   # on the Q/B diagonal, inside the N/D jug handle
    ('F', '57 St', 462, None), ('F', '47-50 Sts-Rockefeller Center', 462, None),
    ('F', '42 St-Bryant Pk/5 Av', 462, None), ('F', '34 St-Herald Sq', 462, None), ('F', '23 St', 462, None),
    ('F', '14 St/6 Av', 462, None), ('F', 'W 4 St-Washington Sq', 462, None),
    # Broadway: 57 St and 49 St in one vertical line; Times Sq down on the diagonal, below the 7
    ('R', '57 St-7 Av', 366, 1090), ('R', '49 St', 366, None),
    ('R', 'Times Sq-42 St/Port Authority Bus Terminal', 410, 1235),
    # streets along the bottom of Central Park: level lines
    ('N', '5 Av/59 St', None, 1058), ('N', 'Lexington Av/59 St', None, 1058), ('N', 'Queensboro Plaza', None, 1058),
    ('F', 'Lexington Av/63 St', None, 1017.5), ('F', 'Roosevelt Island', None, 1017.5), ('F', '21 St-Queensbridge', None, 1017.5),
    # Central Park West and 8 Av: one exact vertical
    ('C', '116-st-bc', 257.5, None), ('C', '14 St/8 Av', 257.5, None),
]

# WALK: (route, from station, to station, [points and stations in between])
# The route's track between the two stations (where it stops or passes) becomes straight links through the
# listed points (named in POINTS, written '@name') and stations. Use it where the topology needs changing:
# a different junction, a branch leaving somewhere else.
WALK = [
    # Queens Blvd express bypass: M/R turn up before 36 St, run straight past Steinway, turn right before
    # 46 St and rejoin; the E carries on along the diagonal and the F joins it there
    ('M', 'Queens Plaza', 'Northern Blvd', ['@qJ1', '36 St', 'Steinway St', '@qC', '46 St', '@qJ3']),
    ('R', 'Queens Plaza', 'Northern Blvd', ['@qJ1', '36 St', 'Steinway St', '@qC', '46 St', '@qJ3']),
    ('E', 'Queens Plaza', 'Northern Blvd', ['@qJ1', '@qJ2', '@qJ3']),
    ('F', '21 St-Queensbridge', 'Northern Blvd', ['@qJ2', '@qJ3']),
    # top of 6 Av: each branch makes one 90-degree turn into the trunk
    ('F', '57 St', '47-50 Sts-Rockefeller Center', ['@sBD', '@sM']),
    ('B', '7-av-bde', '47-50 Sts-Rockefeller Center', ['@sBD', '@sM']),
    ('D', '7-av-bde', '47-50 Sts-Rockefeller Center', ['@sBD', '@sM']),
    ('M', '5 Av/53 St', '47-50 Sts-Rockefeller Center', ['@sE', '@sMc', '@sM']),
    ('E', '5 Av/53 St', '7-av-bde', ['@sE']),
    # bottom of CPW: B/D come west from 7 Av and turn up into the trunk once
    ('A', '59 St-Columbus Circle', '50-st-ce', ['@cpw']),
    ('C', '59 St-Columbus Circle', '50-st-ce', ['@cpw']),
    ('B', '59 St-Columbus Circle', '7-av-bde', ['@cpw']),
    ('D', '59 St-Columbus Circle', '7-av-bde', ['@cpw']),
    # 57 St-7 Av: N/R/W turn east at 60 St, the Q goes on up to 63 St
    ('N', '57 St-7 Av', '5 Av/59 St', ['@b57']),
    ('R', '57 St-7 Av', '5 Av/59 St', ['@b57']),
    ('W', '57 St-7 Av', '5 Av/59 St', ['@b57']),
    ('Q', '57 St-7 Av', 'Lexington Av/63 St', ['@b57', '@b63', '@j63']),
    # 63 St tunnel: F and Q on separate tracks side by side through one Lex/63 St platform
    ('F', '57 St', 'Lexington Av/63 St', ['@j63']),
    ('F', 'Lexington Av/63 St', 'Roosevelt Island', ['@jq']),
    ('Q', 'Lexington Av/63 St', '72 St', ['@jq']),
    # Downtown Brooklyn (see POINTS)
    ('A', 'High St', 'Hoyt-Schermerhorn Sts', ['@jA', 'Jay St-MetroTech', '@j1']),
    ('C', 'High St', 'Hoyt-Schermerhorn Sts', ['@jA', 'Jay St-MetroTech', '@j1']),
    ('F', 'York St', 'Bergen St', ['@jA', 'Jay St-MetroTech', '@j1', '@jFG']),
    ('G', 'Hoyt-Schermerhorn Sts', 'Bergen St', ['@jFG']),
    ('G', 'Hoyt-Schermerhorn Sts', 'Fulton St', ['@dAC', '@jG']),
    # 8 Av: all three through cAC between Canal St and Chambers St
    ('A', 'Canal St', 'Chambers St/WTC/Park Pl/Cortlandt St', ['@cAC']),
    ('C', 'Canal St', 'Chambers St/WTC/Park Pl/Cortlandt St', ['@cAC']),
    ('E', 'Canal St', '@444.5,1682.3', ['@cAC']),
    ('R', 'Court St/Borough Hall', 'Jay St-MetroTech', ['@rCJ']),
    # DeKalb Av (from a point on the Manhattan Bridge approach)
    ('N', '@923.2,1804.8', 'Atlantic Av-Barclays Ctr', ['@dJ1', '@dC', '@dX']),
    ('D', '@923.2,1804.8', 'Atlantic Av-Barclays Ctr', ['@dJ1', '@dC', '@dX']),
    ('Q', '@923.2,1804.8', 'Atlantic Av-Barclays Ctr', ['@dJ1', '@dJR', 'DeKalb Av', '@dX']),
    ('B', '@923.2,1804.8', 'Atlantic Av-Barclays Ctr', ['@dJ1', '@dJR', 'DeKalb Av', '@dX']),
    ('R', 'Jay St-MetroTech', 'Atlantic Av-Barclays Ctr', ['@dJR', 'DeKalb Av', '@dX']),
    ('A', 'Hoyt-Schermerhorn Sts', 'Lafayette Av', ['@dAC', '@jG']),
    ('C', 'Hoyt-Schermerhorn Sts', 'Lafayette Av', ['@dAC', '@jG']),
    # Queens Plaza: R joins E/M on one radial curve just before the station
    ('E', 'Court Sq', 'Queens Plaza', ['@ecq', '@qJR']),
    ('M', 'Court Sq', 'Queens Plaza', ['@ecq', '@qJR']),
    ('R', 'Lexington Av/59 St', 'Queens Plaza', ['@rJ', '@qJR']),
    ('N', 'Lexington Av/59 St', 'Queensboro Plaza', ['@rJ']),
    ('W', 'Lexington Av/59 St', 'Queensboro Plaza', ['@rJ']),
]

# Stations in the rules below can also be '@name' (a POINTS entry) or '@x,y' (the route's node nearest there).

# PLATFORM: (station, route, other route): the other route stops at the same platform as the route, so a
# station with services on separate tracks side by side shows as one wide pill.
PLATFORM = [
    ('Lexington Av/63 St', 'F', 'Q'),
    ('Hoyt-Schermerhorn Sts', 'A', 'G'),
]

# SHARE: (route, follows, from station, to station)
# The route uses the other route's track between the two stations, for services that share a tunnel the
# source happens to draw apart.
SHARE = [
    ('M', 'E', 'Court Sq', 'Queens Plaza'),          # 53 St tunnel: E and M are one pair of tracks
]

# PATHS: (route, from station, to station, [waypoints])
# The route's track between the two stations becomes straight lines through the waypoints: (x, y) points,
# or ('J', x, y) to also place the junction nearest that point exactly there. Every node on the way is moved
# onto the new line (nearest point, in order), and the edges along it follow the line.
PATHS = [
    # Eastern Pkwy: the 4/5 run straight on from Bowling Green; the 2/3 come down from Clark St on one
    # diagonal through Borough Hall and turn once into the trunk, which runs due east past Hoyt and Nevins and
    # turns once onto the diagonal through Atlantic Av
    ('2', 'Wall St', 'Bergen St', [(659.5, 1845), (742, 1845), ('J', 814, 1917), (1066, 1917)]),
    ('4', 'Court St/Borough Hall', 'Nevins St', [('J', 814, 1917)]),
    # 6 Av: one exact vertical from 47-50 Sts to the turn at W 4 St
    ('F', '47-50 Sts-Rockefeller Center', 'Broadway-Lafayette St/Bleecker St', [(462, 1541)]),
    # Fulton St line: level past Lafayette Av, one turn, then one straight line to Broadway Jct
    ('C', 'Lafayette Av', 'Broadway Junction', [(1158, 1954)]),
    # G: one quarter-circle from Hoyt-Schermerhorn down onto the F line
    ('G', 'Hoyt-Schermerhorn Sts', 'Bergen St', arc(872.5, 1981.5, 27.5, 270, 180, 24) + [('J', 845, 1981.5)]),
    # F: one straight line from York St to Carroll St
    ('F', 'York St', 'Carroll St', []),
    # F/G: straight down from Carroll St, one 90-degree turn east into Smith-9 Sts
    ('F', 'Carroll St', 'Smith-9 Sts', [(845, 2086.9)]),
    # top of Lex: each branch turns once between the trunk and its next stop
    ('4', '161 St-Yankee Stadium', '138 St-Grand Concourse', [(518, 486), ('J', 611, 579)]),
    ('5', '3 Av-149 St', '138 St-Grand Concourse', [(611, 561), ('J', 611, 579)]),
    ('6', '3 Av-138 St', '125 St', [('J', 611, 652)]),
    # Eastern Pkwy straight from the Grand Army Plaza corner to Utica Av; the 2/5 leave it at Rogers Junction
    # on one 45-degree line, with President St on that line rather than on the bend
    ('3', 'Grand Army Plaza', 'Crown Heights-Utica Av', [(1217, 2066), ('J', 1404, 1917)]),
    ('2', 'Franklin Av-Medgar Evers College', 'Sterling St', [('J', 1404, 1917)]),
    # 2/3 split: the 2 leaves the Lenox line on one sweeping curve onto 149 St
    ('2', '135 St', '149 St-Grand Concourse', [('J', 454, 621)] + arc(514, 621, 60, 180, 270, 18)),
    # level lines along the bottom of Central Park
    ('N', '5 Av/59 St', 'Queensboro Plaza', []),
    ('F', 'Roosevelt Island', '21 St-Queensbridge', []),
    ('C', '116-st-bc', '14 St/8 Av', []),
    # R/W: one smooth curve from Rector St into South Ferry/Whitehall St; the R carries straight on east and
    # jogs down at 45 degrees halfway to Court St
    ('R', 'Rector St', 'South Ferry/Whitehall St', [(509.5, 1831.4)] + arc(545.5, 1831.4, 36, 180, 90, 18)),
    ('R', 'South Ferry/Whitehall St', 'Court St/Borough Hall', [(683.8, 1867.4), (693.4, 1877.0)]),
    ('R', 'Court St/Borough Hall', 'Jay St-MetroTech', [(809, 1877), (846, 1840)]),
    # Manhattan Bridge to Atlantic Av: one straight diagonal through DeKalb Av
    ('Q', 'Canal St', 'Atlantic Av-Barclays Ctr', [(740.4, 1628)]),
    ('R', 'Atlantic Av-Barclays Ctr', '59 St', []),
    # Broadway: straight down from 57 St past 49 St, one bend onto the diagonal above Times Sq
    ('R', '57 St-7 Av', '34 St-Herald Sq', [(366, 1191)]),
]

# ORDER: (route, from station, to station, lanes left to right when travelling from -> to)
# Fixes the lane order on the route's track between the two stations; branches and neighbouring trunks then
# take their order from it.
ORDER = [
    # Queens Blvd: M and R branch off north, so they ride on the north-west side from Queens Plaza
    # 53 St: the M rides south of the E so it peels off south into 6 Av cleanly; it crosses to the north side
    # where the R joins just before Queens Plaza
    ('M', '5 Av/53 St', 'Court Sq', ['E', 'M']),
    ('R', 'Lexington Av/59 St', 'Queens Plaza', ['R', 'M', 'E']),
    ('M', 'Queens Plaza', '36 St', ['R', 'M', 'F', 'E']),
    # Central Park West, west to east: express pair A D, local pair C B, colours alternating
    ('A', '59 St-Columbus Circle', '125-st-acbd', ['A', 'D', 'C', 'B']),
    # Broadway heading north: the Q on the west side so it carries on north at 57 St while N/R/W turn east
    # without crossing it; 60 St then runs N W R from north to south and the R leaves south at Queens Plaza
    ('R', '34 St-Herald Sq', '57 St-7 Av', ['Q', 'N', 'W', 'R']),
    ('N', '@b57', 'Queensboro Plaza', ['N', 'W', 'R']),
    # Manhattan Bridge heading into Brooklyn: D/N on the north-east side, so they leave together for the jug
    # handle at DeKalb Av, and B/Q on the south-west side (the B crosses D/N where it joins the bridge, the
    # fewest crossings); the R joins B/Q on the south-west side
    ('N', 'Canal St', '@dJ1', ['D', 'N', 'B', 'Q']),
    ('D', 'Broadway-Lafayette St/Bleecker St', 'Grand St', ['M', 'F', 'B', 'D']),   # B/D keep their 6 Av order to the bridge
    ('R', '@dJR', '@dX', ['B', 'Q', 'R']),
    # Fulton St line: A and C keep the order they have on 8 Av, so they don't swap at WTC
    # 8 Av heading south: the E joins from the east at 50 St, locals C/E side by side, A on the west; at
    # Chambers the E ends and the C pairs with the A
    ('E', '50 St', 'Chambers St/WTC/Park Pl/Cortlandt St', ['E', 'C', 'A']),
    ('A', 'Fulton St', 'Euclid Av', ['C', 'A']),
    # Jay St heading south: A/C on the east side, so they turn east to Hoyt without crossing the F; east of
    # Hoyt the G rides on the south side until it turns NE
    ('A', 'Jay St-MetroTech', 'Hoyt-Schermerhorn Sts', ['C', 'A', 'F']),
    ('A', 'Hoyt-Schermerhorn Sts', 'Lafayette Av', ['C', 'A', 'G']),
    # F/G heading south: the G on the east side, where it joins, and on the inside of the turn at Smith-9 Sts
    ('F', 'Bergen St', 'Church Av', ['G', 'F']),
    # Lexington into Brooklyn: 4 and 5 keep their Lex order past Brooklyn Bridge
    ('4', 'Brooklyn Bridge-City Hall/Chambers St', 'Nevins St', ['5', '4']),
    # 2 and 3 keep their 7 Av order all the way to Rogers Junction, where they part
    ('2', 'Chambers St', 'Franklin Av-Medgar Evers College', ['2', '3']),
    # 63 St heading east: the Q on the north side, as it turns north at jq
    ('F', '57 St', 'Roosevelt Island', ['Q', 'F']),
    # 6 Av heading south: M joins from the east, B/D from the west, keeping B east of D as on 53 St
    ('B', '47-50 Sts-Rockefeller Center', 'W 4 St-Washington Sq', ['M', 'F', 'B', 'D']),
]

# GROUP: (route, from station, to station, lanes): on the route's track between the two stations these lanes,
# where they run side by side, form one track group (no gap between them).
GROUP = [
    ('E', '@cAC', '@446.6,1733.2', ['C', 'A']),   # halfway from Canal St to Chambers St, on to WTC: the C pairs with the A
    ('E', 'Queens Plaza', 'Jackson Hts-Roosevelt Av/74 St-Broadway', ['F', 'E']),   # Queens Blvd express pair
    ('N', 'Canal St', '@dJ1', ['D', 'N']),   # Manhattan Bridge: D/N pair (bypassing DeKalb Av) ...
    ('Q', 'Canal St', '@dJ1', ['B', 'Q']),   # ... and B/Q pair
    ('N', '@dJ1', '@dX', ['D', 'N']),   # D/N stay a pair round the DeKalb Av jug handle
    ('Q', '@dJ1', 'Atlantic Av-Barclays Ctr', ['B', 'Q']),   # B/Q stay a pair through DeKalb Av
]

# APART: (route, from station, to station, lanes): on the route's track between the two stations these lanes
# are separate tracks (drawn with a gap), even where they make the same stops. A list inside the list is a set
# of lanes that may still share a track with each other.
APART = [
    ('F', '57 St', 'Roosevelt Island', ['Q', 'F']),   # 63 St tunnel: F and Q on their own tracks
    ('A', 'W 4 St-Washington Sq', '@cAC', [['C', 'E'], 'A']),   # 8 Av: locals | express, to halfway from Canal St to Chambers St
    ('F', 'York St', 'Bergen St', ['F', ['A', 'C']]),   # Jay St: F beside the A/C pair
    ('A', 'Hoyt-Schermerhorn Sts', '@dAC', [['A', 'C'], 'G']),   # out of Hoyt: A/C pair, G beside them
    ('A', '@dAC', 'Euclid Av', ['A', 'C']),   # from just past Hoyt: A express, C local
    ('A', '@dAC', '@jG', ['C', 'A', 'G']),
    ('R', '@dJR', '@dX', ['R', ['Q', 'B']]),   # DeKalb Av: the R on its own track beside Q/B   # ... with the G on its own track beside them until it turns away
    ('E', '@cAC', '@446.6,1733.2', ['E', ['C', 'A']]),   # from halfway to Chambers St on to WTC: E apart, C pairs with A
]

# ANCHOR: (route, from station, to station, lane): on the route's track between the two stations the lane keeps
# the position it has where it arrives at the first station, and the other lanes are placed around it, so it is
# drawn perfectly straight through junctions where lanes join and leave.
ANCHOR = [
    ('F', 'East Broadway', 'Church Av', 'F'),   # F straight from York St past Jay St and Bergen St to Carroll St
    ('A', 'Hoyt-Schermerhorn Sts', 'Euclid Av', 'A'),   # A straight on from Hoyt: only the C moves out where they part, the G stays put
    ('Q', '@j63', '@jq', 'Q'),   # Q straight along 63 St: the F joins beside it, the Q has its one bend at b63
    ('R', 'South Ferry/Whitehall St', '@rCJ', 'R'),
    # DeKalb Av: D/N keep their bridge position round the jug handle and down 4 Av, B/Q theirs down Brighton, so
    # neither pair moves where the other leaves or the R joins
    ('N', '@dJ1', '59 St', 'N'),
    ('Q', '@dJ1', 'Brighton Beach', 'Q'),   # R straight on out of Whitehall St after the W ends, to the Jay St climb
]
