"""Build a woreda (ADM3) -> zone (ADM2) -> region (ADM1) lookup from the geoBoundaries Ethiopia files.

geoBoundaries ADM3 polygons carry only their own name, so the parent zone and region are found by
locating a point inside each woreda polygon within the ADM2 and ADM1 polygons.
Source: geoBoundaries gbOpen, ETH ADM1 / ADM2 / ADM3 (https://www.geoboundaries.org), CC BY 4.0.
Run:  python build_eth_admin_lookup.py        (needs shapely and pandas)
"""
import json
from pathlib import Path
import pandas as pd
from shapely.geometry import shape

folder = Path(__file__).parent

def load(level):
    with open(folder / f'geoBoundaries-ETH-{level}.geojson', encoding='utf-8') as f:
        feats = json.load(f)['features']
    return [(ft['properties']['shapeName'], shape(ft['geometry'])) for ft in feats]

adm1, adm2, adm3 = load('ADM1'), load('ADM2'), load('ADM3')

def parent_of(point, parents):
    for name, geom in parents:
        if geom.contains(point):
            return name
    # point on a border / outside: take the nearest polygon
    return min(parents, key=lambda p: p[1].distance(point))[0]

rows = []
for name, geom in adm3:
    pt = geom.representative_point()
    rows.append({'adm3_woreda': name, 'adm2_zone': parent_of(pt, adm2), 'adm1_region': parent_of(pt, adm1)})

lookup = pd.DataFrame(rows)
lookup.to_csv(folder / 'eth_admin_lookup.csv', index=False)
print(len(lookup), 'woredas ->', lookup['adm2_zone'].nunique(), 'zones,', lookup['adm1_region'].nunique(), 'regions')
