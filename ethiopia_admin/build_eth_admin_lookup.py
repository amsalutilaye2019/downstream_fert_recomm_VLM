"""Build the Ethiopia admin lookup and boundary files used by eda_gold_vs_extracted.ipynb.

Source: OCHA / CSA + BoFED "Ethiopia - Subnational Administrative Boundaries" (COD-AB, 2021 version,
eth_admbnda_adm3_csa_bofedb_2021.*): 1,082 woredas (ADM3) with their zone (ADM2) and region (ADM1).
That version already has the Sidama region (2020) and the South West Ethiopia region (Nov 2021).
SNNP was dissolved in August 2023 into Central Ethiopia and South Ethiopia, so the remaining SNNP zones are
assigned to those two regions from the zone lists of the two new regions:
  Central Ethiopia: Gurage, Hadiya, Kembata Tembaro, Halaba, Silte zones and Yem (special woreda -> zone)
  South Ethiopia:   Wolayita, Gamo, Gofa, Gedeo, South Omo, Konso zones and Ale, Amaro, Basketo, Burji, Derashe
Outputs (same folder): eth_admin_lookup.csv, eth_adm3_regions.geojson (woredas, lon/lat, simplified to ~200 m),
                       eth_regions.geojson (regions = admin 1, dissolved).
Run:  python build_eth_admin_lookup.py        (needs pyshp, pyproj, shapely, pandas)
"""
import json
from pathlib import Path
import pandas as pd
import shapefile
import shapely
from pyproj import CRS, Transformer
from shapely.geometry import MultiPolygon, Polygon, shape, mapping
from shapely.ops import transform, unary_union

folder = Path(__file__).parent
reader = shapefile.Reader(str(folder / 'eth_admbnda_adm3_csa_bofedb_2021'))

# The shapefile is in UTM zone 37N (metres): convert to longitude / latitude (WGS 84)
crs = CRS.from_wkt((folder / 'eth_admbnda_adm3_csa_bofedb_2021.prj').read_text())
to_lonlat = Transformer.from_crs(crs, 'EPSG:4326', always_xy=True).transform

adm1_names = {'Benishangul Gumz': 'Benishangul-Gumuz', 'Gambela': 'Gambella',
              'South West Ethiopia': 'Southwest Ethiopia'}
central_zones = {'Guraghe', 'Hadiya', 'Kembata Tembaro', 'Halaba', 'Siltie', 'Yem Special'}
south_zones = {'Alle', 'Amaro', 'Basketo', 'Burji', 'Derashe', 'Gamo', 'Gedeo', 'Gofa', 'Konso',
               'South Omo', 'Wolayita'}

rows, feats = [], []
for sr in reader.iterShapeRecords():
    rec = sr.record.as_dict()
    woreda, zone, adm1 = rec['ADM3_EN'].strip(), rec['ADM2_EN'].strip(), rec['ADM1_EN'].strip()
    if adm1 == 'SNNP':
        if zone in central_zones:
            region = 'Central Ethiopia'
        elif zone in south_zones:
            region = 'South Ethiopia'
        else:
            raise ValueError(f'SNNP zone not assigned: {zone}')
    else:
        region = adm1_names.get(adm1, adm1)
    geom = transform(to_lonlat, shape(sr.shape.__geo_interface__))
    geom = shapely.set_precision(geom.simplify(0.002, preserve_topology=True), 1e-4)   # ~200 m, keeps the files small
    rows.append({'adm3_woreda': woreda, 'adm2_zone': zone, 'adm1_region_2021': adm1, 'region': region})
    feats.append((region, woreda, zone, adm1, geom))

pd.DataFrame(rows).to_csv(folder / 'eth_admin_lookup.csv', index=False)

with open(folder / 'eth_adm3_regions.geojson', 'w', encoding='utf-8') as f:
    json.dump({'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {'region': r, 'woreda': w, 'zone': z, 'adm1_2021': a},
         'geometry': mapping(g)} for r, w, z, a, g in feats]}, f, separators=(',', ':'))

def region_outline(gs):
    """Dissolve the woredas of a region; close the thin gaps left by the simplification and drop inner holes."""
    u = unary_union([g.buffer(0.003) for g in gs]).buffer(-0.003)
    polys = list(u.geoms) if u.geom_type == 'MultiPolygon' else [u]
    return MultiPolygon([Polygon(p.exterior) for p in polys])

by_region = {}
for r, w, z, a, g in feats:
    by_region.setdefault(r, []).append(g)
with open(folder / 'eth_regions.geojson', 'w', encoding='utf-8') as f:
    json.dump({'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {'region': r}, 'geometry': mapping(region_outline(gs))}
        for r, gs in by_region.items()]}, f, separators=(',', ':'))

print(len(rows), 'woredas;', len(by_region), 'regions:', sorted(by_region))
