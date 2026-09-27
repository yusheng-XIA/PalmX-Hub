#!/usr/bin/env python3
"""Fig. 1a contribution bars, recomputed with a symmetric crop definition.
Source: FAOSTAT QCL bulk file Production_Crops_Livestock_E_All_Data_(Normalized).zip
(file dated 2025-12-23), downloaded 2026-09-24; World (area code 5000), 2023.
Area = element 5312 'Area harvested' of the primary crop; oil = element 5510 'Production'
of the derived crude oil item(s) of the same crop."""
import pandas as pd
d = pd.read_csv('qcl_world_2023_from_bulk.csv')
d = d[(d['Area Code'] == 5000) & (d.Year == 2023)]
def get(item, el):
    r = d[(d.Item == item) & (d['Element Code'] == el)]
    assert len(r) == 1, (item, el); r = r.iloc[0]
    return float(r.Value), int(r['Item Code']), str(r['Item Code (CPC)']).strip("'"), r.Flag
CROPS = {  # crop: (area item, [oil items])
 'Oil palm': ('Oil palm fruit', ['Palm oil', 'Oil of palm kernel']),
 'Soybean': ('Soya beans', ['Soya bean oil']),
 'Rapeseed': ('Rape or colza seed', ['Rapeseed or canola oil, crude']),
 'Sunflower': ('Sunflower seed', ['Sunflower-seed oil, crude']),
 'Groundnut': ('Groundnuts, excluding shelled', ['Groundnut oil']),
 'Cotton (cottonseed)': ('Seed cotton, unginned', ['Cottonseed oil']),
 'Coconut (copra)': ('Coconuts, in shell', ['Coconut oil']),
 'Olive': ('Olives', ['Olive oil']),
 'Sesame': ('Sesame seed', ['Oil of sesame seed']),
 'Linseed': ('Linseed', ['Oil of linseed']),
 'Safflower': ('Safflower seed', ['Safflower-seed oil, crude']),
 'Maize': ('Maize (corn)', ['Oil of maize']),
}
NINE = ['Oil palm','Soybean','Rapeseed','Sunflower','Groundnut','Cotton (cottonseed)','Coconut (copra)','Olive','Sesame']
SCHEMES = {
 'A_9crops_recommended': NINE,
 'B_USDA8_no_sesame': [c for c in NINE if c != 'Sesame'],
 'C_11crops_all_QCL_oilcrop_oils': NINE + ['Linseed','Safflower'],
 'S1_9crops_no_cotton(sensitivity)': [c for c in NINE if c != 'Cotton (cottonseed)'],
 'S2_9crops_plus_maize_area_and_oil(sensitivity)': NINE + ['Maize'],
}
rows = []
for crop,(ai,oils) in CROPS.items():
    a, ac, acpc, af = get(ai, 5312)
    ov = [get(o, 5510) for o in oils]
    rows.append(dict(crop=crop, area_item=ai, area_item_code=ac, area_item_cpc=acpc, area_element='5312 Area harvested',
        area_flag=af, area_harvested_ha_2023=a, oil_items=' + '.join(oils),
        oil_item_codes=' + '.join(str(x[1]) for x in ov), oil_item_cpc=' + '.join(x[2] for x in ov),
        oil_element='5510 Production', oil_flags=' + '.join(x[3] for x in ov),
        oil_production_t_2023=sum(x[0] for x in ov)))
C = pd.DataFrame(rows).set_index('crop')
out, summ = [], []
for s, crops in SCHEMES.items():
    sub = C.loc[crops].copy(); A, O = sub.area_harvested_ha_2023.sum(), sub.oil_production_t_2023.sum()
    sub['share_of_area_pct'] = 100*sub.area_harvested_ha_2023/A
    sub['share_of_oil_pct'] = 100*sub.oil_production_t_2023/O
    sub.insert(0, 'scheme', s); sub = sub.reset_index()
    tot = dict(scheme=s, crop='TOTAL', area_harvested_ha_2023=A, oil_production_t_2023=O, share_of_area_pct=100.0, share_of_oil_pct=100.0)
    out += [sub, pd.DataFrame([tot])]
    op = sub[sub.crop=='Oil palm'].iloc[0]
    summ.append((s, len(crops), A/1e6, O/1e6, op.share_of_area_pct, op.share_of_oil_pct))
f = pd.concat(out, ignore_index=True)
f['area_item_code'] = f['area_item_code'].astype('Int64')
f = f[['scheme'] + [c for c in f.columns if c != 'scheme']]
f['source'] = 'FAOSTAT QCL (Crops and livestock products), World, 2023; bulk file dated 2025-12-23'
f['downloaded'] = '2026-09-24'
f.to_csv('Fig1a_shares_final.tsv', sep='\t', index=False, float_format='%.4f')
print(C[['area_item_code','area_harvested_ha_2023','oil_item_codes','oil_production_t_2023','area_flag','oil_flags']].to_string())
print('\nscheme  n  area_Mha  oil_Mt  OP_area%  OP_oil%')
for r in summ: print('%-48s %2d %8.2f %7.2f %6.2f %6.2f' % r)
po = C.loc['Oil palm']; print('\nold asymmetric (9 area / 9 oils + maize oil): %.2f / %.2f' % (
 100*po.area_harvested_ha_2023/C.loc[NINE].area_harvested_ha_2023.sum(),
 100*po.oil_production_t_2023/(C.loc[NINE].oil_production_t_2023.sum()+C.loc['Maize'].oil_production_t_2023)))
