#!/usr/bin/env python3
"""Recompute the Fig. 1a contribution bars from FAOSTAT QCL (World, 2023).
Source: FAOSTAT Crops and livestock products, bulk file
Production_Crops_Livestock_E_All_Data_(Normalized).zip (release file dated
2025-12-23), downloaded 2026-09-24. The authors' original calculation table was
not found; the crop list below is the combination that reproduces both 8.6%
and 39.5% with this release and must be confirmed by the authors."""
import pandas as pd
d = pd.read_csv('src/faostat_QCL_world_2022_2024.csv'); d = d[d.Year == 2023]
ah = d[d['Element Code'] == 5312].set_index('Item')['Value']
pr = d[d['Element Code'] == 5510].set_index('Item')['Value']
rows = [  # crop (area item), oil item(s)
 ('Oil palm', 'Oil palm fruit', ['Palm oil', 'Oil of palm kernel']),
 ('Soybean', 'Soya beans', ['Soya bean oil']),
 ('Rapeseed', 'Rape or colza seed', ['Rapeseed or canola oil, crude']),
 ('Sunflower', 'Sunflower seed', ['Sunflower-seed oil, crude']),
 ('Groundnut', 'Groundnuts, excluding shelled', ['Groundnut oil']),
 ('Cotton', 'Seed cotton, unginned', ['Cottonseed oil']),
 ('Coconut', 'Coconuts, in shell', ['Coconut oil']),
 ('Olive', 'Olives', ['Olive oil']),
 ('Sesame', 'Sesame seed', ['Oil of sesame seed']),
 ('Maize (oil only)', None, ['Oil of maize']),
]
out = []
for crop, a, oils in rows:
    out.append(dict(crop=crop, faostat_area_item=a or '(area not included)',
                    area_harvested_ha_2023=ah[a] if a else 0.0,
                    faostat_oil_items=' + '.join(oils),
                    oil_production_t_2023=sum(pr[o] for o in oils)))
df = pd.DataFrame(out)
A, O = df.area_harvested_ha_2023.sum(), df.oil_production_t_2023.sum()
df['share_of_area_pct'] = 100 * df.area_harvested_ha_2023 / A
df['share_of_oil_pct'] = 100 * df.oil_production_t_2023 / O
tot = dict(crop='TOTAL', faostat_area_item='', area_harvested_ha_2023=A, faostat_oil_items='',
           oil_production_t_2023=O, share_of_area_pct=100.0, share_of_oil_pct=100.0)
df = pd.concat([df, pd.DataFrame([tot])])
df.to_csv('Fig1a_shares.tsv', sep='\t', index=False, float_format='%.4f')
print(df.to_string())
# sensitivity
po, pk = pr['Palm oil'], pr['Oil of palm kernel']
print('\nSensitivity (FAOSTAT 2023, current release):')
print(' oil share without maize oil (8 oils + PKO): %.2f%%' % (100*(po+pk)/(O-pr['Oil of maize'])))
print(' palm oil only / same denominator: %.2f%%' % (100*po/O))
print(' palm fruit area / "Oilcrops, Oil Equivalent" area: %.2f%%' % (100*ah['Oil palm fruit']/ah['Oilcrops, Oil Equivalent']))
print(' palm oil / all 13 vegetable oils in QCL: %.2f%%' % (100*po/sum(pr[i] for i in ['Palm oil','Oil of palm kernel','Soya bean oil','Rapeseed or canola oil, crude','Sunflower-seed oil, crude','Groundnut oil','Cottonseed oil','Coconut oil','Oil of maize','Olive oil','Oil of sesame seed','Oil of linseed','Safflower-seed oil, crude'])))
