#!/usr/bin/env python3
"""Refit original V10 source crops to a larger borderless paired-photo slot."""
from pathlib import Path
from PIL import Image, ImageOps
import json
root=Path(__file__).resolve().parents[1]
pkg=root/'OilPalm_14Haplotype_QC_V4'
src=pkg/'source-data/V10_phenotypes'
dst=pkg/'source-data/V11_phenotypes'
dst.mkdir(exist_ok=False)
for name in ('Oleifera','Seedless','Dura','Pisifera','Nigerian','Tenera'):
    tile=Image.new('RGB',(621,298),'white')
    for view,x,w in [('fruit',0,378),('bunch',392,229)]:
        photo=Image.open(src/f'{name}_{view}.png').convert('RGB')
        fitted=ImageOps.contain(photo,(w,294),Image.Resampling.LANCZOS)
        tile.paste(fitted,(x+(w-fitted.width)//2,(298-fitted.height)//2))
    tile.save(dst/f'{name}.png',dpi=(600,600))
(dst/'Provenance.json').write_text(json.dumps({'original_crops':str(src),'tile_mm':[26.3,12.6],'tile_pixels':[621,298],'changes':'Aspect-preserving refit from original high-resolution crops; no biological image editing; no borders','source_pdf_provenance':str(src/'Provenance.json')},indent=2)+'\n')
print('Six V11 phenotype tiles prepared.')
