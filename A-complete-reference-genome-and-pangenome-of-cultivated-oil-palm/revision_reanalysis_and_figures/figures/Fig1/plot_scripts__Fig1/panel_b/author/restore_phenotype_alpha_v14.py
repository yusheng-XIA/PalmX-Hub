#!/usr/bin/env python3
"""Recover original PDF transparency in the existing scientific figure assets."""
from pathlib import Path
import csv,json,hashlib
from PIL import Image,ImageOps
root=Path(__file__).resolve().parents[1];pkg=root/'OilPalm_14Haplotype_QC_V4'
source=pkg/'表型_白底_名称标注.pdf'
page=Image.open(root/'tmp/fig1e_v14/phenotype_transparent.png')
assert page.mode=='RGBA' and page.getchannel('A').getextrema()==(0,255)
dest=pkg/'source-data/V14_phenotypes';dest.mkdir(exist_ok=False)
with (pkg/'source-data/V10_phenotypes/Crop_Manifest.tsv').open() as f: records=list(csv.DictReader(f,delimiter='\t'))
boxes={(r['Material'],r['View']):tuple(float(v) for v in r['Crop_Points_Top_Left'].split(',')) for r in records}
checks=[]
for name in ('Oleifera','Seedless','Dura','Pisifera','Nigerian','Tenera'):
    tile=Image.new('RGBA',(621,298),(0,0,0,0))
    for view,x,w in [('fruit',0,378),('bunch',392,229)]:
        crop=page.crop(tuple(round(v*600/72) for v in boxes[name,view]))
        crop.save(dest/f'{name}_{view}.png',dpi=(600,600))
        fitted=ImageOps.contain(crop,(w,294),Image.Resampling.LANCZOS)
        tile.alpha_composite(fitted,(x+(w-fitted.width)//2,(298-fitted.height)//2))
    alpha=tile.getchannel('A'); hist=alpha.histogram()
    assert hist[0]>0 and hist[255]>0
    assert all(tile.getpixel(pt)[3]==0 for pt in ((0,0),(620,0),(0,297),(620,297)))
    tile.save(dest/f'{name}.png',dpi=(600,600))
    checks.append({'Material':name,'Mode':tile.mode,'Transparent_fraction':hist[0]/(621*298),'Opaque_fraction':hist[255]/(621*298),'Corner_alpha':0})
(dest/'Provenance.json').write_text(json.dumps({'source_pdf':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'render_command':'pdftocairo -png -transp -r 600 -singlefile','operation':'Recover source PDF native transparency; reuse V10 crop bounds and V11 aspect-preserving layout; no color-threshold background removal or generative modification','assets':checks},indent=2)+'\n')
print(json.dumps(checks,indent=2))
