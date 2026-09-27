#!/usr/bin/env python3
"""Extract unchanged photographs from the supplied labelled PDF render."""
from pathlib import Path
import csv,json,hashlib
from PIL import Image,ImageOps,ImageDraw

root=Path(__file__).resolve().parents[1]
pkg=root/'OilPalm_14Haplotype_QC_V4'
source=pkg/'表型_白底_名称标注.pdf'
render=root/'tmp/fig1e_v10/phenotype_source_600dpi.png'
dest=pkg/'source-data/V10_phenotypes'
dest.mkdir(exist_ok=False)
# Coordinates are PDF points measured from the top-left of the source page.
boxes={
 'Dura':((24,23,109,94),(276,23,346,112)),
 'Tenera':((110,23,194,94),(351,28,420,112)),
 'Seedless':((195,34,258,94),(427,28,486,112)),
 'Pisifera':((28,140,104,207),(276,132,346,219)),
 'Nigerian':((108,137,194,207),(350,133,421,215)),
 'Oleifera':((193,140,262,207),(428,136,486,214)),
}
im=Image.open(render).convert('RGB')
factor=600/72
records=[];tiles=[]
for group,(fruit,bunch) in boxes.items():
 parts=[]
 for kind,box in [('fruit',fruit),('bunch',bunch)]:
  pixels=tuple(round(v*factor) for v in box)
  crop=im.crop(pixels)
  crop.save(dest/f'{group}_{kind}.png',dpi=(600,600))
  parts.append(crop)
  records.append({'Material':group,'View':kind,'Source_PDF':str(source),'Crop_Points_Top_Left':','.join(map(str,box)),'Pixel_Width':crop.width,'Pixel_Height':crop.height})
 # One 20.5 x 12.1 mm image per material: two source views, no repetition by haplotype.
 tile=Image.new('RGB',(484,286),'white')
 for part,(x,w) in zip(parts,[(0,279),(287,197)]):
  fitted=ImageOps.contain(part,(w,270),Image.Resampling.LANCZOS)
  tile.paste(fitted,(x+(w-fitted.width)//2,(286-fitted.height)//2))
 tile.save(dest/f'{group}.png',dpi=(600,600))
 tiles.append((group,tile))
with (dest/'Crop_Manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
(dest/'Provenance.json').write_text(json.dumps({'source_pdf':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'render_dpi':600,'operation':'Crop labels/white margins and fit photos without changing aspect ratio or biological content; no AI image generation','display_scale':'Views fitted independently; do not infer comparative fruit size from displayed dimensions'},indent=2)+'\n')
sheet=Image.new('RGB',(968,948),'white');draw=ImageDraw.Draw(sheet)
for i,(group,tile) in enumerate(tiles):
 x=(i%2)*484;y=(i//2)*316
 draw.text((x+8,y+5),group,fill='black');sheet.paste(tile,(x,y+25))
sheet.save(root/'tmp/fig1e_v10/phenotype_crop_preview.png')
print('Extracted six labelled material tiles; Seedless mapping remains explicit, not inferred.')
