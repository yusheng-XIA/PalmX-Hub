#!/usr/bin/env Rscript
# Native 180-mm mixed-encoding layout with one phenotype tile per haplotype pair.
args <- commandArgs(TRUE);stopifnot(length(args)%in%c(1L,2L));replace_output<-length(args)==2L&&identical(args[[2]],'--replace')
project<-normalizePath(args[[1]])
pkg<-file.path(project,'OilPalm_14Haplotype_QC_V4')
stem<-'OilPalm_14Haplotype_QC_V19-color'
source('paperplot_helpers.R')
suppressPackageStartupMessages(library(ape))
datafile<-file.path(pkg,'source-data/OilPalm_14Haplotype_QC_V8_plotting_data.tsv')
d<-read.delim(datafile,check.names=FALSE,na.strings='NA',stringsAsFactors=FALSE)
map<-read.delim(file.path(pkg,'source-data/OilPalm_14Haplotype_QC_V7_tree_mapping.tsv'),check.names=FALSE)
tr<-read.tree(file.path(pkg,'source-data/OilPalm_14Haplotype_QC_V7_rooted.nwk'))
stopifnot(nrow(d)==14,nrow(map)==12,setequal(map$Tree_ID,tr$tip.label))
rows_y<-79.5-(seq_len(nrow(d))-1)*5.2
row_by_tip<-match(map$QC_ID[match(tr$tip.label,map$Tree_ID)],d$Sample_ID)
stopifnot(!anyNA(row_by_tip),!anyDuplicated(row_by_tip))
nt<-length(tr$tip.label);nn<-nt+tr$Nnode
x<-y<-rep(NA_real_,nn);x[seq_len(nt)]<-20.8;y[seq_len(nt)]<-rows_y[row_by_tip]
desc<-function(node) if(node<=nt) node else unlist(lapply(tr$edge[tr$edge[,1]==node,2],desc))
depth<-function(node) if(node<=nt) 0 else 1+max(vapply(tr$edge[tr$edge[,1]==node,2],depth,numeric(1)))
root<-setdiff(tr$edge[,1],tr$edge[,2])[1];maxdepth<-depth(root)
layout_node<-function(node) {
 if(node<=nt)return(invisible(NULL))
 kids<-tr$edge[tr$edge[,1]==node,2];invisible(lapply(kids,layout_node))
 x[node]<<-3+17.8*(1-depth(node)/maxdepth);y[node]<<-mean(range(y[kids]))
}
layout_node(root)
group_cols<-c(Oleifera='#E15759','FL haplotypes'='#4E79A7',Dura='#579D49',Pisifera='#AE789F',Nigerian='#D2A522',TN='#925749')
ink<-'#222222'
metric_cols<-c('#4E79A7','#EF8A24','#579D49','#37A3B1','#B59A20','#A25D9F','#AE8C7B','#2C6EAF')
fields<-c('Assembly_Size_Gb','Scaffold_N50_Mb','Corrected_LAI','Kmer_Completeness_pct','Merqury_QV','BUSCO_C_pct','Gap_Runs_chr16','Chr_Telomere_Ends_n')
headers<-c('Assembly\nsize (Gb)','Scaffold\nN50 (Mb)','LAI','k-mer\ncompleteness (%)','QV','BUSCO (%)','Gaps','Telomere')
caps<-c(2.1,150,25,100,80)
axis_left<-c(69.5,85.5,101.6,117.2,133.8)
axis_width<-c(9.1,10.4,8.4,10.3,8.1)
symbol_x<-c(149.3,160.2,172.5)
header_x<-c(axis_left+axis_width/2,symbol_x)
phenotype_map<-c(Oleifera='Oleifera','FL haplotypes'='Seedless',Dura='Dura',Pisifera='Pisifera',Nigerian='Nigerian',TN='Tenera')
slots<-do.call(rbind,lapply(names(phenotype_map),function(group) {
 idx<-which(d$Figure_Group==group);stopifnot(length(idx)==2,diff(idx)==1)
 do.call(rbind,lapply(c('Fruit','Bunch'),function(view) {
  image<-file.path(pkg,'source-data/V14_phenotypes',paste0(phenotype_map[[group]],'_',tolower(view),'.png'))
  if(group=='Pisifera'&&view=='Bunch')image<-file.path(pkg,'source-data/V17_phenotypes/Pisifera_bunch.png')
  data.frame(Material=group,Source_Label=phenotype_map[[group]],View=view,Rows=paste(idx,collapse=';'),X_mm=if(view=='Fruit')40 else 54.3,Y_mm=mean(rows_y[idx])-4.8,Width_mm=if(view=='Fruit')13 else 10.9,Height_mm=9.6,Image=image)
 }))
}))
rownames(slots)<-NULL;stopifnot(nrow(slots)==12,all(file.exists(slots$Image)))
photos<-lapply(slots$Image,png::readPNG);placements<-slots
for(i in seq_len(nrow(slots))) {
 ar<-dim(photos[[i]])[2]/dim(photos[[i]])[1]
 w<-min(slots$Width_mm[i],slots$Height_mm[i]*ar);h<-w/ar
 placements$Width_mm[i]<-w;placements$Height_mm[i]<-h
 placements$X_mm[i]<-slots$X_mm[i]+(slots$Width_mm[i]-w)/2
 placements$Y_mm[i]<-slots$Y_mm[i]+(slots$Height_mm[i]-h)/2
}
scale_meta<-jsonlite::fromJSON(file.path(pkg,'source-data/V17_phenotypes/Scale_Provenance.json'))
pis<-placements[placements$Material=='Pisifera'&placements$View=='Bunch',]
scale_length_mm<-pis$Width_mm*scale_meta$scale_fraction_of_image_width
scale_x<-65.7 # Move annotation into the clear gutter; calibrated length is unchanged.
scale_y<-pis$Y_mm+(mean(c(scale_meta$bar_start_pt[2],scale_meta$bar_end_pt[2]))-scale_meta$source_bbox_pt[2])/scale_meta$image_height_pt*pis$Height_mm

build<-function(spec) {
 els<-list(grid::rectGrob(gp=grid::gpar(fill='white',col=NA)))
 add<-function(g)els[[length(els)+1L]]<<-g
 mm<-function(v)grid::unit(v,'mm')
 line<-function(x0,y0,x1,y1,col=ink,lwd=.4) add(grid::segmentsGrob(mm(x0),mm(y0),mm(x1),mm(y1),gp=grid::gpar(col=col,lwd=lwd)))
 text<-function(label,x,y,size=6,col=ink,just='left',face='plain') add(grid::textGrob(label,mm(x),mm(y),just=just,gp=grid::gpar(fontfamily='Arial',fontsize=size,fontface=face,col=col,lineheight=.95)))
 circle<-function(x,y,r,fill,col=NA,lwd=.35) add(grid::circleGrob(mm(x),mm(y),mm(r),gp=grid::gpar(fill=fill,col=col,lwd=lwd)))
 rect<-function(x,y,w,h,fill) add(grid::rectGrob(mm(x),mm(y),mm(w),mm(h),just=c('left','bottom'),gp=grid::gpar(fill=fill,col=NA)))
 symbol<-function(x,y,pch,size,col=ink) add(grid::pointsGrob(mm(x),mm(y),pch=pch,size=mm(size),gp=grid::gpar(col=col,fill=col,lwd=.5)))
 roundrect<-function(x,y,w,h,fill,r=1) add(grid::roundrectGrob(mm(x),mm(y),mm(w),mm(h),r=mm(r),just=c('left','bottom'),gp=grid::gpar(fill=fill,col=NA)))
 roundrect(1.5,6.6,37.8,86.8,'#F0F7FC')
 roundrect(39.9,6.6,28.5,86.8,'#EEF8F2')
 roundrect(69,6.6,112.5,86.8,'#F5F9FC')
 roundrect(1.5,89.4,37.8,4,'#DFEDF8');rect(1.5,89.4,37.8,2,'#DFEDF8')
 roundrect(39.9,89.4,28.5,4,'#E0F1E8');rect(39.9,89.4,28.5,2,'#E0F1E8')
 roundrect(69,89.4,112.5,4,'#E2EFF7');rect(69,89.4,112.5,2,'#E2EFF7')
 band_cols<-c(Oleifera='#F7EAF0','FL haplotypes'='#E7F1FC',Dura='#E6F3E8',Pisifera='#F0E9F9',Nigerian='#FFF5DA',TN='#EEEAE7')
 for(group in names(band_cols)) {
  rr<-which(d$Figure_Group==group);yy<-min(rows_y[rr])-2.45;hh<-diff(range(rows_y[rr]))+4.9
  roundrect(2,yy,36.5,hh,band_cols[[group]],r=.8)
  roundrect(69.3,yy,111.9,hh,grDevices::adjustcolor(band_cols[[group]],alpha.f=.85),r=.8)
 }
 text('e',3,95.5,8,face='bold')
 text('Haplotype',21.5,87.4,5.5,just='centre',face='bold')
 symbol(5,84.2,16,1.8);text('hap1',7,84.2,5.5)
 symbol(14.7,84.2,17,2.1);text('hap2',16.7,84.2,5.5)
 symbol(25,84.2,5,1.8,col='#777777');text('Reference',27,84.2,5.5)
 text('Material',84.5,95.5,5.5,face='bold')
 lx<-c(95.5,111,131.5,144.5,158.5,176)
 for(i in seq_along(group_cols)){circle(lx[i],95.5,.7,group_cols[i]);text(names(group_cols)[i],lx[i]+1.5,95.5,5.5)}
 text('Phylogeny',21.5,91.5,7,just='centre',face='bold')
 text('Phenotype',52.25,91.5,7,just='centre',face='bold')
 text('Assembly metrics',125,91.5,7,just='centre',face='bold')
 text('Fruit',46.5,85.5,6,just='centre');text('Bunch',59.75,85.5,6,just='centre')
 for(i in seq_len(nrow(placements)))add(grid::rasterGrob(photos[[i]],mm(placements$X_mm[i]),mm(placements$Y_mm[i]),width=mm(placements$Width_mm[i]),height=mm(placements$Height_mm[i]),just=c('left','bottom'),interpolate=TRUE))
 # Preserve the actual source calibration; annotations are independent vectors.
 line(scale_x,scale_y,scale_x+scale_length_mm,scale_y,lwd=.6)
 line(scale_x,scale_y,scale_x,scale_y+.35,lwd=.6)
 line(scale_x+scale_length_mm,scale_y,scale_x+scale_length_mm,scale_y+.35,lwd=.6)
 text('5 cm',scale_x+scale_length_mm/2,scale_y+1.4,5.5,just='centre')
 for(j in seq_along(headers))text(headers[j],header_x[j],86.7,6,col=metric_cols[j],just='centre',face='bold')
 ticks<-list(c(0,1,2),c(0,75,150),c(0,10,25),c(0,50,100),c(0,40,80))
 for(j in 1:5) {
  left<-axis_left[j];right<-left+axis_width[j]
  line(left,8.8,left,81.4,col='#777777',lwd=.35);line(left,81.4,right,81.4,col='#777777',lwd=.4)
  for(k in seq_along(ticks[[j]])) {
   v<-ticks[[j]][k];xx<-left+axis_width[j]*v/caps[j]
   line(xx,81.4,xx,81.9,col='#777777',lwd=.4)
   text(format(v,trim=TRUE),xx,83.0,5.5,just=if(k==1)'left' else if(k==length(ticks[[j]]))'right' else 'centre')
  }
 }
 for(node in unique(tr$edge[,1])) {
  kids<-tr$edge[tr$edge[,1]==node,2]
  line(x[node],min(y[kids]),x[node],max(y[kids]),lwd=.55)
  for(k in kids)line(x[node],y[k],x[k],y[k],lwd=.55)
 }
 for(i in seq_len(nt)) {
  row<-row_by_tip[i];hap2<-grepl('hap2$',tr$tip.label[i])
  symbol(x[i],y[i],if(hap2)17 else 16,if(hap2)2.1 else 1.8,group_cols[[d$Figure_Group[row]]])
 }
 for(i in which(d$Sample_ID%in%c('EO12','EG11')))symbol(20.8,rows_y[i],5,1.8,col='#777777')
 for(i in seq_len(nrow(d))) {
  yy<-rows_y[i];text(d$Figure_Label[i],23.2,yy,6.5)
  for(j in 1:5) {
   v<-d[[fields[j]]][i];left<-axis_left[j]
   if(is.na(v)){text('NA',left+axis_width[j]/2,yy,5.5,just='centre');next}
   stopifnot(v>=0,v<=caps[j]);end<-left+axis_width[j]*v/caps[j]
   value<-if(j==1)sprintf('%.2f',v) else sprintf('%.1f',v)
   if(j<=3){rect(left,yy-1,end-left,2,metric_cols[j]);text(value,end+.65,yy,5.5)}
   else {line(left,yy,end,yy,col=grDevices::adjustcolor(metric_cols[j],alpha.f=.22),lwd=.55);circle(end,yy,.65,metric_cols[j]);text(value,end+.95,yy,5.5)}
  }
  for(j in 6:8) {
   v<-d[[fields[j]]][i];xx<-symbol_x[j-5];tx<-xx+2.2
   if(is.na(v)){circle(xx,yy,.8,'white','#888888');text('NA',tx,yy,6,col='#777777');next}
   if(j==6) {
    circle(xx,yy,1.3,'#F1E5EF')
    theta<-seq(pi/2,pi/2+2*pi*v/100,length.out=80)
    add(grid::polygonGrob(mm(c(xx,xx+1.3*cos(theta))),mm(c(yy,yy+1.3*sin(theta))),gp=grid::gpar(fill=metric_cols[j],col=NA)))
   } else if(j==7)circle(xx,yy,1.05*sqrt((8+48*log10(v+1)/log10(61374))/56),metric_cols[j])
   else {
    fill<-c('29'='#D4E6E8','30'='#A7CDD1','31'='#6BAAB1','32'='#328C99')[[as.character(v)]]
    add(grid::rectGrob(mm(xx-.8),mm(yy-.8),mm(1.6),mm(1.6),just=c('left','bottom'),gp=grid::gpar(fill=fill,col='#777777',lwd=.25)))
   }
   value<-if(j==6)sprintf('%.1f',v) else if(j==7)format(v,big.mark=',',scientific=FALSE,trim=TRUE) else paste0(as.integer(v),'/32')
   text(value,tx,yy,6,just='left')
  }
 }
 do.call(grid::grobTree,els)
}
spec<-pp_render_spec(width_mm=183,height_mm=99,panel_tags=FALSE,text_pt=list(body=6.5,axis_title=6,tick=5.5,legend=5.5,annotation=5.5,panel_title=7,caption=5.5,panel_tag=8))
plot<-build(spec);attr(plot,'pp_vector_builder')<-build
attr(plot,'pp_recipe_evidence')<-list(mode='production',backend='native-grid-simplified-assembly-overview',data=d,tree=tr,images=placements,scale_bar=scale_meta)
attr(plot,'pp_panel_evidence')<-list(data=d,tip_rows=row_by_tip,row_y_mm=rows_y,phenotype_mapping=phenotype_map,style='User-requested Arial monochrome metric simplification',reference_tree_not_copied=TRUE)
render<-file.path(pkg,'checks/V19_color_render');dir.create(render,showWarnings=FALSE)
files<-pp_save_all_with_qa_loop(plot,file.path(render,stem),render_spec=spec,max_iterations=0)
for(ext in c('pdf','svg','png')) {
 dest<-file.path(pkg,paste0(stem,'.',ext));if(file.exists(dest)&&!replace_output)stop('Refusing existing output: ',dest)
 stopifnot(file.copy(files[[ext]],dest,overwrite=replace_output))
}
if(!file.exists(file.path(pkg,'source-data',paste0(stem,'_plotting_data.tsv'))))stopifnot(file.copy(datafile,file.path(pkg,'source-data',paste0(stem,'_plotting_data.tsv')),overwrite=FALSE))
write.table(slots,file.path(pkg,'source-data',paste0(stem,'_phenotype_slots.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(placements,file.path(pkg,'source-data',paste0(stem,'_image_placements.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(Tree_ID=tr$tip.label,Sample_ID=d$Sample_ID[row_by_tip],Row=row_by_tip,Y_mm=y[seq_len(nt)]),file.path(pkg,'source-data',paste0(stem,'_tree_mapping.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(pp_to_json(list(version='V19-color',haplotype_legend=list(position='below Phylogeny and above EO12',title_y_mm=87.4,key_y_mm=84.2),header_compaction=list(legend_y_mm=95.5,previous_y_mm=99,canvas_height_mm=99),background_variant='colored panels and paired-material bands',status='candidate; user review pending',dimensions_mm=c(183,99),font='Arial',panel_label=list(text='e',point_size=8),body_point_sizes=c(5.5,6,6.5,7),quality_data_unchanged=TRUE,tree_tips=12,reference_marker='open diamond; unconnected',floating_material_labels=FALSE,metric_colors=as.list(stats::setNames(metric_cols,fields)),color_revision='Restored V16 per-metric colors and BUSCO/gap/telomere symbols; only shortened headers retained',symbol_reference='V16',symbol_centers_mm=as.list(symbol_x),numeric_alignment_last_three='left',metric_encodings=c('3 zero-based horizontal bars','2 point tracks','BUSCO pie + value; gap bubble + value; telomere square + count, as V16'),row_pitch_mm=5.2,photographs='Fruit and Bunch as separate native-image views, common pair center, no generated images',scale_bar=list(label='5 cm',sample='Pisifera bunch',length_mm=scale_length_mm,x_mm=scale_x,y_mm=scale_y,source=scale_meta),physical_export=attr(files,'qa_export_audit'))),file.path(pkg,'checks',paste0(stem,'_metadata.json')))
writeLines('**Fig. 1e | Phylogeny, representative phenotypes and assembly metrics of oil-palm haplotypes.** The original midpoint-rooted 12-tip cladogram is retained without branch-length scaling; EO12 and EG11 were not included in that tree and are shown as unconnected open diamonds. Circles and triangles identify hap1 and hap2. Fruit and bunch views are independently fitted within fixed image regions and centered between the two haplotype rows of each material; Seedless represents FL and Tenera represents TN, as confirmed by the user. Images are representative, not a common physical-size comparison. The editable 5-cm bar applies only to the Pisifera bunch and is recalibrated from the original PDF bar and image transform. Assembly size, 16-chromosome Scaffold N50 and LAI are zero-baseline bars; k-mer completeness and Merqury QV are point tracks; complete genome BUSCO uses fraction disks, chromosome N-run counts use the original logarithmic-area circles, and telomere-positive ends/32 use squares with the original fill scale; each retains its numeric label. Headers are BUSCO (%), Gaps and Telomere. All source values are unchanged, with displayed rounding only. FL k-mer completeness comes from separate single-assembly rows against the same all-read k-mer denominator, not a combined assembly result; read sets, k and historical evaluation versions are not fully harmonized across materials. BUSCO databases also differ, so these values do not establish a uniform-method ranking.',file.path(pkg,'checks',paste0(stem,'_legend.md')))
cat('V19 relocated haplotype legend, Arial 183 x 99 mm, three bars + two point tracks + V16 terminal symbols restored.\n')
