#!/usr/bin/env Rscript
# F006 v12: integrated 80 x 86 mm portrait with per-haplotype rings; native vector grid; source calls stay unchanged.
source(Sys.getenv('PAPERPLOT_HELPER'))
library(grid)
root <- '${DATA_DIR2}/projects/1-oil_palm/08-Nature_review/7-Fig1'
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else file.path(root, 'results/10-ancestry/F006_Kmer_Pedigree_Linear_v2')
read_tsv <- function(name) read.delim(file.path(out, 'source-data', name), check.names=FALSE)
runs <- read_tsv('Track_Runs_v2.tsv')
lens <- read_tsv('Chromosome_Lengths_v2.tsv')
comp <- read_tsv('Track_Composition_v2.tsv')
local <- read_tsv('Local_Trend_30Mb_v9.tsv')
area_max <- if(length(args)>=2) as.numeric(args[2]) else 1
stopifnot(area_max %in% c(.5,1))
states <- c('Dura-like','Pisifera-like','Oleifera-like','Mixed','Insufficient')
pal <- setNames(c('#E68A2E','#7955A3','#289D8F','#5C7FA3','#DEE2E7'), states)
rows <- unique(comp$Sample)
group_names <- comp$Group[match(rows, comp$Sample)]
yc <- numeric(14); yc[1] <- 18
for(i in 2:14) yc[i] <- yc[i-1] + 3.4 + ifelse(group_names[i]!=group_names[i-1],.3,0)
stopifnot(nrow(lens)==16,length(rows)==14,area_max==1,
          all(abs(tapply(comp$Fraction,comp$Sample,sum)-1)<1e-12))
width_mm <- 80; height_mm <- 86
spec <- pp_render_spec(n_panels=1,width_mm=width_mm,height_mm=height_mm,
                       panel_tags=FALSE,ocr='off',human_review='pending',
                       text_pt=list(body=6,panel_title=6,species=6,tick=6,
                                    legend=6,caption=6,annotation=6,axis_title=6))
spec$width_exception <- 'User explicitly requested 80 x 86 mm; main track height expanded and local area integrated directly below it.'
spec$layout_authorization <- 'User requested all per-haplotype bars replaced by rings, five bottom material rings removed, local area restored to 0-1, and a cohesive 80 x 86 mm layout. Uniform unbold 6-pt text retained.'
build_panel <- function(inputs, context) {
  spec <- context$render_spec
  grobs <- list()
  add <- function(g) grobs[[length(grobs)+1L]] <<- g
  ux <- function(x) unit(x/width_mm,'npc')
  uy <- function(y) unit(1-y/height_mm,'npc')
  txt <- function(label,x,y,size=spec$text_pt$body,face='plain',just='left',color='#1E2530') {
    add(textGrob(label,x=ux(x),y=uy(y),just=just,
                 gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=color)))
  }
  rect <- function(x,y,w,h,fill,border=NA,lwd=.3) {
    add(rectGrob(x=ux(x),y=uy(y),width=ux(w),height=unit(h/height_mm,'npc'),
                 just=c('left','top'),gp=gpar(fill=fill,col=border,lwd=lwd)))
  }
  line <- function(x0,y0,x1,y1,color='#AAB2BC',lwd=.35,lty=1) {
    add(segmentsGrob(ux(x0),uy(y0),ux(x1),uy(y1),gp=gpar(col=color,lwd=lwd,lty=lty)))
  }
  rounded <- function(x,y,w,h,fill,border=NA,r=.8) {
    add(roundrectGrob(ux(x+w/2),uy(y+h/2),width=ux(w),height=unit(h/height_mm,'npc'),
                     r=unit(r,'mm'),gp=gpar(fill=fill,col=border,lwd=.3)))
  }
  donut <- function(values,x,y,r=3.1) {
    stopifnot(abs(sum(values)-1)<1e-10)
    starts <- c(0,head(cumsum(values),-1))
    for(k in seq_along(states)) if(values[k]>0) {
      theta <- seq(pi/2-2*pi*starts[k],pi/2-2*pi*(starts[k]+values[k]),length.out=max(3,ceiling(values[k]*160)))
      xx <- c(x+r*cos(theta),rev(x+r*.55*cos(theta)))
      yy <- c(y-r*sin(theta),rev(y-r*.55*sin(theta)))
      add(polygonGrob(ux(xx),uy(yy),gp=gpar(fill=pal[k],col='white',lwd=.13)))
    }
    add(circleGrob(ux(x),uy(y),r=unit(r,'mm'),gp=gpar(fill=NA,col='#8A929B',lwd=.25)))
  }
  rect(0,0,width_mm,height_mm,'white')
  # Compact labels preserve the category map; full Latin names are in the caption.
  legx <- c(1.5,26.5,51.5,1.5,26.5)
  legy <- c(2.5,2.5,2.5,6,6)
  labels <- c('Dura-like','Pisifera-like','Oleifera-like','Mixed','Insufficient evidence')
  for(k in 1:5) {
    rect(legx[k],legy[k]-.75,1.8,1.5,pal[k],NA)
    txt(labels[k],legx[k]+2.6,legy[k],size=6)
  }
  txt('Uncalibrated',78,6,size=6,just='right',color='#64707D')
  x0 <- 17; main_w <- 54; cw <- main_w/16
  txt('Chromosome (EG11)',x0,10,size=6,face='plain')
  txt('Share',75.5,13.2,size=6,just='centre')
  txt('Group',1.5,13.2,size=6,face='plain')
  txt('Hap',14,13.2,size=6,face='plain',just='centre')
  for(i in 1:16) {
    x <- x0+(i-.5)*cw
    txt(as.character(i),x,13.2,size=6,just='centre')
    line(x,15.1,x,15.7,'#68717B',.25)
  }
  line(x0,15.1,x0+main_w,15.1,'#68717B',.3)
  for(i in 1:15) line(x0+i*cw,16.5,x0+i*cw,max(yc)+1.35,'#CBD0D6',.16,2)
  fills <- c('Dura'='#FBEDE0','Pisifera'='#EFE8F7','E. oleifera'='#E2F3EF',
             'EO12'='#E8EEF5','EG11'='#EEEAF5','Nigerian'='#FBF2DF','TN'='#E3F0FA','FL'='#E7F4E9')
  for(g in unique(group_names)) {
    ys <- yc[group_names==g]
    rounded(1.3,min(ys)-1.35,10.5,max(ys)-min(ys)+2.7,fills[g],r=.5)
    txt(g,2,mean(ys),size=6,face=if(g %in% c('EO12','EG11')) 'plain' else 'italic')
  }
  ci <- match(runs$Chromosome,lens$Chromosome); ri <- match(runs$Sample,rows)
  x <- x0+(ci-1)*cw+runs$Start0/lens$Length_bp[ci]*cw
  ww <- (runs$End0-runs$Start0)/lens$Length_bp[ci]*cw
  rect(x,yc[ri]-1.2,ww,2.4,pal[runs$Class])
  for(i in seq_along(rows)) {
    yy <- yc[i]
    rounded(x0,yy-1.2,main_w,2.4,NA,'#CDD3DB',r=.35)
    hap <- if(rows[i] %in% c('EO12','EG11')) '\u2013' else if(grepl('h1$|HapA$',rows[i])) 'h1' else 'h2'
    txt(hap,14.4,yy,size=6,just='centre')
    d <- comp[comp$Sample==rows[i],]
    values <- d$Fraction[match(states,d$Class)]
    donut(values,75.5,yy,1.35)
  }
  # The physical slot is rebuilt, rather than raster-resizing the v9 panel.
  txt('Local',1.5,72,size=6)
  txt('support',1.5,74.6,size=6)
  txt('30 Mb',1.5,78,size=6,color='#64707D')
  bx <- x0; bw <- main_w; by <- 69; bh <- 11; bcw <- bw/16
  for(ch in lens$Chromosome) {
    ci <- match(ch,lens$Chromosome)
    dc <- local[local$Chromosome==ch,]
    positions <- unique(dc$Center_bp)
    sums <- numeric(length(positions))
    for(k in c(3,4,2,1,5)) {
      d <- dc[dc$Class==states[k],]
      stopifnot(identical(d$Center_bp,positions))
      xx <- bx+(ci-1)*bcw+positions/lens$Length_bp[ci]*bcw
      upper <- sums+d$Fraction
      add(polygonGrob(ux(c(xx,rev(xx))),uy(c(by+bh-pmin(upper,area_max)/area_max*bh,rev(by+bh-pmin(sums,area_max)/area_max*bh))),
                      gp=gpar(fill=pal[k],col=NA)))
      if(k!=5) add(linesGrob(ux(xx),uy(by+bh-upper/area_max*bh),gp=gpar(col='#FFFFFF70',lwd=.15)))
      sums <- upper
    }
    stopifnot(all(abs(sums-1)<1e-12),all(1-dc$Fraction[dc$Class=='Insufficient']<=area_max))
    txt(as.character(ci),bx+(ci-.5)*bcw,82.5,size=6,just='centre')
    if(ci>1) line(bx+(ci-1)*bcw,max(yc)+1.35,bx+(ci-1)*bcw,80.4,'#B9C1CB',.16,2)
  }
  line(bx,by,bx,by+bh,'#5E6873',.3)
  line(bx,by+bh,bx+bw,by+bh,'#5E6873',.3)
  for(v in c(0,.5,1)) {
    yy <- by+bh-v/area_max*bh
    txt(format(v,trim=TRUE),bx-1.2,yy,size=6,just='right')
    line(bx-.6,yy,bx,yy,'#5E6873',.25)
  }
  gTree(children=do.call(gList,grobs))
}
inputs <- list(runs=runs, chromosome_lengths=lens, track_composition=comp,
               local_composition=local)
builder <- function(s) build_panel(inputs,list(render_spec=s))
plot <- builder(spec)
evidence <- list(mode='production', backend='native_grid',
                 recipe='genome_track_reference + stacked_fraction_composition',
                 template='multi-panel-template.R', data=inputs, colors=pal,
                 order=rows, x_mapping='equal chromosome widths; linear native bp within each column',
                 calibrated=FALSE, smoothing=FALSE,
                 local_display=list(window_bp=3e7,step_bp=2e6,interpolation="linear between exact moving-window summaries within each chromosome",edge_rule="windows clipped to chromosome bounds",main_calls_changed=FALSE,y_limits=c(0,area_max),above_limit_clipped=area_max<1))
attr(plot,'pp_vector_builder') <- builder
attr(plot,'pp_recipe_evidence') <- evidence
attr(plot,'pp_render_spec') <- spec

# Final physical allocations; user explicitly approved the portrait layout.
layout <- data.frame(Panel=c('Two-row legend','Groups and haplotypes','Fourteen tracks',
                             'Haplotype share rings','Local area'),
                     X_mm=c(1.5,1.3,17,73.5,17),Y_mm=c(1.5,12,16.5,12,68.5),
                     Width_mm=c(77,14.5,54,5,54),Height_mm=c(6,54,49.5,54,15))
write.table(layout,file.path(out,'checks/Layout_v12.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
jsonlite::write_json(spec,file.path(out,'checks/Render_Spec_v12.json'),auto_unbox=TRUE,pretty=TRUE)
stem <- file.path(out,'Kmer_Pedigree_Linear_v12')
exports <- pp_save_all_with_qa_loop(plot,stem,formats=c('pdf','svg','png'),
                                   render_spec=spec,max_iterations=0,overwrite=FALSE,
                                   qa_out_dir=file.path(out,'checks/QA_v12'),
                                   qa_context=list(figure_type='genome_track',expected_panels=1))
print(exports)
cat('\nEXPORT_COMPLETE\n')
