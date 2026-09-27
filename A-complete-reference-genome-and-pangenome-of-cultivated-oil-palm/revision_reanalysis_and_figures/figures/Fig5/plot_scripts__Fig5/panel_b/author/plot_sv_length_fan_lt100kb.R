#!/usr/bin/env Rscript
# User-requested 50 bp <= SV length <100 kb subset. One locus per observation; angles encode within-type length-bin proportions.
suppressPackageStartupMessages(library(ggplot2))
source('paperplot_helpers.R')
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[1]); out <- normalizePath(args[2])
stem <- file.path(out,'checks','F001_SV_Length_fan_lt100kb')
stopifnot(dir.exists(file.path(out,'checks')),dir.exists(file.path(out,'source-data')))
visible <- file.path(out,paste0('F001_SV_Length_fan_lt100kb.',c('pdf','svg','png')))
if(any(file.exists(visible))) stop('Refusing to replace existing figures')
if(length(Sys.glob(paste0(stem,'*')))) stop('Refusing to overwrite an existing candidate')
d <- read.delim(input,check.names=FALSE)
types <- c('DEL','INS','DUP','INV','TRA')
expected <- c(73536L,86268L,11234L,608L,6668L)
stopifnot(nrow(d)==178314L,!anyDuplicated(d$SV_ID),all(is.finite(d$Length_bp)),all(d$Length_bp>=50),identical(as.integer(table(factor(d$SV_Type,levels=types))),expected))
original_counts <- expected
d <- d[d$Length_bp < 100000, , drop=FALSE]
expected <- c(73536L,86268L,11142L,396L,6519L)
stopifnot(nrow(d)==177861L,identical(as.integer(table(factor(d$SV_Type,levels=types))),expected))
filter_summary <- data.frame(SV_Type=types,Original_N=original_counts,Retained_N=expected,Excluded_GE100kb=original_counts-expected)
write.table(filter_summary,file.path(out,'source-data','F001_SV_Length_fan_lt100kb_Filter.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
breaks <- c(50,100,500,1000,5000,10000,50000,100000)
bins <- c('50–100 bp','100–500 bp','0.5–1 kb','1–5 kb','5–10 kb','10–50 kb','50–100 kb')
cols <- c('#F8EDDB','#EED2AA','#E4B47C','#D58A49','#C16532','#A84424','#833022')
d$Bin <- cut(d$Length_bp,breaks=breaks,right=FALSE,labels=bins)
stopifnot(!anyNA(d$Bin))
tab <- as.data.frame(table(SV_Type=factor(d$SV_Type,levels=types),Length_Bin=d$Bin),stringsAsFactors=FALSE)
names(tab)[3] <- 'Count'
tab$Total <- expected[match(tab$SV_Type,types)]
tab$Fraction <- tab$Count/tab$Total
tab$Percent <- 100*tab$Fraction
tab$Angle_Deg <- 72*tab$Fraction
stopifnot(sum(tab$Count)==177861L,all(abs(tapply(tab$Fraction,tab$SV_Type,sum)-1)<1e-12))
centres <- 90-72*(0:4); polygons <- list(); k <- 0
for(i in seq_along(types)) {
 s <- tab[tab$SV_Type==types[i],]; s <- s[match(bins,s$Length_Bin),]
 mid <- centres[i]*pi/180
 cx <- 60+2*cos(mid); cy <- 51+2*sin(mid)
 start <- centres[i]+36
 for(j in seq_along(bins)) {
  sweep <- s$Angle_Deg[j]
  if(sweep>0) {
   angles <- seq(start,start-sweep,length.out=max(3,ceiling(sweep/.2)+1))*pi/180
   k <- k+1
   polygons[[k]] <- data.frame(x=c(cx,cx+31*cos(angles),cx),y=c(cy,cy+31*sin(angles),cy),id=k,Bin=bins[j])
  }
  start <- start-sweep
 }
 stopifnot(abs(start-(centres[i]-36))<1e-10)
}
poly <- do.call(rbind,polygons)
labels <- data.frame(x=60+44*cos(centres*pi/180),y=51+44*sin(centres*pi/180),Name=c('Deletions','Insertions','Duplications','Inversions','Translocations'),N=paste0('n = ',format(expected,big.mark=',',trim=TRUE)))
legend_edges <- seq(20,100,length.out=8)
legend <- data.frame(xmin=head(legend_edges,-1),xmax=tail(legend_edges,-1),Bin=bins)
boundary <- data.frame(x=legend_edges,Label=c('0.05','0.1','0.5','1','5','10','50','100'))
p <- ggplot()+
 geom_polygon(data=poly,aes(x,y,group=id,fill=Bin),colour='#69533E',linewidth=.14)+
 geom_text(data=labels,aes(x,y+1.8,label=Name),family='Arial',size=8/ggplot2::.pt)+
 geom_text(data=labels,aes(x,y-2.2,label=N),family='Arial',size=7/ggplot2::.pt)+
 geom_rect(data=legend,aes(xmin=xmin,xmax=xmax,ymin=109,ymax=111,fill=Bin),colour=NA)+
 geom_text(data=boundary,aes(x,y=105.8,label=Label),family='Arial',size=7/ggplot2::.pt)+
 annotate('text',x=60,y=116,label='SV length (kb)',family='Arial',size=8/ggplot2::.pt)+
 annotate('text',x=60,y=4,label='SVs <100 kb; proportions within each type',family='Arial',size=7/ggplot2::.pt)+
 scale_fill_manual(values=setNames(cols,bins),limits=bins,drop=FALSE,guide='none')+
 coord_fixed(xlim=c(0,120),ylim=c(0,120),expand=FALSE)+labs(x=NULL,y=NULL)+
 theme_void(base_family='Arial',base_size=8)+theme(plot.margin=margin(0,0,0,0),plot.background=element_rect(fill='white',colour=NA))
roles <- c('panel_title','axis_title','species','tick','legend','caption','annotation','body')
spec <- pp_render_spec(n_panels=1L,width_mm=120,height_mm=120,mode='preview',panel_tags=FALSE,text_pt=as.list(setNames(rep(8,length(roles)),roles)),human_review='pending',expected_labels=labels$Name)
sha <- strsplit(system2('sha256sum',shQuote(input),stdout=TRUE),' +')[[1]][1]
stopifnot(sha=='39b3a1f8e26025737c777ba1e7f021fc6c901609605655373c3544224c66f2c4')
attr(p,'pp_panel_evidence') <- list(source=input,sha256=sha,summary=tab,unit='SV locus',denominator='within SV type after restricting to 50 <= Length_bp < 100000',filter_summary=filter_summary,bin_rule='[lower, upper)',palette=setNames(cols,bins),sector_deg=72,radial_translation_mm=2,radius_mm=31,type_order_clockwise=types,limits=c('TRA length is translocated reference segment span','Different evidence layers have different detection limits'))
write.table(tab,file.path(out,'source-data','F001_SV_Length_fan_lt100kb_Bins.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
base_theme <- pp_production_theme
pp_production_theme <- function(spec) base_theme(spec)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank(),plot.margin=margin(0,0,0,0))
stopifnot(inherits(pp_normalize_production(p,spec)$theme$axis.line,'element_blank'))
files <- pp_save_all_with_qa_loop(p,stem,preset='nature',formats=c('pdf','svg','png'),max_iterations=0L,render_spec=spec,qa_out_dir=file.path(out,'checks','fan_lt100kb_visual'),qa_context=list(family='composition',allow_grid='off'))
writeLines(c('# SV length fan preview','', 'Five equal 72-degree sectors; fixed radius; each sector translated 2 mm radially to create white gaps. Clockwise from top: DEL, INS, DUP, INV, TRA. Within each type, clockwise slices encode ascending length bins. Angles and areas within a sector are proportional to locus counts; each type is independently normalized to 100%. Type totals are annotated.', '', 'Seven shared bins: [50,100), [100,500), [500,1000), [1000,5000), [5000,10000), [10000,50000), [50000,100000) bp. User-requested restriction to Length_bp <100000 removes 453 loci (DUP 92, INV 212, TRA 149) and retains 177861; all denominators and annotated n refer to retained loci; zero-count bins retained in Source Data.', '', paste('Source:',input),paste('SHA256:',sha),'', 'Length is the existing SVLEN_Median_bp per catalog locus. TRA uses translocated reference-segment span, not breakpoint distance. DEL/INS use read-and-assembly evidence; DUP/INV/TRA use SyRI rearrangement evidence. Length distributions also reflect detection differences. Sample-carrier frequency is not used.', '', '120 x 120 mm. Candidate preview; final human acceptance pending. Custom fan geometry implements the composition-count contract; no exact fan recipe/template exists. No significance test; no biological reanalysis.'),paste0(stem,'_README.md'))
stopifnot(all(file.copy(unname(files),file.path(out,basename(unname(files))),overwrite=FALSE)))
cat('Export status:',attr(files,'qa_contract')$status,'\n')
