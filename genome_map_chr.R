library(circlize)
library(stringr)
library(patchwork)
library(gcookbook)
library(grid)
library(ggplot2)
library(data.table)
#####chr
#cnvr_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\cnvr.csv"
#cen_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\cen_position.txt"
#recurrent_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\nstd\\ClinGen recurrent CNV .bed file V1.1-hg19.bed"
#refseq_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\ncbiRefSeq.txt\\ncbiRefSeq.txt"
#omim_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\OMIM Genes.aed"
#hi_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\ClinGen_haploinsufficiency_gene_GRCh37.bed"
#ti_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\data\\ClinGen_triplosensitivity_gene_GRCh37.bed"
#output_dir <- "C:\\Users\\huliang\\Desktop\\data and script\\plot_gviz\\chrs\\"

cnvr_dir <- "/root/project/nipt2022/cnv_data/cnv_data/cnvr.csv"
cen_dir <- "/root/project/nipt2022/cnv_data/cnv_data/cen_position.txt"
recurrent_dir <- "/root/project/nipt2022/cnv_data/cnv_data/ClinGen recurrent CNV .bed file V1.1-hg19.bed"
refseq_dir <- "/root/project/nipt2022/cnv_data/cnv_data/ncbiRefSeq.txt"
omim_dir <- "/root/project/nipt2022/cnv_data/cnv_data/OMIM Genes.aed"
hi_dir <- "/root/project/nipt2022/cnv_data/cnv_data/ClinGen_haploinsufficiency_gene_GRCh37.bed"
ti_dir <- "/root/project/nipt2022/cnv_data/cnv_data/ClinGen_triplosensitivity_gene_GRCh37.bed"
output_dir <- "/root/project/nipt2022/cnv_data/cnv_data/gviz/"

draw_one_chr <- function(chr){
  print(chr)
cnvr <- read.csv(cnvr_dir)
cnvr$chr <- gsub(23,'X',cnvr$chr)
cnvr <- cnvr[which(cnvr$chr == chr ), ]
chrband <- read.delim(cen_dir,stringsAsFactors = F)
chrband <- chrband[which(chrband$X.chrom == paste0('chr',chr)), ]
chrband$band <- chrband$gieStain
chrband$gieStain <- gsub('gneg',5,chrband$gieStain)
chrband$gieStain <- gsub('gpos','',chrband$gieStain)
chrband$gieStain <- gsub('acen',50,chrband$gieStain)
chrband$gieStain <- gsub('gvar',50,chrband$gieStain)
chrband$gieStain <- gsub('stalk',50,chrband$gieStain)
chrband$gieStain <- as.integer(chrband$gieStain)
recurrent <- read.table(recurrent_dir,stringsAsFactors = F)
recurrent <- recurrent[which(recurrent$V1 == paste0('chr',chr) ), ]
recurrent$V2 <- as.integer(recurrent$V2)
recurrent$V3 <- as.integer(recurrent$V3)
recurrent <- recurrent[nrow(recurrent):1, ]
acen_start <- chrband[which(chrband$band == 'acen'), ][1,2]
acen_end <- chrband[which(chrband$band == 'acen'), ][2,3]
acen <- chrband[which(chrband$band == 'acen'), ][1,3]
last_tick <- floor(chrband[nrow(chrband),3]/1e7)

#known_gene <- fread("C:\\Users\\huliang\\Desktop\\data and script\\data\\knownGene.txt\\knownGene.txt",data.table=F)
#known_gene <- known_gene[known_gene$V2 == paste0('chr',chr),4:5]
#known_gene <- known_gene[!duplicated(known_gene[ ,1]), ]
refseq_gene <- fread(refseq_dir,data.table=F)
refseq_gene <- refseq_gene[refseq_gene$V3 == paste0('chr',chr),]
refseq_gene$size <- refseq_gene$V6 - refseq_gene$V5
out <- data.frame()
for(i in unique(refseq_gene$V13)){
  new <- refseq_gene[which(refseq_gene$V13 == i), ]
  new <- new[which(new$size == max(new$size)), ][1, ]
  out <- rbind(out,new)
}
refseq_gene <- out[,c(5:6,2) ]
#refseq_gene <- refseq_gene[which(substr(refseq_gene$V2,1,2) == 'NM'), ]
#refseq_gene <- refseq_gene[which(substr(refseq_gene$V2,1,2) == 'NM'),]
#refseq_gene <- refseq_gene[,c(5:6,2) ]
omim_gene <- fread(omim_dir,data.table=F)
omim_gene <- omim_gene[-c(1:9),c(1:3,9,7)]
omim_gene <- omim_gene[which(omim_gene$`omim:phenoMapKey(aed:Integer)` != 0), ]
omim_gene <- omim_gene[which(omim_gene$`bio:sequence(aed:String)` == paste0('chr',chr)), ]
omim_gene$dominant <- ifelse(grepl('dominant',omim_gene[ ,5])==T,2,1)
haplo <- read.delim(hi_dir,header = F)
haplo <- haplo[-1, ]
haplo <- haplo[which(haplo$V1 == paste0('chr',chr)), ]
haplo$V5 <- ifelse(haplo$V5 >= 30,0,haplo$V5)
haplo <- haplo[which(haplo$V5 %in% 1:3), ]
triplo <- read.delim(ti_dir,header = F)
triplo <- triplo[-1, ]
triplo <- triplo[which(triplo$V5 %in% 1:3), ]
triplo <- triplo[which(triplo$V1 == paste0('chr',chr)), ]
haplo$type <- 1
if(nrow(triplo) > 0){
  triplo$type <- 2
  haplo <- rbind(haplo,triplo)
}
#if(chr %in% c(1:5)){
  refseq <- refseq_gene[1, ]
  for(i in 2:nrow(refseq_gene)){
    if(refseq_gene[i,1] - refseq[nrow(refseq),2] >= 1000){
      refseq <- rbind(refseq, refseq_gene[i, ])
    }
    else{
      if(refseq_gene[i,2] - refseq[nrow(refseq),2] >= 0){
        refseq[nrow(refseq),2] <- refseq_gene[i,2]
        if(substr(refseq_gene[i,3],1,2)=='NM' ){
          refseq[nrow(refseq),3] <- 'NM'
        }
      }
    }
  }
#}
  omim <- omim_gene[1, ]
  for(i in 2:nrow(omim_gene)){
    if(omim_gene[i,2] - omim[nrow(omim),3] >= 1000){
      omim <- rbind(omim, omim_gene[i, ])
    }
    else{
      #if(omim_gene[i,2] - omim[nrow(omim),2] >= 0){
      omim[nrow(omim),3] <- omim_gene[i,3]
      #}
    }
  }
#}
#else{
#  refseq <- refseq_gene
#}
remove <- c()
for(i in 1: nrow(cnvr)){
  if((cnvr[i,2] >= acen_start & cnvr[i,2] <= acen_end) |
     (cnvr[i,3] >= acen_start & cnvr[i,2] <= acen_end)){
    remove <- c(remove,i)
  }
  if(cnvr[i,7] == '6' ){
    if((cnvr[i,2] >= 29794756 & cnvr[i,2] <= 33086926) |
       (cnvr[i,3] >= 29794756 & cnvr[i,2] <= 33086926)){
      remove <- c(remove,i)
    }
  }
} 
if(length(remove) > 0 ){
  cnvr <- cnvr[-remove, ]
}
cnvr_df <- data.frame()
for(k in c('dup','del')){
  a <- cnvr[which(cnvr$Freq >= 2 & cnvr$indel == k ), ]
  a$id <- 0
  a <- a[order(a$start,-a$size), ]
  ids  <- 1:nrow(a)
  track_id <- 1
  while(length(ids) >0 ){
    track1 <- min(ids)
    for(i in 2:nrow(a)){
      if(a[i,2] > a[track1[length(track1)],3] & i %in% ids){
        track1 <- c(track1,i)
        next
      }
    }
    a[track1,6] <- track_id
    ids <- ids[ids %in% track1 == F]
    track_id <- track_id + 1
  }
  cnvr_df <- rbind(cnvr_df,a)
}
cnvr_df[ ,6] <- 0.4 * cnvr_df[ ,6] 
cnvr_df[which(cnvr_df$indel == 'del'),6] <- -1 *  cnvr_df[which(cnvr_df$indel == 'del'),6] - 2.1
#f_dup = colorRamp2(breaks = c(12,113,1130), colors = c("#87CEFF",  "#1C86EE",'#0000CD'))
#f_del = colorRamp2(breaks = c(12,113,1130), colors = c("#FFAEB9",  "#FF3E96",'#8B2252'))
f_dup = colorRamp2(breaks = c(2,12,113,1130), colors = c("#87CEFF","#1E90FF",  "#1874CD",'#104E8B'))
f_del = colorRamp2(breaks = c(2,12,113,1130), colors = c("#FFAEB9",  "#FF3E96",'#8B2252',"#7D26CD"))
f_grey = colorRamp2(breaks = c(10,25,50,75,100), colors = c("grey90","grey75",'grey50','grey75','grey90'))
f_omim = colorRamp2(breaks = c(1,2), colors = c("#C1FFC1","#2E8B57"))
f_haplo = colorRamp2(breaks = c(0,1,2,3,30,40), colors = c("grey80","#008B8B",'#8B008B','#8B008B','#0000FF','#CD00CD'))#DarkCyan

p <- ggplot() + theme_classic() 
for(i in 1:nrow(cnvr_df)){
  p <- p + annotate('segment',x=cnvr_df[i,2],y=cnvr_df[i,6],
                    xend = cnvr_df[i,3],yend = cnvr_df[i,6], 
                    linewidth = 2,colour=ifelse(cnvr_df[i,8] == 'dup',f_dup(cnvr_df[i,4]),f_del(cnvr_df[i,4]))
                    )
}
for(i in 1:nrow(chrband)){
  p <- p + annotate('segment',x=chrband[i,2],y=0,xend = chrband[i,3],yend = 0, 
                    linewidth = ifelse(chrband[i,6] == 'acen',1,ifelse(chrband[i,6] == 'stalk',2,3)),
                    colour=f_grey(chrband[i,5]))
}
p <- p + annotate('segment',x= acen - 200000,y=0,xend = acen + 200000,yend = 0, linewidth = 1,colour='black')
for(i in 1:nrow(recurrent)){
  p <- p + annotate('segment',x=recurrent[i,2],y=-0.5,xend = recurrent[i,3],yend = -0.5, linewidth = 3,
                    colour=ifelse(recurrent[i,9]=='0,0,0','grey50','#FF6F00'))
}
for(i in 1:nrow(refseq)){
  p <- p + annotate('segment',x=refseq[i,1],y=-1,
                    xend = refseq[i,2],yend = -1, linewidth = 3,
                    colour=ifelse(substr(refseq[i,3],1,2) == 'NM','#4169E1','grey90'))
}
for(i in 1:nrow(omim)){
  p <- p + annotate('segment',x=omim[i,2],y=-1.5,
                    xend = omim[i,3],yend = -1.5, linewidth = 3,
                    #colour=ifelse(omim[i,3] == 0,'grey90','#006400'))
                    colour=f_omim(omim[i,6]))
}
for(i in 1:nrow(haplo)){
  p <- p + annotate('segment',x=haplo[i,2],y=-2,
                    xend = haplo[i,3],yend = -2, linewidth = 3,
                    #colour=ifelse(omim[i,3] == 0,'grey90','#006400'))
                    colour=f_haplo(haplo[i,6]))
}
p <- p + theme(axis.line.y = element_blank(),
               axis.title = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.major.x = element_line(colour='grey95'),
               panel.border = element_rect(fill = 'transparent',colour='black'))
if(chr %in% c(1,2)){
  p <- p + scale_x_continuous(breaks = c(0:floor(last_tick/2))*2e7,labels = paste0(c(0:floor(last_tick/2))*20,'Mb'))
}
else{
  p <- p + scale_x_continuous(breaks = 0:last_tick*1e7,labels = paste0(0:last_tick*10,'Mb'))
}

p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,y=0, label =paste0('chr',chr),colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,y=-0.5, label ='Recurrent regions',colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,y=-1, label ='RefSeq genes',colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,y=-1.5, label ='OMIM genes',colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,y=-2, label ='HI/TI genes',colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,
                  y= max(cnvr_df[which(cnvr_df$indel == 'dup'), ]$id)/2, 
                  label ='Duplications',colour='black',size=3)
p <- p + annotate('text',x=-1 * chrband[nrow(chrband),2]/20,
                  y= (min(cnvr_df[which(cnvr_df$indel == 'del'), ]$id) + max(cnvr_df[which(cnvr_df$indel == 'del'), ]$id))/2,  
                  label ='Deletions',colour='black',size=3)
print(chr)
return(p)
}

p <- draw_one_chr(1)
pdf(paste0(output_dir,'chr1_cnvr_browser_2.pdf'),width=13,height=3)
p
dev.off()

p <- draw_one_chr(1)
for(k in c(2:5)){
  p <- p / draw_one_chr(k)
}
p <- p + plot_annotation(tag_levels = 'I')
pdf(paste0(output_dir,'chr1-5_cnvr_browser_2.pdf'),width=12,height=15)
p 
dev.off()

p <- draw_one_chr(6)
for(k in c(7:10)){
  p <- p / draw_one_chr(k)
}
p <- p + plot_annotation(tag_levels = 'I')
pdf(paste0(output_dir,'chr6-10_cnvr_browser_2.pdf'),width=12,height=15)
p 
dev.off()

p <- draw_one_chr(11)
for(k in c(12:15)){
  p <- p / draw_one_chr(k)
}
p <- p + plot_annotation(tag_levels = 'I')
pdf(paste0(output_dir,'chr11-15_cnvr_browser_2.pdf'),width=12,height=15)
p 
dev.off()

p <- draw_one_chr(16)
for(k in c(17:20)){
  p <- p / draw_one_chr(k)
}
p <- p + plot_annotation(tag_levels = 'I')
pdf(paste0(output_dir,'chr16-20_cnvr_browser_2.pdf'),width=12,height=15)
p
dev.off()

p <- draw_one_chr(21)
for(k in c(22,'X')){
  p <- p / draw_one_chr(k)
}
#p <- p + plot_spacer() + plot_spacer()
p <- p + plot_annotation(tag_levels = 'I')
pdf(paste0(output_dir,'chr21-X_cnvr_browser_2.pdf'),width=12,height=9)
p 
dev.off()


#p <- draw_one_chr(1)
#for(k in c(2:12)){
#  p <- p / draw_one_chr(k)
#}
#q <- draw_one_chr(13)
#for(k in c(14:22,'X')){
#  q <- q / draw_one_chr(k)
#}
#q <- q + plot_spacer()
#pdf(paste0('C:\\Users\\huliang\\Desktop\\data and script\\plot_gviz\\chr_cnvr_browser_2.pdf'),width=24,height=30)
#p | q
#dev.off()


