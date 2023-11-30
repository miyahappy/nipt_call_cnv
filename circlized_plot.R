library(circlize)
library(data.table)

cen_position <- read.delim("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cytoBandIdeo.txt",header = F)
cen_position <- cen_position[which(cen_position$V5 %in% c('acen','gvar','stalk')), ]
masked <- cen_position
masked$type <- 'chr_acen'
masked[nrow(masked)+1, ] <- c('chr6',29794756,33086926,'p','hla','hla')
masked <- masked[ ,c(1:3,6)]
names(masked) <- c('chr','start','end','type')
common_1kg <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\common_1000g.bed",data.table = F)
common_1kg <- common_1kg[which(common_1kg$V14 %in%c('0.05 to 0.1','0.1 to 0.2','0.2 to 0.5','Over 0.5')), ]
common_1kg <- common_1kg[which(common_1kg$V11 %in% c("alu insertion","insertion","line1 insertion")==F), ]
common_1kg <- common_1kg[ ,1:3]
common_1kg[ ,4] <- 'common'
names(common_1kg) <- c('chr','start','end','type')
masked <- rbind(masked,common_1kg)
masked$start <- as.integer(masked$start)
masked$end <- as.integer(masked$end)

all_cnv_freq <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnv_bins_freq.csv")
all_cnv_freq <- all_cnv_freq[ ,-1]
dup_freq <- all_cnv_freq[which(all_cnv_freq$indel == 'dup'), ]
del_freq <- all_cnv_freq[which(all_cnv_freq$indel == 'del'), ]
bed_dup <- data.frame(chr = paste0('chr',dup_freq$chr),
                      start= dup_freq$start,end = dup_freq$start + 10000,
                      value1 = dup_freq$Freq)
bed_dup <- bed_dup[order(bed_dup$chr,bed_dup$start), ]
bed_dup$id <- paste(bed_dup$chr,bed_dup$start,sep='_')
bed_del <- data.frame(chr = paste0('chr',del_freq$chr),
                      start= del_freq$start,
                      end = del_freq$start + 10000,
                      value1 = del_freq$Freq)
bed_del <- bed_del[order(bed_del$chr,bed_del$start), ]
bed_del$id <- paste(bed_del$chr,bed_del$start,sep='_')
#bed_del$value1 <- -1 * bed_del$value1
bed_all <- merge(bed_del,bed_dup,by='id',all = T)
#bed_all <- rbind(bed_dup,bed_del)
#
bed_all[ ,2][which(is.na(bed_all[ ,2]))] <- bed_all[ ,6][which(is.na(bed_all[ ,2]))]
bed_all[ ,3][which(is.na(bed_all[ ,3]))] <- bed_all[ ,7][which(is.na(bed_all[ ,3]))]
bed_all[ ,4][which(is.na(bed_all[ ,4]))] <- bed_all[ ,8][which(is.na(bed_all[ ,4]))]
bed_all <- bed_all[ ,c(1,3:5,9)]
names(bed_all) <- c('chr','start','end','value1','value2')
bed_all$chr <- substr(bed_all$chr,1,regexpr('_',bed_all$chr)-1)
bed_all$value1 <- bed_all$value1/113000
bed_all$value2 <- bed_all$value2/113000
bed_all$value1 <- ifelse(bed_all$value1 > 0.01,0.01,bed_all$value1)
bed_all$value2 <- ifelse(bed_all$value2 > 0.01,0.01,bed_all$value2)
bed_all$value1 <- -1*bed_all$value1
bed_all$chr <- gsub('chr23','chrX',bed_all$chr)
masked_row <- c()
#bed_all <- bed_all[which(bed_all$chr =='chr22'), ]
for(i in 1:nrow(bed_all)){
  select <- masked[which(masked$chr == bed_all[i,1]&
                         masked$start <= bed_all[i,2]&
                         masked$end >= bed_all[i,3]), ]
  if(nrow(select)>0){
    masked_row <- c(masked_row,i)
    #print(bed_all[i,2])
    #print(nrow(select))
  }
  print(i)
}
bed_all <- bed_all[-masked_row,]
write.csv(bed_all,"E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnv_bins_freq_masked.csv")

#bed_all[is.na(bed_all)] <- 0 
#bed_all_min <- bed_all[which((is.na(bed_all$value1) | bed_all$value1 >= -0.001) &
#                               (is.na(bed_all$value2) | bed_all$value2 <= 0.001) ), ]
#bed_all_max <- bed_all[which(bed_all$value1 <= -0.001 | bed_all$value2 >= 0.001), ]

#bed_all[which(bed_all$chr == 'chr16' & bed_all$start )]
#bed_all_zoomed <- bed_all[which(bed_all$chr == 'chr16' & 
  #                                bed_all$start >=1e7 & bed_all$end <= 1.1e7 ), ]
#bed_all_zoomed$name <- 'chr16_zoom'
#zoom_data <- bed_all_zoomed[ ,c(6,2,1)]
#zoom_sector <- data.frame(name = 'chr16_zom',start = 1e7,end = 1.1e7,cate ='chr16_zoom' )
bed_del <- bed_all[ ,1:4]
bed_dup <- bed_all[ ,c(1:3,5)]
bed_dup_low <- bed_dup[which(bed_dup$value2 < 0.001), ]
bed_dup_high <- bed_dup[which(bed_dup$value2 > 0.001), ]
bed_del_low <- bed_del[which(bed_del$value1 > -0.001), ]
bed_del_high <- bed_del[which(bed_del$value1 < -0.001), ]
#bed_dup_na.omit <- na.omit(bed_dup)
#bed_del_na.omit <- na.omit(bed_del)
#load(system.file(package = "circlize", "extdata", "DMR.RData"))
#head(DMR_hyper)
cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
#cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
cnvr <- cnvr[ ,-1]
cnvr$chr <- gsub('23','X',cnvr$chr)
bed_cnvr <- data.frame(chr = paste0('chr',cnvr$chr),start = cnvr$start,end = cnvr$end,value1=cnvr$Freq)
cnv <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\all_cnv_events_10kbins_50kb_rmdup.csv",data.table = F)
#cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
cnv <- cnv[ ,-1]
cnv$chr <- paste0('chr',cnv$chr)
cnv$chr <- gsub('chr23','chrX',cnv$chr)
bed_cnv <- data.frame(chr = cnv$chr,start = cnv$start,end = cnv$end,indel=cnv$indel)
bed_cnv_del <- bed_cnv[which(bed_cnv$indel == 'del'),1:3]


files <- list.files('C:\\Users\\DELL\\Desktop\\count_10k\\')
#count <- data.frame()
path <- paste0('C:\\Users\\DELL\\Desktop\\count_10k\\',files[1])
count <- fread(path,data.table = F)[ ,1:4]
for(i in files[-1]){
  path <- paste0('C:\\Users\\DELL\\Desktop\\count_10k\\',i)
  new <- fread(path,data.table = F)[ ,4]
  count <- cbind(count,new)
  print(i)
}
total <- apply(count[ ,4:ncol(count)],2,sum) 
count[ ,4:ncol(count)] <- count[ ,4:ncol(count)]/total
#covarage <- data.frame(chr= count$V1,start = count$V2 ,end=count$V3,
#                       value1 = apply(count[ ,4:ncol(count)],1,mean))
#covarage <- covarage[which(covarage$value1 > 0), ]
covarage <- data.frame(chr= count$V1,start = count$V2 ,end=count$V3,
                       value1 = apply(count[ ,4:ncol(count)],1,mean))
covarage <- covarage[which(covarage$value1 > 0), ]
covarage$value1 <- (covarage$value1 * 35)/10000
covarage$value1 <- ifelse(covarage$value1 > 0.21,0.21,covarage$value1)


circos.clear()
f1 <- function(){
  circos.par("start.degree" = 90,gap.after = c(rep(0.5, 22), 10),cell.padding = c(0.001, 1e7, 0.001, 1e7))
  circos.initializeWithIdeogram(track.height=0.3,species = "hg19",
                                chromosome.index = c(paste0('chr',1:22),'chrX'),
                                ideogram.height = convert_height(3.5, "mm"))
  circos.genomicTrack(bed_dup_high,  track.height = 0.08,ylim = c(0.001,0.011),
                      bg.lwd=0.3,
                      panel.fun = function(region, value, ...) {
                        for(h in c(0.005,0.01)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "#AAAAAA")
                        }    
                        circos.genomicPoints(region, value, pch = 16, cex = 0.1,col=c('blue'))
                      },  track.margin = c(0.0, 0))
  circos.yaxis(side = "left", at = c(0.001,0.005,0.01), lwd =0.6,
               labels = c('0.1%','0.5%','>1%'),
               tick.length = convert_x(1, "mm"), 
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  circos.genomicTrack(bed_dup_low,track.height = 0.08,ylim = c(-0.00005,0.001),
                      bg.lwd=0.3,
                      panel.fun = function(region, value, ...) {
                        for(h in c(0.0005)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "#AAAAAA")
                        }    
                        circos.genomicPoints(region, value, pch = 16, cex = 0.1,col=c('blue'))
                      },  track.margin = c(0.005, 0.0))
  circos.yaxis(side = "left", at = c(0,0.0005), lwd =0.6,
               labels = c('0','0.05%'),
               tick.length = convert_x(1, "mm"), 
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  ##############################
  circos.genomicTrack(bed_del_low,  track.height = 0.08,ylim = c(-0.001,0.00005),
                      bg.lwd=0.3,
                      panel.fun = function(region, value, ...) {
                        for(h in c(-0.0005)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "#AAAAAA")
                        }    
                        circos.genomicPoints(region, value, pch = 16, cex = 0.1,col=c('red'))
                      },  track.margin = c(0.0, 0))
  circos.yaxis(side = "left", at = seq(-0.001, 0, by = 0.0005),lwd =0.6,
               labels = c('0.1%','0.05%','0'),
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  circos.genomicTrack(bed_del_high,  track.height = 0.08,ylim = c(-0.011,-0.001),
                      bg.lwd=0.3,
                      panel.fun = function(region, value, ...) {
                        for(h in seq(-0.01, -0.001, by = 0.005)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "#AAAAAA")
                        }    
                        circos.genomicPoints(region, value, pch = 16, cex = 0.1,col=c('red'))
                      },  track.margin = c(0.005, 0))
  circos.yaxis(side = "left", at = seq(-0.01, -0.001, by = 0.005),lwd =0.6, 
               labels = c('>1%','0.5%','0.1%'),
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  circos.genomicTrack(covarage, track.height = 0.08,ylim = c(0,0.21),
                      bg.lwd=0.3,
                      panel.fun = function(region, value, ...) {
                        circos.genomicPoints(region, value, pch = 16, cex = 0.05,col=c('orange'))
                      },  track.margin = c(0.0, 0))
  circos.yaxis(side = "left", at = seq(0,0.2, by = 0.1),lwd =0.6,
               labels = c('0','0.1','0.2'),
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  circos.genomicTrack(masked[ ,1:3], panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "red", border = NA, ...)
    
  })
}
f2 <- function(){
  #circos.par("start.degree" = 90)
  #circos.par("start.degree" = 90)
  circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0))
  circos.initializeWithIdeogram(track.height=0.3,species = "hg19",
                                chromosome.index = 'chr16',
                                ideogram.height = convert_height(3.5, "mm"))
  #text(0, 0, 'freq', cex = 1)
  circos.track(bed_all_zoomed,  track.height = 0.15,ylim = c(-0.01,0.01),
               panel.fun = function(region, value, ...) {
                 for(h in seq(-0.01, 0.01, by = 0.005)) {
                   circos.lines(CELL_META$cell.xlim, c(h, h), lty = 2, col = "#AAAAAA")
                 }    
                 circos.points(region, value, pch = 16, cex = 0.2,col=c('red','blue'))
               },  track.margin = c(0.0, 0))
  circos.yaxis(side = "left", at = seq(-0.01, 0.010, by = 0.005), 
               labels = c('>1%','0.5%','0','0.5%','>1%'),
               sector.index = get.all.sector.index()[1], labels.cex = 0.5)
  
}
pdf(paste0('E:\\bioinformatics\\project\\NIPT\\cnv\\circlized\\circlized_20230404.pdf'),width=10,height=10)
correspondance <- data.frame(cate =  'chr16',
                             start = 1e7,
                             end   = 1.1e7,
                             name = 'chr16',
                             start.1 = 1e7,
                             end.1 = 1.1e7)
#circos.nested(f1(), f2(), correspondance)
f1()
dev.off()




