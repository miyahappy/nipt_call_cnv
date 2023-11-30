# select cnvrs which covered known pathogenic regions.
source("C:\\Users\\DELL\\Desktop\\cnv_script\\cnv_lib.R",encoding ='utf-8')

omim_patho <- arrange_omim()
cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
#cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
cnvr <- cnvr[ ,-1]
cnvr$chr <- paste0('chr',cnvr$chr)
cnvr$chr <- gsub('chr23','chrX',cnvr$chr)
n <- 0
cnvr_omim <- data.frame()
for(i in 1:nrow(cnvr)){
  cnvr_id <- paste0(cnvr[i,6],'-',cnvr[i,4],'-',cnvr[i,7])
  chr <- cnvr[i,6]
  start <- cnvr[i,1]
  end <- cnvr[i,2]
  indel <- cnvr[i,7]
  size <- cnvr[i,8]
  freq <- cnvr[i,3]
  #cnvr_range <- start:end
  search <- omim_patho[which(omim_patho$chr == chr &
                            ((omim_patho[ ,2]<= start &omim_patho[ ,3]>= start )|
                            (omim_patho[ ,2]<= end &omim_patho[ ,3]>= end )|
                            (omim_patho[ ,2]>= start &omim_patho[ ,3] <= end )|
                            (omim_patho[ ,2]<= start &omim_patho[ ,3] >= end ))
  ), ]
  if(nrow(search) > 0){
    n <- n +1
    print(n/i)
    omim_gene <- paste(search$gene_name,sep=';') 
    omim_inheritance <- paste(search$inheritance,sep=';') 
    new <- data.frame(cnvr_id =cnvr_id,chr=chr,start=start,end=end,indel=indel,
                      size=size,freq=freq,
                      omim_gene= omim_gene,
                      omim_inheritance=omim_inheritance)
    cnvr_omim <- rbind(cnvr_omim,new)
  }
}
cnvr_exon <- cnvr_exon()


cnvr_overlap_omim <- cnvr_exon[ ,c(1:7,12)]
cnvr_overlap_omim <- cnvr_overlap_omim[!duplicated(cnvr_overlap_omim), ]
length(unique(cnvr_overlap_omim$cnvr_id))/nrow(cnvr)
table(cnvr_overlap_omim$orf)

cnvr_omim_xlr <- cnvr_exon[cnvr_exon$omim_inheritance %in% c('XLR'), ]
xlr_table <- cnvr_omim_xlr[ ,c(5,7,8)]
xlr_table <- xlr_table[which(xlr_table$freq > 0),]
omim_xlr_overlapped <- data.frame()
for(i in c('del','dup')){
  for(k in unique(xlr_table$omim_gene)){
    df <- xlr_table[which(xlr_table$omim_gene == k & xlr_table$indel == i), ]
    new <-  data.frame(omim_gene = k,indel = i,freq = sum(df$freq)) 
    omim_xlr_overlapped <- rbind(omim_xlr_overlapped,new)
  }
}
omim_xlr_overlapped <- omim_xlr_overlapped[which(omim_xlr_overlapped$freq > 0), ]
omim_patho_chrx <- omim_patho[which(omim_patho$chr =='chrX' ), ]
# percentage of xlr omim genes covered by cnvrs.
length(unique(omim_xlr_overlapped$omim_gene))/length(unique(omim_patho_chrx$gene_name))
length(unique(cnvr_omim_xlr$cnvr_id))/nrow(cnvr)

cnvr_overlap_dominant <- cnvr_exon[cnvr_exon$omim_inheritance %in% c('AD','AD/AR','XLD')&
                                   cnvr_exon$indel %in% c('del','dup'), ]
table(cnvr_overlap_dominant[which(cnvr_overlap_dominant$indel =='dup'), ]$orf)
table(cnvr_overlap_dominant[which(cnvr_overlap_dominant$indel =='del'), ]$orf)

length(unique(cnvr_overlap_dominant$cnvr_id))
length(unique(cnvr_overlap_dominant$cnvr_id))/nrow(cnvr)
cnvr_overlap_dominant_del <- cnvr_overlap_dominant[which(cnvr_overlap_dominant$indel =='del'),]
length(unique(cnvr_overlap_dominant_del$cnvr_id))/nrow(cnvr[which(cnvr$indel =='del' ), ])
cnvr_overlap_dominant_dup <- cnvr_overlap_dominant[which(cnvr_overlap_dominant$indel =='dup'),]
length(unique(cnvr_overlap_dominant_dup$cnvr_id))/nrow(cnvr[which(cnvr$indel =='dup' ), ])

cnvr_overlap_dominant1 <- cnvr_overlap_dominant[ ,1:7]
cnvr_overlap_dominant1 <- cnvr_overlap_dominant1[!duplicated(cnvr_overlap_dominant1), ]
sum(cnvr_overlap_dominant1$freq)/nrow(dnacopy_output)

cnv <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\all_cnv_events_10kbins_50kb_rmdup.csv",data.table = F)
#cnvr <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr.csv")
cnv <- cnv[ ,-1]
cnv$chr <- paste0('chr',cnv$chr)
cnv$chr <- gsub('chr23','chrX',cnv$chr)
n <- 0
cnv_omim <- data.frame()
for(i in 1:nrow(cnv)){
  cnv_id <- paste0(cnv[i,2],'-',cnv[i,3],'-',cnv[i,4],'-',cnv[i,9])
  chr <- cnv[i,2]
  start <- cnv[i,3]
  end <- cnv[i,4]
  indel <- cnv[i,9]
  size <- end - start + 10000
  sample_id <- cnv[i,1]
  #cnvr_range <- start:end
  search <- omim_patho[which(omim_patho$chr == chr &
                            ((omim_patho[ ,2]<= start &omim_patho[ ,3]>= start )|
                             (omim_patho[ ,2]<= end &omim_patho[ ,3]>= end )|
                             (omim_patho[ ,2]>= start &omim_patho[ ,3] <= end )|
                             (omim_patho[ ,2]<= start &omim_patho[ ,3] >= end ))
  ), ]
  if(nrow(search) > 0){
    n <- n +1
    print(n/i)
    omim_gene <- paste(search$gene_name,sep=';') 
    omim_inheritance <- paste(search$inheritance,sep=';') 
    new <- data.frame(cnv_id =cnv_id,chr=chr,start=start,end=end,indel=indel,
                      size=size,omim_gene= omim_gene,sample_id = sample_id,
                      omim_inheritance=omim_inheritance)
    cnv_omim <- rbind(cnv_omim,new)
  }
}
write.csv(cnv_omim,"E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnv_omim.csv")
cnv_exon <- cnv_exon()
cnv_overlap_omim <- cnv_exon
cnv_overlap_omim <- cnv_overlap_omim[!duplicated(cnv_overlap_omim), ]
nrow(cnv_overlap_omim)/nrow(cnv)
length(unique(cnv_overlap_omim$sample_id))/113000


cnv_omim_xlr <- cnv_exon[cnv_exon$omim_inheritance %in% c('XLR'), ]
nrow(cnv_omim_xlr)/nrow(cnv)
length(unique(cnv_omim_xlr$sample_id))/113000

cnv_overlap_dominant <- cnv_exon[cnv_exon$omim_inheritance %in% c('AD','AD/AR','XLD')&
                                     cnv_exon$indel %in% c('del','dup'), ]
table(cnv_overlap_dominant$exon)
nrow(cnv_overlap_dominant)/nrow(cnv)
length(unique(cnv_overlap_dominant$sample_id))/113000
length(unique(cnv_overlap_dominant[which(cnv_overlap_dominant$exon==T & 
                                           cnv_overlap_dominant$indel=='dup'),8]))
length(unique(cnv_overlap_dominant[which(cnv_overlap_dominant$exon==T & 
                                           cnv_overlap_dominant$indel=='del'),8]))
length(unique(cnv_overlap_dominant[which(cnv_overlap_dominant$exon==F & 
                                           cnv_overlap_dominant$indel=='dup'),8]))
length(unique(cnv_overlap_dominant[which(cnv_overlap_dominant$exon==F & 
                                           cnv_overlap_dominant$indel=='del'),8]))

table(cnv_overlap_dominant[which(cnv_overlap_dominant$indel =='dup'), ]$exon)
table(cnv_overlap_dominant[which(cnv_overlap_dominant$indel =='del'), ]$exon)


ad_genes <- as.data.frame(table(cnv_overlap_dominant$gene))
write.csv(ad_genes,"E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\AD_genes.csv")

dosage_gene <- read.csv("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\Clingen-Dosage-Sensitivity-2023-02-09.csv",header = F)
dosage_gene <- dosage_gene[ ,c(1,3,4)]
names(dosage_gene) <- c('gene','haplo','triplo')
hs_gene <- dosage_gene[which(grepl('(Sufficient)|(Little)|(Emerging)',dosage_gene$haplo)),c(1,2)]
ts_gene <- dosage_gene[which(grepl('(Sufficient)|(Little)|(Emerging)',dosage_gene$triplo)), c(1,3)]
cnv_overlap_dosage <- merge(cnv_overlap_dominant,hs_gene,by='gene',all.x = T,all.y = F)
cnv_overlap_dosage <- merge(cnv_overlap_dosage,ts_gene,by='gene',all.x = T,all.y = F)
cnv_overlap_dosage <- cnv_overlap_dosage[which(cnv_overlap_dosage$exon == F), ]
cnv_overlap_dosage <- cnv_overlap_dosage[which((cnv_overlap_dosage$indel=='del' &
                                               cnv_overlap_dosage$haplo != '')|
                                               (cnv_overlap_dosage$indel=='dup' &
                                               cnv_overlap_dosage$triplo != '')), ]
nrow(unique(cnv_overlap_dosage[cnv_overlap_dosage$indel=='del',c(2,9)]))
nrow(unique(cnv_overlap_dosage[cnv_overlap_dosage$indel=='dup',c(2,9)]))
length(unique(cnv_overlap_dosage[cnv_overlap_dosage$indel=='del',2]))
length(unique(cnv_overlap_dosage[cnv_overlap_dosage$indel=='dup',2]))


rhd_start <- 25598981
rhd_end <- 25656936
rhd <- cnv[which(cnv$chr == 'chr1' & cnv$indel == 'del' &
                 ((cnv$start <= rhd_start & cnv$end >= rhd_end)|
                 (cnv$start >= rhd_start & cnv$end <= rhd_end)|
                 (cnv$start >= rhd_start & cnv$start <= rhd_end & cnv$end >= rhd_end)|
                 (cnv$start <= rhd_start & cnv$end <= rhd_end & cnv$end >= rhd_start)) ), ]


names(rhd) <- c('nipt_id',names(rhd)[-1])
rhd_id <- gsub('.csv','',rhd$sample_id)
library(ggplot2)
ggplot(rhd,aes(x=value)) + geom_histogram(binwidth = 0.05)

nipt <- access_niptdatabase()
nipt <- nipt[ ,c(1,42,48)]
nipt <- na.omit(nipt)
nipt_rdh <- nipt[which(nipt$nipt_id %in% rhd_id),]

write.table(nipt_rdh[ ,3],'E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\rdh_id_card.txt',
            row.names = F,col.names = F,quote = F)

all_rhd <- fread("E:\\bioinformatics\\project\\NIPT\\2021_nipt_clinic\\lis\\others\\rhd_2017-01-01.csv",data.table = F)
all_rhd <- all_rhd[all_rhd$SFZH != '', ]
all_rhd <- all_rhd[ ,c(5,7)]
all_rhd <- all_rhd[!duplicated(all_rhd), ]
rhd_cnv <- all_rhd[which(all_rhd$SFZH %in% nipt_rdh[ ,3] ), ]
names(rhd_cnv) <- c('rhd','id_cards')
rhd_cnv <- merge(rhd_cnv,nipt_rdh,by='id_cards',all=F)
rhd_cnv$nipt_id <- paste0(rhd_cnv$nipt_id,'.csv')
rhd_cnv <- merge(rhd_cnv,rhd,by='nipt_id',all = F)
ggplot(rhd_cnv,aes(x=rhd,y=value)) + geom_boxplot() +geom_jitter()

table(rhd_cnv$TESTRESULT)

all_rhd <- fread("E:\\bioinformatics\\project\\NIPT\\2021_nipt_clinic\\lis\\others\\rhd_2017-01-01.csv",data.table = F)
all_rhd <- all_rhd[all_rhd$SFZH != '', ]
all_rhd <- all_rhd[ ,c(5,7)]
all_rhd <- all_rhd[!duplicated(all_rhd), ]
neg_rhd <- all_rhd[which(grepl('阴性',all_rhd$TESTRESULT) ), ]
nipt <- access_niptdatabase()
nipt <- nipt[ ,c(1,42,48)]
nipt <- na.omit(nipt)
nipt_rdh <- nipt[which(nipt$id_cards %in% unique(neg_rhd$SFZH)),]
nipt_rdh$nipt_id <- paste0(nipt_rdh$nipt_id,'.csv')
names(cnv) <- c('nipt_id',names(cnv)[-1])
rhd_cnv <- merge(cnv,nipt_rdh,by='nipt_id',all=F)
rhd_cnv <- rhd_cnv[which(rhd_cnv$indel == 'del'), ]
length(unique(rhd_cnv$nipt_id))
rhd_cnv <- rhd_cnv[which(rhd_cnv$chr == 'chr1' & 
                           rhd_cnv$start >= 25000000 & rhd_cnv$start <= 26000000 &
                           rhd_cnv$end >= 25000000 & rhd_cnv$end <= 26000000 ), ]
length(unique(rhd_cnv$nipt_id))


