# access mysql database in master.
library(RMySQL)
library(data.table)

bind_cbs_result <- function(){
  cbs_path <- 'C:\\Users\\DELL\\Desktop\\cnv_10k\\cnv_cbs\\'
  cbs_files <- list.files(cbs_path)
  cbs_result <- data.frame()
  for(i in cbs_files){
    new <- read.csv(paste0(cbs_path,i))
    cbs_result <- rbind(cbs_result,new)
    print(i)
  }
  cbs_result$padj <- p.adjust(cbs_result$pval,'BH')
  cbs_result <- cbs_result[which(cbs_result$padj <= 1e-3), ]
  cbs_result$size <- cbs_result$loc.end - cbs_result$loc.start + 10000
  write.csv(cbs_result,'E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cbs_result.csv',row.names = F)
}

access_niptdatabase <- function(){
  con <- dbConnect(MySQL(), host='192.168.5.246', dbname="nipt2", user="win7", password="huliang")
  dbSendQuery(con,'SET NAMES gbk')
  nipt <- dbReadTable(con,'table_analysis')
  nipt <- nipt[which(nipt$reads_num >= 3.5 & nipt$GC < 42), ]
  nipt <- nipt[which(nipt$doctor_advice == '合格' & nipt$report_status == '已发送'), ]
  info <- dbReadTable(con,'table_sample')
  nipt <- merge(nipt,info,by='sample_id',all=F)
  return(nipt)
}

arrange_omim <- function(){
  omim <- read.delim("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\OMIM Genes.aed")
  omim <- omim[-c(1:9),]
  omim <- omim[ ,-c(5:6)]
  names(omim) <- c('chr','start','end','gene_name','disorder',
                   'gene_title','key','symbol','symbol_list','omim_id')
  omim_patho <- omim[which(omim$disorder != ''), ]
  omim_patho$inheritance <- '1'
  omim_patho$disorder <- gsub('\\?.*(;|$)','',omim_patho$disorder)
  omim_patho$disorder <- gsub('\\{.*(;|$)','',omim_patho$disorder)
  omim_patho$disorder <- gsub('\\[.*(;|$)','',omim_patho$disorder)
  omim_patho <- omim_patho[which(omim_patho$disorder %in% c('','[') == F), ]
  for(i in 1:nrow(omim_patho)){
    if(grepl('Autosomal recessive', omim_patho[i,5])==T & grepl('Autosomal dominant', omim_patho[i,5])==F){
      omim_patho[i,11] <- 'AR'
    }
    if(grepl('Autosomal dominant', omim_patho[i,5])==T & grepl('Autosomal recessive', omim_patho[i,5])==F){
      omim_patho[i,11] <- 'AD'
    }
    if(grepl('Autosomal dominant', omim_patho[i,5])==T & grepl('Autosomal recessive', omim_patho[i,5])==T){
      omim_patho[i,11] <- 'AD/AR'
    }
    if(grepl('Y-linked', omim_patho[i,5])==T ){
      omim_patho[i,11] <- 'YL'
    }
    if(grepl('X-linked recessive', omim_patho[i,5])==T ){
      omim_patho[i,11] <- 'XLR'
    }
    if(grepl('X-linked dominant', omim_patho[i,5])==T ){
      omim_patho[i,11] <- 'XLD'
    }
  }
  omim_patho$inheritance <- gsub('1','other',omim_patho$inheritance)
  return(omim_patho)
  write.csv(omim_patho, "E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\omim_patho.csv")
}

cnvr_exon <- function(){
  cnvr_omim$gene <- substr(cnvr_omim$omim_gene,1,regexpr('\\(',cnvr_omim$omim_gene)-2)
  known_gene <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\knownGene.txt\\knownGene.txt",data.table = F)
  ensembl <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\knownGene.txt\\Protein Coding Ensembl Genes.aed",data.table = F)
  ensembl <- ensembl[-c(1:11), ]
  ensembl_id <- ensembl[ ,c(4,10)]
  names(ensembl_id) <- c('gene_name','ensembl_id')
  ensembl_id <- ensembl_id[which(ensembl_id$gene_name %in% cnvr_omim$gene ), ]
  known_gene$ensembl_id <- substr(known_gene$V1,1,15)
  known_gene <- merge(ensembl_id,known_gene,by='ensembl_id',all.x = F)
  
  cnvr_exon <- data.frame()
  cnvr_omim$exon <- ''
  cnvr_omim$orf <- ''
  for(i in 1:nrow(cnvr_omim)){
    start <- cnvr_omim[i,3]
    end <- cnvr_omim[i,4]
    gene <- cnvr_omim[i,10]
    gene_df <- known_gene[which(known_gene$gene_name == gene), ]
    if(nrow(gene_df) > 0 ){
      gene_df <- gene_df[which(gene_df$V5 - gene_df$V4 == max(gene_df$V5 - gene_df$V4)), ]
      exon_start <- strsplit(gene_df[1,11],',')[[1]]
      exon_end <- strsplit(gene_df[1,12],',')[[1]]
      exon_df <- data.frame(start = exon_start,end = exon_end)
      if( start < exon_df[1,1] &  exon_df[nrow(exon_df),2] < end ){
        cnvr_omim[i,ncol(cnvr_omim)] <- F
      }
      else{
        cnvr_omim[i,ncol(cnvr_omim)] <- T
      }
      for(k in 1:nrow(exon_df)){
        if((exon_df[k,1] < start & start < exon_df[k,1]) |
           (exon_df[k,1] < end & end < exon_df[k,1])     |
           (exon_df[k,1] < start &  end < exon_df[k,1])  |
           (start < exon_df[k,1] &  exon_df[k,1] < end )
        ){
          print(i)
          new <- cnvr_omim[i,]
          #new$exon = paste(new$exon,k,';')
          cnvr_exon <- rbind(cnvr_exon,new)
          break
        }
        else{
          print('intron')
        }
      }
    }
  }
  return(cnvr_exon)
  write.csv(cnvr_exon,"E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnvr_exon.csv")
}

cnv_exon <- function(){
  cnv_omim$gene <- substr(cnv_omim$omim_gene,1,regexpr('\\(',cnv_omim$omim_gene)-2)
  known_gene <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\knownGene.txt\\knownGene.txt",data.table = F)
  ensembl <- fread("E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\knownGene.txt\\Protein Coding Ensembl Genes.aed",data.table = F)
  ensembl <- ensembl[-c(1:11), ]
  ensembl_id <- ensembl[ ,c(4,10)]
  names(ensembl_id) <- c('gene_name','ensembl_id')
  ensembl_id <- ensembl_id[which(ensembl_id$gene_name %in% cnv_omim$gene ), ]
  known_gene$ensembl_id <- substr(known_gene$V1,1,15)
  known_gene <- merge(ensembl_id,known_gene,by='ensembl_id',all.x = F)
  
  cnv_exon <- data.frame()
  cnv_omim$exon <- ''
  for(i in 1:nrow(cnv_omim)){
    start <- cnv_omim[i,3]
    end <- cnv_omim[i,4]
    gene <- cnv_omim[i,10]
    gene_df <- known_gene[which(known_gene$gene_name == gene), ]
    if(nrow(gene_df) > 0 ){
      gene_df <- gene_df[which(gene_df$V5 - gene_df$V4 == max(gene_df$V5 - gene_df$V4)), ]
      exon_start <- strsplit(gene_df[1,11],',')[[1]]
      exon_end <- strsplit(gene_df[1,12],',')[[1]]
      exon_df <- data.frame(start = exon_start,end = exon_end)
      if( start < exon_df[1,1] &  exon_df[nrow(exon_df),2] < end ){
        cnv_omim[i,ncol(cnv_omim)] <- F
      }
      else{
        cnv_omim[i,ncol(cnv_omim)] <- T
      }
      for(k in 1:nrow(exon_df)){
        if((exon_df[k,1] < start & start < exon_df[k,1]) |
           (exon_df[k,1] < end & end < exon_df[k,1])     |
           (exon_df[k,1] < start &  end < exon_df[k,1])  |
           (start < exon_df[k,1] &  exon_df[k,1] < end )
        ){
          print(i)
          new <- cnv_omim[i,]
          #new$exon = paste(new$exon,k,';')
          cnv_exon <- rbind(cnv_exon,new)
          break
        }
        else{
          print('intron')
        }
      }
    }
  }
  return(cnv_exon)
  write.csv(cnv_exon,"E:\\bioinformatics\\project\\NIPT\\cnv\\tables\\cnv_exon.csv")
}