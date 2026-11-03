library(tidyverse)
library(readxl)
library(openxlsx)

coldata_349<-read.csv("__PATH_001__",header = T,row.names = 1)

info_149<-read_excel("__PATH_002__",sheet = 1) %>%
  dplyr::rename_with(~c("sampleid","Total_reads","Total_reads_million","filter",
                        "clean_reads","Sec_flt_reads","clean_reads_final","clean_rate",
                        "clean_bp","rRNA","rRNA_rate","miRNA","miRNA_rate","tpiRNA",
                        "tpiRNA_rate","mRNA","mRNA_rate","other_RNA","other_RNA_rate",
                        "Uniq_Genome","Uniq_Genome_rate","Multi_Genome","Multi_Genome_rate","overall")) %>%
  dplyr::select(sampleid,rRNA_rate)
info_200<-read_excel("__PATH_003__",sheet = 1) %>%
  dplyr::rename_with(~c("sampleid","Total_reads","filter",
                        "clean_reads","Sec_flt_reads","clean_reads_final",
                        "clean_bp","rRNA","rRNA_rate","miRNA","miRNA_rate","tpiRNA",
                        "tpiRNA_rate","mRNA","mRNA_rate","other_RNA","other_RNA_rate",
                        "Uniq_Genome","Uniq_Genome_rate","Multi_Genome","Multi_Genome_rate","overall"),.cols = 1:22) %>%
  dplyr::select(sampleid,rRNA_rate)



info_rRNA<-rbind(info_149,info_200) %>%
  dplyr::filter(sampleid %in% coldata_349$sample) %>%
  mutate(rRNA_rate=as.numeric(rRNA_rate)*100)

p<-ggplot(data = info_rRNA,aes(x=rRNA_rate))+
  geom_histogram(bins = 30,fill="#00BFC4",color="black")+
  theme_bw()+
  ggtitle("rRNA,YRNA,vtRNA mapping rate")+
  xlab("Mapping rate (%)")+
  ylab("Number of samples")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        axis.text = element_text(face="bold"),
        panel.border = element_blank())

ggsave(p,filename="投稿期刊/cell reports/rRNAYRNAvtRNA_mapping rate.pdf",width = 6,height = 4)

write.xlsx(info_rRNA,"投稿期刊/cell reports/rRNA_mapping rate.xlsx")

