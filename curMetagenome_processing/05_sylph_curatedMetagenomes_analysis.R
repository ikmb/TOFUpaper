library(data.table)
library(tidyverse)
library(viridis)
library(Matrix)
library(irlba)
# Read and bind all sylph profile files into one data object
file_table <- data.table()
for (i in seq(1,18)){
  all_profiles <- list.files(path=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_',i,"/sylph"), pattern = '\\_profile.tbl$', ignore.case = T)
  file_table <- rbind(file_table, 
                      data.table(files=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_',i,"/sylph/",all_profiles), id=gsub("_profile.tbl","",all_profiles))
  ) %>% 
    filter(id %in% samples_to_process$Run) %>%
    unique()
}

all_profiles_secondrun <- list.files(path=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_secondrun_1/sylph'), pattern = '\\_profile.tbl$', ignore.case = T)
file_table <- rbind(file_table, 
                    data.table(files=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_secondrun_1/sylph/',all_profiles_secondrun), id=gsub("_profile.tbl","",all_profiles_secondrun))
)

file_table <- file_table %>% 
  filter(id %in% samples_to_process$Run) %>%
  unique()

temp <- file_table$files %>% 
  lapply(., fread, sep="\t")
data <- rbindlist( temp )

data_ann <- data %>%
  select(Sample=Sample_file, Taxonomic_abundance, Contig_name) %>% 
  mutate(Sample=gsub("_[0-9]_raw.fastq.gz","", Sample)) %>%
  pivot_wider(id_cols=Sample,names_from=Contig_name, values_from = Taxonomic_abundance)

## Retrieve contig to taxonomic name table for gtdb database
#wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz
gtdb <- fread("/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/ikmb_repository/databases/GTDB-TK/release220/taxonomy/gtdb_taxonomy.tsv",sep="\t", header=F)


# Add taxonomic info to abundance table
data_withtax <- data %>% 
  mutate(V1=gsub("gtdb_genomes_reps_r220/database/[A-Z]+/[0-9]+/[0-9]+/[0-9]+/","",Genome_file) %>% 
           gsub("_genomic.fna.gz","",.)) %>% 
  left_join(., gtdb %>% mutate(V1=gsub("^[A-Z]+_","",V1))) %>%
  separate(V2, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = ";")


fwrite(data_withtax, file="sylph_data_with_tax.csv")



## Aitchison PCA after Maltes fast CLR and PCA script:

pos_ind <- data_withtax %>%
  select(Sample=Sample_file, Taxonomic_abundance, Contig_name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  mutate(Sample=gsub("_[0-9]_raw.fastq.gz","", Sample)) %>% 
  mutate(Sample=gsub("_unpaired_raw.fastq.gz","", Sample)) %>%
  select(Sample, Tax=Genus, Taxonomic_abundance) %>% mutate(Tax=gsub("_[A-Z]+$","",Tax)) %>%
  group_by(Sample,Tax) %>% 
  summarise(total_tax_abundance_for_taxlevel= sum(Taxonomic_abundance, na.rm=T)) %>% 
  ungroup() %>%
  pivot_wider(id_cols=Sample, names_from=Tax, values_from = total_tax_abundance_for_taxlevel ) %>% #colnames() %>% .[grepl("Bacter",.)]
  select(Sample, g__Bacteroides, g__Prevotella) %>%
  mutate_if(is.numeric, ~replace_na(., 1*10**-6)) %>% 
  mutate(logbactoprev=log10(g__Bacteroides / g__Prevotella) ) %>%
  select(Sample, logbactoprev) %>%
  ungroup() %>%
  mutate(pos=seq(1:16462)) %>%
  arrange(desc(abs(logbactoprev)))

all_genus <- data_withtax %>%
  select(Sample=Sample_file, Tax=Genus, Taxonomic_abundance) %>%
  mutate(Sample=gsub("_[0-9]_raw.fastq.gz","", Sample)) %>%
  mutate(Tax=gsub("_[A-Z]+$","",Tax)) %>%
  dplyr::rename(sample=Sample,clade=Tax) %>% 
  group_by(sample, clade) %>% 
  summarise(Taxonomic_abundance=sum(Taxonomic_abundance)) %>%
  ungroup()
totgen=all_genus$clade %>% unique() %>% length()

all_genus_long_clr = all_genus %>%
  mutate(val_ln = log(10000*Taxonomic_abundance+1))
all_genus_long_clr_smry = all_genus_long_clr %>%
  group_by(sample) %>%
  summarise(ln_cum = sum(val_ln), n=totgen) %>%
  mutate(geomean=exp(ln_cum/n))
all_genus_long_clr = all_genus_long_clr %>%
  left_join(all_genus_long_clr_smry) %>%
  mutate(clr = log(((10000*Taxonomic_abundance)+1)/geomean)) %>%
  mutate(clr = ifelse(clr<0,0,clr)) %>%
  select(sample, clade, clr)#batch,


all_genus_long_clr_for_SM = all_genus_long_clr %>%
  mutate(samplef = as.factor(sample), cladef=as.factor(clade)) 

M <- with(all_genus_long_clr_for_SM, sparseMatrix(i=as.numeric(samplef),
                                                  j=as.numeric(cladef),
                                                  x=as.numeric(clr),
                                                  dimnames=list(levels(samplef), levels(cladef))))


pc <- M %*% irlba(M, nv=5, nu=0, center=colMeans(M), right_only=TRUE)$v
pc2 = irlba::prcomp_irlba(M)


# pca_plot <- pc %>%
#   as.matrix() %>%
#   as.data.frame() %>%
#   rownames_to_column("Sample") %>%
#   mutate(Sample=gsub("_unpaired_raw.fastq.gz","",Sample)) %>%
#   left_join(pos_ind %>% mutate(Sample=gsub("_unpaired_raw.fastq.gz","",Sample)), by="Sample") %>%
#   arrange((abs(logbactoprev))) %>%
#   ggplot()+
#   geom_point(aes(x=V1,y=V2, color=logbactoprev))+
#   scale_color_viridis_c(name = expression(paste("log"["10"],bgroup("(",over("Bacteroides","Prevotella"),")"))))+
#   theme_bw() +
#   theme(text = element_text(size=14),
#         legend.direction = "horizontal",
#         axis.text = element_text(color="black", size=14),
#         panel.grid=element_blank(),
#         legend.position = 'bottom')

pca_plot <- pc2$x %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("pos") %>%
  # mutate(Sample=gsub("_unpaired_raw.fastq.gz","",Sample)) %>%
  mutate(pos=as.numeric(pos)) %>%
  left_join(pos_ind %>% mutate(Sample=gsub("_unpaired_raw.fastq.gz","",Sample)), by="pos") %>%
  arrange((abs(logbactoprev))) %>%
  #filter(abs(logbactoprev)>=3) %>%
  ggplot()+
  geom_point(aes(x=PC1,y=PC2, color=logbactoprev))+
  scale_color_viridis_c(name = expression(paste("log"["10"],bgroup("(",over("Bacteroides","Prevotella"),")"))))+
  theme_bw() +
  #guides(color=guide_legend(title="Log10(Bacteroides/Prevotella)"))+
  theme(text = element_text(size=14),
        legend.direction = "horizontal",
        axis.text = element_text(color="black", size=14),
        panel.grid=element_blank(),
        legend.position = 'bottom')+
  labs(x=paste0("PCA1 (",( summary(pc2)$importance[2,1] * 100)%>% round(.,digits=1),"%)"), 
       y=paste0("PCA2 (",( summary(pc2)$importance[2,2] * 100) %>% round(.,digits=1),"%)"))


pca_plot
ggsave(filename="sylph/pcaplot_17k_v2.svg",pca_plot, width=6, height = 6 )               
ggsave(filename="sylph/pcaplot_17k_v2.png",pca_plot, width=6, height = 6 )  
ggsave(filename="sylph/pcaplot_17k_v2.pdf",pca_plot, width=6, height = 6 )  
