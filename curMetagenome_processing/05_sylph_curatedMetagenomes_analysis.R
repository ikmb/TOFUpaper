library(data.table)
library(tidyverse)
library(viridis)
library(Matrix)
library(irlba)
library(curatedMetagenomicData)
samples_to_process <-  sampleMetadata  %>% select(Run=NCBI_accession)  %>% 
  mutate(Run =  gsub(";.*", "", Run)) %>% unique()

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

#write out sample metadata for profiled samples:
sampleMetadata %>% dplyr::rename(id=NCBI_accession) %>% left_join(data_ann %>% select(id=Sample), .,by="id") %>% select(id, PMID, gender, age) %>% mutate(id=gsub("_unpaired_raw.fastq.gz","",id)) %>% fwrite("curatedmetagenomics_sample_list_sex_age.csv")
sampleMetadata %>% dplyr::rename(id=NCBI_accession) %>% left_join(data_ann %>% select(id=Sample), .,by="id") %>% select(gender, age) %>% group_by(gender) %>% summarise(age=mean(age, na.rm = T))


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
ggsave(filename="pcaplot_17k_v2.svg",pca_plot, width=6, height = 6 )               
ggsave(filename="pcaplot_17k_v2.png",pca_plot, width=6, height = 6 )  
ggsave(filename="pcaplot_17k_v2.pdf",pca_plot, width=6, height = 6 )  


sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 24.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Berlin
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] curatedMetagenomicData_3.8.0   TreeSummarizedExperiment_2.8.0 Biostrings_2.68.1              XVector_0.40.0                 SingleCellExperiment_1.22.0    SummarizedExperiment_1.30.2   
# [7] Biobase_2.60.0                 GenomicRanges_1.52.0           GenomeInfoDb_1.36.2            IRanges_2.34.1                 S4Vectors_0.38.1               BiocGenerics_0.46.0           
# [13] MatrixGenerics_1.12.3          matrixStats_1.3.0              irlba_2.3.5.1                  Matrix_1.6-5                   viridis_0.6.5                  viridisLite_0.4.2             
# [19] lubridate_1.9.4                forcats_1.0.0                  stringr_1.5.1                  dplyr_1.1.4                    purrr_1.0.2                    readr_2.1.5                   
# [25] tidyr_1.3.1                    tibble_3.2.1                   ggplot2_3.5.1                  tidyverse_2.0.0                data.table_1.15.4             
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.16.0             jsonlite_1.8.8                MultiAssayExperiment_1.26.0   magrittr_2.0.3                ggbeeswarm_0.7.2              fs_1.6.4                     
# [7] zlibbioc_1.46.0               vctrs_0.6.5                   memoise_2.0.1                 DelayedMatrixStats_1.22.6     RCurl_1.98-1.14               htmltools_0.5.8.1            
# [13] S4Arrays_1.0.6                AnnotationHub_3.8.0           curl_5.2.1                    BiocNeighbors_1.18.0          plyr_1.8.9                    DECIPHER_2.28.0              
# [19] cachem_1.1.0                  mime_0.12                     lifecycle_1.0.4               pkgconfig_2.0.3               rsvd_1.0.5                    R6_2.5.1                     
# [25] fastmap_1.2.0                 GenomeInfoDbData_1.2.10       shiny_1.8.1.1                 digest_0.6.36                 colorspace_2.1-0              AnnotationDbi_1.62.2         
# [31] scater_1.28.0                 ExperimentHub_2.8.1           RSQLite_2.3.7                 vegan_2.6-6.1                 beachmat_2.16.0               filelock_1.0.3               
# [37] fansi_1.0.6                   timechange_0.3.0              mgcv_1.9-1                    httr_1.4.7                    abind_1.4-5                   compiler_4.3.3               
# [43] bit64_4.0.5                   withr_3.0.0                   BiocParallel_1.34.2           DBI_1.2.3                     MASS_7.3-60                   rappdirs_0.3.3               
# [49] DelayedArray_0.26.7           permute_0.9-7                 tools_4.3.3                   vipor_0.4.7                   beeswarm_0.4.0                ape_5.8                      
# [55] interactiveDisplayBase_1.38.0 httpuv_1.6.15                 glue_1.8.0                    nlme_3.1-165                  promises_1.3.0                grid_4.3.3                   
# [61] mia_1.8.0                     cluster_2.1.6                 reshape2_1.4.4                generics_0.1.3                gtable_0.3.5                  tzdb_0.4.0                   
# [67] hms_1.1.3                     BiocSingular_1.16.0           ScaledMatrix_1.8.1            utf8_1.2.4                    ggrepel_0.9.6                 BiocVersion_3.17.1           
# [73] pillar_1.9.0                  yulab.utils_0.1.7             later_1.3.2                   splines_4.3.3                 BiocFileCache_2.8.0           treeio_1.24.3                
# [79] lattice_0.22-6                bit_4.0.5                     tidyselect_1.2.1              DirichletMultinomial_1.42.0   scuttle_1.10.2                gridExtra_2.3                
# [85] stringi_1.8.4                 lazyeval_0.2.2                yaml_2.3.8                    codetools_0.2-20              BiocManager_1.30.23           cli_3.6.3                    
# [91] xtable_1.8-4                  munsell_0.5.1                 Rcpp_1.0.12                   dbplyr_2.3.0                  png_0.1-8                     parallel_4.3.3               
# [97] assertthat_0.2.1              blob_1.2.4                    sparseMatrixStats_1.12.2      bitops_1.0-7                  decontam_1.20.0               tidytree_0.4.6               
# [103] scales_1.3.0                  crayon_1.5.3                  rlang_1.1.4                   KEGGREST_1.40.0 