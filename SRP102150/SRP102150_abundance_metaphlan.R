library(data.table)
library(tidyverse)
library(phyloseq)
library(microbiome) 
library(DESeq2)
library(vegan)
library(cowplot)
library(ggpubr)
library(Maaslin2)
aggregate_top_taxa2 <- function(x, top, level) {
  x <- aggregate_taxa(x, level)
  
  tops <- top_taxa(x, top)
  tax <- tax_table(x)
  
  inds <- which(!rownames(tax) %in% tops)
  
  tax[inds, level] <- "Other"
  
  tax_table(x) <- tax
  
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  
  aggregate_taxa(x, level)
}

setwd("/home/sukmb465/Documents/Projects/TOFU-MAaPO/SRP102150")
# mag_dt_metaphlan <- rbindlist(lapply(list.files("~/Documents/Projects/TOFU-MAaPO/SRP102150"), fread))

SRP102150_metadata <- fread("SRP102150_metadata.csv")


mag_dt_metaphlan <- fread("metaphlan_abs_abundances.txt")


# mag_abudance_wide <- pivot_wider(mag_dt_metaphlan[,-"bin"],names_from = "SampleID", values_from = "scaled_weighted_contigs_TPM_per_bin", values_fill = 0)

# fwrite(mag_abudance_wide, file="merged_MAG_tpm_abundance.tbl")

mag_dt_metaphlan <- mag_dt_metaphlan %>% select(-clade_taxid)

taxsplit = strsplit(mag_dt_metaphlan$clade_name, split="|", fixed=TRUE)
taxunitmatrix = matrix(NA, ncol=max(sapply(taxsplit, length)), nrow=length(taxsplit))
colnames(taxunitmatrix) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxunitmatrix)]
rownames(taxunitmatrix) = mag_dt_metaphlan$taxa
for (i in 1:nrow(taxunitmatrix)){
  taxunitmatrix[i, 1:length(taxsplit[[i]])] <- taxsplit[[i]]
}
taxunitmatrix = gsub("[a-z]__", "", taxunitmatrix)
# taxunitmatrix = phyloseq::tax_table(taxunitmatrix)

taxunitmatrix[taxunitmatrix[,"Genus"]=='',"Genus"] <- NA_real_
taxunitmatrix[taxunitmatrix[,"Species"]=='',"Species"] <- NA_real_


colnames(mag_dt_metaphlan) <- colnames(mag_dt_metaphlan) %>% gsub("_metaphlan","",.)

ps <- phyloseq(tax_table(taxunitmatrix), 
               otu_table(data.frame(row.names = rownames(taxunitmatrix),mag_dt_metaphlan[,!"clade_name"]), taxa_are_rows = T),
               sample_data(data.frame(row.names = SRP102150_metadata$Run, Sample = SRP102150_metadata$Run, SRP102150_metadata))
)

ps %>% aggregate_taxa(level = "Phylum") %>% # aggregate all ASVs into the level phyloum
  rarefy_even_depth() %>% # make all samples with the same sequencing depth using rarefaction
  plot_bar(x="Run", fill="Phylum") + 
  facet_wrap(~ Desc, scales = "free")


ps %>% aggregate_top_taxa2(top = 11, level = "Order") %>% #tax_table()
  #rarefy_even_depth() %>%
  plot_bar(x="Sample", fill="Order") +
  facet_wrap(~ Desc, scales = "free") +
  scale_fill_brewer(palette="Paired") +# change colors of bars
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=8,angle=90)
  )

rich <-
  ps %>% 
  rarefy_even_depth() %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  estimate_richness(measures = "Shannon") %>% # get alpha diversity
  rownames_to_column("Sample") %>% #format
  left_join(sample_data(ps) %>% data.frame(), by = "Sample") # here I am combining the output of estimate richness with the sample data of the phyloseq object. I need to transform it to data.frame first

# Check output
rich %>% head()

rich %>% 
  arrange(Desc) %>% 
  wilcox.test(Shannon ~ Desc, data = .)

#Fancy plot
ggplot(rich, aes(x = Desc, y = Shannon)) +
  geom_boxplot() +
  geom_jitter()

# check how many species can be detected:
ps %>% aggregate_taxa(level = "Species") %>% aggregate_taxa(level = "Species")


ps.ord <-
  ps %>% 
  rarefy_even_depth() %>% # let's rarefy the data to make sure that sequencing depth is similar for all samples
  aggregate_taxa(level = "Species") %>%
  ordinate("NMDS", "bray")
p1 = plot_ordination(ps, ps.ord, type="samples", color="Desc")

print(p1)

adonis2(phyloseq::distance(ps, method = "bray") ~ sample_data(ps)$Desc)

# rbind(data.table(
#         divergence=divergence(subset_samples(ps, Desc == "Case"), apply(abundances(subset_samples(ps, Desc == "Case")), 1, median), method="bray"),
#         diagnose=rep("Case",50)),
#       data.table(
#         divergence=divergence(subset_samples(ps, Desc == "Control"), apply(abundances(subset_samples(ps, Desc == "Control")), 1, median), method="bray"),
#         diagnose=rep("Control",50))
#       ) %>% ggplot(., aes(x=diagnose,y=divergence)  )+geom_boxplot()



ps.to.dseq <-
  ps %>%
  aggregate_taxa(level = "Phylum")

otu_table(ps.to.dseq) <- ps.to.dseq %>% otu_table()+1

# Create DESeq2 object from 
dseq <-
  ps.to.dseq %>% 
  phyloseq_to_deseq2(design = ~ Desc)

res <-
  DESeq(dseq)


res %>% colData %>% head()

# Extract the result table
res.df <-
  res %>% 
  results(contrast=c("Desc","Case","Control"),tidy = T)
View(res.df)
fwrite(res.df, file="deseq2_results_genus.csv")

#Visualize what we got out of it
res.df %>% head()

res.df.to.plot <-
  res.df %>% 
  filter(padj < 0.05) %>% # keep only results with adjusted P value less than 0.05
  mutate(Genus = row) %>% # Create Genus column
  left_join(tax_table(ps.to.dseq )@.Data %>% data.frame(), by = "Genus") %>% # Add taxonomy information from the phyloseq object.
  # Arrange the data for a prettier plot
  arrange(log2FoldChange) %>% 
  mutate(Genus = factor(Genus, levels = Genus %>% unique()))

head(res.df.to.plot)


#Plot

ggplot(res.df.to.plot, aes(x = log2FoldChange, y = Genus)) +
  geom_jitter(aes(col = Genus, size = baseMean))  +
  geom_vline(xintercept = 0)


firstplot <- ps %>% aggregate_taxa(level = "Phylum") %>% # aggregate all ASVs into the level phyloum
  rarefy_even_depth() %>% # make all samples with the same sequencing depth using rarefaction
  plot_bar(x="Run", fill="Phylum") + 
  facet_wrap(~ Desc, scales = "free_x")+theme_bw()+theme(legend.position = "bottom")


plot_tosave <- plot_grid(
  firstplot+theme(axis.title.x = element_blank()),
  plot_grid(plot_ordination(ps, ps.ord, type="samples", color="Desc")+ 
              # stat_ellipse(type = "norm", linetype = 2) +
              stat_ellipse(type = "t") +
              labs(color="Diagnose")+
              theme_bw(),
            ggplot(res.df.to.plot, aes(x = log2FoldChange, y = Genus)) +
              geom_point(aes(col = Genus, size = baseMean))  +
              geom_vline(xintercept = 0)+guides(col="none"),
            nrow = 1,
            ncol = 2, 
            rel_widths = c(1,1),
            labels = c("b)","c)")
  ),
  nrow=2,
  ncol = 1,
  rel_heights = c(1,1),
  labels=c("a)")
)

ggsave(filename="arranged_plot_SRP102150.svg",plot_tosave, width=12, height = 11 )

# Maaslin2 on the MetaPhlAn results

# maaslin_metaphlan_input <- mag_dt_metaphlan %>% as.data.frame() %>% column_to_rownames("clade_name")

maaslin_metaphlan_input <- ps %>%
  aggregate_taxa(level = "Genus") %>% otu_table()

maaslin_metaphlan <- Maaslin2(
  maaslin_metaphlan_input, maaslin_input_metadata, 'maaslin2_metaphlan',
  fixed_effects = c('Desc', 'sex'),
  normalization="CLR",
  #random_effects = c('Bases', 'AvgSpotLen'),
  reference = "Desc,Control",
  standardize = FALSE)

maaslin_metaphlan_input2 <- ps %>%
  aggregate_taxa(level = "Phylum") %>% otu_table()

maaslin_metaphlan <- Maaslin2(
  maaslin_metaphlan_input2, maaslin_input_metadata, 'maaslin2_metaphlan_phylum',
  fixed_effects = c('Desc', 'sex'),
  normalization="CLR",
  #random_effects = c('Bases', 'AvgSpotLen'),
  reference = "Desc,Control",
  standardize = FALSE)


# Maaslin2 on the HUMAnN Results ####

maaslin_input_data <- fread("humann_merged_pathabundance.tsv")
colnames(maaslin_input_data) <- gsub("_Abundance","",colnames(maaslin_input_data))
colnames(maaslin_input_data)[1] <- "Pathway"
# maaslin_input_data <- data.frame(row.names = maaslin_input_data$`# Pathway`,maaslin_input_data[,-1])
maaslin_input_data <- maaslin_input_data %>%
  #rename(Pathway=`# Pathway`) %>%
  mutate(Pathway=gsub("\\|.*","",Pathway)) %>%
  group_by(`Pathway`) %>%
  summarise_all(.funs = sum, na.rm = TRUE)
maaslin_input_data <-   data.frame(row.names = maaslin_input_data$Pathway, maaslin_input_data %>% select(-Pathway))
rownames(maaslin_input_data) <- gsub(" ","_",rownames(maaslin_input_data)) %>% gsub(":","",.)
  

  
maaslin_input_metadata <- data.frame(row.names=SRP102150_metadata$Run,SRP102150_metadata)


maaslin_fit_data_SRA <- Maaslin2(
  maaslin_input_data, maaslin_input_metadata, 'maaslin2_humann2',
  fixed_effects = c('Desc', 'sex'),
  #random_effects = c('Bases', 'AvgSpotLen'),
  reference = "Desc,Control",
  standardize = FALSE)



maaslin_humann_plot_SRA <- maaslin_fit_data_SRA$results %>% 
  arrange(qval) %>%
  head(4) %>% 
  pull(feature) %>%
  paste(.,collapse = "|") %>%
  grepl(.,rownames(maaslin_input_data)) %>%
  maaslin_input_data[.,] %>%
  t() %>%
  cbind(., maaslin_input_metadata  %>% mutate(Diagnose=Desc) %>% select(Diagnose) ) %>%
  pivot_longer(cols = -c(Diagnose)) %>%
  ggplot(.,aes(x=Diagnose,y=value))+
  geom_boxplot()+
  geom_jitter(alpha=0.5)+
  facet_wrap(~ name, scales = "free", labeller = label_wrap_gen(multi_line = TRUE, width = 20))+
  theme_minimal() +
  ylab("Abundance") +
  xlab(element_blank())+
  stat_pvalue_manual(data.frame(name=maaslin_fit_data$results %>% head(4) %>% pull(feature) %>% gsub("\\.","-",.) %>% gsub(".plants_and_fungi.","(plants_and_fungi)",.),#%>% gsub(" ","_",.) %>% gsub(":","",.) %>% gsub("\\.","-",.)
                                pval=paste0("q-val: ",maaslin_fit_data$results %>% head(4) %>% pull(qval) %>% signif(., 3) %>% as.character()),
                                group1 = "Case",
                                group2 = "Control",
                                y.position = c( 1300, 120000, 120000, 17000)
  ), label = "pval")
  # geom_text(data= data.frame(name=maaslin_fit_data$results %>% head(4) %>% pull(feature) %>% gsub("\\.","-",.) %>% gsub(".plants_and_fungi.","(plants_and_fungi)",.),#%>% gsub(" ","_",.) %>% gsub(":","",.) %>% gsub("\\.","-",.)
  #                            pval=maaslin_fit_data$results %>% head(4) %>% pull(pval) %>% round(., 5) %>% as.character()
  #                            ,Diagnose = "Case"
  #                            ,value = 1000
  #                            ), 
  #           aes(label = pval), nudge_y = 1, vjust = -1)

# maaslin_humann_plot

plot_to_save2 <- plot_grid(
  firstplot+
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90), legend.text = element_text(size = 9), legend.title = element_text(size = 10)
                  #, axis.ticks.y = element_blank(), axis.text.y = element_blank()
                  )+
    guides(fill=guide_legend(nrow=4,byrow=T)),
  plot_grid(plot_ordination(ps, ps.ord, type="samples", color="Desc")+ 
              stat_ellipse(type = "t") +
              labs(color="Diagnose")+
              theme_bw()+
              theme(legend.position = "bottom", legend.text = element_text(size = 11), legend.title = element_text(size = 12))+
              guides(color=guide_legend(nrow=1,byrow=T)),
            maaslin_humann_plot_SRA+
              theme_bw()+
              theme(axis.text.x = element_text(size = 12))+
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))),#+theme(plo = unit(0.1, "lines") + unit(0.5, "lines")),
            nrow = 1,
            ncol = 2, 
            rel_widths = c(0.5,1),
            labels = c("b)","c)"),
            label_size = 16
  ),
  nrow=2,
  ncol = 1,
  rel_heights = c(1,1),
  labels=c("a)"),
  label_size = 16
)

plot_to_save2

ggsave(filename="SRP102150/arranged_plot_SRP102150_2.svg",plot_to_save2, width=10, height = 10 )               


