library(ggthemes)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(lubridate)
library(ggpattern)
library(cowplot)


if (Sys.info()[["nodename"]] == "ukshikmb-nl056") {
  dir_to_cluster <- "/run/user/1000/gvfs/sftp:host=caucluster.rz.uni-kiel.de,user=sukmb465/"
} else if (Sys.info()[["nodename"]] == "MBP-von-ewacker.i-kmb.de") {
  dir_to_cluster <- "~/caucluster_mount/"
} else {
  dir_to_cluster <- "~/caucluster_mount/"
}

comparison_sample_ids <- list.files(path=paste0(dir_to_cluster,"work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/TOFUDEV160/maxbin2"),include.dirs = T) 

# Bin Scoring Plot: ####
checkm_paths <- vector()
for(i in seq(1:length(comparison_sample_ids))){
  checkm_paths[i] <- (paste0(dir_to_cluster,"work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/TOFUDEV160/checkm/",comparison_sample_ids[i],"/",comparison_sample_ids[i],"_checkm_table.tsv"))
}

checkm_tables <- sapply(checkm_paths, function(x) fread(x, sep = "\t", header = T) , simplify=F)  %>% rbindlist(use.names = T, idcol="Filename") %>% 
  mutate(SampleID=gsub("_checkm_table.tsv","",basename(Filename))) %>% 
  select(-Filename) %>%
  mutate(score=Completeness - 0.5*(Contamination) ) %>%
  arrange(-score) %>%
  mutate(bin_rank = row_number()) 


nfmag_path <- vector()
nfmag_path_unbinned <- vector()

nfmag_tables <- fread(paste0(dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/nfcore_mag/checkm_postprocessing/nfcore_mag_checkm_table.tsv"), sep = "\t", header=T) %>%
  mutate(SampleID=gsub("MEGAHIT-.*-","",`Bin Id`) %>% gsub("\\..*$","",.)) %>% 
  mutate(score=Completeness - 0.5*(Contamination) ) %>%
  arrange(-score) %>%
  mutate(bin_rank = row_number()) 


atlas_table <- fread(paste0(dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/ATLAS/batch1batch2_postprocessing/atlas_checkm_table.tsv")) %>%
  mutate(score=Completeness - 0.5*(Contamination) ) %>%
  arrange(-score) %>%
  mutate(SampleID=gsub("_.*$","",`Bin Id`)) %>%
  mutate(bin_rank = row_number())


metafun_table <- fread(paste0(dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/metafun/results/metagenome/bin_results/metafun_checkm_table.tsv"))  %>% 
  mutate(score=Completeness - 0.5*(Contamination) ) %>%
  arrange(-score) %>%
  mutate(SampleID=gsub("_.*$","",`Bin Id`)) %>%
  mutate(bin_rank = row_number())


color_df <- data.frame(labels=c("TOFU-MAaPO",
                                "nf-core/mag",
                                "ATLAS",
                                "metaFun"
                                ), 
                       values=c("#BEAED4",
                                "#7FC97F",
                                "#FDC086",
                                "#386CB0"
                                #brewer.pal(5, "Accent")
                                )
)

tofu_nfcore_comparison <- list(tofumaapo=checkm_tables,
                               nfmag=nfmag_tables,
                               atlas=atlas_table,
                               metaFun=metafun_table
) %>% 
  rbindlist(use.names =T, idcol="origin" ) %>%
  mutate(origin=gsub("tofumaapo","TOFU-MAaPO",origin),
         origin=gsub("nfmag","nf-core/mag",origin),
         origin=gsub("atlas","ATLAS",origin),
         origin=gsub("metafun","metaFun",origin),
         origin=factor(origin, 
                       levels=(c(
                         "TOFU-MAaPO",
                         "nf-core/mag",
                         "ATLAS",
                         "metaFun")
                       ))) %>%
  #filter out negative scores
  filter(score>0)%>%
  ggplot(.)+
  aes(x=bin_rank,y=score, colour=origin)+
  theme_bw()+
  geom_line(linewidth=1.5, alpha=0.8)+
  xlab("Bin Rank") + 
  ylab("Score") +
  labs(color="") +
  coord_cartesian(xlim = c(-1,3200),
                  ylim = c(-1, 101), expand = F) +
  scale_colour_manual("",labels=color_df$labels,
                      values = setNames(color_df$values, color_df$labels), 
                      
  )+
  theme(legend.position="inside", 
        legend.justification = c("left", "bottom"),
        legend.position.inside = c(0.005, 0.01),
        legend.box.just = "left",
        legend.box.margin = margin(0,0,5,0),
        legend.margin = margin(-10,0,0,0),
        text = element_text(size=22), 
        panel.grid = element_blank(),
        axis.text = element_text(color="black", size=22),
        plot.margin = unit(c(5.5,5.5,0,5.5),"points"),
        #panel.grid.major = element_blank()
        panel.grid.minor = element_blank()
  )+
  guides(color=guide_legend(nrow=5,byrow=TRUE,label.theme = element_text(size=20),reverse = F))+
  labs(y="Score [%]")
tofu_nfcore_comparison 
# ggsave(tofu_nfcore_comparison, file="SRP102150/100sample_sra_nfcore_tofumaapo_bin_scoring.png", width=10, height=6)
# ggsave(tofu_nfcore_comparison, file="SRP102150/100sample_sra_nfcore_tofumaapo_bin_scoring.svg", width=10, height=6)


#Runtime benchmark: ####




#### NFCORE needs the second run added to it!!!! ####
runtime_nfcore <- system(command = paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/nfcore_mag/cmdbench.log_firstrun"), intern=T) %>% 
    gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% 
    trimws() %>% 
    as.double() +
  system(command = paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/nfcore_mag/cmdbench.log"), intern=T) %>%
    gsub("Runtime: ","",.) %>% 
    gsub(" second\\(s\\)","",.) %>% 
    trimws() %>% 
    as.double()

runtime_tofu <- system(command =paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/TOFU-MAaPO/SRP102150/cmdbench_fivebinners.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()
runtime_atlas <- system(command = paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/ATLAS/batch1/cmdbench.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()+
  system(command = paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/ATLAS/batch2/cmdbench_firstrun.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()+
  system(command = paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/ATLAS/batch2/cmdbench_secondrun.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()
# runtime_atlas_cobinning <- system(command = "grep \"Runtime: \" /run/user/1000/gvfs/sftp:host=134.245.63.226,user=ewacker/dpool/ewacker/metagenomics/sra_SRP102150_100samples/atlas_cobinning/cmdbench.log", intern=T) %>% gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()


runtime_metafun <- system(command =paste0("grep \"Runtime: \" ",dir_to_cluster,"work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/metafun/cmdbench_metafun.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()+
  system(command =paste0("grep \"Runtime: \" ",dir_to_cluster,"/work_ikmb/sukmb465/projects/TOFUpaper/benchmarks/metagenomics/metafun/cmdbench_metafun2.log"), intern=T) %>% 
  gsub("Runtime: ","",.) %>% gsub(" second\\(s\\)","",.) %>% trimws() %>% as.double()

runtimes <- data.frame(Tool=factor(c("tofumaapo",
                                     "nfcore_mag",
                                     "runtime_atlas",
                                     "metafun"
                                     ),
                                   level=c("tofumaapo",
                                           "nfcore_mag",
                                           "runtime_atlas",
                                           "metafun")),#[c(2,1,5,4,3)]),
                       Runtime=c(lubridate::as.duration(runtime_tofu),
                                 lubridate::as.duration(runtime_nfcore),
                                 lubridate::as.duration(runtime_atlas),
                                 lubridate::as.duration(runtime_metafun)))


runtime_plot <- ggplot(runtimes, aes(x=Tool, y=Runtime %>% as.numeric("hours"), pattern=Tool, fill=Tool))+
  geom_bar_pattern(stat="identity",pattern_color= "#FDC086",pattern_fill="#FDC086",pattern_density = .45, pattern_angle = 45)+
  scale_fill_manual(values = color_df$values)+#[c(3,2,1,4,5)]) +
  scale_pattern_manual(values = c(nfcore_mag = "none", 
                                  tofumaapo = "none", 
                                  runtime_atlas="none",
                                  metafun="none"))+
  coord_cartesian(ylim = c(-2, 505), expand = F) +
  geom_text(aes(label=Runtime %>% as.numeric("hours") %>% round(digits=2), size=24),position = position_dodge(width = 0.8), vjust = -0.6, size=8)+
  xlab("") + 
  ylab("Runtime [h]") +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=22), 
        axis.text = element_text(color="black", size=20, angle = 00),
        panel.grid = element_blank(),
        axis.text.x = element_text(color="black", size=18, angle = 90, hjust = 0.5, vjust = 0.5))+
  scale_x_discrete(breaks = c("nfcore_mag", 
                              "tofumaapo",
                              "runtime_atlas",
                              "metafun"), 
                   label = c("nf-core/mag",
                             "TOFU-MAaPO",
                             "ATLAS",
                             "metaFun"
                   ), position = "bottom")

runtime_plot

# Combined plot ####

plot_magbenchmark <- plot_grid(
  tofu_nfcore_comparison+theme(legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.position = "bottom"),
  runtime_plot,
  nrow=1,
  ncol = 2,
  rel_heights = c(1,1),
  rel_widths = c(1,0.75),
  label_size = 24,
  labels=c("a)","b)"),
  scale = 0.90,
  align = "v",
  axis = "b",
  greedy = F
  
)


plot_magbenchmark

ggsave(filename="assembly_benchmarks/arranged_benchmark_100samples_v2.svg",plot_magbenchmark, width=12, height = 6, bg = "white" )               
ggsave(filename="assembly_benchmarks/arranged_benchmark_100samples_v2.jpg",plot_magbenchmark, width=12, height = 6, bg = "white" )  


# Return Bins per Pipeline table:
list(tofumaapo=checkm_tables,
           nfmag=nfmag_tables,
           atlas=atlas_table,
           metafun=metafun_table
      ) %>% 
      rbindlist(use.names =T, idcol="origin" ) %>%
       mutate(origin=gsub("tofumaapo","TOFU-MAaPO",origin),
              origin=gsub("nfmag","nf-core/mag",origin),
              origin=gsub("atlas","ATLAS",origin),
              origin=gsub("metafun","metaFun",origin),
                             origin=factor(origin, 
                                          levels=(c(
                                                "TOFU-MAaPO",
                                                  "nf-core/mag",
                                                "ATLAS",
                                                "metaFun")
                                            ))) %>%
       filter(score>0) %>% #filter out negative scores
  group_by(origin) %>% 
  # filter(score>0.5) %>%   #check for only medium and high quality bins
  summarise(binspertool=n())

