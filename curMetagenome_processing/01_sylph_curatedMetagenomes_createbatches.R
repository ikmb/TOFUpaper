library(tidyverse)
library(curatedMetagenomicData)
# samples_to_process <- sampleMetadata  %>% select(Run=NCBI_accession) %>% unique()


samples_to_process <-  sampleMetadata  %>% select(Run=NCBI_accession)  %>% 
  mutate(Run =  gsub(";.*", "", Run)) %>% unique()

batch_sizes <- c(rep(1000, ceiling(((nrow(samples_to_process)) - 1000) / 1000)))


remaining_rows <- nrow(samples_to_process) - sum(batch_sizes)
# batch_sizes <- seq(100, 100, by = 100)
while (remaining_rows > 0) {
  if (remaining_rows >= 1000) {
    batch_sizes <- c(batch_sizes, 1000)
    remaining_rows <- remaining_rows - 1000
  } else {
    batch_sizes <- c(batch_sizes, remaining_rows)
    remaining_rows <- 0
  }
}

# Initialize an index for the starting row of each batch
start_indices <- c(1, cumsum(batch_sizes[-length(batch_sizes)]) + 1)



# Split the dataframe into batches
batches <- lapply(start_indices, function(start) {
  samples_to_process %>%
    dplyr::slice(start:(start + batch_sizes[match(start, start_indices)] - 1))
})




# write_batch_file <- function(batch_df, batch_index) {
#   # Extract the first column and collapse into a comma-separated string
#   column_content <- paste0(
# "#!/bin/bash
# NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra \'[",paste(batch_df$Run, collapse = ","),"]\' --exact-matches --outdir batch_",batch_index," --apikey \"${NCBI_API_KEY}\" --sylph_db /dpool/ewacker/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing")
#   
#   # Write the content to a separate file named after the batch index
#   file_name <- paste0("/run/user/1000/gvfs/sftp:host=172.27.1.16,user=ewacker/dpool/ewacker/sylph_processing/batch_", batch_index, ".sh")
#   writeLines(column_content, file_name)
# }
write_batch_file <- function(batch_df, batch_index) {
  # Extract the first column and collapse into a comma-separated string
  column_content <- paste0(
    "#!/bin/bash
NXF_VER=23.10.1 ./nextflow run ikmb/tofu-maapo -r 6887e5c -profile custom -c tofu.config --sra \'[",paste(batch_df$Run, collapse = ","),"]\' --exact-matches --outdir batch_",batch_index," --apikey \"${NCBI_API_KEY}\" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing")
  
  # Write the content to a separate file named after the batch index
  file_name <- paste0("/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_sra/batch_", batch_index, ".sh")
  writeLines(column_content, file_name)
}


# Apply the function to each batch
for (i in seq_along(batches)) {
  write_batch_file(batches[[i]], i)
}



# After previous batches ran through, some IDs were still missing, will queue them up once again ####
file_table <- data.table()
for (i in seq(1,18)){
  all_profiles <- list.files(path=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_',i,"/sylph"), pattern = '\\_profile.tbl$', ignore.case = T)
  file_table <- rbind(file_table, 
                      data.table(files=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_',i,"/sylph/",all_profiles), id=gsub("_profile.tbl","",all_profiles))
  ) %>% 
    filter(id %in% samples_to_process$Run) %>%
    unique()
}

secondrun_runs <- samples_to_process$Run[!(samples_to_process$Run %in%(rbind(file_table, data.table(files=paste0('/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_',i,"/sylph/",all_profiles), id=gsub("_profile.tbl","",all_profiles))
) %>% pull(id) ))] %>% .[!is.na(.)]
  
secondru_call <- paste0("nextflow run /work_beegfs/sukmb465/TOFU-MAaPO --sra \'[",paste(secondrun_runs, collapse = ","),"]\' --outdir /work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_","secondrun_1"," --apikey \"${NCBI_API_KEY}\" --sylph_db /work_beegfs/sukmb465/projects/TOFUpaper/sylph_db/gtdb-r220-c200-dbv1.syldb --sylph_processing --exact_matches \n")
writeLines(secondru_call, paste0("/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/batch_", "secondrun_1", ".sh") )

