library(lubridate)
library(data.table)
library(tidyverse)

bench_times <- fread("/run/user/1000/gvfs/sftp:host=medcluster.medfdm.uni-kiel.de/work_beegfs/sukmb465/projects/TOFUpaper/sylph_processing/duration_run.csv",header=F)
time_vec <- bench_times %>% pull(V2) %>% .[-1] %>% .[!(. == "")] # %>% lubridate::hm() #gsub("H ","",.) %>% gsub("M 
total_duration <- sum(
  sapply(time_vec, function(x) {
    if (is.na(x)) return(dminutes(0))  # Handle NA values
    
    h <- as.numeric(str_extract(x, "\\d+(?=h)"))
    m <- as.numeric(str_extract(x, "\\d+(?=m)"))
    
    # Set NA values to 0
    h <- ifelse(is.na(h), 0, h)
    m <- ifelse(is.na(m), 0, m)
    return(h * 60 + m)  # Convert to duration
  })
)

total_hours <- as.numeric(total_duration) %/% 60
total_minutes <- as.numeric(total_duration) %% 60

# Output result
sprintf("%dH %dM", total_hours, total_minutes)

sum(dminutes(
as.numeric(str_extract(time_vec[10], "\\d+(?=h)"))*60 + as.numeric(str_extract(time_vec[10], "\\d+(?=m)")) )
)
