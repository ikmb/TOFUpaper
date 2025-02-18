#!/bin/bash

log_file="batch_execution.log"

# Loop through the batch scripts
for file in batch_*.sh; do
    [ -e "$file" ] || continue  # Skip if no matching files
#    echo "Processing: $file"
    script="$file"
    echo "Executing $script..."
    
    # Attempt to run the script
    bash "$script" 2>&1
    exit_status=$?
    
    # Check if the script exited successfully
    if [ $exit_status -eq 0 ]; then
        echo "$script executed successfully." >> "$log_file"
        rsync -a --delete /work_beegfs/sukmb465/emptyfolder/ ./work/
    else
        echo "$script failed." >> "$log_file"
        i=$(echo "$script" | grep -o -E '[0-9]+')
        cp .nextflow.log "$i.log.nextflow.log"
    fi
done