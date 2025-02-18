# Create the destination directory if it doesn't exist
mkdir -p combined_bin_results

# Loop over all directories starting with "SRR"
for dir in ../batch2/SRR*/; do
  # Check if the binning/tool directory exists inside the current SRR directory
  if [ -d "${dir}binning/tool" ]; then
    # Copy all fasta files from binning/tool to the combined_results directory
    cp "${dir}binning/tool"/*.fasta combined_bin_results/
  fi
done

# Loop over all directories starting with "SRR"
for dir in ../batch1/SRR*/; do
  # Check if the binning/tool directory exists inside the current SRR directory
  if [ -d "${dir}binning/DASTool/bins" ]; then
    # Copy all fasta files from binning/tool to the combined_results directory
    cp "${dir}binning/DASTool/bins"/*.fasta combined_bin_results/
  fi
done