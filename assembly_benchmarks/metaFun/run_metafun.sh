conda activate metafun
metafun -module RAWREAD_QC --inputDir /dpool/ewacker/metagenomics/TOFU-MAaPO/SRP102150/SRP102150_raw_reads -p 62 --filter human
metafun -module ASSEMBLY_BINNING -p 62 --semibin2_mode human_gut
metafun -module BIN_ASSESSMENT -p 62 -m /dpool/ewacker/metagenomics/sample_list_raw_nfcore.csv -c 1