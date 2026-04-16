# Create the destination directory if it doesn't exist
mkdir -p /dpool/ewacker/metagenomics/TOFU-MAaPO/SRP102150/metafun/results/metagenome/bin_results
cd /dpool/ewacker/metagenomics/TOFU-MAaPO/SRP102150/metafun/results/metagenome/bin_results

cp /dpool/ewacker/metagenomics/TOFU-MAaPO/SRP102150/metafun/results/metagenome/BIN_ASSESSMENT/results/metagenome/BIN_ASSESSMENT/bins_quality_passed/*.fa .

singularity shell -B /dpool /dpool/ewacker/metagenomics/TOFU-MAaPO/quickstart/singularity_cache/quay.io-biocontainers-checkm-genome-1.1.3--py_1.img

checkm lineage_wf -t 60 -x fa /dpool/ewacker/metagenomics/TOFU-MAaPO/SRP102150/metafun/results/metagenome/bin_results . --tab_table --file metafun_checkm_table.tsv
