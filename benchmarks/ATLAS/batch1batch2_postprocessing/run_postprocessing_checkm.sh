singularity shell -B /dpool /dpool/ewacker/metagenomics/TOFU-MAaPO/quickstart/singularity_cache/quay.io-biocontainers-checkm-genome-1.1.3--py_1.img

checkm lineage_wf -t 60 -x fasta /dpool/ewacker/metagenomics/ATLAS/batch1batch2_postprocessing/combined_bin_results . --tab_table --file atlas_checkm_table.tsv
