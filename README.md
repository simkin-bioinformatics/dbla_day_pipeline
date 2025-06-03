# Build a sif file if you don't have one yet
'''bash
sudo singularity build DBLa_clean_cluster_classify_combine.sif DBLa_clean_cluster_classify_combine.def
'''

# Edit the config file
Open snakemake_prototype.yaml and edit it.  For initial testing purposes none of the values need to be changed

# Run the bash script
'''sh
bash run_DBLa_clean_cluster_classify_combine.sh
'''