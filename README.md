# Clone this repository into a folder of your choosing on a linux machine that has singularity installed.
```sh
git clone https://github.com/simkin-bioinformatics/dbla_day_pipeline.git
```

# Build a sif file if you don't have one yet
```sh
sudo singularity build DBLa_clean_cluster_classify_combine.sif DBLa_clean_cluster_classify_combine.def
```

# Edit the config file
Open snakemake_prototype.yaml and edit it.  For initial testing purposes none of the values need to be changed

Note that you can't currently choose an input folder with more than 153 samples

# Run the bash script
```sh
bash run_DBLa_clean_cluster_classify_combine.sh
```
