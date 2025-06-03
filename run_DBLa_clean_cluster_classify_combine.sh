# change directory to the absolute path of the current working directory to
# avoid binding issues with symlinks
newhome=$(pwd -P)
cd $newhome

# import variables from yaml
yml (){
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# remove whitespace from variables
rmwt () {
   no_hash=$(echo -e $1 | sed -e 's/\#.*$//')
   echo -e $no_hash | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

# create singularity bindings
establish_binds () {
    eval $(yml snakemake_prototype.yaml)
    input_folder=$(rmwt $input_folder)
    output_folder=$(rmwt $output_folder)
    sif_file=$(rmwt $sif_file)

    singularity_bindings="-B $input_folder:/opt/input -B $output_folder:/opt/output"
}


establish_binds
mkdir -p $output_folder
singularity run $singularity_bindings $sif_file snakemake -s snakemake_prototype.smk -c $clean_cpus
