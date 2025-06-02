configfile: 'snakemake_prototype.yaml'

rule all:
	'''
	this rule currently only runs the make_desc_file and concatenate_files steps
	'''
	input:
		catted_r1=expand('/opt/output'+'/{experiment}_read1_concatenated.fastq', experiment=config['experiment_name']),

rule make_desc_file:
	'''
	this pipeline currently assumes fewer than 153 samples in the input folder.
	This assumption could probably be relaxed (e.g. by creating sample barcodes
	equal to the number of samples and writing these to MIDs.tsv, or by pairing
	forward and reverse reads together to make 153 x 153 unique combinations)
	but more testing is probably needed before implementing a solution like
	this.
	'''
	params:
		input_folder="/opt/input",
		R1_suffix=config['R1_suffix'],
		R2_suffix=config['R2_suffix']
	output:
		desc_file='/opt/output/MID_sample_mappings.desc'
	script:
		'scripts/create_desc_file.py'

rule concatenate_files:
	'''
	currently outputs un-compressed catted fastq files - not sure yet if the Day
	pipeline accepts compressed fastq files or not. Also, Day pipeline is
	currently hard-coded for specific MIDs so don't change these yet.
	'''
	input:
		mid_file=config['MID_file'],
		seekdeep_fastq_directory='/opt/input',
		desc_file='/opt/output/MID_sample_mappings.desc'
	params:
		R1_suffix=config['R1_suffix'],
		R2_suffix=config['R2_suffix']
	output:
		temporary_prepend_folder=directory('/opt/output/{experiment}_temporary_prepended_files'),
		catted_r1='/opt/output/{experiment}_read1_concatenated.fastq',
		catted_r2='/opt/output/{experiment}_read2_concatenated.fastq'
	script:
		'scripts/prepend_barcodes.py'

rule clean_dbla:
	'''
	'''
	input:
		catted_r1='/opt/output/{experiment}_read1_concatenated.fastq',
		catted_r2='/opt/output/{experiment}_read2_concatenated.fastq',
		desc_file='/opt/output/MID_sample_mappings.desc'
	params:
		output_folder='/opt/output',
		clean_percent_identity=config['clean_percent_identity'],
		clean_cpus=config['clean_cpus']
	output:
		cleaned='/opt/output'+'/{experiment}_concatenated_cleaned.fasta' #not sure if this name formatting exactly matches what clean step will output
	shell:
		'/opt/conda/envs/new_env/bin/python /opt/cleanDBLalpha.py -o {params.output_folder} -r {input.catted_r1} -R {input.catted_r2} -d {input.desc_file} --perID {params.clean_percent_identity} --cpu {params.clean_cpus} --verbose'

rule cluster_dbla:
	'''
	'''
	input:
		cleaned='/opt/output'+'/{experiment}_concatenated_cleaned.fasta' #not sure if this name formatting exactly matches what clean step will output
	params:
		output_folder='/opt/output',
		cluster_percent_identity=config['cluster_percent_identity'],
		cluster_cpus=config['cluster_cpus']
	output:
		cluster_file='/opt/output'+'/{experiment}_concatenated_cleaned_renamed_centroids.fasta'
	shell:
		'/opt/conda/envs/new_env/bin/python /opt/clusterDBLa.py -o {params.output_folder} -r {input.cleaned} --perID {params.cluster_percent_identity} --cpu {params.cluster_cpus} --verbose'

rule classify_dbla_first:
	'''
	'''
	input:
		cluster_file='/opt/output'+'/{experiment}_concatenated_cleaned_renamed_centroids.fasta'
	params:
		classify_threshold=config['classify_threshold'],
		output_folder='/opt/output'
	output:
		classify_first='/opt/output'+'/{experiment}_concatenated_cleaned_renamed_centroids_nhmmOut.txt'
	shell:
		'/opt/conda/envs/new_env/bin/python /classifyDBLalpha/reads_to_domains/allocate_reads.py -r {input.cluster_file} -E {params.classify_threshold} -o {params.output_folder} --noUproc --splitIsolates'

rule classify_dbla_second:
	'''
	'''
	input:
		classify_first='/opt/output'+'/{experiment}_concatenated_cleaned_renamed_centroids_nhmmOut.txt'
	params:
	output:
		classify_second='{experiment}_reads_to_domains.csv'
	shell:
		'/opt/conda/envs/new_env/bin/python /classifyDBLalpha/reads_to_domains/reads_to_domains.py --hmm {input.classify_first} --out {output.classify_second}'

rule combine_dbla:
	'''
	this script needs extensive modification - the underlying Rscript currently
	has hardcoded inputs and outputs that need to accept snakemake arguments instead
	'''
	input:
	params:
	output:
	shell:
		'/opt/conda/envs/new_env/bin/Rscript combine_step.R'