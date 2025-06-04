configfile: '/opt/working_dir/snakemake_prototype.yaml'

rule all:
	'''
	this rule currently only runs the make_desc_file and concatenate_files steps
	'''
	input:
		#cleaned=expand('/opt/output/{experiment}_DBLa_cleaned.fasta', experiment=config['experiment_name']) #not sure if this name formatting exactly matches what clean step will output
		#cluster_file=expand('/opt/output/{experiment}_DBLa_cleaned_renamed_centroids.fasta', experiment=config['experiment_name'])
		#classify_first=expand('/opt/output/{experiment}_DBLa_cleaned_renamed_centroids_nhmmOut.txt', experiment=config['experiment_name'])
		#classify_second=expand('/opt/output/{experiment}_reads_to_domains.csv', experiment=config['experiment_name'])
		combined_dbla_file=expand('/opt/output/{experiment}_DBLa_binary_ups.csv', experiment=config['experiment_name'])

rule make_desc_file:
	'''
	this pipeline currently assumes fewer than 153 samples in the input folder.
	This assumption could probably be relaxed (e.g. by creating sample barcodes
	equal to the number of samples and writing these to MIDs.tsv, or by pairing
	forward and reverse reads together to make 153 x 153 unique combinations)
	but more testing is probably needed before implementing a solution like
	this.
	'''
	input:
		mid_file=config['MID_file']
	params:
		input_folder='/opt/input',
		R1_suffix=config['R1_suffix'],
		R2_suffix=config['R2_suffix']
	output:
		desc_file='/opt/output/{experiment}_MID_sample_mappings.desc'
	script:
		'/opt/working_dir/scripts/create_desc_file.py'

rule concatenate_files:
	'''
	currently outputs un-compressed catted fastq files - not sure yet if the Day
	pipeline accepts compressed fastq files or not. Also, Day pipeline is
	currently hard-coded for specific MIDs so don't change these yet.
	'''
	input:
		mid_file=config['MID_file'],
		seekdeep_fastq_directory='/opt/input',
		desc_file='/opt/output/{experiment}_MID_sample_mappings.desc'
	params:
		R1_suffix=config['R1_suffix'],
		R2_suffix=config['R2_suffix']
	output:
		temporary_prepend_folder=temp(directory('/opt/output/{experiment}_temporary_prepended_files')),
		catted_r1='/opt/output/{experiment}_concatenated_R1.fastq',
		catted_r2='/opt/output/{experiment}_concatenated_R2.fastq',
	script:
		'/opt/working_dir/scripts/prepend_barcodes.py'

rule clean_dbla:
	'''
	'''
	input:
		catted_r1='/opt/output/{experiment}_concatenated_R1.fastq',
		catted_r2='/opt/output/{experiment}_concatenated_R2.fastq',
		desc_file='/opt/output/{experiment}_MID_sample_mappings.desc'
	params:
		output_folder='/opt/output',
		clean_percent_identity=config['clean_percent_identity'],
		clean_cpus=config['clean_cpus'],
		experiment='{experiment}'
	output:
		cleaned='/opt/output/{experiment}_DBLa_cleaned.fasta'
	shell:
		'''
		python /opt/cleanDBLalpha.py -o {params.output_folder} -r {input.catted_r1} -R {input.catted_r2} -d {input.desc_file} --perID {params.clean_percent_identity} --cpu {params.clean_cpus} --verbose
		mv {params.output_folder}/{params.experiment}_concatenated_R_DBLa_cleaned.fasta {params.output_folder}/{params.experiment}_DBLa_cleaned.fasta
		'''

rule cluster_dbla:
	'''
	'''
	input:
		cleaned='/opt/output/{experiment}_DBLa_cleaned.fasta' #not sure if this name formatting exactly matches what clean step will output
	params:
		output_folder='/opt/output',
		cluster_percent_identity=config['cluster_percent_identity'],
		cluster_cpus=config['cluster_cpus']
	output:
		cluster_file='/opt/output/{experiment}_DBLa_cleaned_renamed_centroids.fasta',
		binary_otu='/opt/output/{experiment}_DBLa_cleaned_renamed_otuTable_binary.txt'
	shell:
		'python /opt/clusterDBLa.py -o {params.output_folder} -r {input.cleaned} --perID {params.cluster_percent_identity} --cpu {params.cluster_cpus} --verbose'

rule classify_dbla_first:
	'''
	This program is messy - it requires a trailing / in the output folder name
	in order to send files to the desired output folder location.
	'''
	input:
		cluster_file='/opt/output/{experiment}_DBLa_cleaned_renamed_centroids.fasta'
	params:
		classify_threshold=config['classify_threshold'],
		output_folder='/opt/output'
	output:
		classify_first='/opt/output/{experiment}_DBLa_cleaned_renamed_centroids_nhmmOut.txt'
	shell:
		'python /opt/programs/classifyDBLalpha/reads_to_domains/allocate_reads.py -r {input.cluster_file} -E {params.classify_threshold} -o {params.output_folder}/ --noUproc --splitIsolates'

rule classify_dbla_second:
	'''
	'''
	input:
		classify_first='/opt/output/{experiment}_DBLa_cleaned_renamed_centroids_nhmmOut.txt'
	output:
		classify_second='/opt/output/{experiment}_reads_to_domains.csv'
	shell:
		'python /opt/programs/classifyDBLalpha/reads_to_domains/reads_to_domains.py --hmm {input.classify_first} --out {output.classify_second}'

rule combine_dbla:
	'''
	binary_file is the classify_second output of rule classify_dbla_second.
	ups_file is the binary_otu output of rule cluster_dbla
	'''
	input:
		binary_file='/opt/output/{experiment}_DBLa_cleaned_renamed_otuTable_binary.txt',
		ups_file='/opt/output/{experiment}_reads_to_domains.csv'
	output:
		combined_dbla_file='/opt/output/{experiment}_DBLa_binary_ups.csv'
	script:
		'/opt/working_dir/scripts/combine_step.R'