'''
This is an 'un-demultiplexer' or maybe a 're-multiplexer'. It takes
demultiplexed data as input and makes it look like data that hasn't been
demultiplexed yet (for input into Karen Day's pipeline). It uses a desc file to
prepend barcodes to sample reads and concatenate them all into a single file.

Make sure your MIDs.tsv file contains at least as many MIDs as the number of
samples in your desc file. The MIDs.tsv file also needs to have the same MIDs as
are used by the cleanDBLa script of the Day pipeline.
'''

MID_file=snakemake.input.mid_file
fastq_directory = snakemake.input.seekdeep_fastq_directory
desc_mapping_file = snakemake.input.desc_file
edited_fastq_directory = snakemake.output.temporary_prepend_folder
catted_r1 = snakemake.output.catted_r1
catted_r2 = snakemake.output.catted_r2
r1_suffix=snakemake.params.R1_suffix
r2_suffix=snakemake.params.R2_suffix

MIDIC = dict([line.strip().split('\t') for line in open(MID_file)])

def list_files_in_directory(directory_path):
	import os
	"""Lists all files in the specified directory.

	Args:
		directory_path: The path to the directory.
	
	Returns:
		A list of strings, where each string is a file name in the directory.
		Returns an empty list if the directory does not exist or if an error occurs.
	"""
	try:
		files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
		return files
	except FileNotFoundError:
		print(f"Error: Directory not found: {directory_path}")
		return []
	except Exception as e:
		print(f"An error occurred: {e}")
		return []

def prepend_barcodes(input_fastq_file, output_fastq_file, barcode):
	import gzip
	original_fastq = gzip.open(input_fastq_file, mode='rt')
	edited_fastq = open(output_fastq_file, 'w')
	for line_number, line in enumerate(original_fastq):
		if line_number%4==1:
			edited_fastq.write(barcode)
		elif line_number%4==3:
			edited_fastq.write('C'*len(barcode))
		edited_fastq.write(line)

def match_barcodes(desc_mapping_file, fastq_directory, edited_fastq_directory):
	import subprocess
	file_list = list_files_in_directory(fastq_directory)
	subprocess.call(['mkdir', edited_fastq_directory])
	for line in open(desc_mapping_file):
		if '#' not in line:
			line = line.strip().split()
			sample_id = line[0]
			for file in file_list:
				file_id = file.replace(r1_suffix+'.fastq.gz', '')
				file_id = file_id.replace(r2_suffix+'.fastq.gz', '')
				if sample_id==file_id:
					input_fastq_file=f'{fastq_directory}/{file}'
					output_fastq_file=f"{edited_fastq_directory}/{file.replace('.fastq.gz', '_edited.fastq')}"
					if file.endswith(r1_suffix+'.fastq.gz'):
						barcode = MIDIC[line[1]]
					elif file.endswith(r2_suffix+'.fastq.gz'):
						barcode = MIDIC[line[2]]
					else:
						print(file, 'does not end with', r1_suffix, 'or', r2_suffix)
					prepend_barcodes(input_fastq_file, output_fastq_file, barcode)

def concatenate_files(edited_fastq_directory, catted_r1, catted_r2, r1_suffix, r2_suffix):
	edited_file_list = list_files_in_directory(edited_fastq_directory)
	edited_file_list.sort()
	concatenated_read1 = open(catted_r1,'w')
	concatenated_read2 = open(catted_r2, 'w')
	for file in edited_file_list:
		if r1_suffix in file:
			for line in open(f"{edited_fastq_directory}/{file}",'r'):
				concatenated_read1.write(line)
		elif r2_suffix in file:
			for line in open(f"{edited_fastq_directory}/{file}",'r'):
				concatenated_read2.write(line)

match_barcodes(desc_mapping_file, fastq_directory, edited_fastq_directory)
concatenate_files(edited_fastq_directory, catted_r1, catted_r2, r1_suffix, r2_suffix)