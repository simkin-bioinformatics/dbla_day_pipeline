'''
randomly assigns MID numbers to the samples in a folder of demultiplexed fastq
files. Useful for Karen Day's pipeline that expects a desc file and
un-demultiplexed input data.
'''
import os
MID_file=snakemake.input.mid_file
fastq_directory = snakemake.params.input_folder
shared_suffix_R1=snakemake.params.R1_suffix
shared_suffix_R2=snakemake.params.R2_suffix
desc_file = open(snakemake.output.desc_file, 'w')
desc_file.write("#ID-Number  AF-MID  BR-MID\n")
MID_dict = dict([line.strip().split('\t') for line in open('MIDs.tsv')])

def list_files_in_directory(directory_path):
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

file_list = list_files_in_directory(fastq_directory)
file_list.sort()
sample_name_list = []
for file in file_list:
    sample_name = file.replace(shared_suffix_R1+'.fastq.gz', '')
    sample_name = sample_name.replace(shared_suffix_R2+'.fastq.gz', '')
    if sample_name not in sample_name_list:
        sample_name_list.append(sample_name)

barcodes=len(MID_dict)
for sample_number, sample in enumerate(sample_name_list):
    first_barcode=(sample_number//barcodes)+1
    second_barcode=(sample_number%barcodes)+1
    desc_file.write(f"{sample}\t{first_barcode}\t{second_barcode}\n")
# print(file_list)