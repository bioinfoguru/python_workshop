from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaMemCommandline
from Bio import SeqIO
import pandas as pd
import sys
import os
from multiprocessing.pool import ThreadPool
import subprocess
import pysam
import pybedtools


N_CORES = 1
dir_path_win = "C:\\fusion_finder\\"
dir_path = "/mnt/c/Users/fusion_finder/Input/"

bwa_exec = "/mnt/c/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
reference_transcriptome = dir_path + "dummy_transcriptome.fna"#"transcriptome.fa"

sam_out1 = "bwa_aligned.sam"
read_file_1 = dir_path + "dummy_1_val_1.fq"#"SRR1575220_1.trim.fastq"
read_file_2 = dir_path + "dummy_2_val_2.fq"#"SRR1575220_2.trim.fastq"

#### fuction to get sequences of the reads from the fastq file. ####
def get_seqs(fastq_file, list_read_IDs):
    print("Reading fastq file", fastq_file)
    with open(fastq_file, "r") as f:
        records = list(SeqIO.parse(f, "fastq"))
        read_IDs = [record.id for record in records]
        sequences = [record.seq for record in records]
        quality_id = [record.name for record in records]
        quality_scores = [ "".join([chr(score + 33) for score in record.letter_annotations["phred_quality"]]) for record in records]
        # print(quality_scores)
        df = pd.DataFrame({"ReadID": read_IDs, "sequence": sequences, "QualityId":quality_id, "Quality": quality_scores})

    print("Filtering reads for", fastq_file)
    df = df[df["ReadID"].isin(list_read_IDs)]
    return df

# function for indexing and alignment
def index_align(fastq_file, prefix1, read1, read2, out_sam):

    subprocess.run([bwa_exec, "index", fastq_file, \
                "-p", prefix1], cwd=dir_path)

    subprocess.run([bwa_exec, "mem", prefix1,\
                read1, read2, "-o", out_sam], cwd=dir_path)
                 

#index_align(reference_transcriptome, "I2", read_file_1, read_file_2, sam_out1)

print("Alignment Done.")

# convert SAM file to BAM file and sort the BAM file by read name using samtools
temp1 = subprocess.run(["samtools", "view", "-b", "-o", dir_path + sam_out1[:-4] + ".bam", dir_path + sam_out1], capture_output=True, cwd=dir_path)

temp1 = subprocess.run(["samtools", "view", "-f", "4", dir_path + sam_out1[:-4] + ".bam"], capture_output=True, cwd=dir_path)
print(type(temp1))
temp1_list = temp1.stdout.decode("utf-8").split("\n")
UR_list = [x.split("\t")[0] for x in temp1_list]
print(len(UR_list))


# sort the unmapped reads by read name
pysam.sort("-n", "-o", dir_path + sam_out1[:-4] + "_unmapped_sorted.bam", dir_path + sam_out1[:-4] + "_unmapped.bam")

print("Converting BAM to BEDPE")
bam_file = pybedtools.BedTool(dir_path + sam_out1[:-4] + "_unmapped_sorted.bam")
bam_file.bam_to_bed(bedpe=True).saveas(dir_path + sam_out1[:-4] + "_sorted.bedpe")

df_bedpe = pd.read_csv(dir_path + sam_out1[:-4] + "_sorted.bedpe", sep="\t", header=None)
df_bedpe.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"]
print(df_bedpe.head())
print(df_bedpe.shape)


df_bedpe = df_bedpe[df_bedpe["chr1"] != "."]

df_bedpe[["chr1", "chr1_2"]] = df_bedpe["chr1"].str.split(".", expand=True)
df_bedpe[["chr2", "chr2_2"]] = df_bedpe["chr2"].str.split(".", expand=True)
print(df_bedpe.head())
df_bedpe_pt = df_bedpe[df_bedpe["chr1"] == "Pt"]
df_bedpe_mt = df_bedpe[df_bedpe["chr1"] == "Mt"]

df_bedpe = df_bedpe[df_bedpe["chr1"] != df_bedpe["chr2"]]

df_bedpe = pd.concat([df_bedpe, df_bedpe_pt, df_bedpe_mt], axis=0)

print(df_bedpe.head())

df_bedpe["chr1"] = df_bedpe["chr1"] + "." + df_bedpe["chr1_2"]
df_bedpe["chr2"] = df_bedpe["chr2"] + "." + df_bedpe["chr2_2"]
df_bedpe = df_bedpe.drop(["chr1_2", "chr2_2"], axis=1)

print(df_bedpe.head())
print(df_bedpe.shape)

DR_list = df_bedpe["name"].unique().tolist()
print(len(DR_list))

# create artifical transcriptome using the bedpe file

df_bedpe_right = df_bedpe.iloc[:, [3,4,5,6,7,9]]
df_bedpe_left = df_bedpe.iloc[:, [0,1,2,6,7,8]]

df_bedpe_right.to_csv(dir_path + "rightside.bed", sep="\t", header=False, index=False)
df_bedpe_left.to_csv(dir_path + "leftside.bed", sep="\t", header=False, index=False)


pybedtools.bedtool.BedTool(dir_path + "rightside.bed").sequence(fi=reference_transcriptome, s=True, name=True, fo=dir_path + "rightside.fa")
pybedtools.bedtool.BedTool(dir_path + "leftside.bed").sequence(fi=reference_transcriptome, s=True, name=True, fo=dir_path + "leftside.fa")

df_right = pd.read_csv(dir_path + "rightside.fa", sep="\t", header=None)
df_left = pd.read_csv(dir_path + "leftside.fa", sep="\t", header=None)

# concatenate the two dataframes
df_right_left = pd.concat([df_right, df_left], axis=1)
print(df_right_left.head())
print(df_right_left.shape)

UR_DR_list = list(set(UR_list + DR_list))

###### get sequences of UR_DR_list from both the fastq files. ######
print("Getting sequences of UR_DR_list from both the fastq files.")
print(read_file_1)
cmd = "split -l 4000000 "+ "'"+read_file_1+"'" +" --additional-suffix=_splitfile_"
subprocess.run(cmd, shell=True, cwd=dir_path)
# get the list of split files
split_files = [f for f in os.listdir(dir_path) if f.endswith("_splitfile_")]

pool = ThreadPool(processes=N_CORES)
jobs=[]

for x in split_files:
    async_result = pool.apply_async(get_seqs, args=(dir_path+x,UR_DR_list))
    jobs.append(async_result)

df_final = pd.DataFrame()
for t in jobs:
    out = t.get()
    df_final = pd.concat([df_final, out])

    # Write dataframe to a file.
df_final["ReadID"] = ">" + df_final["ReadID"].astype(str)
df_final.to_csv(dir_path+"reads_1.fasta", index=False, sep="\n", header=False)
os.system("rm *_splitfile_")

# read_file_2
cmd = "split -l 4000000 "+ read_file_2 +" --additional-suffix=_splitfile_"
subprocess.run(cmd, shell=True, cwd=dir_path)
# get the list of split files
split_files = [f for f in os.listdir(dir_path) if f.endswith("_splitfile_")]

pool = ThreadPool(processes=N_CORES)
jobs=[]
for x in split_files:
    async_result = pool.apply_async(get_seqs, args=(dir_path+x,UR_DR_list))
    jobs.append(async_result)

df_final = pd.DataFrame()
for t in jobs:
    out = t.get()
    df_final = pd.concat([df_final, out])

df_final["ReadID"] = ">" + df_final["ReadID"].astype(str)
df_final.to_csv(dir_path+"reads_2.fasta", index=False, sep="\n", header=False)

os.system("rm *_splitfile_")
##### Done getting sequences of UR_DR_list#####

print("Alignemt 2")
index_align(reference_transcriptome, "I2", dir_path + "reads_1.fasta", dir_path + "reads_2.fasta", "bwa_aligned_2.sam")

pysam.view("-b", "-o", dir_path + "bwa_aligned_2.bam", dir_path + "bwa_aligned_2.sam")
pysam.sort("-n", "-o", dir_path + "bwa_aligned_2_sorted.bam", dir_path + "bwa_aligned_2.bam")

pysam.view("-b", "-f", "2", "-F", "4", "-o", dir_path + "CM.bam", dir_path + "bwa_aligned_2_sorted.bam")

temp2 = subprocess.run(["samtools", "view", dir_path + "CM.bam"], capture_output=True, cwd=dir_path)
temp2_list = temp2.stdout.decode("utf-8").split("\n")
CM_list = [x.split("\t")[0] for x in temp2_list]
CM_list = list(set(CM_list))
print("Length of CM_list is ", len(CM_list))

pybedtools.bedtool.BedTool(dir_path + "bwa_aligned_2.bam").bam_to_bed(bedpe=True).saveas(dir_path + "PE.bedpe")

df_PE = pd.read_csv(dir_path + "PE.bedpe", sep="\t", header=None)
df_PE.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"]
print(df_PE.head())
print(df_PE.shape)

df_PE = df_PE[df_PE["chr1"] != "."]

df_PE[["chr1", "chr1_2"]] = df_PE["chr1"].str.split(".", expand=True)
df_PE[["chr2", "chr2_2"]] = df_PE["chr2"].str.split(".", expand=True)
print(df_PE.head())
df_PE = df_PE[df_PE["chr1"] != df_PE["chr2"]]
print(df_PE.head())
print(df_PE.shape)


df_PE["chr1"] = df_PE["chr1"] + "." + df_PE["chr1_2"]
df_PE["chr2"] = df_PE["chr2"] + "." + df_PE["chr2_2"]
df_PE = df_PE.drop(["chr1_2", "chr2_2"], axis=1)

print(df_PE.head())
print(df_PE.shape)

PE_list = df_PE["name"].unique().tolist()
print(len(PE_list))

df_PE.to_csv(dir_path + "PE.bedpe", sep="\t", header=True, index=False)
df_PE.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"]

print("Readng GTF file...")
df_gtf = pd.read_csv(dir_path + "dummy.gtf", sep="\t", header=None)
df_gtf.columns = ["chr", "source", "type", "start", "end", "score", "strand", "frame", "attribute"]
df_gtf = df_gtf[df_gtf["type"] == "transcript"]

df_gtf[["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]] = df_gtf["attribute"].str.split(";", expand=True)
df_gtf[["c1_name", "c1_value"]] = df_gtf["c1"].str.split(" ", expand=True)
df_gtf[["c2_name", "c2_value"]] = df_gtf["c2"].str.split(" ", expand=True)
df_gtf[["c3_name", "c3_value"]] = df_gtf["c3"].str.split(" ", expand=True)
df_gtf[["c4_name", "c4_value"]] = df_gtf["c4"].str.split(" ", expand=True)
df_gtf[["c5_name", "c5_value"]] = df_gtf["c5"].str.split(" ", expand=True)
df_gtf[["c6_name", "c6_value"]] = df_gtf["c6"].str.split(" ", expand=True)
df_gtf[["c7_name", "c7_value"]] = df_gtf["c7"].str.split(" ", expand=True)
df_gtf[["c8_name", "c8_value"]] = df_gtf["c8"].str.split(" ", expand=True)

df_gtf = df_gtf.drop(["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"], axis=1)

df_gtf["c1_value"] = df_gtf["c1_value"].str.replace('"', '')
df_gtf["c2_value"] = df_gtf["c2_value"].str.replace('"', '')
df_gtf["c3_value"] = df_gtf["c3_value"].str.replace('"', '')
df_gtf["c4_value"] = df_gtf["c4_value"].str.replace('"', '')
df_gtf["c5_value"] = df_gtf["c5_value"].str.replace('"', '')
df_gtf["c6_value"] = df_gtf["c6_value"].str.replace('"', '')
df_gtf["c7_value"] = df_gtf["c7_value"].str.replace('"', '')
df_gtf["c8_value"] = df_gtf["c8_value"].str.replace('"', '')

df_gtf = df_gtf.rename(columns={"c1_value": "gene_id", "c2_value": "transcript_id", "c3_value": "gene_name", "c4_value": "gene_source", "c5_value": "gene_biotype", "c6_value": "transcript_name", "c7_value": "transcript_source", "c8_value": "transcript_biotype"})
df_gtf = df_gtf.drop(["attribute"], axis=1)
df_gtf = df_gtf.drop(["c1_name", "c2_name", "c3_name", "c4_name", "c5_name", "c6_name", "c7_name", "c8_name"], axis=1)

print("GTF after processing")
print(df_gtf.head())
print(df_gtf.shape)


print("Processing PE file...")
df_PE[["chr1", "chr1_2"]] = df_PE["chr1"].str.split(">", expand=True)
df_PE[["chr1_gene", "chr1_3"]] = df_PE["chr1_2"].str.split("::", expand=True)
df_PE[["chr1_read", "chr1_range", "chr1_4"]] = df_PE["chr1_3"].str.split(":", expand=True)
df_PE = df_PE.drop(["chr1_2", "chr1_3"], axis=1)

df_PE[["chr2", "chr2_2"]] = df_PE["chr2"].str.split(">", expand=True)
df_PE[["chr2_gene", "chr2_3"]] = df_PE["chr2_2"].str.split("::", expand=True)
df_PE[["chr2_read", "chr2_range", "chr2_4"]] = df_PE["chr2_3"].str.split(":", expand=True)
df_PE = df_PE.drop(["chr2_2", "chr2_3"], axis=1)

print("df_PE after processing")
print(df_PE.head())
