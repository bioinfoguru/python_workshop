import pandas as pd
import sys
import os
from multiprocessing.pool import ThreadPool
import subprocess
import pysam
import pybedtools

dir_path = "/mnt/c/Users/bioinfo guru/Downloads/fusion_finder/Input/"

# read the bedpe file into a dataframe
df_PE = pd.read_csv(dir_path + "CM.bedpe", sep="\t", header=None)
df_PE.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"]
print(df_PE[["chr1","chr2"]].head())
print(df_PE.shape)



print("Reading GTF file")
# read GTF file into a dataframe
df_gtf = pd.read_csv(dir_path + "Genome_transcripts.gtf", sep="\t", header=None)
df_gtf.columns = ["chr", "source", "type", "start", "end", "score", "strand", "frame", "attribute"]

df_gtf_attrib = pd.DataFrame(df_gtf["attribute"].str.split("; ").tolist())
df_gtf_attrib = df_gtf_attrib.replace('"', '', regex=True)

for x in range(0,len(df_gtf_attrib.columns)):
    df_gtf_attrib[x] = df_gtf_attrib[x].str.split(" ").str[1]


df_gtf_attrib.columns = ["gene_id", "transcript_id", "db_xref", "gbkey", "gene_name", "transcript_name", "transcript_source", "transcript_biotype", "c9"]


print(df_gtf_attrib.head())

df_gtf = pd.concat([df_gtf, df_gtf_attrib], axis=1)
df_gtf = df_gtf.drop(["attribute"], axis=1)

print(df_gtf.head())


df_PE_left = df_PE.copy()
df_PE_right = df_PE.copy()


print("Processing PE file...")
df_PE_left[["chr1_1", "chr1_2"]] = df_PE_left["chr1"].str.split(">", expand=True)
df_PE_left[["chr1_gene", "chr1_3"]] = df_PE_left["chr1_1"].str.split("::", expand=True)
df_PE_left[["chr1_read", "chr1_range_strand"]] = df_PE_left["chr1_3"].str.split(":", expand=True)
df_PE_left[["chr1_range", "chr1_strand"]] = df_PE_left["chr1_range_strand"].str.split("(", expand=True)
df_PE_left[["chr1_start", "chr1_end"]] = df_PE_left["chr1_range"].str.split("-", expand=True)
df_PE_left["chr1_strand"] = df_PE_left["chr1_strand"].str.replace(")", "")

df_PE_left = df_PE_left.drop(["chr1_2", "chr1_3", "chr1_range_strand"], axis=1)


df_PE_right[["chr2_1", "chr2_2"]] = df_PE_right["chr2"].str.split(">", expand=True)
df_PE_right[["chr2_gene", "chr2_3"]] = df_PE_right["chr2_2"].str.split("::", expand=True)
df_PE_right[["chr2_read", "chr2_range_strand"]] = df_PE_right["chr2_3"].str.split(":", expand=True)
df_PE_right[["chr2_range", "chr2_strand"]] = df_PE_right["chr2_range_strand"].str.split("(", expand=True)
df_PE_right[["chr2_start", "chr2_end"]] = df_PE_right["chr2_range"].str.split("-", expand=True)
df_PE_right["chr2_strand"] = df_PE_right["chr2_strand"].str.replace(")", "")

df_PE_right = df_PE_right.drop(["chr2_2", "chr2_3", "chr2_range_strand"], axis=1)

print(df_PE_right[["chr2_range"]].head())

df_PE_left_right = pd.concat([df_PE_left, df_PE_right], axis=1)
df_PE_left_right = df_PE_left_right.drop(["chr1", "chr2"], axis=1)

df_PE_left_right.to_csv(dir_path + "PE_left_right.txt", sep="\t", index=False)



print("df_PE right after processing")
print(df_PE_right[["chr1_read","chr2_read"]].head())

df_gtf_left = df_gtf[df_gtf["transcript_id"].isin(df_PE_left["chr1_read"])]
df_gtf_final = pd.merge(df_gtf, df_PE_left_right, left_on="transcript_id", right_on="chr1_read", how="left")
df_gtf_final = df_gtf_final[df_gtf_final["chr1_read"].notnull()]
df_gtf_final = df_gtf_final.drop(["chr1_1", "chr2_1"], axis=1)

df_gtf_final.to_csv(dir_path + "gtf_final.txt", sep="\t", index=False)
