import pandas as pd

resultFilepath="/data/mitsuki/out/altorf/evolve/result/GCF_000022125.1_ASM2212v1.tab2"
names_lst=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
result_df=pd.read_csv(resultFilepath,delimiter='\t',names=names_lst)
thres=1e-5
filtered_df=result_df[result_df["evalue"]<thres]
for qseqid in filtered_df["qseqid"]:
    print(qseqid)