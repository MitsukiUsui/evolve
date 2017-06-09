import sys
import pandas as pd
import random
from Bio import SeqIO

def create_query(inFilepath, outFilepath, addDecoy=False):
	def is_typical(seq_record):
		if len(seq_record)>=6 and len(seq_record)%3==0:
			return True
		else:
			return False
	
	print("\tIN :{}".format(inFilepath))
	print("\tOUT:{}".format(outFilepath))
	
	with open(outFilepath,'w') as f:
		totalCDS=0
		processedCDS=0
		for seq_record in SeqIO.parse(inFilepath, "fasta"):
			totalCDS+=1
			if is_typical(seq_record):
				processedCDS+=1
				
				target=[(2,(1,len(seq_record)-2),False),#Frame Number, (start, end), is reverse complement 
						(3,(2,len(seq_record)-1),False),
						(4,(0,len(seq_record)),  True),
						(5,(1,len(seq_record)-2),True),
						(6,(2,len(seq_record)-1),True)]

				for frameNum, (start,end), revComp in target:
					seq=seq_record.seq[start:end]
					if revComp:
						subSeq=seq.reverse_complement()
					
					products=[]
					headers=[]					  
					product=str(seq.translate(table=11)).replace('*','X')		  
					products.append(product)
					headers.append(">F{0}|{1}".format(frameNum,seq_record.id.strip()))
					
					if addDecoy:
						dr=product[::-1]#decoy reverse
						products.append(dr)
						headers.append(">DRF{0}|{1}".format(frameNum,seq_record.id.strip()))
						ds=''.join(random.sample(product,len(product)))#decoy shuffle
						products.append(ds)
						headers.append(">DSF{0}|{1}".format(frameNum,seq_record.id.strip()))
						
					for product,header in zip(products, headers):					 
						f.write(header+'\n')
						f.write(product+'\n')
						
		print("\tDONE: extraction from {0}/{1} CDSs".format(processedCDS, totalCDS))
		
def main(targetFilepath):
	print("START processing {}".format(targetFilepath))
	target_df=pd.read_csv(targetFilepath)
	for basename in target_df["ftp_basename"]:
		inFilepath="/data/mitsuki/data/refseq/cds_from_genomic/{}_cds_from_genomic.fna".format(basename)
		outFilepath="/data/mitsuki/out/altorf/evolve/query/{}.query".format(basename)
		create_query(inFilepath, outFilepath, addDecoy=True)
		
if __name__=="__main__":
	main(sys.argv[1])
