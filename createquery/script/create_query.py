import sys
import random
from collections import Counter
import pandas as pd
import numpy as np
from scipy import stats
import Bio
from Bio import SeqIO


class AaBias():
	""" remember amino acid bias
	
	Attributes:
		self.aa_lst:	  list of 1 char representation of 20 Amino Acids
		self.aaCount_dct: dictionary of amino acid (+ X) count
		self.aaFreq_arr:  numpy array of culmitive frequency of each Amino Acids
	"""
	
	def __init__(self, filepath=None):
		"""initializes attributes and updates if file is given"""
		
		self.aa_lst=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
		self.all_lst=self.aa_lst+["X"]
		self.aaCount_dct={}
		for aa in self.all_lst:
			self.aaCount_dct[aa]=0
		self.aaFreq_arr=np.zeros(len(self.all_lst))
		
		if filepath is not None:
			for seq_record in SeqIO.parse(filepath, "fasta"):
				self.update_count(seq_record.seq)
			self.calc_freq()
				
	
	def update_count(self,aaSeq):
		"""get amino acid sequence and update count dictionary"""
		
		counter=Counter(aaSeq)
		for k,v in counter.items():
			if k in self.aa_lst:
				self.aaCount_dct[k]+=v
			elif k in ["X"]:
				self.aaCount_dct[k]+=v
			else:
				#sys.stderr("CAUGHT unexpected amino acid: {}".format(k))
				#sys.exit(1)
				if k!="U":	
					print("\tUNDIFINED amino acid: {0},{1}".format(k,v))
				self.aaCount_dct["X"]+=v	
			

	def calc_freq(self):
		"""calculate frequency based on count dictionary"""
		
		totalCount=0
		for _,v in self.aaCount_dct.items():
			totalCount+=v
		
		culmFreq=0.
		for i,aa in enumerate(self.all_lst):
			culmFreq+=self.aaCount_dct[aa]/totalCount
			self.aaFreq_arr[i]=culmFreq
		self.aaFreq_arr[-1]=1 #correct rounding error
			 
	def generate_random(self,seqLen):
		"""generate random sequence of length seqLen, according to frequency"""
		
		retVal=""
		rand_arr=np.random.rand(seqLen)
		for rand in rand_arr:
			for i,culmFreq in enumerate(self.aaFreq_arr):
				if culmFreq>=rand:
					break
			retVal+=self.all_lst[i]
			
		return retVal
		
		
	def kl(self, other):
		"""calculate kullback leibler distance between two AaBias instance"""
		
		assert(self.aa_lst==other.aa_lst)
		
		#convert count dictionary to numpy array
		a1=np.zeros(len(self.aa_lst))
		a2=np.zeros(len(other.aa_lst))
		for i,aa in enumerate(self.aa_lst):
			a1[i]=self.aaCount_dct[aa]+1
			a2[i]=other.aaCount_dct[aa]+1#to avoid 0 division
		return stats.entropy(a1,a2)
				
	def _show(self):
		"""[DEBUG] show """
		for i,aa in enumerate(self.all_lst):
			print(aa, self.aaCount_dct[aa], self.aaFreq_arr[i])


def is_typical(seq_record):
	"""define prerequisite for CDS to be base of queries"""
	if len(seq_record)>=6 and len(seq_record)%3==0:
		return True
	else:
		return False
	   
		
def get_ab_lst(filepath):
	""" get list of AaBias instance for all 6 reading frame
	
	Arguments:
		filepath: cds_from_genomic file
	"""
	
	ab_lst=[AaBias() for _ in range(6)]
	for seq_record in SeqIO.parse(filepath, "fasta"):
		if is_typical(seq_record):
			#Frame Number, (start, end), is reverse complement 
			target=[(1,(0,len(seq_record)-3),False),# -3 to drop last stop codon
					(2,(1,len(seq_record)-2),False),
					(3,(2,len(seq_record)-1),False),
					(4,(0,len(seq_record)),  True),
					(5,(1,len(seq_record)-2),True),
					(6,(2,len(seq_record)-1),True)]

			for frameNum, (start,end), revComp in target:
				seq=seq_record.seq[start:end]
				if revComp:
					seq=seq.reverse_complement()
				aaSeq=str(seq.translate(table=11)).replace('*','X')
				ab_lst[frameNum-1].update_count(aaSeq)

	for i in range(6):
		ab_lst[i].calc_freq()
	
	return ab_lst

	
			
def create_query(inFilepath, outFilepath, addDecoy=False, ab_lst=None):
	""" create query
	
	create query based on the translational product of 5 alternative reading frame 
	
	Arguments:
		inFilepath:  cds_from_genomic file
		outFilepath: fasta file to be written query
		addDecoy:	 if True, add decoy
		ab_lst:		 list of AaBias for each reading frame
	"""
	

	print("\t\tIN  :{}".format(inFilepath))

	query_lst=[]#each element is (header, query)
	decoy_lst=[]#devides list between query & decoy in order to split location of query & decoy in the output file
	
	totalCDS=0
	processedCDS=0
	
	#compose query and decoy(DR/DS/DTS) if required
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
					seq=seq.reverse_complement()
				
				queryHeader=">F{0}|{1}".format(frameNum,seq_record.description.strip())
				querySeq=str(seq.translate(table=11)).replace('*','X')
				query_lst.append((queryHeader,querySeq))
				
				if addDecoy:
					#decoy reverse
					#drHeader=">DRF{0}|{1}".format(frameNum,seq_record.id.strip())
					#drSeq=querySeq[::-1]
					#decoy_lst.append((drHeader,drSeq))
				   
					#decoy shuffle
					dsHeader=">DSF{0}|{1}".format(frameNum,seq_record.id.strip())
					dsSeq=''.join(random.sample(querySeq,len(querySeq)))
					decoy_lst.append((dsHeader,dsSeq))
					
					#decoy total shuffle
					dsHeader=">DTSF{0}|{1}".format(frameNum,seq_record.id.strip())
					dsSeq=ab_lst[frameNum-1].generate_random(len(querySeq))
					decoy_lst.append((dsHeader,dsSeq))
					
			 
	#output query_lst and decoy_lst into outFilepath
	with open(outFilepath,'w') as f:
		for header,seq in query_lst:
			f.write(header+'\n')
			f.write(seq+'\n')
		if addDecoy:
			for header,seq in decoy_lst:
				f.write(header+'\n')
				f.write(seq+'\n')
		
	print("\t\tcreated from {0}/{1} CDSs".format(processedCDS, totalCDS))
	print("\t\tOUT :{}".format(outFilepath))

		
def main(genusNum):
	print("START processing genus {}".format(genusNum))
	
	catalogFilepath="/home/mitsuki/altorf/evolve/createdatabase/out/genus_{}.txt".format(genusNum)
	catalog_df=pd.read_csv(catalogFilepath)
	dbFilepath="/data/mitsuki/data/blast/db/genus_{}.fasta".format(genusNum)
	abDb=AaBias(dbFilepath)

	print("DONE calculating Amino Acid Bias of {}".format(dbFilepath))
	
	for i,basename in enumerate(catalog_df["ftp_basename"]):
		print("PROCESSING {}/{} ({})...".format(i+1,catalog_df.shape[0],basename))
		
		inFilepath="/data/mitsuki/data/refseq/cds_from_genomic/{}_cds_from_genomic.fna".format(basename)
		outFilepath="/data/mitsuki/out/altorf/evolve/query/{}.query".format(basename)
		
		ab_lst=get_ab_lst(inFilepath)
		print("\tDONE calculating Amino Acid Bias for 6 frame")
		
		create_query(inFilepath, outFilepath, addDecoy=True, ab_lst=ab_lst)
		print("\tDONE creating query") 
		
if __name__=="__main__":
	main(sys.argv[1])
