import pandas as pd
import numpy as np

#soak up current promoters.bed file
prom = pd.read_csv('genelines.bed',header=None,index_col=None,sep='\t').values
annot = pd.read_csv('annot.gff3',header=None,index_col=None,sep='\t',comment='#').values

#write out 5' UTR BED
utr5 = open('utr5.bed','w')
utr3 = open('utr3.bed','w')
#and also write CDS span gene lines
cds = open('genelines2.bed','w')

#the order is the same in the universe.txt as it is in annot.gff3
#after all, annot.gff3 was used to craft universe.txt
#however, promoters.bed might have lost some along the way, so err on the side of caution
geneinds = []
for i in np.arange(annot.shape[0]):
	if 'gene' in annot[i,:]:
		geneinds.append(i)
#and need an "end of file" geneinds as well
geneinds.append(annot.shape[0]+1)

#loop over the genes
for i in np.arange(prom.shape[0]):
	#dig out the subannot and filter it to cds
	subannot = annot[geneinds[i]:geneinds[i+1],:]
	noncds = []
	for j in np.arange(subannot.shape[0]):
		if 'CDS' not in subannot[j,:]:
			noncds.append(j)
	temp = np.delete(subannot,noncds,axis=0)
	#if there's anything found, continue
	if temp.size:
		if prom[i,5]=='+':
			#"left to right" transcript
			newpos = np.min(temp[:,3])
			holdlist = [prom[i,0],prom[i,1],newpos,prom[i,3],prom[i,4],prom[i,5]]
			utr5.write('\t'.join([str(x) for x in holdlist])+'\n')
			newpos = np.max(temp[:,4])
			holdlist = [prom[i,0],newpos,prom[i,2],prom[i,3],prom[i,4],prom[i,5]]
			utr3.write('\t'.join([str(x) for x in holdlist])+'\n')
		else:
			#"right to left" transcript
			newpos = np.max(temp[:,4])
			holdlist = [prom[i,0],newpos,prom[i,2],prom[i,3],prom[i,4],prom[i,5]]
			utr5.write('\t'.join([str(x) for x in holdlist])+'\n')
			newpos = np.min(temp[:,3])
			holdlist = [prom[i,0],prom[i,1],newpos,prom[i,3],prom[i,4],prom[i,5]]
			utr3.write('\t'.join([str(x) for x in holdlist])+'\n')
		#regardless of right to leftness or left to rightness, picking up the CDS span
		holdlist = [prom[i,0],np.min(temp[:,3]),np.max(temp[:,4]),prom[i,3],prom[i,4],prom[i,5]]
		cds.write('\t'.join([str(x) for x in holdlist])+'\n')

utr5.close()
utr3.close()
cds.close()