import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--Input', dest='input', required=True, type=argparse.FileType('r'))
parser.add_argument('--Remove', dest='remove', required=True, type=argparse.FileType('r'))
parser.add_argument('--Output', dest='output', required=True, type=argparse.FileType('w'))
parser.add_argument('--Kick', dest='kick', action='store_true')
args = parser.parse_args()

#soak up files
input = pd.read_csv(args.input,header=None,index_col=None,sep='\t').values
remove = pd.read_csv(args.remove,header=None,index_col=None,sep='\t').values

for i in range(input.shape[0]):
	rmsub = remove[remove[:,0]==input[i,0],:]
	if input[i,-1]=='+':
		#"left to right" gene
		#start of remove needs to be before the end of the promoter region
		#and end of remove needs to be after the start of the promoter region
		rmsub = rmsub[rmsub[:,1]<input[i,2],:]
		if rmsub.size:
			input[i,1] = min([input[i,2],max([input[i,1],np.max(rmsub[:,2])])])
	if input[i,-1]=='-':
		#"right to left" gene
		rmsub = rmsub[rmsub[:,2]>input[i,1],:]
		if rmsub.size:
			input[i,2] = max([input[i,1],min([input[i,2],np.min(rmsub[:,1])])])
	#so, now that we've evaluated all the possible overlaps, do we delete this gene?
	if input[i,1]==input[i,2] and args.kick:
		continue
	#if not, write it out
	args.output.write('\t'.join([str(x) for x in input[i,:]])+'\n')