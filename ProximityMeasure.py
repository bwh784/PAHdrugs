#!/usr/bin/env python

import getopt, sys, urllib, time
import random as rd
import networkx as nx
from scipy.stats import norm
import numpy as np

def ds(G,node,geneset):
	D=[]
	for i in range(len(geneset)):
		if(nx.has_path(G,node,target=geneset[i])):
			Lt=nx.shortest_path_length(G,node,target=geneset[i])
			D.append(Lt)
	if(len(D)>0):
		Ls=min(D)
		return Ls
	else:
		return -1
	

def distance(G,geneset1,geneset2):
	D=0
	idx=0
	for i in range(len(geneset1)):
		dt=ds(G,geneset1[i],geneset2)
		if(dt>=0):
			idx=idx+1
			D=D+dt
	if(idx>0):
		return float(D)/idx
	else:
		return -1



def main():

	file0=open("PH_genes_in_PPI201806.txt","r")
	file1=open("Human_Interactome_proteins.txt","r")	
	file2=open("Human_Interactome.txt","r")
	file3=open("PH_all_drugs_moduledrugs_targets.txt","r")
	file4=open("PH_all_drugs_moduledrugs_target_PH_proximity.txt","w")
	
	
 	
	subprot=[]
	for line in file0:
		tmp=line.split()
		subprot.append(tmp[0])
	file0.close()

	
		
	prot=[]
	for line in file1:
		tmp=line.split()
		prot.append(tmp[0])
	file1.close()	

 	G=nx.Graph()
	for line in file2:
		tmp=line.split()
		G.add_node(tmp[1])
		G.add_node(tmp[3])
		G.add_edge(tmp[1],tmp[3])	
	file2.close()
	
	#p=nx.shortest_path_length(G)
	
	idx=0
	for line in file3:
		tmp=line.split("\t")
		DT=[]
		idx=idx+1
		print idx,tmp[0]
		file4.write("%s\t"%tmp[0])
		for i in range(1,len(tmp)):
			if(G.has_node(tmp[i].strip())):
				DT.append(tmp[i].strip())
		#print DG
		d=distance(G,DT,subprot)
		
		file4.write("%s\t"%round(d,2))
					
# 		
# 			
		RR=1000
		Rdd=[]
		for k in range(RR):
			i=0
			#print k
			randDT=[]
			dsize=len(DT)
			while (i<dsize):
				r=rd.randint(0,len(prot)-1)
				randDT.append(prot[r])
				i=i+1	
			dd=distance(G,randDT,subprot)
			if(dd>0):
				Rdd.append(dd)
			
		u=np.mean(Rdd)
		o=np.std(Rdd)
		if(u>d):
			p=norm(u,o).cdf(d)
		else:
			p=1-norm(u,o).cdf(d)		
		file4.write("%s\t%s\t%s\t%s\n"%(u,o,p,dsize))		
		

	
if __name__ == "__main__":
    main()