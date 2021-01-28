# ! //usr/bin/env python
import sys, os, string, math, networkx
import numpy
import community

    
def main():
	#if (len(sys.argv)<>3):
		#print 'Please provide two file'
		#sys.exit()
			
			
	file1 = open("PHAlldrugtarget-disease_BipartiteNetwork.txt","r")
	file2 = open("PHAlldrugtarget-disease_BipartiteNetwork_Modularity.txt","w")
	
 	G=networkx.Graph()
 	
	
	for line in file1:
		tmp=line.split()
		G.add_edge(tmp[0].strip(),tmp[1].strip())		
	file1.close()
	

 	C_best=community.best_partition(G)

 	for key in C_best:
 		file2.write("%s\t%s\n"%(key,C_best[key]))	
 	file2.close()
 	
 	M=community.modularity(C_best, G)
 	print M
 	


  	

if __name__=='__main__':
	main()
