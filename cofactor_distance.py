# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:50:02 2018

@author: hraanan
"""

''' 

This script calculates the distance between cofactors on a PDB file.
input: directory containing PDB files
output: file containing a list of cofactors and the distance between them

1. extract list of cofactors from file 
2. open PDBs and get the coordinats of each cofactor that in the list.
3. if there are more than 1 cofactor in the PDB make list of the cofactor pairs that in distance of less than min_dis 

important variables:
proteins_list= list of PDBs
lig= list of cofactor deatils (strings) 
ligands = list of lig's 
dis_min=20 # the minimal distance between cofactor to include in chain
nebrs= list of pers of cofactor in distance less than min_dis (i,j)
nebrsdis=list of pers of cofactor in distance less than min_dis include the distance (i,j,dis)
'''

import sys

sys.path.append("F:\programs\span")

import center # functions get_atom_list and get_center
import math
import os
from Bio import PDB
pdbl=PDB.PDBList()



    




dis_max=14 # the max and minimal distance between cofactor to include in chain
dis_min=4


AA=['PHE','TRP','TYR','ALA','CYS','ASP','GLU','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL']
CF=[' DA',' DC',' DG',' DT','  A','  C','  G','  U','HOH','UNK','DOD'] # list of cofactor to ignore
Metals=['FE','MN','CU','CO','NI','W','MO','V'] 
chl=['BCB','BCL','CHL','CL0','CLA','PHO','PMR','BPB','BPH']
heme=['B12','CNC','COH','DHE','F43','HAS','HDD','HDE','HEA','HEC','HEM','SRM']
pyr_atom_list=['C3A','C3B','C3C','C3D']
NAD=['NAD','ADJ','ENP','NAP','NDP','NJP','NZQ','XNP']

Distance=[]
in_flie_name=sys.argv[1]
in_file=open(in_flie_name,'r')

Alledges_file=open('edges_'+in_flie_name,'w')#"chain_from_"+sys.argv[1]+'_to_'+sys.argv[2]+'.txt',"w") # Output file
#Alledges_file=open('edges_test','w')#"chain_from_"+sys.argv[1]+'_to_'+sys.argv[2]+'.txt',"w") # Output file
Alledges_file.write('source\ttarget\tdistance\tprotein\n')

edges={}
#%%
# make dict of  microen andlist of  pdbs
microen_dict={}
pdbs_list=[]
cofactors_list=[]
with open(in_flie_name,'r') as microen_file:
    for line in microen_file:
        if line[:-2].lower() not in pdbs_list:
            pdbs_list.append(line[:-1].lower())
        
with open('groups_align_cof_and_nonredundant_2.12.19_ratio_ca_25_rmsd_4_ratio_0.1_md_2_with_sam.txt','r') as microen_file:
    microen_file.readline()
    for line in microen_file:
        #print(line)
        line=line.split('\t')
        microen_dict[line[0]]=line[1][:-1]
#        if line[0].split('.')[0] not in pdbs_list:
#            pdbs_list.append(line[0].split('.')[0])
        if line[0].split('.')[1].split('_')[0] not in cofactors_list:
            cofactors_list.append(line[0].split('.')[1].split('_')[0])
    microen_file.close()

#%%

pdbdir='F:\pdb_download_8_2018\pdb/'
for prot_index,protein in enumerate(pdbs_list):
    print('prot_index:'+str(prot_index))
    parser = PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(protein,pdbdir+protein[1:3]+'/pdb'+protein+'.ent')
#%%    
#geting the coordinates of all cofactors in the pdb
    for model in structure:
        if model.id==0:
         Ligands=[]
         lig=[]
         for chain in model:
              for residue in chain:
                 res_id=protein+'.'+residue.resname+'_'+chain.id+'_'+str(residue.id[1])
                 res_id=res_id.replace(" ","")                    
                 if residue.resname not in AA and residue.resname not in CF:
                     print ('res_id:'+res_id)                    
                     if res_id in microen_dict.keys():
                        print(chain.id,residue.id)
                        coord = []
                        cof_atom_list=['all']                    
                        if residue.resname in chl or residue.resname in heme:
                            cof_atom_list=pyr_atom_list
                        elif center.is_in_list(residue.resname)==True:
                            cof_atom_list=center.get_atom_list(residue.resname)
                  
                        coord=center.get_atom_coord_list(residue,cof_atom_list)   
                        x=0
                        y=0
                        z=0
                        lig=protein,chain.id,residue.id[1],microen_dict.get(res_id),x,y,z,coord,residue.resname  # lig= list of cofactor deatils (strings) 
                        print(lig)
                        Ligands.append(lig)
#%%     
   
    x=[]
    y=[]
    z=[]
    nebrs=[]
    nebrsdis=[]
    allperdis=[]
    allper=[]  
#make list of x's y's z's of all cofactors in the PDB
    for cluster in Ligands:
        x.append(cluster[4])
        y.append(cluster[5])
        z.append(cluster[6])
        
# if there are more than 1 cofactor in the PDB, make list of the cofactor pairs that in distance of less than min_dis 
    if len(x)>1 :
        
        for i in range(0,len(x)):
           #print('********************************************************',Ligands[i][8]) 
           atoms_i=Ligands[i][7]    
  #         print(i)    
           for j in range(0,len(x)):
               
                  #  print(i,j)
                   # Dis=math.sqrt((pow((x[i]-x[j]),2))+(pow((y[i]-y[j]),2))+(pow((z[i]-z[j]),2)))
                                        
                    if i==j:
                        continue
                    atoms_j=Ligands[j][7]    
                    Dis=100
                    for at_i in atoms_i:
                        for at_j in atoms_j:
                            if at_i==at_j:
                                continue
                            #print (at_i[0],at_j[0])
                            at_dis=math.sqrt((pow((at_i[0]-at_j[0]),2))+(pow((at_i[1]-at_j[1]),2))+(pow((at_i[2]-at_j[2]),2)))    
                            if at_dis<Dis:
                                Dis=at_dis
                                                               
                                #print('Dis:')
                                #print(at_dis)
                    n=[i,j]
                    allper.append(n)
                    n=[i,j,Dis]
                    allperdis.append(n)
                                        
                    if Ligands[i][8]==Ligands[j][8]:
                        dis_min=4
                    else:
                        dis_min=0
                    print('dis_min',dis_min)
                    if Dis>dis_min and Dis<dis_max:
#                        print(i,j)                         
#                        print('Dis:')
#                        print(Dis)
                        n=[i,j]
                        nebrs.append(n)
                        n=[i,j,Dis]
                        nebrsdis.append(n)
                        

    for x in range (0,len(Ligands)):
        Ligands[x]=Ligands[x]+(x,)
      #  print(Ligands[x])
    onesidenebrs=[]
    cng=0
    for j in range(len(Ligands)):
                i=Ligands[j]
               # print('--',i[0],i[7],i[3],i[1],i[2])
                

    for i in nebrs:
        if i not in onesidenebrs and i[::-1] not in onesidenebrs and int(i[0])!=int(i[1]):
            onesidenebrs.append(i)
    onesidenebrsdis=[]
    for i in onesidenebrs:
        for j in nebrsdis:
            if i[0]==j[0] and i[1]==j[1]:
                onesidenebrsdis.append(j)
    
    for i in nebrs:
        #print(i)
        if int(i[0])==int(i[1]):
            nebrs.remove(i)
    print('Done neighbors procsses...............')
    print('nebrs:',nebrs)
    

    chain=[]
    groups=[]
    lignebrs=[]
##    for i in Ligands:
##      print(i)
    for i in range(0,len(Ligands)):
        lignebrs.append([]) 
        for j in onesidenebrs:
            if j[0]==i and j[1] not in lignebrs[i]:
                lignebrs[i].append(j[1])
            if j[1]==i and j[0] not in lignebrs[i]:
                lignebrs[i].append(j[0])
 #   print ('onesidenebrs:',onesidenebrsdis)
    for i in onesidenebrsdis:
        print(i)
        print(Ligands[i[0]])
        fst=Ligands[i[0]][3]#+';'+Ligands[i[0]][8]
        scd=Ligands[i[1]][3]#+';'+Ligands[i[1]][8]
#        x=[Ligands[i[0]][4],Ligands[i[1]][4]]
#        y=[Ligands[i[0]][5],Ligands[i[1]][5]]
#        z=[Ligands[i[0]][6],Ligands[i[1]][6]]        
#        Dis=math.sqrt((pow((x[0]-x[1]),2))+(pow((y[0]-y[1]),2))+(pow((z[0]-z[1]),2)))        
        dis=i[2]        
        key=[fst,scd]
        if Ligands[i[0]][3] != Ligands[i[1]][3]:
            key=sorted(key)
        #key[0]=key[0].split(';')
        #key[1]=key[1].split(';')
        #Alledges_file.write(key[0][0]+'\t'+key[1][0]+'\t'+key[0][1]+'\t'+key[1][1]+'\t'+str(dis)+'\t'+protein+'\n')          
        Alledges_file.write(str(int(key[0]))+'\t'+str(int(key[1]))+'\t'+str(dis)+'\t'+protein+'\n')          
        edge=key[0][0]+'-'+key[1][0]
        key=edge
        if key in edges:
            edges[key] += 1
        else:
            edges[key] = 1

Alledges_file.close()


print("end")
