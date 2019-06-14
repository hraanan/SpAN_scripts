# This script extracts 15A diameter microenvironments from PDB
# input: directory containing PDB files and file contaning the atoms for each cofactor (manual_cofactor_atoms_list.txt)
# output: 15A microenvironments for each cofactor of each protein in a PDB format
import Bio
import os
import sys
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
import math
import numpy
from collections import Counter
import random 
from Bio.PDB import *
import gzip


    
def get_center(residue,al): # get the geometric centar of list of atom coordinatins 
    coord = []
    
    for atom in residue:
       # print(al)
        if atom.name in al or al==['all']:
        #   print(atom.coord)
           at=atom.coord
           x=at[0]
           y=at[1]
           z=at[2]
           atcord=[x,y,z]
           coord.append(atcord)
    x=0
    y=0
    z=0
    i=0
    for point in coord:
        i=i+1
        x=x+point[0]
        y=y+point[1]
        z=z+point[2]
    x=x/i
    y=y/i
    z=z/i
    center=numpy.array([x,y,z])    
    return center;






######initial variabels############# 
rootdir ='f:/pdb_download_8_2018/pdb/' # folder of the pdb original filesdownladed form the pdb
i=0
EC=""
cng=0     
fd=['1fdn'] # used to test script on knowen pdb file
AA=['PHE','TRP','TYR','ALA','CYS','ASP','GLU','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL']
CF=[' DA',' DC',' DG',' DT','  A','  C','  G','  U','HOH','UNK','UNX'] # list of cofactor names to ignore
Metals=['FE','MN','CU','CO','NI','W','MO','V'] 
chl=['BCB','CLA','CHL','BCL','CL0','PMR','PHO'] #chlorophyll cofactors list
heme=['HEA','HAS','2FH','522','89R','DDH','DHE','HES','HDD','HDE','HDM','HEB','HEC','HEM','HEO','HEV','HP5','MH0','N7H','NTE','OBV','SRM','VER']
pyr_atom_list=['C3A','C3B','C3C','C3D']
NAD=['NAD','ADJ','ENP','NAP','NDP','NJP','NZQ','XNP']

#####################
organic_cofactors_list=[]
cofactors_dict={}
organic_cofactors_pdb_file=open('F:\programs\span\data/manual_cofactor_atoms_list.txt','r')
for line in organic_cofactors_pdb_file:
    line=line.split('\t')
    organic_cofactors_list.append(line[1])
    if line[1] not in cofactors_dict:
        #print(len(line[4]))        
        if len(line[4])>1:
            cofactors_dict[line[1]]=[line[3],line[4][:-1]]
        else:
            cofactors_dict[line[1]]=[line[3]]

pdbl=PDB.PDBList()
Error_out=open("microfolds_out.txt","w")


for subdir, dirs, files in os.walk(rootdir):
    for file in files:

    
        try:
            line=file
            protein=line[3:7]
            #print ('pdb_code:'+protein)
            protein=protein.lower()
            Error_out.write('pdb_code:'+protein+'\n')
        
            parser = PDB.PDBParser(PERMISSIVE=1,get_header=1,QUIET=1)
            curdir=os.getcwd()
            filename=rootdir+protein[1:3]+'/'+file
            #print(filename)
            
            final_file=rootdir+protein[1:3]+'/pdb'+protein+'.ent'
            #print ('unziping')
            # unzipping gz file 
            gz = gzip.open(filename, 'rb') 
            with open(final_file, 'wb') as out: 
                out.writelines(gz) 
            gz.close()
            #print ('unziping done')
            #os.remove(filename)
            # openning pdb file 
            structure = parser.get_structure(protein,rootdir+protein[1:3]+'/pdb'+protein+'.ent')

             
            for model in structure:
                if model.id==0:
                 atom_list = Selection.unfold_entities(model, 'A') # A for atoms
                 ns = NeighborSearch(atom_list)
                 lig=[]
                 for chain in model:
                     for residue in chain:
                          if residue.resname not in AA and residue.resname not in CF:
                              #print(chain.id,residue.resname)
                         
                              if residue.resname in cofactors_dict:
                                nuc=0
                                for al in cofactors_dict[residue.resname]:
                                   # print(al) 
                                    #print(chain.id,residue.resname)
                                    atom_in_res=[]
                                    #for atom in residue:
                                     #   atom_in_res.append(atom.element)
                                       
                                    #if any(x in Metals for x in atom_in_res)==False:
                                        #print ('not metal')
                                     #   continue
                                     
                                    al=al.split(';')
                                    center = get_center(residue,al)
                                    #print ('center',center)
                                    lig=protein,chain.id,residue.id[1],residue.resname,center
                                    #print(lig)
                                    all_neighbors = ns.search(center, 15.0,"R") # 15.0 for distance in angstrom
                                    microfold_name=protein+'.'+residue.resname+'_'+ chain.id +'_'+str(residue.id[1])
                                    microfold_name=microfold_name.replace(' ','')
                                    microfold_name=microfold_name.replace('/','_')
                                    microfold_dir=residue.resname
                                    microfold_dir=microfold_dir.replace(' ','')                                                            
                                    if nuc==1:
                                        cof='ADE'
                                       # if residue.resname in NAD:
                                       #     cof=residue.resname
                                        microfold_name=protein+'.'+cof+'_'+ chain.id +'_'+str(residue.id[1])+'_'+residue.resname
                                        microfold_name=microfold_name.replace(' ','')
                                        microfold_name=microfold_name.replace('/','_')
                                        microfold_dir=cof
                                        microfold_dir=microfold_dir.replace(' ','')
                                   # print(microfold_name)
                                    if not os.path.exists('f:/microfolds_8_2018_1/organic/'+microfold_dir):
                                        os.makedirs('f:/microfolds_8_2018_1/organic/'+microfold_dir)
                                    Select = Bio.PDB.Select
                                    class MicroSelect(Select):
                                       def accept_residue(self, residue):
                                           if residue in all_neighbors and residue.resname!='HOH':
                                               return 1
                                           else:
                                               return 0
                                    io=PDBIO()
                                    io.set_structure(structure)
                                            
                                    io.save('f:/microfolds_8_2018_1/organic/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())
                                    
                                    nuc=1
                                continue
                              if residue.resname in chl:
                                
                                    atom_in_res=[]
                                    center = get_center(residue,pyr_atom_list)
                                    #print ('center',center)
                                    lig=protein,chain.id,residue.id[1],residue.resname,center
                                    #print(lig)
                                    all_neighbors = ns.search(center, 15.0,"R") # 15.0 for distance in angstrom
                                    microfold_name=protein+'.'+residue.resname+'_'+ chain.id +'_'+str(residue.id[1])
                                    microfold_name=microfold_name.replace(' ','')
                                    microfold_name=microfold_name.replace('/','_')
                                    microfold_dir=residue.resname
                                    microfold_dir=microfold_dir.replace(' ','')                                                            
                                    if not os.path.exists('f:/microfolds_8_2018_1/chl/'+microfold_dir):
                                        os.makedirs('f:/microfolds_8_2018_1/chl/'+microfold_dir)
                                    Select = Bio.PDB.Select
                                    class MicroSelect(Select):
                                       def accept_residue(self, residue):
                                           if residue in all_neighbors and residue.resname!='HOH':
                                               return 1
                                           else:
                                               return 0
                                    io=PDBIO()
                                    io.set_structure(structure)
                             
                                    io.save('f:/microfolds_8_2018_1/chl/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())
                                    
                                    continue
                              if residue.resname in heme:
                                
                                    atom_in_res=[]
                                    center = get_center(residue,pyr_atom_list)
                                    #print ('center',center)
                                    lig=protein,chain.id,residue.id[1],residue.resname,center
                                    #print(lig)
                                    all_neighbors = ns.search(center, 15.0,"R") # 15.0 for distance in angstrom
                                    microfold_name=protein+'.'+residue.resname+'_'+ chain.id +'_'+str(residue.id[1])
                                    microfold_name=microfold_name.replace(' ','')
                                    microfold_name=microfold_name.replace('/','_')
                                    microfold_dir=residue.resname
                                    microfold_dir=microfold_dir.replace(' ','')                                                            
                                    if not os.path.exists('f:/microfolds_8_2018_1/heme/'+microfold_dir):
                                        os.makedirs('f:/microfolds_8_2018_1/heme/'+microfold_dir)
                                    Select = Bio.PDB.Select
                                    class MicroSelect(Select):
                                       def accept_residue(self, residue):
                                           if residue in all_neighbors and residue.resname!='HOH':
                                               return 1
                                           else:
                                               return 0
                                    io=PDBIO()
                                    io.set_structure(structure)
                                       
                                    io.save('f:/microfolds_8_2018_1/heme/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())
                                   
                                    continue
                                
                              atom_in_res=[]
                              for atom in residue:
                                    atom_in_res.append(atom.element)
                                   
                              if any(x in Metals for x in atom_in_res)==False:
                                    #print ('not metal')
                                  continue
                                    
                              al=['all']
                              center = get_center(residue,al)
                                #print ('center',center)
                              lig=protein,chain.id,residue.id[1],residue.resname,center
                                #print(lig)
                              all_neighbors = ns.search(center, 15.0,"R") # 15.0 for distance in angstrom
                              microfold_name=protein+'.'+residue.resname+'_'+ chain.id +'_'+str(residue.id[1])
                              microfold_name=microfold_name.replace(' ','')
                              microfold_dir=residue.resname
                              microfold_dir=microfold_dir.replace(' ','')
                             # print(microfold_name)
                              if not os.path.exists('f:/microfolds_8_2018_1/metals/'+microfold_dir):
                                    os.makedirs('f:/microfolds_8_2018_1/metals/'+microfold_dir)
                              Select = Bio.PDB.Select
                              class MicroSelect(Select):
                                   def accept_residue(self, residue):
                                       if residue in all_neighbors and residue.resname!='HOH':
                                           return 1
                                       else:
                                           return 0
                              io=PDBIO()
                              io.set_structure(structure)
                                #print('/home/hraanan/MicrofoldsPDBs/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())                        
                              io.save('f:/microfolds_8_2018_1/metals/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())
    
                                
        except:
            
            Error_out.write( protein )
            continue
                               
    
Error_out.close()
#prot.close()
print("end")
