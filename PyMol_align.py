# -*- coding: utf-8 -*-
"""
This script dose pairwise alignment of microenvironmets using PyMol align package. 
input: a list of pairs of microenvironments and the microenvironments PDB files.
output: a file containing list of alignment results. 
@author: hraanan
"""
print ('welcom to pymol align')
import os
import time
import numpy
import math
import __main__
import sys
import center 

__main__.pymol_argv = ['pymol','-c'] # Pymol: quiet and no GUI
#from time import sleep
import pymol
pymol.finish_launching()

def file_len(count_file_name):
    lines=0
    count_file=open(count_file_name,'r')
    for line in count_file:
        lines=lines+1
    
    count_file.close()
    return lines


dirpath = os.getcwd()

in_file_name=sys.argv[1]
#in_file_name='/home/hraanan/span/align_lists/align_cofactors.txt'
in_file=open(in_file_name,'r')
#print(in_file_name)
microen_dir='/scratch/hr253/span/new_all/'


count_file_name=dirpath+'/cnt/cnt_'+in_file_name #_'+in_file_name.split('.')[1]+'.txt'
count_file_len=0
if os.path.exists(count_file_name)==False:
    count_file=open(count_file_name,'w')
else:
    count_file_len=file_len(count_file_name)
    count_file=open(count_file_name,'a')

#out_file.write('Sourec\tTarget\tSl\tTl\tLigand\tRMSD\tAlign CA\tRaw alignment score\tAligned Residues\tLigand center distance'+'\t'+'Structural Distance'+'\n')


out_file_name=dirpath+'/out/out_'+in_file_name+str(count_file_len) #_'+in_file_name.split('.')[1]+'.txt'
out_file=open(out_file_name,'a')              
err_file_name=dirpath+'/err/err_'+in_file_name+str(count_file_len) #_'+in_file_name.split('.')[1]+'.txt'
err_file=open(err_file_name,'a')              

run_file_name=dirpath+'/run/run_'+in_file_name #+str(count_file_len) #_'+in_file_name.split('.')[1]+'.txt'
run_file=open(run_file_name,'w')              
run_file.close()

for i in range(count_file_len):
        in_file.next()
t=0
start = time.time()    

print('start align from:'+str(count_file_len))    
for line in in_file:
    
    pymol.cmd.reinitialize()
    try:
        #print(line)
            
        t=t+1
        count_file.write(str(t+count_file_len)+'\n')
        #if t>1000:
         #   break
        line=line.split('\t')
        microen1=line[0][:-4]
        microen2=line[1][:-6]
        chain1=line[0].split('_')[-2]
        chain2=line[1].split('_')[-2]
        #print(chain1,chain2)
        #PDB2='3cw9.01A_B_991'
        lig=line[0].split('_')[0]
        lig1=lig.split('.')[1]
        Fldr1=microen_dir+lig1+'/'
        if lig1=='ADE':
            lig1=lig1+'_'+line[0].split('_')[1]
        
        #PDB1=lig1.split('.')[0]
        
        lig=line[1].split('_')[0]
        lig2=lig.split('.')[1]
        Fldr2=microen_dir+lig2+'/'
        if lig2=='ADE':
            lig2=lig2+'_'+line[0].split('_')[1]
        #PDB2=lig2.split('.')[0]
        
        
        
        pdbFile_1=str(Fldr1)+microen1+".pdb"
        pdbFile_2=str(Fldr2)+microen2+".pdb"
        pymol.cmd.load(pdbFile_1,'PDB1')	
        pymol.cmd.load(pdbFile_2,'PDB2')
        Ql=pymol.cmd.count_atoms('PDB1 and name CA')
        Tl=pymol.cmd.count_atoms('PDB2 and name CA')
        if Ql<10 or Tl<10:
            #print('too short microen')
            continue
        #print('Ql,Tl:'+str(Ql)+str(Tl))
        res_num1=microen1.split('_')
        res_num2=microen2.split('_')
       
        res_num1=res_num1[-1]
        res_num2=res_num2[-1]             
        #print('res_num1:',res_num1)                    
        #print('res_num2:',res_num2)
       
        x=pymol.cmd.align('PDB1 and name CA','PDB2 and name CA',quiet=1)
        
        align_out=[x[0],x[1],x[5],x[6]]
        #print(align_out)
        atomlist=[]                    
        #print(lig1,lig2) 
        cof_atom_list=center.get_atom_list(lig1)
        if cof_atom_list==['all']:
            atomlist=pymol.cmd.get_model('PDB1 and chain '+chain1+' and resi '+res_num1, 1).get_coord_list()
        else:
            for atom in cof_atom_list:
                #print('PDB1 and chain '+chain1+' and resi '+res_num1+' and name '+atom)
                atomlist=atomlist+(pymol.cmd.get_model('PDB1 and chain '+chain1+' and resi '+res_num1+' and name '+atom, 1).get_coord_list())
        #print(atomlist)
        center1=center.get_center(atomlist)    
        #print(center1)       
        if center1=='NA':
            #print('center1=NA')    
            continue
        atomlist=[]                    
                                
        cof_atom_list=center.get_atom_list(lig2)
        if cof_atom_list==['all']:
            atomlist=pymol.cmd.get_model('PDB2 and chain '+chain2+' and resi '+res_num2, 1).get_coord_list()
        else:
            for atom in cof_atom_list:
                atomlist=atomlist+(pymol.cmd.get_model('PDB2 and chain '+chain2+' and resi '+res_num2+' and name '+atom, 1).get_coord_list())
        center2=center.get_center(atomlist)    
        #print(center2)         
        if center2=='NA':
            #print('center2=NA')        
            continue
        atomlist=[]             
        #print(center1,center2)
        Dis=math.sqrt((pow((center1[0]-center2[0]),2))+(pow((center1[1]-center2[1]),2))+(pow((center1[2]-center2[2]),2)))
        #print(Dis)
        if Dis>15:
           # print('large distance')
            continue
        #print(x)    
        D= Ql+Tl-(2*x[6])  
        align_out=str(Ql)+'\t'+str(Tl)+'\t'+lig1+'_'+lig2+'\t'+str(x[0])+'\t'+str(x[1])+'\t'+str(x[5])+'\t'+str(x[6])+'\t'+str(Dis)+'\t'+str(D)
        #print(align_out)
        out_file.write(microen1+'\t'+microen2+'\t'+align_out+'\n')	
    except Exception as e:
        err_file.write('\t'.join(line))
        #err_file.write(e)
        continue 
pymol.cmd.quit()
end = time.time()
print ('Run time is:'+str(end - start))            
in_file.close()
err_file.close() 
out_file.close()
os.remove(run_file_name)
print('run_file_removed')


	

	

