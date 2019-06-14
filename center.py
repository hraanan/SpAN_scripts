# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 11:03:05 2018
This script returns information about cofactors
input:file contaning the atoms for each cofactor (manual_cofactor_atoms_list.txt)
output: returns information about cofactors
@author: hraanan
"""
import numpy
import time
#import sys

cof_dict={}
atoms_file=open('f:/programs/span/data/manual_cofactor_atoms_list.txt','r')
atoms_file.readline()
al_dict={}
for line in atoms_file:
    line=line.split('\t')
    cof_dict[line[1]]=line[2]
    if line[3]=='all':
        continue
    al_dict[line[1]]=[line[3],line[4][:-1]]

#time.sleep(5)

def get_atom_list(cof):
    cof=cof.upper()
    if cof not in al_dict and cof[0:3]!='ADE':
        return ['all']
    ade=0 #if not ade read the first list of atoms 
    cof=cof.split('_')
    #print(cof)
    if cof[0]=='ADE':
        ade=1 #if not ade read the second list of atoms 
        cof=cof[1]
    else:
        cof=cof[0]
    
    #print(ade)
    atom_list=al_dict.get(cof)[ade]
    return atom_list.split(';')
     
def get_atom_coord_list(residue,al): # get the geometric centar of list of atom coordinatins 
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
    return coord
     
def get_center(atom_coord_list):
    x=0
    y=0
    z=0    
        
    for atom in atom_coord_list:
        
        x=x+atom[0]
        y=y+atom[1]
        z=z+atom[2]
        
    try:
            x=x/len(atom_coord_list)
            y=y/len(atom_coord_list)
            z=z/len(atom_coord_list)
            center=numpy.array([x,y,z])
    except:
            return 'NA'   
    center=numpy.array([x,y,z])
    return center;

def has_ade(cof):
    cof=cof.upper()
    if al_dict.get(cof)[1]!='na':
        return True
    else:
        return False
def is_in_list(cof):
    cof=cof.upper()
    if cof in al_dict:
        return True
    else:
        return False
def cof_type(cof):
    cof=cof.upper()
    if is_in_list(cof):
        return cof_dict.get(cof)
    else:
        return 'NA'

print('center.py imported successfully')

#test get_atom_list function

#print(is_in_list('sam'))
#print(get_atom_list('MTE'))
#print(get_atom_list('NAD'))
#print(get_atom_list('ADE_FAD'))
#
#print(has_ade('NAD'))
#print(has_ade('CLA'))
#print(cof_type('cla'))
