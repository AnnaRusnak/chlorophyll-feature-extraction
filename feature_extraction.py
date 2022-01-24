#importing 
from Bio.PDB import NeighborSearch, PDBParser, Selection, StructureBuilder, xpdb, PDBIO
import pandas as pd
import csv

io=PDBIO()

import os
files=os.listdir("/Users/annarusnak/Dropbox (Brown)/Grad School classes/Fall 2021/project/pdb_files/CLA:CHL") #all the files in a folder 
#path=os.getcwd()

tail_atoms=[]
for i in range(1,21):
    tail_atoms.append("<Atom C" + str(i)+">")

files.sort()

P = ['C', 'N', 'O', 'S']
L = ['C', 'N', 'O', 'MG']
    
PL=[]

for i in range(0, len(L)):
    for j in range(0, len(P)):
        PL.append(L[i]+P[j])
        
columns=PL
columns.append('label')
df = pd.DataFrame()
csv_file = 'features.txt'

dict_data=[]

for filename in files:
    if filename.endswith(".ent"):
        #print(filename) #-> if you wanna check what pdbs you're working on 

    
         ####################################################################################################
        #  if pdb is too large and residue numbers repeat themselves, and sloppy parser gives you a warning  #
        #  and renumbers them                                                                                #
        #  meaning the residue numbers in original pdb and residue numbers in output dataframe are gonna     #
        #  be different                                                                                      #
         ####################################################################################################
            
    # parsing a pdb in     
        sloppyparser = PDBParser(PERMISSIVE=True, structure_builder=xpdb.SloppyStructureBuilder())
        structure_id = "pdb"
        print('doing filename ', filename)
        structure = sloppyparser.get_structure(structure_id, filename)
    
        
        #constructing a cKDTree
        all_atoms=Selection.unfold_entities(structure, 'A')
        atoms_tree = NeighborSearch(all_atoms)
    
        chl_kinds=['CHL', 'CLA']

     # looping through pdb 
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in chl_kinds:
                        dict={}
                        for val in PL:
                            dict[val]=0
                        #define a full pdb name of chlorophyll
                        residuename=residue.get_resname()
                        residuenumber=residue.id[1]
                        chlorophyll_id=residuename+" "+str(chain.id)+" "+str(residuenumber)
                        
                        chlorophyll_atoms=Selection.unfold_entities(residue, 'A') # all atoms in a chlorophyll
                        
                        for atom in chlorophyll_atoms:
                            if str(atom) in tail_atoms:
                                pass
                            else:
                                #print(atom)
                                nearby_atoms=atoms_tree.search(atom.coord, 5) # all neighbors within 5A
                                
                                for neighbor in nearby_atoms:
                                    parent=neighbor.get_parent()
                                    if parent.get_resname() != "HOH" and parent.get_resname()+" "+str(parent.get_full_id()[2])+" "+str(parent.id[1]) != chlorophyll_id:
                                        if atom.element+neighbor.element in dict:
                                            dict[atom.element+neighbor.element] = dict[atom.element+neighbor.element]+1
                                            dict['label']=residuename
                        #dict_data.append(dict)
                        df=df.append(dict, ignore_index = True)
                        #print(dict)

#print(dict_data[-10])
df = df[df.label != 0]
print(df.head(n=10))

df.to_csv('features.csv')

                   
                                

