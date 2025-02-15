# timber

import os
import sys
import lomap
import pandas as pd
import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter
from rdkit.Chem import PandasTools
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from .ligprep_tools import check_file

##############################################################################

# no radial star map option - separate code
# TO DO: use molecule names in output image
def run_lomap(input_sdf,dir_name='lomap_dir',output_name='mapping'):

    print('Running LOMAP\n')

    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    # load molecules
    suppl=Chem.SDMolSupplier(input_sdf,removeHs=False)

    for m in suppl:
        name=m.GetProp('_Name')
        writer=SDWriter(dir_name+'/'+name+'.sdf')
        writer.write(m)
        writer.flush()

    # allow charge changing with ecrscore
    # NB removed stereo in : /usr/prog/cadd/amber_tools/alchemistry2/lib/python3.9/site-packages/lomap/graphgen.py (search cjd edit)
    db_mol = lomap.DBMolecules(dir_name, output=True, ecrscore=0.0, name=output_name)

    strict, loose = db_mol.build_matrices()

    nx_graph = db_mol.build_graph()

    # move the lomap files
    os.system('mv %s.* ./%s' % (output_name,dir_name))
    os.system('mv %s_score_with_connection.txt ./%s' % (output_name,dir_name))

    # input CSV file for timber
    df=pd.read_csv(dir_name+'/'+output_name+'_score_with_connection.txt',skiprows=1,header=None,usecols=(0,1,2,3,4,5,6,7),names=['Index_1','Index_2','Name1','Name2','Str_sim','Eff_sim','Loose_sim','Connect'])

    # fix the molecule names
    for index,row in df.iterrows():
        df.at[index,'Name1']=str(row['Name1']).strip().replace('.sdf','').replace('.mol2','')
        df.at[index,'Name2']=str(row['Name2']).strip().replace('.sdf','').replace('.mol2','')
        df.at[index,'Connect']=str(row['Connect']).strip()

    df=df[df['Connect']=='Yes']
    df=df.reset_index(drop=True)

    df[['Name1','Name2']].to_csv(output_name+'.csv',index=False)

def connect_rbfe_map(mapping_csv,input_sdf):

    if not check_file(mapping_csv):
        raise Exception('Error: cannot find %s!\n' % (mapping_csv))

    if not check_file(input_sdf):
        raise Exception('Error: cannot find %s!\n' % (input_sdf))

    df=pd.read_csv(mapping_csv)

    if len(list(df.columns))!=2:
        raise Exception('Error: RBFE map requires two columns!\n')

    name1_col=list(df.columns)[0]
    name2_col=list(df.columns)[1]

    # create a graph of the initial mapping
    G = nx.Graph()

    for index,row in df.iterrows():
        name1=row[name1_col]
        name2=row[name2_col]

        if name1 not in list(G.nodes()):
            G.add_node(name1)
        if name2 not in list(G.nodes()):
            G.add_node(name2)

    for index,row in df.iterrows():
        G.add_edge(row[name1_col],row[name2_col])

    if nx.is_connected(G):
        return None
    else:
        print('Warning: network is disconnected. Adding nodes.\n')
        ligand_df=PandasTools.LoadSDF(input_sdf)
        updated_G=add_nodes_rbfe_map(G.copy(),ligand_df)
        print('Over-writing %s\n' % (mapping_csv))
        write_mapping_file(mapping_csv,updated_G)

def add_nodes_rbfe_map(G,ligand_df):

    def rdkit_smi(mol):
        return AllChem.MolToSmiles(Chem.RemoveHs(mol))

    ligand_df['Smiles']=ligand_df['ROMol'].apply(rdkit_smi)
    PandasTools.AddMoleculeColumnToFrame(ligand_df,'Smiles','ROMol')
    ligand_df['Fingerprint']=None
    for index,row in ligand_df.iterrows():
        ligand_df.at[index,'Fingerprint']=FingerprintMols.FingerprintMol(row['ROMol'], minPath=1, maxPath=7, fpSize=2048,bitsPerHash=2, useHs=True, tgtDensity=0.0,minSize=128)

    # disconnected components, sorted by graph size
    components = sorted([G.subgraph(c).copy() for c in nx.connected_components(G)],key = lambda x: len(x.nodes()),reverse=True)

    # get the largest cycle
    core_G_names=[]
    core_G_fp=[]
    for node in list(components[0].nodes()):
        core_G_names.append(node)
        idx=ligand_df[ligand_df['ID']==node].index.tolist()[0]
        core_G_fp.append(ligand_df.loc[idx,'Fingerprint'])

    # for each other cycle, connect first ligand to closest similar core_G cycle
    for cyc in components[1:]:
        name=list(cyc.nodes())[0]
        idx=ligand_df[ligand_df['ID']==name].index.tolist()[0]
        fp=ligand_df.loc[idx,'Fingerprint']
        sim_idx=np.argmax(DataStructs.BulkTanimotoSimilarity(fp,core_G_fp))
        print('Connecting %s to %s\n' % (name,core_G_names[sim_idx]))
        G.add_edge(name,core_G_names[sim_idx])

    if not nx.is_connected(G):
        raise Exception('Error: could not fully connect network!\n')

    return G

def write_mapping_file(output_file,Tree,reverse=False):

    with open(output_file,'w') as f:
        f.write('Name1,Name2\n')
        for val in nx.generate_edgelist(Tree, data=False):
            name1=val.split()[0].strip()
            name2=val.split()[1].strip()
            line=name1+','+name2+'\n'
            f.write(line)

            if reverse:
                line_rev=name2+','+name1+'\n'
                f.write(line_rev)

