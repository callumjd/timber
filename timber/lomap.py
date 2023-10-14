# timber

import os
import lomap
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDWriter

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

