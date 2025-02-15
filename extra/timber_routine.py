#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
import argparse
import json
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter,rdmolops
import timber

##############################################################################

def run(mode):

    # input settings
    input_sdf='T4LysozymeL99A_ligands.sdf'
    prot='prot.pdb'
    extra=[]
    protocol='one-step'
    ff='gaff2'
    ti_repeats=1
    equil_ns=1
    prod_ns=2
    monte_wat=0
    hmass=True
    protein_ff='ff19SB'
    water_ff='tip3p'

    exp_data=None
    col=4
    units='uM'

    res_seq=0 # abfe run - single residue number or list of 3 protein atoms
    lig_boresch_smarts=None # optional smarts of 3 ligand atoms

    # this is a dict, setting lambda windows
    #with open('schedule.json','r') as f:
        #schedule=json.load(f)

    schedule={'complex_ligands':9,  # one-step, three-step vdw, absolute
          'solvent_ligands':9,
          'complex_decharge':5,  # three-step, absolute-three-step
          'solvent_decharge':5,
          'complex_recharge':5,  # three-step
          'solvent_recharge':5,
          'complex_restraint':5} # absolute

    ### MODE: SUBMIT
    if mode=='submit':

        # make a lomap mapping file
        timber.run_lomap(input_sdf,dir_name='lomap_dir',output_name='mapping')
        timber.connect_rbfe_map('mapping.csv',input_sdf)

        # load molecules
        suppl=Chem.SDMolSupplier(input_sdf,removeHs=False)
        all_mols=[m for m in suppl]
        all_names=[str(m.GetProp('_Name')) for m in suppl]

        suppl_noH=Chem.SDMolSupplier(input_sdf,removeHs=True)
        all_mols_noH=[m for m in suppl_noH]

        # MCSS for ligand alignment
        full_mcs=timber.get_mcs(all_mols_noH,strict=True)

        # Align all molecules to first one in the set
        for m in all_mols:
            if m.GetProp('_Name')!=all_names[0]:
                timber.rms_fit(all_mols[0],m,mcss=full_mcs.smartsString,mcss_exclusion=None,bak_seed=None,tolerance=2.0,ene_cutoff=500,initial_rms_tol=1.0)

        # read in the lomap network
        df=pd.read_csv('mapping.csv')
        df=timber.clean_mapping_file(df)

        # prepare all dir lig1->lig2
        for index,row in df.iterrows():
            name1=str(row['Name1'])
            if protocol in ['one-step','three-step']:
                name2=str(row['Name2'])
                timber.run_rbfe_setup(all_mols[all_names.index(name1)],all_mols[all_names.index(name2)],ff=ff,full_mcs=full_mcs.smartsString,align=True)
            elif protocol in ['absolute','absolute-three-step']:
                # for absolute runs
                timber.run_abfe_setup(all_mols[all_names.index(name1)],ff=ff)

        # get ion numbers for salt concentration. water number is approximated as prot_res*50
        prot_chg,prot_res=timber.protein_charge(prot)
        pos_ion,neg_ion=timber.return_salt(nwat=prot_res*50,conc=0.15,charge=prot_chg+int(rdmolops.GetFormalCharge(all_mols[0])))

        # make build.leap file and build the prmtop
        prep_files,pdb_files=timber.parse_extra(extra,prot)
        timber.build_ti_leap(prot,prep_files=prep_files,pdb_files=pdb_files,protocol=protocol,pos_ion=pos_ion,neg_ion=neg_ion,ff=ff,protein_ff=protein_ff,water_ff=water_ff,ion_ff='ionsjc_tip3p')

        timber.run_build(df,protocol=protocol,hmass=hmass,use_openff=(ff=='openff'))

        if protocol in ['absolute','absolute-three-step']:
            # absolute run: write DISANG file - requires an amber numbered res_seq to set residue index for Boresch restraints
            timber.write_abfe_disang(df,protocol=protocol,res_seq=res_seq,lig_boresch_smarts=None)

        # submit windows
        timber.run_prod(df,protocol=protocol,ti_repeats=ti_repeats,schedule=schedule,hmass=hmass,equil_ns=equil_ns,prod_ns=prod_ns,monte_water=monte_wat)

        # example to extend the prod sampling (using same prod.in file)
        #timber.extend_prod(df)

    elif mode=='analysis':
        # read in the lomap network
        df=pd.read_csv('mapping.csv')

        # run_analysis once runs are complete
        timber.run_analysis(df)

        # cinnabar
        if exp_data:
            timber.run_cinnabar(df,data=exp_data,col=col,units=units)
        else:
            timber.run_cinnabar(df)

        # endpoint analysis - currently only for RBFE
        #timber.run_endpoint_mdanalysis(df,prot)

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='Full TIMBER routines (RBFE/ABFE)\n')
    parser.add_argument('-m','--mode',dest='m',help='Mode: submit or analysis',type=str,required=True,choices=['submit','analysis'])

    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################
    run(mode=args['m'])

# MAIN
if __name__ == '__main__':
    main()

