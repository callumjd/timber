#!/usr/bin/env python
import argparse
import glob
import sys
import os
import copy
import pandas as pd
import numpy as np
import pytraj as pt

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

# Run pmemd.decomp to get atomic free energies
# This is only intended for protein-ligand RBFE runs
# There is a mis-match in ddG due to ti_lam_sch etc not supported

###############################################################################
## Routines ##
###############################################################################

def check_file(file_in):
    output=False
    if os.path.exists(file_in) and os.path.getsize(file_in)>0:
        output=True
    return output

def mapping_tuples(csv_file):

    map_list=[]

    # assumes first line are comments
    with open(csv_file,'r') as f:
        f.readline()
        for line in f:
            if line[0]!='#':
                lig1=line.split(',')[0].strip()
                if len(line.split(','))==2:
                    lig2=line.split(',')[1].strip()
                    map_list.append((lig1,lig2))
                else:
                    map_list.append((lig1))

    return map_list

def write_cluster(file_name,media):

    text1='''
#!/bin/bash
#$ -N afep_cpu
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=604800,m_mem_free=4G

dir_list='afep_list.dat'
MY_DIR=$( sed -n "${SGE_TASK_ID}p" < $dir_list )
cd $MY_DIR

'''
    with open(file_name,'w') as f:
        f.write(text1)
        f.write('/home/dicksca3/AMBER_files/Decomp/amber22/bin/pmemd.decomp -O -i decomp_run.in -o mdout -c prod.nc -x newsnap.nc -ref prod.rst7 -p ../../../%s_ligands.prmtop -decomp decomp.log -e mden -r restart\n' % (media))
        f.write('rm newsnap.nc\n')
        f.write('\n')
        f.write('cd ../\n')

def write_afep_data(name,output_list):

    with open('complex/complex_ligands.pdb','r') as f:
        with open(name+'.pdb','w') as f_out:
            f_out.write('HEADER spectrum b, green_white_red\n')
            for line in f:
                if line.split()[0]=='ATOM':
                    at_idx=int(line.split()[1])
                    f_out.write(line[0:60]+'%6.2f\n' % (output_list[at_idx-1]))
                elif line.split()[0]=='TER':
                    f_out.write(line)

    with open(name+'.dat','w') as f_out:
        for i in range(0,len(output_list)):
            f_out.write('%d %lf\n' % (i+1,output_list[i]))

def write_decomp_input(prot_resi_range):

    with open('prod.in','r') as f:
        input_file=f.readlines()

    for i in range(0,len(input_file)):
        if 'nmropt = 0,' in input_file[i]:
            start_idx=i

    with open('decomp_run.in','w') as f:
        for i in range(0,len(input_file)):
            if ('gti_add_sc' not in input_file[i] and 'gti_lam_sch' not in input_file[i] and 'gti_scale_beta' not in input_file[i] and 'gti_vdw_exp' not in input_file[i] and 'gti_ele_exp' not in input_file[i]):
                if i==start_idx+1:
                    f.write('\n')
                    f.write('ntwd=1,reweight=1,\n')
                    f.write('ligmask=\":1,2\",\n')
                    f.write('proteinmask=\":3-%d\",\n' % (prot_resi_range+2))
                    f.write('decompmask=\"*\",\n')
                    f.write('\n')
                else:
                    f.write(input_file[i].replace('scalpha = 0.5','scalpha = 0.2').replace('scbeta = 1.0','scbeta = 12.0'))

def submit_afep(prot_resi_range,ti_rep):

    for media in ['complex','solvent']:
        os.chdir(media)
        for rep in range(1,ti_rep+1):
            os.chdir('one-step_rep%d/%s_ligands' % (rep,media))

            lambda_list=glob.glob('lambda_*')
            lambda_list.sort(key=lambda x: int(x.split('_')[-1]))

            for i in range(0,len(lambda_list)):
                os.chdir('lambda_%d' % (i))
                write_decomp_input(prot_resi_range)
                os.chdir('../')

            with open('afep_list.dat','w') as f:
                for val in lambda_list:
                    f.write('%s\n' % (val))

            write_cluster('afep_task.sh',media)
            os.system('qsub -t 1-%d afep_task.sh' % (len(lambda_list)))

            os.chdir('../../')
        os.chdir('../')

def get_lambda(filename):

    with open(filename,'r') as f:
        for line in f:
            if 'clambda' in line:
                lambda_value=float(line.split('clambda =')[1].split()[0].strip(','))

    return round(lambda_value,4)

def parse_decomp(fname):
    
    with open(fname,'r') as f:
        startlines=f.readlines()[0:5]
    
    n_steps=int(startlines[0].strip().split()[2])
    
    groups={}
    for i in range(1,len(startlines)):
        if 'Ligand-Ligand' in startlines[i]:
            groups.update({'Ligand-Ligand DV/DL':float(startlines[i].strip().split()[2])/n_steps})
        elif 'Ligand-Protein' in startlines[i]:
            groups.update({'Ligand-Protein DV/DL':float(startlines[i].strip().split()[2])/n_steps})
        elif 'Ligand-Cofactor' in startlines[i]:
            groups.update({'Ligand-Cofactor DV/DL':float(startlines[i].strip().split()[2])/n_steps})
        elif 'Ligand-Water/Bulk' in startlines[i]:
            groups.update({'Ligand-Water/Bulk DV/DL':float(startlines[i].strip().split()[2])/n_steps})
    
    df=pd.read_csv(fname,skiprows=5,delim_whitespace=True)
    
    for col in list(df.columns):
        if col!='Atom':
            df[col]=df[col]/n_steps
    
    return df,groups

def FE_atom_DVDL(atom_dvdl_array,lambda_values):

    weights=None
    # gaussian quadrature weights
    if lambda_values==[0.1127,0.5,0.88729]:
        weights=np.array([0.27777,0.44444,0.27777])
    elif lambda_values==[0.04691,0.23076,0.5,0.76923,0.95308]:
        weights=np.array([0.11846,0.23931,0.28444,0.23931,0.11846])
    elif lambda_values==[0.02544,0.12923,0.29707,0.5,0.70292,0.87076,0.97455]:
        weights=np.array([0.06474,0.13985,0.19091,0.20897,0.19091,0.13985,0.06474])
    elif lambda_values==[0.01592,0.08198,0.19331,0.33787,0.5,0.66213,0.80669,0.91802,0.98408]:
        weights=np.array([0.04064, 0.09032, 0.13031, 0.15617, 0.16512, 0.15617, 0.13031, 0.09032, 0.04064])
    elif lambda_values==[0.00922,0.04794,0.11505,0.20634,0.31608,0.43738,0.56262,0.68392,0.79366,0.88495,0.95206,0.99078]:
        weights=np.array([0.02359,0.05347,0.08004,0.10158,0.11675,0.12457,0.12457,0.11675,0.10158,0.08004,0.05347,0.02359])

    atom_ddg=[]

    for i in range(np.shape(atom_dvdl_array)[0]):
        dvdl_slice=atom_dvdl_array[i,:]
        if isinstance(weights,np.ndarray):
            run_sum=0
            for gs in range(0,len(lambda_values)):
                run_sum+=(dvdl_slice[gs]*weights[gs])
        else:
            run_sum=np.trapz(dvdl_slice,x=lambda_values)

        atom_ddg.append(run_sum)

    return atom_ddg

def process_afep(ti_rep):

    # FIXME: Also output the lig-lig ddG breakdowns

    # ddG per atom, avg over ti_rep ...

    store_ddg=[]
    for media in ['complex','solvent']:
        os.chdir(media)

        system=pt.load('%s_ligands.inpcrd' % (media),'%s_ligands.prmtop' % (media))
        atom_range=len(pt.select_atoms(system.top,'*'))

        # should match in either media
        if media=='complex':
            ligand_range=len(pt.select_atoms(system.top,':LIG,MOD'))

        media_ddg=np.zeros(atom_range,dtype=np.float64)

        for rep in range(1,ti_rep+1):
            atom_dvdl=[]

            os.chdir('one-step_rep%d/%s_ligands' % (rep,media))

            lambda_list=glob.glob('lambda_*')
            lambda_list.sort(key=lambda x: int(x.split('_')[-1]))

            lambda_values=[get_lambda(x+'/prod.out') for x in lambda_list]

            for i in range(0,len(lambda_list)):
                local_df,groups=parse_decomp('%s/decomp.log' % (lambda_list[i]))
                atom_dvdl.append(np.array(list(local_df['DV/DL'])))

            atom_ddg=FE_atom_DVDL(np.column_stack(atom_dvdl),lambda_values)
            media_ddg+=np.array(atom_ddg)

            os.chdir('../../')
        os.chdir('../')

        media_ddg=media_ddg/ti_rep
        store_ddg.append(media_ddg)

    print('AFEP ddG: %lf\n' % (np.sum(store_ddg[0])-np.sum(store_ddg[1])))

    final_ddg=copy.deepcopy(store_ddg[0])

    # for just the ligands, correct complex-solv
    for i in range(0,len(final_ddg)):
        if i<ligand_range:
            final_ddg[i]=final_ddg[i]-store_ddg[1][i]

    return final_ddg

###############################################################################
# RUN
###############################################################################

def run(input_csv,overwrite,name='decomp_afep'):

    if not check_file(input_csv):
        print('Error: cannot find %s!\n' % (input_csv))
        sys.exit()

    map_list=mapping_tuples(input_csv)

    for pair in map_list:
        dir_name=pair[0]+'~'+pair[1]
        os.chdir(dir_name)

        # get number of residues in the protein
        com=pt.load('complex/complex_ligands.inpcrd','complex/complex_ligands.prmtop')
        prot_resi_range=len(pt.select_atoms(com.top,'@N'))

        # get ti repeats
        ti_rep=len(glob.glob('complex/one-step_rep*'))

        # check if we have run decomp
        myfile='complex/one-step_rep1/complex_ligands/lambda_0/decomp.log'

        if not os.path.exists(myfile):
            # Submit the afep analysis jobs
            print('\nSubmit AFEP jobs: %s' % (dir_name))
            submit_afep(prot_resi_range,ti_rep)
        elif overwrite:
            # Submit the afep analysis jobs
            print('\nOverwrite AFEP jobs: %s' % (dir_name))
            submit_afep(prot_resi_range,ti_rep)
        else:
            # Process the outputs
            print('\nProcess AFEP jobs: %s' % (dir_name))
            atom_ddg_output=process_afep(ti_rep)
            write_afep_data(name,atom_ddg_output)            

        os.chdir('../')

###############################################################################
## Parse arguments ##
###############################################################################

def main(argv=None):
    ## Command line arguments ##
    parser = argparse.ArgumentParser(description='Wrapper for atomic contribution analysis (AFEP)\n')

    parser.add_argument('-i','--input',dest='i',help='CSV file with ligand mappings',type=str,required=True)
    parser.add_argument('-o','--overwrite',dest='o',help='Overwrite existing analysis',action='store_true')
    parser.add_argument('-n','--name',dest='n',help='Output name',type=str,required=False,default='afep_results')

    args=vars(parser.parse_args())

###############################################################################
## Run selected protocols ##
###############################################################################

    run(input_csv=args['i'],overwrite=args['o'],name=args['n'])

# MAIN
if __name__ == '__main__':
    main()

