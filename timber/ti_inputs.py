# timber

import os
import glob

##############################################################################
# Modify the cluster submit script and MD inputs as needed
##############################################################################

def clean_mapping_file(df):

    todrop=[]
    name1_col=list(df.columns)[0]
    for index,row in df.iterrows():
        if '#' in str(row[name1_col]):
            todrop.append(index)

    df=df.drop(index=todrop)
    df=df.reset_index(drop=True)
    return df

def write_ti_cluster_script(prmtop,nodel=False):
    # h_rt=259100 gives 3 days / 72 hrs

    text1='''#!/bin/bash

#$ -N TI
#$ -cwd
#$ -S /bin/bash
#$ -l cuda=1
#$ -l gpu_card=1
#$ -l h_rt=259100,m_mem_free=2G
#$ -j y
#$ -q default.q

module purge
module load Amber/22.4-foss-2021b-AmberTools-22.5-CUDA-11.4.1
source /usr/prog/Amber/22.4-foss-2021b-AmberTools-22.5-CUDA-11.4.1/amber.sh

# SET CUDA_VISIBLE_DEVICES
export CUDA_VISIBLE_DEVICES=`/cm/shared/apps/nibri/sciComp/get_gpu_map.sh`

'''

    text2='''

dir_list='dir_list.dat'
MY_DIR=$( sed -n "${SGE_TASK_ID}p" < $dir_list )

cd $MY_DIR

pmemd -O -i min.in -p $prmtop -c $inpcrd -ref $inpcrd -o min.out -e min.en -inf min.info -r min.rst7

pmemd.cuda -O -i heat.in -p $prmtop -c min.rst7 -ref min.rst7 -o heat.out -e heat.en -inf heat.info -r heat.rst7 -x heat.nc

pmemd -O -i cpu_press.in -p $prmtop -c heat.rst7 -ref heat.rst7 -o cpu_press.out -e cpu_press.en -inf cpu_press.info -r cpu_press.rst7 -x cpu_press.nc

pmemd.cuda -O -i press.in -p $prmtop -c cpu_press.rst7 -ref cpu_press.rst7 -o press.out -e press.en -inf press.info -r press.rst7 -x press.nc

pmemd.cuda -O -i equil.in -p $prmtop -c press.rst7 -ref press.rst7 -o equil.out -e equil.en -inf equil.info -r equil.rst7 -x equil.nc

pmemd.cuda -O -i prod.in -p $prmtop -c equil.rst7 -ref equil.rst7 -o prod.out -e prod.en -inf prod.info -r prod.rst7 -x prod.nc

'''

    text3='''
for f in heat.nc cpu_press.nc press.nc equil.nc
do
if test -f ${f};
then
rm ${f}
fi
done

# If KEEPTRAJ file is present, keep the "prod.nc" file
if ! test -f KEEPTRAJ;
then
if test -f prod.nc;
then
rm prod.nc
fi
fi

cd ../

'''

    with open('run_prod.sh','w') as f:
        f.write(text1)
        f.write('export prmtop="%s"\n' % (prmtop))
        f.write('export inpcrd="%s"\n' % (prmtop.replace('prmtop','inpcrd')))
        f.write(text2)
        if not nodel:
            f.write(text3)

def write_ti_cluster_extend():

    extended=len(glob.glob('lambda_0/extend*en'))

    with open('run_prod.sh','r') as f:
        data=f.readlines()

    with open('run_extend.sh','w') as f:
        for line in data:
            if 'pmemd' not in line and 'cd ../' not in line:
                f.write(line)

        if extended==0:
            f.write('pmemd.cuda -O -i prod.in -p $prmtop -c prod.rst7 -ref prod.rst7 -o extend0.out -e extend0.en -inf extend0.info -r extend0.rst7 -x extend0.nc\n')
        else:
            f.write('pmemd.cuda -O -i prod.in -p $prmtop -c extend%d.rst7 -ref extend%d.rst7 -o extend%d.out -e extend%d.en -inf extend%d.info -r extend%d.rst7 -x extend%d.nc\n' % (extended,extended,extended+1,extended+1,extended+1,extended+1,extended+1,extended+1))

        f.write('\n')
        f.write('cd ../')
        f.write('\n')

def write_ti_inputs(protocol,schedule,ti_masks,ti_mask_len,lambda_val,all_lambda_list,hmass=True,equil_ns=1,prod_ns=5,monte_water=0,n_atoms1=0,n_atoms2=0):

    equil_ns=float(equil_ns)
    prod_ns=float(prod_ns)
   
    # monte carlo water sampling settings
    if hmass:
        nmc = 500
        nmd = 500
    else:
        nmc = 1000
        nmd = 1000

    # for MBAR
    all_lambda=''
    if len(all_lambda_list)>0:
        for val in all_lambda_list:
            all_lambda=all_lambda+str(val)+' '

    # prepare settings depending on protocol
    if hmass:
        dt=0.004
        equil_nstlim=int((equil_ns/4)*1e6)
        prod_nstlim=int((prod_ns/4)*1e6)
    else:
        dt=0.002
        equil_nstlim=int((equil_ns/2)*1e6)
        prod_nstlim=int((prod_ns/2)*1e6)

    icfe=1 # turn on TI
    disang=''

    # one-step 
    if protocol=='one-step':
        ti1=':1'
        ti2=':2'
        charge_mask=''
        ifsc=1
        nmropt=0
        noshakemask=':1,2'

    # three-step
    elif protocol=='three-step':
        ti1=':1'
        ti2=':2'
        nmropt=0
        noshakemask=':1,2'

        if schedule in ['complex_ligands','solvent_ligands']:
            ifsc=1
            chg1='@'+str(n_atoms1-ti_mask_len[0]+1)+'-'+str(n_atoms1)
            chg2='@'+str(n_atoms1+n_atoms2-ti_mask_len[1]+1)+'-'+str(n_atoms1+n_atoms2)
            charge_mask=chg1+'|'+chg2
        elif schedule in ['complex_decharge','solvent_decharge']:
            ifsc=0
            charge_mask=':2@'+ti_masks[0].split('@')[1]
            ti_masks=['',''] # over-ride ti_masks
        elif schedule in ['complex_recharge','solvent_recharge']:
            ifsc=0
            charge_mask=':1@'+ti_masks[1].split('@')[1]
            ti_masks=['',''] # over-ride ti_masks

    # absolute
    elif protocol in ['absolute','absolute-three-step']:
        nmropt=0

        if schedule=='complex_restraint':
            ti1=':1'
            ti2=':2'
            charge_mask=''
            ifsc=0
            noshakemask=':1,2'
            nmropt=1
            disang='../../../../disang.shift.RST'

        elif schedule in ['complex_decharge','solvent_decharge']:
            ti1=':1'
            ti2=':2'
            charge_mask=':2'
            ifsc=0
            noshakemask=':1,2'
            if schedule=='complex_decharge':
                nmropt=1
                disang='../../../../disang.pair.RST'

        elif schedule in ['complex_ligands','solvent_ligands']:
            if schedule=='complex_ligands':
                nmropt=1
                disang='../../../../disang.RST'

            ti_masks=[':1',''] # over-ride ti_masks
            ti1=':1'
            ti2=''
            ifsc=1
            noshakemask=':1'

            if protocol=='absolute':
                charge_mask=''
            elif protocol=='absolute-three-step':
                charge_mask=':1'

    # input file templates
    minimise='minimisation\n' \
    '&cntrl\n' \
    'imin = 1, ntmin = 2,\n' \
    'maxcyc = 1000, ncyc= 500,\n' \
    'ntpr = 20, ntwe = 20,\n' \
    'ntb = 1,\n' \
    'ntr = 1, restraint_wt = 25.00,\n' \
    'restraintmask=\"!:WAT,Na+,K+,Cl-\",\n' \
    '\n' \
    'gti_add_sc = 1,\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.2, scbeta = 17.3,\n' \
    'logdvdl = 0,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
    '/\n' \
    '\n' \
    % (icfe,lambda_val,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1])

    heat='heating\n' \
    '&cntrl\n' \
    'imin = 0, nstlim = 20000, irest = 0, ntx = 1, dt = 0.001,\n' \
    'nmropt = 1,\n' \
    'ntt = 1, temp0 = 300.0, tempi = 5.0, tautp = 1.0,\n' \
    'ntb = 1,\n' \
    'ntc = 2, ntf = 1,\n' \
    'ioutfm = 1, iwrap = 1,\n' \
    'ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,\n' \
    'cut = 9,\n' \
    '\n' \
    'ntr = 1, restraint_wt = 25.00,\n' \
    'restraintmask=\"!:WAT,Na+,K+,Cl-\",\n' \
    '\n' \
    'gti_add_sc = 1,\n' \
    'tishake = 1, noshakemask=\"%s\",\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.2, scbeta = 17.3,\n' \
    'logdvdl = 0,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
    '/\n' \
    '&ewald\n' \
    '/\n' \
    '&wt\n' \
    'type=\"TEMP0\",\n' \
    'istep1 = 0, istep2 = 8000,\n' \
    'value1 = 5.0, value2 = 300.0,\n' \
    '/\n' \
    '&wt type = \"END\"\n' \
    '/\n' \
    '&wt\n' \
    'type=\"END\",\n' \
    '&end\n' \
    'DISANG=%s\n' \
    '/\n' \
    '\n' \
    % (noshakemask,icfe,lambda_val,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1],disang)

    press='pressurising\n' \
    '&cntrl\n' \
    'imin = 0, nstlim = 100000, irest = 1, ntx = 5, dt = 0.001,\n' \
    'temp0 = 300.0, tautp = 1.0,\n' \
    'ntp = 1, pres0 = 1.0, taup = 2.0, ntt = 3, gamma_ln = 2.0,\n' \
    'ntb = 2,\n' \
    'ntc = 2, ntf = 1,\n' \
    'ioutfm = 1, iwrap = 1,\n' \
    'ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,\n' \
    'cut = 9, barostat = 2,\n' \
    '\n' \
    'ntr = 1, restraint_wt = 25.00,\n' \
    'restraintmask=\"!:WAT,Na+,K+,Cl-\",\n' \
    'nmropt=%d,\n' \
    '\n' \
    'gti_add_sc = 1,\n' \
    'tishake = 1, noshakemask=\"%s\",\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.2, scbeta = 17.3,\n' \
    'logdvdl = 0,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
    '/\n' \
    '&ewald\n' \
    'skinnb=5.0,\n' \
    '/\n' \
    '&wt\n' \
    'type=\"END\",\n' \
    '&end\n' \
    'DISANG=%s\n' \
    '/\n' \
    '\n' \
    % (nmropt,noshakemask,icfe,lambda_val,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1],disang)

    # older SC potentials to prevent crash
    equil='Equil %.3f ns NPT\n' \
    '&cntrl\n' \
    'imin = 0, ntx = 5, irest = 1,\n' \
    'nstlim = %d, dt = %.3f,\n' \
    'ntwe = %d, ntwx = %d, ntpr = %d, ntwr = %d,\n' \
    'ntc = 2, ntf = 1,\n' \
    'temp0 = 300.0,\n' \
    'ntp = 1, pres0 = 1.0, taup = 2.0, tautp = 1.0, ntt = 3, gamma_ln = 2.0, ntb = 2,\n' \
    'ioutfm = 1,ntxo = 2,iwrap = 0,ig = -1,\n' \
    'cut = 9, barostat = 2,\n' \
    'mcwat = %d, mcresstr = \"WAT\", nmc = %d, nmd = %d,\n' \
    'nmropt = %d,\n' \
    '\n' \
    'gti_add_sc = 1,\n' \
    'gti_lam_sch = 1,\n' \
    'tishake = 2, gti_syn_mass = 1,\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.2, scbeta = 50.0,\n' \
    'logdvdl = 1,\n' \
    'ifmbar = 1,\n' \
    'mbar_states = %d,\n' \
    'mbar_lambda = %s,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
    '/\n' \
    '&ewald\n' \
    '/\n' \
    '&wt type=\"REST\", istep1=0,istep2=10000,value1=0.0001,value2=0.1,  /\n' \
    '&wt type=\"REST\", istep1=10000,istep2=%d,value1=0.1,value2=1.0,  /\n' \
    '&wt type=\"END\" /\n' \
    'DISANG=%s\n' \
    '/\n' \
    '\n' \
    % (equil_ns,equil_nstlim,dt,int(1000/(dt*100)),int(1000/(dt*100)),int(1000/(dt*100)),int(5000/(dt*100)),monte_water,nmc,nmd,nmropt,icfe,lambda_val,int(len(all_lambda_list)),all_lambda,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1],equil_nstlim,disang)

    prod='Prod %.3f ns NPT\n' \
    '&cntrl\n' \
    'imin = 0, ntx = 5, irest = 1,\n' \
    'nstlim = %d, dt = %.3f,\n' \
    'ntwe = %d, ntwx = %d, ntpr = %d, ntwr = %d,\n' \
    'ntc = 2, ntf = 1,\n' \
    'temp0 = 300.0,\n' \
    'ntp = 1, pres0 = 1.0, taup = 2.0, tautp = 1.0, ntt = 3, gamma_ln = 2.0, ntb = 2,\n' \
    'ioutfm = 1,ntxo = 2,iwrap = 0,ig = -1,\n' \
    'cut = 9, barostat = 2,\n' \
    'mcwat = %d, mcresstr = \"WAT\", nmc = %d, nmd = %d,\n' \
    'nmropt = %d,\n' \
    '\n' \
    'gti_add_sc = 1,\n' \
    'gti_lam_sch = 1,\n' \
    'gti_scale_beta = 1,\n' \
    'gti_vdw_exp = 2,\n' \
    'gti_ele_exp = 2,\n' \
    'tishake = 2, gti_syn_mass = 1,\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.5, scbeta = 1.0,\n' \
    'logdvdl = 1,\n' \
    'ifmbar = 1,\n' \
    'mbar_states = %d,\n' \
    'mbar_lambda = %s,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
    '/\n' \
    '&ewald\n' \
    '/\n' \
    '&wt\n' \
    'type=\"END\",\n' \
    '&end\n' \
    'DISANG=%s\n' \
    '/\n' \
    '\n' \
    % (prod_ns,prod_nstlim,dt,int(1000/(dt*100)),int(1000/(dt*100)),int(1000/(dt*100)),int(5000/(dt*100)),monte_water,nmc,nmd,nmropt,icfe,lambda_val,int(len(all_lambda_list)),all_lambda,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1],disang)

    # write files
    with open('min.in','w') as f:
        f.write(minimise)

    with open('heat.in','w') as f:
        f.write(heat)

    # short CPU pressure run
    with open('cpu_press.in','w') as f:
        f.write(press.replace('nstlim = 100000','nstlim = 500'))

    with open('press.in','w') as f:
        f.write(press)

    with open('equil.in','w') as f:
        f.write(equil)

    with open('prod.in','w') as f:
        f.write(prod)

def write_cmdline_ti_pair(sdf_file,ff,mcss,align,protocol='one-step',hmass=True):

    with open('local_cmdline_pair.py','w') as f:
        f.write('import sys\n')
        f.write('import timber\n')
        f.write('import pandas as pd\n')
        f.write('from rdkit import Chem\n')
        f.write('\n')
        f.write('sdf_file=\"%s\"\n' % (sdf_file))
        f.write('ff=\"%s\"\n' % (ff))
        if mcss:
            f.write('mcss=\"%s\"\n' % (mcss))
        else:
            f.write('mcss=None\n')
        f.write('align=%s\n' % (align))
        f.write('protocol=\"%s\"\n' % (protocol))
        f.write('hmass=%s\n' % (hmass))
        f.write('\n')
        f.write('suppl=Chem.SDMolSupplier(sdf_file,removeHs=False)\n')
        f.write('all_mols=[m for m in suppl]\n')
        f.write('all_names=[str(m.GetProp(\"_Name\")) for m in all_mols]\n')
        f.write('idx=int(list(sys.argv)[1])\n')
        f.write('df=pd.read_csv(\"arrayfiles/pairname%d.csv\" % (idx), dtype={\"Name1\":\"str\",\"Name2\":\"str\"})\n')
        f.write('name1=str(df.loc[0,\"Name1\"])\n')
        f.write('name2=str(df.loc[0,\"Name2\"])\n')
        f.write('\n')
        f.write('timber.run_rbfe_setup(all_mols[all_names.index(name1)],all_mols[all_names.index(name2)],ff=ff,full_mcs=mcss,align=align)\n')
        f.write('timber.run_build(df,protocol=protocol,hmass=hmass,use_openff=(ff==\"openff\"))\n')
        f.write('\n')

def write_cmdline_ti_pair_array(task_number):

    text1='''#!/bin/bash
#$ -N setup
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=14400,m_mem_free=2G
#$ -j y
'''

    text2='''
module purge
module load Amber/22-AmberTools22-CUDA11
module load OpenEye
module load PythonDS
conda activate /usr/prog/cadd/amber_tools/alchemistry2

export PYTHONPATH=$PYTHONPATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/
export PATH=$PATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/bin/

'''

    with open('tisetup.sh','w') as f:
        f.write(text1)
        f.write('#$ -t 1-%d\n' % (task_number))
        f.write(text2)
        f.write('\n')
        f.write('MY_DIR=${SGE_TASK_ID}\n')
        f.write('\n')
        f.write('python local_cmdline_pair.py $MY_DIR\n')

def write_cmdline_ti_prod(mapfile,protocol,repeats,schedule_json,hmass=True,equil_ns=1,prod_ns=5,monte_wat=0,nodel=False):

    with open('local_cmdline_prod.py','w') as f:
        f.write('import sys\n')
        f.write('import json\n')
        f.write('import timber\n')
        f.write('import pandas as pd\n')
        f.write('\n')
        f.write('df=pd.read_csv(\"%s.csv\")\n' % (mapfile))
        f.write('df=timber.clean_mapping_file(df)\n')
        f.write('protocol=\"%s\"\n' % (protocol))
        f.write('ti_repeats=%d\n' % (repeats))
        f.write('hmass=%s\n' % (hmass))
        f.write('equil_ns=%d\n' % (equil_ns))
        f.write('prod_ns=%d\n' % (prod_ns))
        f.write('monte_wat=%d\n' % (monte_wat))
        f.write('nodel=%s\n' % (nodel))
        f.write('with open(\"%s\",\"r\") as f:\n' % (schedule_json))
        f.write('    schedule=json.load(f)\n')
        f.write('\n')
        f.write('timber.run_prod(df,protocol=protocol,ti_repeats=ti_repeats,schedule=schedule,hmass=hmass,equil_ns=equil_ns,prod_ns=prod_ns,monte_water=monte_wat,nodel=nodel)\n')
        f.write('\n')

def write_cmdline_ti_prod_array(hold_id=None):

    text1='''#!/bin/bash
#$ -N prod_submit 
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=14400,m_mem_free=2G
#$ -j y
'''

    text2='''
module purge
module load Amber/22-AmberTools22-CUDA11
module load OpenEye
module load PythonDS
conda activate /usr/prog/cadd/amber_tools/alchemistry2

export PYTHONPATH=$PYTHONPATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/
export PATH=$PATH:/usr/prog/cadd/amber_tools/timber/versions/0.1/bin/

mv local_cmdline_pair.py ./arrayfiles
mv tisetup.sh ./arrayfiles
mv setup.o* ./arrayfiles

python local_cmdline_prod.py

mv local_cmdline_prod.py ./arrayfiles

'''
    with open('tiprod.sh','w') as f:
        f.write(text1)
        if hold_id:
            f.write('#$ -hold_jid %d\n' % (hold_id))
        f.write(text2)

