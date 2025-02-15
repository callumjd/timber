# timber

import os
import glob

##############################################################################
# Modify the cluster submit script and MD inputs as needed
##############################################################################

def write_md_inputs(time,watmd,hmass):

    min1='''minimize
 &cntrl
  imin=1,maxcyc=400,ncyc=50,
  ntb=1,ntp=0,
  ntf=1,ntc=1,
  ntpr=50,
  ntwr=200,
  cut=9.0,
  ntr=1, restraint_wt = 1.00, restraintmask='!:WAT & !:Na+,K+,Cl- & !@H=',
 /
 '''

    min2='''minimize
  &cntrl
   imin=1,maxcyc=10000,ncyc=5000,
   ntb=1,ntp=0,
   ntf=1,ntc=1,
   ntpr=50,
   ntwr=200,
   cut=9.0,
   ntr=1, restraint_wt = 0.10, restraintmask='!:WAT & !:Na+,K+,Cl- & !@H=',
 /
 '''

    min3='''minimize
  &cntrl
   imin=1,maxcyc=10000,ncyc=5000,
   ntb=1,ntp=0,
   ntf=1,ntc=1,
   ntpr=50,
   ntwr=200,
   cut=9.0,
 /
 '''

    heat='''heating 100K
  &cntrl
   imin=0, ntx=1, irest=0,
   ntc=2, ntf=2, tol=0.0000001,
   nstlim=2500, ntt=3, gamma_ln=1.0,
   ig=-1,
   ntpr=100, ntwr=10000,ntwx=100,
   dt=0.002,nmropt=1,
   ntb=1,ntp=0,cut=9.0,ioutfm=1,ntxo=2,
   ntr=1, restraint_wt = 5.00, restraintmask='!:WAT & !:Na+,K+,Cl- & !@H=',
  /
  &wt type='TEMP0', istep1=0, istep2=2500,
                  value1=0.0, value2=100.0 /
  &wt type='END' /
  /
  &ewald
  skinnb=3.0,
 /
 '''

    press='''heating 500 ps 310K
  &cntrl
   imin=0, ntx=5, irest=1,
   ntc=2, ntf=2,tol=0.0000001,
   nstlim=250000, ntt=3, gamma_ln=1.0,
   ig=-1,
   ntpr=100, ntwr=10000,ntwx=100,
   dt=0.002,nmropt=1,
   ntb=2,taup=2.0,cut=9.0,ioutfm=1,ntxo=2,
   ntp=1,
   ntr=1, restraint_wt = 5.00, restraintmask='!:WAT & !:Na+,K+,Cl- & !@H=',
   /
   &wt type='TEMP0', istep1=0, istep2=250000,
                  value1=100.0, value2=310.0 /
   &wt type='END' /
   /
   &ewald
   skinnb=3.0,
 /
 '''

    fix='''pro 500 ps 310K ligand only
  &cntrl
   imin=0, ntx=5, irest=1,
   ntc=2, ntf=2, tol=0.0000001,
   nstlim=250000, ntt=3, gamma_ln=1.0,
   temp0=310.0,
   ntpr=5000, ntwr=500000, ntwx=5000,
   dt=0.002, ig=-1,
   ntb=2, cut=9.0, ioutfm=1, ntxo=2,
   ntp=1, barostat = 1,
   ntr=1, restraint_wt = 2.00, restraintmask=':UNL|@CA',
   /
   &ewald
   skinnb=3.0,
 /
 '''

    watmd_file='''pro 10 ns high frame-rate 
  &cntrl
   imin=0, ntx=5, irest=1,
   ntc=2, ntf=2, tol=0.0000001,
   nstlim=5000000, ntt=3, gamma_ln=1.0,
   temp0=310.0,
   ntpr=125, ntwr=500000, ntwx=125,
   dt=0.002, ig=-1,
   ntb=2, cut=9.0, ioutfm=1, ntxo=2,
   ntp=1, barostat = 2,
   ntr=1, restraint_wt=2.0, restraintmask='@CA'
 /
 '''

    watmd_hmass='''pro 10 ns high frame-rate 
  &cntrl
   imin=0, ntx=5, irest=1,
   ntc=2, ntf=2, tol=0.0000001,
   nstlim=2500000, ntt=3, gamma_ln=1.0,
   temp0=310.0,
   ntpr=62, ntwr=250000, ntwx=62,
   dt=0.004, ig=-1,
   ntb=2, cut=9.0, ioutfm=1, ntxo=2,
   ntp=1, barostat = 2,
   ntr=1, restraint_wt=2.0, restraintmask='@CA'
 /
 '''

    # Write files
    f=open('min1.in','w')
    f.write(min1)
    f.close()

    f=open('min2.in','w')
    f.write(min2)
    f.close()

    f=open('min3.in','w')
    f.write(min3)
    f.close()

    f=open('heat.in','w')
    f.write(heat)
    f.close()

    f=open('press.in','w')
    f.write(press)
    f.close()

    f=open('fix.in','w')
    f.write(fix)
    f.close()

    if watmd:
        f=open('watmd.in','w')
        if hmass:
            f.write(watmd_hmass)
        else:
            f.write(watmd_file)
        f.close()

    with open('prod.in','w') as f:
        f.write('production %.3f ns\n' % (time))
        f.write('&cntrl\n')
        f.write(' imin=0, ntx=5, irest=1,\n')
        f.write(' ntc=2, ntf=2, tol=0.0000001,\n')
        if hmass:
           f.write(' nstlim=%d, ntt=3, gamma_ln=2.0,\n' % (int((time/4)*1e6)))
        else:
            f.write(' nstlim=%d, ntt=3, gamma_ln=2.0,\n' % (int((time/2)*1e6)))
        f.write(' taup=5.0,tautp=2.0,pres0=1.0,\n')
        f.write(' temp0=310.0,\n')
        if hmass:
            f.write(' ntpr=2500, ntwr=250000, ntwx=2500,\n')
            f.write(' dt=0.004, ig=-1,\n')
        else:
            f.write(' ntpr=5000, ntwr=500000, ntwx=5000,\n')
            f.write(' dt=0.002, ig=-1,\n')
        f.write(' ntb=2, cut=9.0, ioutfm=1, ntxo=2,\n')
        f.write(' ntp=1, barostat = 2,\n')
        f.write('/\n')

def write_md_cluster_script(sys_name,watmd,dir_list,task_number,hold_id=None):
    # h_rt=604800 gives 7 days / 168 hrs

    text1='''
#$ -S /bin/bash
#$ -N AmberMD 
#$ -cwd
#$ -l cuda=1
#$ -l gpu_card=1
#$ -l h_rt=604800,m_mem_free=4G
#$ -j y
#$ -q default.q
'''

    text2='''
module unload Amber
module load Amber/22-AmberTools22-CUDA11
source /usr/prog/Amber/22-AmberTools22-CUDA11/amber.sh
'''

    f=open('mdrun.sh','w')
    f.write(text1)
    f.write('#$ -t 1-%d\n' % (task_number))
    if hold_id:
        f.write('#$ -hold_jid %d\n' % (hold_id))
    f.write(text2)
    f.write('dir_list=\'%s\'\n' % os.path.abspath(dir_list))
    f.write('MY_DIR=$( sed -n "${SGE_TASK_ID}p" < $dir_list )\n')
    f.write('\n')
    f.write('cd $MY_DIR\n')
    f.write('\n')
    f.close()

    with open('mdrun.sh','a') as f:
        f.write('export CUDA_VISIBLE_DEVICES=`/cm/shared/apps/nibri/sciComp/get_gpu_map.sh`\n')
        f.write('\n## Minimise ##\n')
        f.write('sander -O -i min1.in -o min1.out -p ../%s.prmtop -c ../%s.inpcrd -r min1.rst -ref ../%s.inpcrd\n' % (sys_name,sys_name,sys_name))
        f.write('\n')
        f.write('pmemd.cuda -O -i min2.in -o min2.out -p ../%s.prmtop -c min1.rst -r min2.rst -ref min1.rst\n' % (sys_name))
        f.write('\n')
        f.write('pmemd.cuda -O -i min3.in -o min3.out -p ../%s.prmtop -c min2.rst -r min3.rst -ref min2.rst\n' % (sys_name))
        f.write('\n')
        f.write('## Equilibration ##\n')
        f.write('pmemd.cuda -O -i heat.in -o heat.out -p ../%s.prmtop -c min3.rst -r heat.rst -x heat.nc -ref min3.rst\n' % (sys_name))
        f.write('\n')
        f.write('pmemd.cuda -O -i press.in -o press.out -p ../%s.prmtop -c heat.rst -r press.rst -x press.nc -ref heat.rst\n' %(sys_name))
        f.write('\n')
        f.write('pmemd.cuda -O -i fix.in -o fix.out -p ../%s.prmtop -c press.rst -r fix.rst -x fix.nc -ref press.rst\n' % (sys_name))
        f.write('\n')
        f.write('## Production ##\n')
        f.write('pmemd.cuda -O -i prod.in -o prod.out -p ../%s.prmtop -c fix.rst -r prod.rst -x prod.nc\n' % (sys_name))
        f.write('\n')
        if watmd:
            f.write('## WATMD ##\n')
            f.write('pmemd.cuda -O -i watmd.in -o watmd.out -p ../%s.prmtop -c prod.rst -r watmd.rst -x watmd.nc -ref prod.rst\n' % (sys_name))
            f.write('\n')
        f.write('cd ../../\n')

def write_pbsa_input(filename='mmpbsa.in',time=2):

    with open(filename,'w') as f:
        f.write('Input file for running PB and GB in serial\n')
        f.write('&general\n')
        f.write(' endframe=%d,interval=2,keep_files=0,strip_mask=":WAT,Na+,K+,Cl-,PA,PC,OL,CHL",\n' % (time*100))
        f.write('/\n')
        f.write('&gb\n')
        f.write(' igb=2, saltcon=0.100,\n')
        f.write('/\n')
        f.write('&pb\n')
        f.write(' istrng=0.100, inp=1, radiopt=0,\n')
        f.write('/\n')

def write_cmdline_leap_bash(prot,sys_name,gbsa,hmass,use_openff,ff):

    with open('local_cmdline_leap.py','w') as f:
        f.write('import os\n')
        f.write('import timber\n')
        f.write('\n')
        f.write('prot=\"%s\"\n' % (prot))
        f.write('sys_name=\"%s\"\n' % (sys_name))
        f.write('gbsa=%s\n' % (gbsa))
        f.write('hmass=%s\n' % (hmass))
        f.write('use_openff=%s\n' % (use_openff))
        f.write('ff=\"%s\"\n' % (ff))
        f.write('\n')
        if use_openff:
            f.write('os.system(\"run_ligprep -i mol.sdf -n UNL -ff gaff\")\n')
        else:
            f.write('os.system(\"run_ligprep -i mol.sdf -n UNL -ff %s\" % (ff))\n')

        f.write('timber.cmdline_leap(prot,sys_name,gbsa,hmass,use_openff)\n')

def write_cmdline_array(dir_list,task_number):

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

    with open('mdsetup.sh','w') as f:
        f.write(text1)
        f.write('#$ -t 1-%d\n' % (task_number))
        f.write(text2)
        f.write('\n')
        f.write('dir_list=\'%s\'\n' % os.path.abspath(dir_list))
        f.write('MY_DIR=$( sed -n "${SGE_TASK_ID}p" < $dir_list )\n')
        f.write('\n')
        f.write('cd $MY_DIR\n')
        f.write('python ../local_cmdline_leap.py\n')
        f.write('cd ../')

def watmd_cpptraj(prot,apo=False):

    with open('wat_image.trajin','w') as f:
        if apo:
            f.write('parm ../prot_APO.prmtop\n')
        else:
            f.write('parm ../prot_UNL.prmtop\n')
        f.write('parm %s\n' % (prot))
        f.write('reference %s parm %s\n' % (prot,prot))
        f.write('loadcrd watmd.nc watmd\n')
        f.write('crdaction watmd autoimage\n')
        f.write('crdaction watmd rms reference mass @C,N\n')
        f.write('crdout watmd wat_image.nc\n')
        f.write('crdaction watmd strip :WAT,Na+,K+,Cl-\n')
        f.write('crdaction watmd average wat_avg.pdb nowrap\n')

def write_analysis_array(dir_list,task_number,nodel,prot=None,hold_id=None,gbsa=False,watmd=False,apo=False):

    text1='''#!/bin/bash
#$ -N analysis 
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=604800,m_mem_free=2G
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
    with open('mdanalysis.sh','w') as f:
        f.write(text1)
        f.write('#$ -t 1-%d\n' % (task_number))
        if hold_id:
            f.write('#$ -hold_jid %d\n' % (hold_id))
        f.write(text2)
        f.write('\n')
        f.write('dir_list=\'%s\'\n' % os.path.abspath(dir_list))
        f.write('MY_DIR=$( sed -n "${SGE_TASK_ID}p" < $dir_list )\n')
        f.write('\n')
        f.write('cd $MY_DIR\n')
        # APO RUN
        if apo:
            if prot:
                f.write('analysis_csv -i prod.nc -p ../prot_APO.resi.prmtop -r %s -lig \"\" -sdf \"\"\n' % (prot))
            else:
                f.write('analysis_csv -i prod.nc -p ../prot_APO.resi.prmtop -lig \"\" -sdf \"\"\n')

        # NOT APO RUN
        else:
            if gbsa:
                if prot:
                    f.write('analysis_csv -i prod.nc -p ../prot_UNL.resi.prmtop -r %s -lig \":UNL\" -sdf ../UNL.sdf -gbsa\n' % (prot))
                else:
                    f.write('analysis_csv -i prod.nc -p ../prot_UNL.resi.prmtop -lig \":UNL\" -sdf ../UNL.sdf -gbsa\n')
            else:
                if prot:
                    f.write('analysis_csv -i prod.nc -p ../prot_UNL.resi.prmtop -r %s -lig \":UNL\" -sdf ../UNL.sdf\n' % (prot))
                else:
                    f.write('analysis_csv -i prod.nc -p ../prot_UNL.resi.prmtop -lig \":UNL\" -sdf ../UNL.sdf\n')
        # NODEL
        if not nodel:
            f.write('rm heat.nc press.nc fix.nc prod.nc\n')

        if watmd:
            f.write('cpptraj -i wat_image.trajin')
            f.write('\n')
            if apo:
                f.write('python /usr/prog/cadd/amber_tools/timber/versions/0.1/extra/cyWATMD/cyWATMD.py -i wat_image.nc -p ../prot_APO.prmtop -nf 40000 -o watmd_out\n')
            else:
                f.write('python /usr/prog/cadd/amber_tools/timber/versions/0.1/extra/cyWATMD/cyWATMD.py -i wat_image.nc -p ../prot_UNL.prmtop -nf 40000 -o watmd_out\n')
            if not nodel:
                f.write('rm watmd.nc wat_image.nc\n')

        f.write('cd ../../\n')

def write_md_collect(hold_id=None):

    text1='''#!/bin/bash
#$ -N final_data
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

collect_md_data

'''
    with open('mdfinal.sh','w') as f:
        f.write(text1)
        if hold_id:
            f.write('#$ -hold_jid %d\n' % (hold_id))
        f.write(text2)

