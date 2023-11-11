# timber

import os

##############################################################################
# Modify the cluster submit script and MD inputs as needed
##############################################################################

def write_ti_cluster_script(prmtop):

    text1='''#!/bin/bash

#$ -N TI
#$ -cwd
#$ -S /bin/bash
#$ -l cuda=1
#$ -l gpu_card=1
#$ -l h_rt=259100,m_mem_free=2G
#$ -j y
#$ -q default.q

module load Amber

'''

    text2='''
$AMBERHOME/bin/pmemd -O -i min.in -p $prmtop -c $inpcrd -ref $inpcrd -o min.out -e min.en -inf min.info -r min.rst7

$AMBERHOME/bin/pmemd.cuda -O -i heat.in -p $prmtop -c min.rst7 -ref min.rst7 -o heat.out -e heat.en -inf heat.info -r heat.rst7 -x heat.nc

$AMBERHOME/bin/pmemd -O -i cpu_press.in -p $prmtop -c heat.rst7 -ref heat.rst7 -o cpu_press.out -e cpu_press.en -inf cpu_press.info -r cpu_press.rst7 -x cpu_press.nc

$AMBERHOME/bin/pmemd.cuda -O -i press.in -p $prmtop -c cpu_press.rst7 -ref cpu_press.rst7 -o press.out -e press.en -inf press.info -r press.rst7 -x press.nc

$AMBERHOME/bin/pmemd.cuda -O -i equil.in -p $prmtop -c press.rst7 -ref press.rst7 -o equil.out -e equil.en -inf equil.info -r equil.rst7 -x equil.nc

$AMBERHOME/bin/pmemd.cuda -O -i prod.in -p $prmtop -c equil.rst7 -ref equil.rst7 -o prod.out -e prod.en -inf prod.info -r prod.rst7 -x prod.nc
'''

    with open('run_prod.sh','w') as f:
        f.write(text1)
        f.write('export prmtop=%s\n' % (prmtop))
        f.write('export inpcrd=%s\n' % (prmtop.replace('prmtop','inpcrd')))
        f.write(text2)

def write_ti_inputs(protocol,schedule,ti_masks,ti_mask_len,lambda_val,all_lambda_list,hmass=True,equil_ns=1,prod_ns=5,monte_water=0,n_atoms1=0,n_atoms2=0):

    equil_ns=float(equil_ns)
    prod_ns=float(prod_ns)
    
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
        nmropt=1

        if schedule=='complex_restraint':
            ti1=':1'
            ti2=':2'
            charge_mask=''
            ifsc=0
            noshakemask=':1,2'

        elif schedule in ['complex_decharge','solvent_decharge']:
            ti1=':1'
            ti2=':2'
            charge_mask=':2'
            ifsc=0
            noshakemask=':1,2'

        elif schedule in ['complex_ligands','solvent_ligands']:
            if protocol=='absolute':
                ti_masks=[':1',''] # over-ride ti_masks
                ti1=':1'
                ti2=''
                charge_mask=''
                ifsc=1
                noshakemask=':1'
            elif protocol=='absolute-three-step':
                ti_masks=[':1',':2'] # over-ride ti_masks
                ti1=':1'
                ti2=':2'
                charge_mask=':1,2'
                ifsc=1
                noshakemask=':1,2'

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
    'gti_lam_sch = 1,\n' \
    'gti_scale_beta = 1,\n' \
    'gti_vdw_exp = 2,\n' \
    'gti_ele_exp = 2,\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.5, scbeta = 1.0,\n' \
    'logdvdl = 0,\n' \
    'ifsc = %d,\n' \
    'crgmask = \"%s\",\n' \
    'timask1 = \"%s\", timask2 = \"%s\",\n' \
    'scmask1 = \"%s\", scmask2 = \"%s\",\n' \
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
    'gti_lam_sch = 1,\n' \
    'gti_scale_beta = 1,\n' \
    'gti_vdw_exp = 2,\n' \
    'gti_ele_exp = 2,\n' \
    'tishake = 1, noshakemask=\"%s\",\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.5, scbeta = 1.0,\n' \
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
    'DISANG=disang.RST\n' \
    % (noshakemask,icfe,lambda_val,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1])

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
    'gti_lam_sch = 1,\n' \
    'gti_scale_beta = 1,\n' \
    'gti_vdw_exp = 2,\n' \
    'gti_ele_exp = 2,\n' \
    'tishake = 1, noshakemask=\"%s\",\n' \
    'icfe = %d, clambda = %lf, scalpha = 0.5, scbeta = 1.0,\n' \
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
    'DISANG=disang.RST\n' \
    % (nmropt,noshakemask,icfe,lambda_val,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1])

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
    'mcwat = %d, mcresstr = \"WAT\", nmc = 1000, nmd = 1000,\n' \
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
    'DISANG=disang.RST\n' \
    '\n' \
    % (equil_ns,equil_nstlim,dt,int(1000/(dt*1000)),int(1000/(dt*1000)),int(1000/(dt*1000)),int(5000/(dt*1000)),monte_water,nmropt,icfe,lambda_val,int(len(all_lambda_list)),all_lambda,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1])

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
    'mcwat = %d, mcresstr = \"WAT\", nmc = 1000, nmd = 1000,\n' \
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
    'DISANG=disang.RST\n' \
    '\n' \
    % (prod_ns,prod_nstlim,dt,int(1000/(dt*1000)),int(1000/(dt*1000)),int(1000/(dt*1000)),int(5000/(dt*1000)),monte_water,nmropt,icfe,lambda_val,int(len(all_lambda_list)),all_lambda,ifsc,charge_mask,ti1,ti2,ti_masks[0],ti_masks[1])

    # write files
    with open('min.in','w') as f:
        f.write(minimise)

    with open('heat.in','w') as f:
        f.write(heat)

    # short CPU pressure run
    with open('cpu_press.in','w') as f:
        f.write(press.replace('nstlim = 100000','nstlim = 10000'))

    with open('press.in','w') as f:
        f.write(press)

    with open('equil.in','w') as f:
        f.write(equil)

    with open('prod.in','w') as f:
        f.write(prod)

