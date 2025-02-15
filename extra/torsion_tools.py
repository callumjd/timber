#!/usr/bin/env python
import sys
import os
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.PropertyMol import *
from rdkit.Chem import SDWriter
from rdkit.Chem import rdmolops

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Mutators
from pyevolve import Initializators
from pyevolve import GAllele

##############################################################################

def convert_sdf(rdmol):
    with open('tmpAni.xyz','w') as f_out:
        f_out.write('%d\n' % (len(rdmol.GetAtoms())))
        f_out.write('comment\n')
        for i in range(0,len(rdmol.GetAtoms())):
            pos=rdmol.GetConformer().GetAtomPosition(i)
            f_out.write('%s %lf %lf %lf\n' % (rdmol.GetAtomWithIdx(i).GetSymbol(),pos.x,pos.y,pos.z))

def xtb_torsion_scan(rd_mol,mol_xyz,torsion_list,freeze_list,n_steps,step_size,gbsa_flag,formal_charge,output_name=None):
    
    assert isinstance(torsion_list,tuple), 'Provide tuple of atom indices!'
    if freeze_list!=None:
        assert isinstance(freeze_list,list), 'Provide list of freeze atom indices!'

    ene_values=[]
    dihed_values=[]
    sd_mols=[]    

    start_angle=rdMolTransforms.GetDihedralDeg(rd_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1)

    with open('cntrl_xtb.in','w') as f:
        f.write('$constrain\n')
        f.write('force constant=5\n')
        f.write('dihedral: %d,%d,%d,%d,auto\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))

        if freeze_list is not None:
            for freeze in freeze_list:
                f.write('dihedral: %d,%d,%d,%d,auto\n' % (freeze[0],freeze[1],freeze[2],freeze[3]))

        f.write('$scan\n')
        f.write('1: %lf,%lf,%d\n' % (start_angle,start_angle+(n_steps*step_size),n_steps))

    if gbsa_flag==None:
        os.system('xtb %s --opt -c %d -I cntrl_xtb.in > xtb_raw_out' % (mol_xyz,formal_charge))
    else:
        os.system('xtb %s --opt -c %d -I cntrl_xtb.in --alpb %s > xtb_raw_out' % (mol_xyz,formal_charge,gbsa_flag))

    # get the energies and convert to kcal/mol
    with open('xtbscan.log','r') as f:
        for line in f:
            if 'SCF done' in line:
                ene_values.append(float(line.split()[2]))
            elif 'energy:' in line:
                ene_values.append(float(line.split()[1]))
        
    zero_point=min(ene_values)
    for i in range(0,len(ene_values)):
        ene_values[i]=(ene_values[i]-zero_point)*627.503

    # convert xyz structures to sdf
    with open('xtbscan.log','r') as f:
        data=f.readlines()

    for i in range(0,n_steps):
        local_mol=Chem.Mol(rd_mol)

        local_data=data[(i*len(rd_mol.GetAtoms()))+2+(i*2):(i*len(rd_mol.GetAtoms()))+2+(i*2)+len(rd_mol.GetAtoms())]

        for j in range(0,len(rd_mol.GetAtoms())):
            x=float(local_data[j].split()[1])
            y=float(local_data[j].split()[2])
            z=float(local_data[j].split()[3])

            local_mol.GetConformer().SetAtomPosition(j,(x,y,z))

        sd_mols.append(local_mol)
        dihed_values.append(rdMolTransforms.GetDihedralDeg(local_mol.GetConformer(),torsion_list[0]-1,torsion_list[1]-1,torsion_list[2]-1,torsion_list[3]-1))

    output_profile=None
    if output_name is not None:
        output_name=str(output_name)
        output_profile=str(output_name)
    else:
        output_name='torsion_mols'
        output_profile='torsion_profile'

    writer=SDWriter(output_name+'_xtb.sdf')
    counter=0
    for mol in sd_mols:
        this_scan=('%d-%d-%d-%d' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))

        pm=PropertyMol(mol)
        pm.SetProp('xtb-kcal',ene_values[counter])
        pm.SetProp('torsion-angle',dihed_values[counter])
        pm.SetProp('torsion-scanned',this_scan)
        writer.write(pm)
        counter+=1
    writer.flush()

    # output the energy profile
    with open(output_profile+'_xtb.dat','w') as f_out:
        f_out.write('# scan results %d %d %d %d\n' % (torsion_list[0],torsion_list[1],torsion_list[2],torsion_list[3]))
        for i in range(0,len(dihed_values)):
            f_out.write('%lf %lf\n' % (float(dihed_values[i]),float(ene_values[i])))

    os.system('rm %s cntrl_xtb.in charges wbo xtbopt.log xtbscan.log xtbopt.xyz xtbrestart xtb_raw_out .xtboptok' % (mol_xyz))

def write_dihe(crdfile,a,b,c,d):
    with open('get_dihe.trajin','w') as f:
        f.write('trajin %s\n' % (crdfile))
        f.write('dihedral @%d @%d @%d @%d out ang_dihe.dat\n' % (a,b,c,d))

def write_min(maxcyc,ncyc,shake,igb):
    with open('min_dihe.in','w') as f:
        f.write('Basic minimization with a restraint\n')
        f.write('\n')
        f.write('&cntrl\n')
        f.write('imin=1,\n')
        f.write('cut=1000,\n')
        f.write('maxcyc=%d,\n' % (maxcyc))
        f.write('ncyc=%d,\n' % (ncyc))
        f.write('igb=%d,\n' % (igb))
        f.write('nmropt=1,\n')
        f.write('ntb=0,\n')
        f.write('ntc=%d,\n' % (shake))
        f.write('ntf=%d,\n' % (shake))
        f.write('/\n')
        f.write('&wt type=\'DUMPFREQ\', istep1=10 /\n')
        f.write('&wt type=\'END\'   /\n')
        f.write('DISANG=dihedral.RST\n')
        f.write('LISTIN=POUT\n')
        f.write('LISTOUT=POUT\n')
        f.write('DUMPAVE=jar.log\n')

def write_RST(ang,a,b,c,d):
    with open('dihedral.RST','w') as f:
        f.write('torsion restraint\n')
        f.write('&rst iat=%d,%d,%d,%d r1=%lf, r2=%lf, r3=%lf, r4=%lf, rk2=5000., rk3=5000., /\n' % (a,b,c,d,ang-100,ang,ang,ang+100))

def write_extract(crdfile,frame):
    with open('extract_start.trajin','w') as f:
        f.write('trajin %s %d %d 1\n' % (crdfile,frame,frame))
        f.write('trajout start.inpcrd restart\n')

def amber_scan(prmtop,crdfile,torsion_list,output_file,gbsa_flag=False):

    a=torsion_list[0]
    b=torsion_list[1]
    c=torsion_list[2]
    d=torsion_list[3]

    maxcyc=1000
    ncyc=50
    shake=1

    if gbsa_flag==True:
        igb=5
    else:
        igb=0

    final_ene=[]

    write_dihe(crdfile,a,b,c,d)
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/cpptraj %s<get_dihe.trajin>out' % (prmtop))
    ang_dihe=np.loadtxt('ang_dihe.dat',usecols=(1,),skiprows=1)
    os.system('rm out ang_dihe.dat get_dihe.trajin')

    # Run scan on each frame of input crdfile
    for i in range(0,np.shape(ang_dihe)[0]):

        # Write the sander input minimization file
        write_min(maxcyc,ncyc,shake,igb)

        # Pull starting frame
        write_extract(crdfile,i+1)
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/cpptraj %s<extract_start.trajin>out' % (prmtop))

        # Write RST file
        write_RST(ang_dihe[i],a,b,c,d)

        # Run sander
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/sander -O -i min_dihe.in -o min_dihe.out -p %s -c start.inpcrd -r min_dihe.rst' % (prmtop))

        # Store result
        dat=[line for line in open('mdinfo') if 'EAMBER' in line]
        if dat:
            final_ene.append(float(dat[0].split()[2]))
        else:
            final_ene.append('empty')

        # Clean up
        os.system('rm min_dihe.rst min_dihe.out start.inpcrd mdinfo dihedral.RST extract_start.trajin out jar.log min_dihe.in')

    # Debase energy results to zero and print out
    min_ene=min(final_ene)

    with open(output_file,'w') as f:
        f.write('# Amber profile %d %d %d %d\n' % (a,b,c,d))
        for i in range(0,np.shape(ang_dihe)[0]):
            if final_ene[i]!='empty':
                f.write('%lf %lf\n' % (ang_dihe[i],float(final_ene[i])-min_ene))

def create_mdcrd(mol,mol_ff,suppl,prmtop,output_file):
    counter=1
    for frame in suppl:
        with open('mol'+str(counter)+'.pdb','w') as f:
            for i in range(0,len(mol.GetAtoms())):
                pos=frame.GetConformer().GetAtomPosition(i)
                f.write('ATOM     %2d  %3s %s   1      %7.3f %7.3f %7.3f\n' % (i+1,mol_ff.atoms[i].name,mol_ff.name,pos.x,pos.y,pos.z))
        counter+=1

    with open('trajin','w') as f:
        for i in range(0,len(suppl)):
            f.write('trajin mol%s.pdb\n' % (i+1))
        f.write('trajout %s mdcrd\n' % (output_file))

    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/cpptraj %s <trajin>out' % (prmtop))
    os.system('rm out mol*pdb trajin')

def fit_GA_torsion(gaussian,sander,angles):

    def eval_func_2graphs(chromosome):
        score=0.0
        rmsd=calc_fitted_2graphs(chromosome[0],chromosome[1],chromosome[2],chromosome[3],chromosome[4],chromosome[5],chromosome[6],chromosome[7],chromosome[8],chromosome[9])
        score=1/rmsd
        return score

    def calc_fitted_2graphs(PKIDVF1,PHASE1,PKIDVF2,PHASE2,PKIDVF3,PHASE3,PKIDVF4,PHASE4,PKIDVF5,PHASE5):
        values=[]
        for i,ANGLE in enumerate(angles):
            EDIHEDRAL1=PKIDVF1*(1+math.cos(math.radians(1*ANGLE-PHASE1)))
            EDIHEDRAL2=PKIDVF2*(1+math.cos(math.radians(2*ANGLE-PHASE2)))
            EDIHEDRAL3=PKIDVF3*(1+math.cos(math.radians(3*ANGLE-PHASE3)))
            EDIHEDRAL4=PKIDVF4*(1+math.cos(math.radians(4*ANGLE-PHASE4)))
            EDIHEDRAL5=PKIDVF5*(1+math.cos(math.radians(5*ANGLE-PHASE5)))
            values.append(EDIHEDRAL1+EDIHEDRAL2+EDIHEDRAL3+EDIHEDRAL4+EDIHEDRAL5+sander[i])
        values=np.array(values)
        values=values-min(values)
        err=[]
        for i,value in enumerate(values):
            err.append(value-gaussian[i])
        err=np.array(err)
        rmsd=np.sqrt(np.sum(err**2)/len(err))
        return rmsd

    #initialize alleles
    setOfAlleles=GAllele.GAlleles()

    #variables
    pkidvf=np.arange(-2,5.1,0.1)
    phase=np.arange(0,360.1,60)

    # 1
    setOfAlleles.add(GAllele.GAlleleList(pkidvf))
    setOfAlleles.add(GAllele.GAlleleList(phase))
    # 2
    setOfAlleles.add(GAllele.GAlleleList(pkidvf))
    setOfAlleles.add(GAllele.GAlleleList(phase))
    # 3
    setOfAlleles.add(GAllele.GAlleleList(pkidvf))
    setOfAlleles.add(GAllele.GAlleleList(phase))
    # 4
    setOfAlleles.add(GAllele.GAlleleList(pkidvf))
    setOfAlleles.add(GAllele.GAlleleList(phase))
    # 5
    setOfAlleles.add(GAllele.GAlleleList(pkidvf))
    setOfAlleles.add(GAllele.GAlleleList(phase))

    #initialize genome with defined alleles
    genome = G1DList.G1DList(10)
    genome.setParams(allele=setOfAlleles)

    #define evaluator function
    genome.evaluator.set(eval_func_2graphs)
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    genome.initializator.set(Initializators.G1DListInitializatorAllele)

    ga = GSimpleGA.GSimpleGA(genome)
    ga.selector.set(Selectors.GRouletteWheel)
    ga.setCrossoverRate(0.85)
    ga.setElitism(True)
    ga.setGenerations(5000)
    ga.evolve(freq_stats=1000)

    best=ga.bestIndividual()

    return best

##############################################################################
