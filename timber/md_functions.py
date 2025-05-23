# timber

import os
import sys
import math 
import glob
from rdkit import Chem
from rdkit.Chem import AllChem,SDWriter
from .ligprep_tools import check_file,setup_pdb,setup_hmass
from .openff_ligprep import off_prmtop_converter

##############################################################################

def initial_prot_check(prot):

    with open('check.leap','w') as f:
        f.write('source leaprc.protein.ff19SB\n')
        f.write('source leaprc.phosaa19SB\n')
        f.write('source leaprc.DNA.bsc1\n')
        f.write('source leaprc.water.tip3p\n')
        f.write('loadamberparams frcmod.ionsjc_tip3p\n')
        f.write('mol=loadpdb %s\n' % (prot))
        f.write('saveamberparm mol prmtop inpcrd\n')
        f.write('quit\n')

    os.system('tleap -f check.leap>out')

    status=True

    with open('out','r') as f:
        for line in f:
            if 'Parameter file was not saved' in line:
                status=False
            elif 'Failed to generate parameters' in line:
                status=False
            elif '!FATAL ERROR' in line:
                status=False

    if status:
        os.system('rm check.leap out prmtop inpcrd leap.log\n')

    return status

def build_md_leap(prot,lig=None,sys_name='prot_UNL',prep_files=None,pdb_files=None,pos_ion=0,neg_ion=0,ff='gaff2',protein_ff='ff19SB',water_ff='tip3p',ion_ff='ionsjc_tip3p',gbsa=False,file_name=None):

    convert_prep={'off':'loadoff','lib':'loadoff','prep':'loadamberprep','frcmod':'loadamberparams','mol2':'loadmol2','zinc':'source','add':'loadamberparams'}

    if not check_file(prot):
        raise Exception('Error: cannot find %s\n' % (prot))

    if not file_name:
        file_name='build.leap'

    with open(file_name,'w') as f:
        if ff=='openff':
            if protein_ff=='ff19SB':
                f.write('source leaprc.protein.ff14SB\n')
                f.write('source leaprc.phosaa14SB\n')
            else:
                # TEST required - CMAPs in protein ffs?
                f.write('source leaprc.protein.%s\n' % (protein_ff))
                f.write('source leaprc.phosaa19SB\n')
        else:
            f.write('source leaprc.protein.%s\n' % (protein_ff))
            f.write('source leaprc.phosaa19SB\n')
        f.write('source leaprc.water.%s\n' % (water_ff))
        # DNA force field
        f.write('source leaprc.DNA.bsc1\n')
        # RNA force field
        f.write('source leaprc.RNA.OL3\n')
        f.write('source leaprc.modrna08\n')
        f.write('loadamberparams frcmod.%s\n' % (ion_ff))
        # this is for MG ions etc
        f.write('loadamberparams frcmod.ions234lm_126_tip3p\n')

        if lig:
            if ff=='gaff' or ff=='openff':
                f.write('source leaprc.gaff\n')
                f.write('\n')
                f.write('loadamberparams missing_gaff.frcmod\n')
            elif ff=='gaff2':
                f.write('source leaprc.gaff2\n')
                f.write('\n')
                f.write('loadamberparams missing_gaff2.frcmod\n')

        if prep_files and len(prep_files)>0:
            for val in prep_files:
                my_type=val.split('.')[-1]
                f.write('%s %s\n' % (convert_prep[my_type],os.path.abspath(val)))

        if gbsa:
            f.write('set default PBRadii mbondi2\n')

        if lig:
            f.write('loadoff %s.off\n' % (lig))

        f.write('prot=loadpdb %s\n' % (prot))

        if lig:
            f.write('lig=loadpdb %s.pdb\n' % (lig))

        pdb_str=''
        not_wat=''
        for i in range(0,len(pdb_files)):
            f.write('mol%d=loadpdb %s\n' % (i,pdb_files[i]))
            pdb_str=pdb_str+'mol'+str(i)+' '
            if not is_water(pdb_files[i]):
                not_wat=not_wat+'mol'+str(i)+' '

        if lig:
            f.write('system=combine{prot lig %s}\n' % (pdb_str))
        else:
            f.write('system=combine{prot %s}\n' % (pdb_str))
        f.write('solvatebox system TIP3PBOX 12\n')
        if pos_ion>0 or neg_ion>0:
            f.write('addIonsRand system Na+ %d\n' % (pos_ion))
            f.write('addIonsRand system Cl- %d\n' % (neg_ion))
        else:
            f.write('addIonsRand system Na+ 0\n')
            f.write('addIonsRand system Cl- 0\n')

        f.write('savepdb system %s.pdb\n' % (sys_name))
        f.write('saveamberparm system %s.prmtop %s.inpcrd\n' % (sys_name,sys_name))

        if gbsa:
            f.write('com=combine {prot lig %s}\n'% (not_wat))
            f.write('rec=combine {prot %s}\n'% (not_wat))
            f.write('saveamberparm com com.prmtop com.inpcrd\n')
            f.write('saveamberparm rec rec.prmtop rec.inpcrd\n')
            f.write('saveamberparm lig lig.prmtop lig.inpcrd\n')

        f.write('quit\n')
        f.write('\n')

def protein_charge(prot):

    # simple determination of protein net charge state
    # also return number of residues

    chg_dict={'ARG':1,'LYS':1,'HIP':1,'GLU':-1,'ASP':-1,'ASH':0,'GLH':0,'LYN':0,'TPO':-2}

    chg=0
    res_n=0
    with open(prot,'r') as f:
        for line in f:
            if len(line.split())>0:
                if line.split()[0]=='ATOM' or line.split()[0]=='HETATM':
                    if line[13:15]=='CA':
                        res=line[17:20]
                        res_n+=1
                        if res in chg_dict:
                            chg+=int(chg_dict[res])

    return chg,res_n

def return_salt(nwat,conc,charge):

    # SPLIT method to determine counter-ions
    # Machado, Pantano JCTC 2020, 16, 3
    # https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953

    N0=(nwat*conc)/56.0

    Npos=int(math.ceil(N0-(charge/2)))
    Nneg=int(math.ceil(N0+(charge/2)))

    return Npos,Nneg

def is_water(pdb_file):
    output=False
    with open(pdb_file,'r') as f:
        for line in f:
            if 'HOH' in line:
                 output=True
                 break
            elif 'WAT' in line:
                 output=True
                 break

    return output

def parse_extra(extra,prot):

    # make sure 'pdb' is at the end of this list
    allowed=['frcmod','lib','off','prep','mol2','zinc','add','pdb']

    init_files=[]
    for file_name in extra:
        file_list=glob.glob(file_name)
        if len(file_list)>0:
            for val in file_list:
                my_type=val.split('.')[-1]
                if my_type in allowed:
                    if (check_file(val) and val!=prot):
                        init_files.append(val)
                    else:
                        raise Exception('Error: cannot find file %s\n' % (val))
                else:
                    raise Exception('Error: tleap cannot parse %s file: %s\n' % (my_type,val))
        else:
            raise Exception('Error: cannot find file %s\n' % (file_name))

    prep_files=[]
    pdb_files=[]

    for file_name in init_files:
        my_type=file_name.split('.')[-1]
        if my_type in allowed[0:-1]:
            prep_files.append(file_name)

    for file_name in init_files:
        my_type=file_name.split('.')[-1]
        if my_type=='pdb':
            pdb_files.append(file_name)

    pdb_files.sort(key=lambda x: os.path.getsize(x),reverse=True)

    for val in pdb_files:
        if is_water(val):
            pdb_files.append(pdb_files.pop(pdb_files.index(val)))

    return prep_files,pdb_files

def run_md_build(prot,lig,lig_chg,sys_name,extra,gbsa,ff,file_name=None):

    # parse extra files 
    if extra!=None:
        prep_files,pdb_files=parse_extra(extra,prot)
    else:
        prep_files=[]
        pdb_files=[]

    # protein charge, calculated when no co-factors
    if len(prep_files)==0 and len(pdb_files)==0:
        prot_chg,prot_res=protein_charge(prot)
        pos_ion,neg_ion=return_salt(nwat=prot_res*50,conc=0.15,charge=prot_chg+lig_chg)
    elif len(prep_files)==0 and len(pdb_files)==1 and is_water(pdb_files[0]):
        prot_chg,prot_res=protein_charge(prot)
        pos_ion,neg_ion=return_salt(nwat=prot_res*50,conc=0.15,charge=prot_chg+lig_chg)
    else:
        pos_ion=0
        neg_ion=0

    # build.leap file
    if lig:
        build_md_leap(prot,'UNL',sys_name=sys_name,prep_files=prep_files,pdb_files=pdb_files,pos_ion=pos_ion,neg_ion=neg_ion,gbsa=gbsa,ff=ff,file_name=file_name)
    else:
        build_md_leap(prot,lig=None,sys_name=sys_name,prep_files=prep_files,pdb_files=pdb_files,pos_ion=pos_ion,neg_ion=neg_ion,gbsa=gbsa,ff=ff,file_name=file_name)

def cmdline_leap(prot,sys_name,gbsa,hmass,use_openff):

    os.system('tleap -f build.leap>out')

    if check_file(sys_name+'.prmtop'):
        print('System topology built %s\n' % (sys_name))
    else:
        print('Error: failed to build system topology\n')
        sys.exit()

    # add residue information to the prmtop
    setup_pdb(sys_name+'.prmtop',prot,sys_name+'.resi.prmtop')

    # hmass prmtop
    if hmass:
        setup_hmass(sys_name+'.prmtop')

    # openff parameters
    if use_openff:
        print('Applying openff parameters\n')
        off_prmtop_converter(sys_name,'UNL',xml='openff_unconstrained-2.0.0.offxml',gbsa=gbsa)

def setup_multi_ligand(all_ligands,prot,sys_name,extra,gbsa,ff,build_file=None):

    all_name=[]
    for m in all_ligands:
        name=m.GetProp('_Name')

        if os.path.isdir(name):
            print('Warning: directory %s exists! Skipping.\n' % (name))
        else:
            os.mkdir(name)
            writer=SDWriter(name+'/mol.sdf')
            writer.write(m)
            writer.flush()

            if not build_file:
                run_md_build(prot,'UNL',int(Chem.rdmolops.GetFormalCharge(m)),sys_name,extra,gbsa,ff,file_name=name+'/build.leap')
            else:
                os.system('cp %s %s/build.leap' % (build_file,name))

            all_name.append(name+'/')

    return all_name

def check_ligand_naming(name_list):

    output=True
    for name in name_list:
        if ' ' in name:
            output=False
            break
        elif '~' in name:
            output=False
            break
        elif len(name_list)!=len(set(name_list)):
            output=False
            break

    return output

