#!/usr/bin/env python
import sys
import os
import copy
import numpy as np
import string
import glob as glob
from operator import itemgetter
from rdkit import Chem
from rdkit.Chem import SDWriter
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Descriptors3D
from rdkit.Chem.PropertyMol import *
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import argparse
import distutils.spawn
from ligprep_tools import *
from torsion_tools import *
import logging

# hide the OpenFF warning about OpenEye TK
logging.getLogger().setLevel(logging.ERROR)

from openff.toolkit.topology import Molecule
from openff.fragmenter.fragment import WBOFragmenter,PfizerFragmenter

#
# Split mol into fragments and fit each torsion to XTB profile
#
# https://github.com/openforcefield/openff-fragmenter
#

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

##############################################################################

def fit_criteria(qm_ene_profile,md_ene_profile,criteria=100):
    output=False

    sum_squares=np.sum(np.square(md_ene_profile-qm_ene_profile))
    # criteria kcal/mol
    if sum_squares>criteria:
        output=True
    return output

def generate_extended_conformer(mol,get_extended):
    out_mols=[]
    copy_mol=copy.deepcopy(mol)
    conf_gen=200
    kmeans_clusters=5 # arbitrary. Could set to num rot torsions
    ene_cut=5.0 # prune conformers 5 kcal/mol above lowest energy conformer

    potential = AllChem.ETKDGv3()
    potential.useSmallRingTorsions = True
    potential.useMacrocycleTorsions = True
    AllChem.EmbedMolecule(copy_mol, potential)

    ids = AllChem.EmbedMultipleConfs(copy_mol, numConfs=conf_gen, params=potential)
    min_out=AllChem.MMFFOptimizeMoleculeConfs(copy_mol)

    energies=[]
    for i in range(0,len(min_out)):
        energies.append(float(min_out[i][1]))

    low_ene = copy.deepcopy(copy_mol)
    low_ene.RemoveAllConformers()
    for i in range(0,len(copy_mol.GetConformers())):
        if float(energies[i]-min(energies))<ene_cut:
            conf=copy_mol.GetConformer(i)
            low_ene.AddConformer(conf, assignId=True)

    # just in case we prune more conformers than no. clusters 
    if len(low_ene.GetConformers())<kmeans_clusters:
        kmeans_clusters=len(low_ene.GetConformers())

    all_diheds=enumerateTorsions(low_ene)

    vector_to_cluster=np.zeros((len(all_diheds),conf_gen))

    j=0
    for i in range(0,len(low_ene.GetConformers())):
            mol_dihed=np.zeros(len(all_diheds))
            counter=0
            for val in all_diheds:
                mol_dihed[counter]=float(rdMolTransforms.GetDihedralRad(low_ene.GetConformer(i),val[0],val[1],val[2],val[3]))
                counter+=1
            vector_to_cluster[:,j]=mol_dihed
            j+=1

    X = StandardScaler().fit_transform(vector_to_cluster.T)
    kmeans = KMeans(n_clusters=kmeans_clusters).fit(X)

    cluster_mols=[]
    for i in range(0,kmeans_clusters):
        for j in range(0,len(low_ene.GetConformers())):
            if kmeans.labels_[j]==i:

                local_copy=copy.deepcopy(mol)
                for at in range(0,len(local_copy.GetAtoms())):
                    pos = low_ene.GetConformer(j).GetAtomPosition(at)
                    local_copy.GetConformer().SetAtomPosition(at,pos)

                cluster_mols.append(local_copy)
                break

    pm1=0
    pm2=0
    pm3=0

    pm_prop={'0':'PMI1','1':'PMI2','2':'PMI3'}

    inertia_pm=[]
    for clus_mol in cluster_mols:
        pm=PropertyMol(Chem.Mol(clus_mol))

        pm.SetProp('PMI1',str(Descriptors3D.PMI1(clus_mol)))
        pm1+=Descriptors3D.PMI1(clus_mol)

        pm.SetProp('PMI2',str(Descriptors3D.PMI2(clus_mol)))
        pm2+=Descriptors3D.PMI2(clus_mol)

        pm.SetProp('PMI3',str(Descriptors3D.PMI3(clus_mol)))
        pm3+=Descriptors3D.PMI3(clus_mol)

        inertia_pm.append(pm)

    PM_avg=[pm1/len(cluster_mols),pm2/len(cluster_mols),pm3/len(cluster_mols)]

    inertia_pm.sort(key=lambda x: float(x.GetProp(pm_prop[str(PM_avg.index(max(PM_avg)))])), reverse=True)

    return inertia_pm[0:get_extended]

def get_dihe_mapping(ref_mol,to_fit,idx_1,idx_2,tol=0.05):
 
    r_x_1=ref_mol.GetConformer().GetAtomPosition(idx_1).x
    r_y_1=ref_mol.GetConformer().GetAtomPosition(idx_1).y
    r_z_1=ref_mol.GetConformer().GetAtomPosition(idx_1).z

    r_x_2=ref_mol.GetConformer().GetAtomPosition(idx_2).x
    r_y_2=ref_mol.GetConformer().GetAtomPosition(idx_2).y
    r_z_2=ref_mol.GetConformer().GetAtomPosition(idx_2).z

    for f_at in to_fit.GetAtoms():
        f_x=to_fit.GetConformer().GetAtomPosition(f_at.GetIdx()).x
        f_y=to_fit.GetConformer().GetAtomPosition(f_at.GetIdx()).y
        f_z=to_fit.GetConformer().GetAtomPosition(f_at.GetIdx()).z

        if abs(r_x_1-f_x)<tol and abs(r_y_1-f_y)<tol and abs(r_z_1-f_z)<tol:
            update_1=f_at.GetIdx()

        if abs(r_x_2-f_x)<tol and abs(r_y_2-f_y)<tol and abs(r_z_2-f_z)<tol:
            update_2=f_at.GetIdx()

    dihe_1=[]
    dihe_1_atnum=[]
    for bond in to_fit.GetAtomWithIdx(update_1).GetBonds():
        a1=bond.GetBeginAtomIdx()
        b2=bond.GetEndAtomIdx()

        if a1!=update_1 and a1!=update_2:
            dihe_1.append(a1)
            #dihe_1_atnum.append(to_fit.GetAtomWithIdx(a1).GetAtomicNum())
            dihe_1_atnum.append(len(to_fit.GetAtomWithIdx(a1).GetBonds()))
        if b2!=update_1 and b2!=update_2:
            dihe_1.append(b2)
            #dihe_1_atnum.append(to_fit.GetAtomWithIdx(b2).GetAtomicNum())
            dihe_1_atnum.append(len(to_fit.GetAtomWithIdx(b2).GetBonds()))

    update_0=dihe_1[dihe_1_atnum.index(max(dihe_1_atnum))]

    dihe_2=[]
    dihe_2_atnum=[]
    for bond in to_fit.GetAtomWithIdx(update_2).GetBonds():
        a1=bond.GetBeginAtomIdx()
        b2=bond.GetEndAtomIdx()

        if a1!=update_1 and a1!=update_2:
            dihe_2.append(a1)
            #dihe_2_atnum.append(to_fit.GetAtomWithIdx(a1).GetAtomicNum())
            dihe_2_atnum.append(len(to_fit.GetAtomWithIdx(a1).GetBonds()))
        if b2!=update_1 and b2!=update_2:
            dihe_2.append(b2)
            #dihe_2_atnum.append(to_fit.GetAtomWithIdx(b2).GetAtomicNum())
            dihe_2_atnum.append(len(to_fit.GetAtomWithIdx(b2).GetBonds()))

    update_3=dihe_2[dihe_2_atnum.index(max(dihe_2_atnum))]

    return (update_0,update_1,update_2,update_3)

def tor_lines(file_in):
    with open(file_in,'r') as f:
        start=0
        end=0
        i=0
        for line in f:
            if 'DIHE' in line:
                start=i+1
            elif 'IMPROPER' in line:
                end=i-1
            i+=1

    with open(file_in,'r') as f:
        data=f.readlines()

    data=data[start:end]

    return data

##############################################################################
## MAIN ##
##############################################################################

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Fragment TorFit: fit torsions to XTB scans of fragmented mols')
    parser.add_argument('-i', type=str, help='Input SDF file for fitting',required=True)
    parser.add_argument('-ff',help='Small molecule force field (default: gaff2)',dest='ff',choices=['gaff','gaff2','frosst'],default='gaff2',required=False)
    parser.add_argument('-wbo',help='Use WBO fragmentation scheme rather than Pfizer (slow)',required=False,action='store_true')
    parser.add_argument('-criteria',help='Sum squares fit criteria to determine torsion fitting (default = 100)',required=False,default=100)
    parser.add_argument('-gbsa',dest='gbsa',help='Use implicit water model',action='store_true',required=False,default=False)
    args=vars(parser.parse_args())

    ##############################################################################
    ## MOL SETUP ##
    ##############################################################################

    #if not distutils.spawn.find_executable('parmed'):
    #    print('Error: parmed not loaded\n')
    #    sys.exit()

    # Load the SDF file
    if check_file(args['i']):

        input_type=args['i'].split('.')[-1]
        if input_type=='sdf':
            mol=Chem.SDMolSupplier(args['i'],removeHs=False)[0]
        else:
            print('Error: cannot load %s\n' % (args['i']))
            sys.exit()

    else:
        print('Error: cannot find %s\n' % (args['i']))
        sys.exit()

    # Make a copy of rdkit mol
    rdmol=copy.deepcopy(mol)

    # Make an openff Molecule
    mol = Molecule(args['i'])

    print('\nGenerating fragments ...\n')

    # Get fragments. WBO is slow! Default to Pfizer scheme
    if args['wbo']:
        frag_engine = WBOFragmenter()
    else:
        frag_engine = PfizerFragmenter()
    result = frag_engine.fragment(mol)

    out_frags=[]
    if len(result.fragments)>0:
        for fragment in result.fragments:
            frag=Chem.MolFromSmiles(fragment.smiles)

            map_to_idx={}
            for atom in frag.GetAtoms():
                map_to_idx.update({int(atom.GetProp('molAtomMapNumber')):atom.GetIdx()})

            idx_1=int(map_to_idx[fragment.bond_indices[0]])
            idx_2=int(map_to_idx[fragment.bond_indices[1]])

            frag.SetProp('dihe',str(idx_1)+'-'+str(idx_2))
            frag=AllChem.AddHs(frag)
            AllChem.EmbedMolecule(frag)
            frag_extend=generate_extended_conformer(frag,1)
            out_frags.append(frag_extend[0])

    # No fragments found - use input mol
    else:

        # Parameters and mol_ff for the input
        os.system('/usr/prog/cadd/amber_tools/scripts/run_ligprep.py -i %s -ff %s' % (args['i'],args['ff']))

        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_lig.leap>out')

        # Molecule_ff
        mol_ff=Molecule_ff(name='UNL')

        # output and save all data
        for info in ['Details','Bonds','Angles','Dihedrals']:
            get_info('details.in',info,'prmtop')
            os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/parmed -i details.in > %s_gaff' % (info))
            os.system('rm details.in')

        # get all mol FF data
        mol_ff=save_atom(mol_ff,'Details_gaff',fields=['name','type','charge'])
        mol_ff=save_bond(mol_ff,'Bonds_gaff')
        mol_ff=save_angle(mol_ff,'Angles_gaff')
        mol_ff=save_dihed(mol_ff,'Dihedrals_gaff')
        os.system('rm Details_gaff Bonds_gaff Angles_gaff Dihedrals_gaff')

        # list all torsions present in mol
        torsionList = enumerateTorsions(rdmol)

        # Get unique torsion list by atom type
        unique_mol_tor,unique_mol_id=unique_tor(torsionList,rdmol,mol_ff)

        torsions_to_fit=unique_mol_id

        for mytor in torsions_to_fit:
            local_frag=copy.deepcopy(rdmol)
            local_frag.SetProp('dihe',str(mytor[1])+'-'+str(mytor[2]))
            local_frag.SetProp('local_tor',str(mytor[0])+'-'+str(mytor[1])+'-'+str(mytor[2])+'-'+str(mytor[3]))

            out_frags.append(local_frag)

    # fit each fragment
    for i in range(0,len(out_frags)):
        os.mkdir('piece_'+str(i))
        os.chdir('piece_'+str(i))

        writer=SDWriter('mol.sdf')
        writer.write(out_frags[i])
        writer.flush()

        os.chdir('../')

    print('Running XTB scans ...\n')

    # run this torsion scans
    was_fit=[]
    for i in range(0,len(out_frags)):
        os.chdir('piece_'+str(i))

        os.system('/usr/prog/cadd/amber_tools/scripts/run_ligprep.py -i mol.sdf -ff %s' % (args['ff']))

        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_lig.leap>out')

        print('\npiece_'+str(i))
        template=Chem.SDMolSupplier('mol.sdf',removeHs=False)[0]
        parm_mol=Chem.SDMolSupplier('UNL.sdf',removeHs=False)[0]

        idx_1=int(template.GetProp('dihe').split('-')[0])
        idx_2=int(template.GetProp('dihe').split('-')[1])

        if 'local_tor' in list(out_frags[i].GetPropNames()):
            tor_str=out_frags[i].GetProp('local_tor')
            tor=[int(tor_str.split('-')[0]),int(tor_str.split('-')[1]),int(tor_str.split('-')[2]),int(tor_str.split('-')[3])]
        else:
            tor=get_dihe_mapping(template,parm_mol,idx_1,idx_2)

        # Molecule_ff
        mol_ff=Molecule_ff(name='UNL')

        # output and save all data
        for info in ['Details','Bonds','Angles','Dihedrals']:
            get_info('details.in',info,'prmtop')
            os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/parmed -i details.in > %s_gaff' % (info))
            os.system('rm details.in')

        # get all mol FF data
        mol_ff=save_atom(mol_ff,'Details_gaff',fields=['name','type','charge'])
        mol_ff=save_bond(mol_ff,'Bonds_gaff')
        mol_ff=save_angle(mol_ff,'Angles_gaff')
        mol_ff=save_dihed(mol_ff,'Dihedrals_gaff')
        os.system('rm Details_gaff Bonds_gaff Angles_gaff Dihedrals_gaff')

        # GET ATOM TYPES OF FIT
        type_1=mol_ff.atoms[tor[0]].atom_type
        type_2=mol_ff.atoms[tor[1]].atom_type
        type_3=mol_ff.atoms[tor[2]].atom_type
        type_4=mol_ff.atoms[tor[3]].atom_type
        print('\nFitting torsion: %s %s %s %s\n' % (type_1,type_2,type_3,type_4))

        torIdx=mol_ff.query_dihedral_idx(mol_ff.query_dihedral_atomIdx(tor[0],tor[1],tor[2],tor[3]))

        # RUN XTB SCAN
        convert_sdf(parm_mol)
        if args['gbsa']:
            xtb_torsion_scan(parm_mol,'tmpAni.xyz',(tor[0]+1,tor[1]+1,tor[2]+1,tor[3]+1),None,36,10,'h2o',int(rdmolops.GetFormalCharge(parm_mol)))
        else:
            xtb_torsion_scan(parm_mol,'tmpAni.xyz',(tor[0]+1,tor[1]+1,tor[2]+1,tor[3]+1),None,36,10,None,int(rdmolops.GetFormalCharge(parm_mol)))

        create_mdcrd(parm_mol,mol_ff,Chem.SDMolSupplier('torsion_mols_xtb.sdf',removeHs=False),'prmtop','amber_xtb.mdcrd')

        # INITIAL MM QUALITY 
        amber_scan('prmtop','amber_xtb.mdcrd',(tor[0]+1,tor[1]+1,tor[2]+1,tor[3]+1),'initial_MM.dat',gbsa_flag=args['gbsa'])

        qm_ene_profile=np.genfromtxt('torsion_profile_xtb.dat',skip_header=1,usecols=(1))
        md_ene_profile=np.genfromtxt('initial_MM.dat',skip_header=1,usecols=(1))
        sum_squares=np.sum(np.square(md_ene_profile-qm_ene_profile))

        if fit_criteria(qm_ene_profile,md_ene_profile,float(args['criteria'])):
            print('\nProceeding with fit: sum squares %lf' % (sum_squares))

            # ZERO SCAN
            mol_ff.dihedrals[torIdx].zero_dihed(torIdx)
            write_frcmod('set_zero.frcmod',types=[],bonds=[],angles=[],dihedrals=[mol_ff.dihedrals[torIdx]])

            if len(glob.glob('missing_*frcmod'))>0:
                frcmod_found=glob.glob('missing_*frcmod')[0]
                build_parm('UNL',args['ff'],'make_zero.leap',prmtop_name='zero',frcmod_file=['/usr/prog/cadd/amber_tools/parameters/frmod.add',frcmod_found,'set_zero.frcmod'])
            else:
                build_parm('UNL',args['ff'],'make_zero.leap',prmtop_name='zero',frcmod_file=['/usr/prog/cadd/amber_tools/parameters/frmod.add','set_zero.frcmod'])

            os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_zero.leap>out')
            os.system('rm make_zero.leap out inpcrd set_zero.frcmod')

            amber_scan('zero.prmtop','amber_xtb.mdcrd',(tor[0]+1,tor[1]+1,tor[2]+1,tor[3]+1),'amber_zero.dat',gbsa_flag=args['gbsa'])
            os.system('rm zero.prmtop')

            # PERFORM FITTING
            amber_zero=np.loadtxt('amber_zero.dat',usecols=(1,),skiprows=1)
            amber_angles=np.loadtxt('amber_zero.dat',usecols=(0,),skiprows=1)
            os.system('rm amber_zero.dat')

            fitted=fit_GA_torsion(qm_ene_profile,amber_zero,amber_angles)
            frc_fitted=[round(float(fitted[8]),3),round(float(fitted[6]),3),round(float(fitted[4]),3),round(float(fitted[2]),3),round(float(fitted[0]),3)]
            phase_fitted=[round(float(fitted[9]),3),round(float(fitted[7]),3),round(float(fitted[5]),3),round(float(fitted[3]),3),round(float(fitted[1]),3)]

            mol_ff.dihedrals[torIdx].update_dihed_frc(frc_fitted)
            mol_ff.dihedrals[torIdx].update_dihed_phase(phase_fitted)
            mol_ff.assign_dihedral_idiv(mol_ff.dihedrals[torIdx])

            write_frcmod('fitted_torsion.frcmod',types=[],bonds=[],angles=[],dihedrals=[mol_ff.dihedrals[torIdx]])

            # CHECK FIT QUALITY
            if len(glob.glob('missing_*frcmod'))>0:
                frcmod_found=glob.glob('missing_*frcmod')[0]
                build_parm('UNL',args['ff'],'make_fitted.leap',prmtop_name='fitted',frcmod_file=['/usr/prog/cadd/amber_tools/parameters/frmod.add',frcmod_found,'fitted_torsion.frcmod'])
            else:
                build_parm('UNL',args['ff'],'make_fitted.leap',prmtop_name='fitted',frcmod_file=['/usr/prog/cadd/amber_tools/parameters/frmod.add','fitted_torsion.frcmod'])

            os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_fitted.leap>out')
            os.system('rm make_fitted.leap out inpcrd leap.log')

            amber_scan('fitted.prmtop','amber_xtb.mdcrd',(tor[0]+1,tor[1]+1,tor[2]+1,tor[3]+1),'amber_fitted.dat',gbsa_flag=args['gbsa'])
            os.system('rm fitted.prmtop')

            amber_fitted=np.loadtxt('amber_fitted.dat',usecols=(1,),skiprows=1)

            fitted_sum_squares=np.sum(np.square(amber_fitted-qm_ene_profile))
            if fitted_sum_squares<sum_squares:
                print('\nKeeping fit: %s %s %s %s\n' % (type_1,type_2,type_3,type_4))
                print('Initial %lf vs Fitted %lf\n' % (sum_squares,fitted_sum_squares))
                was_fit.append(i)

        else:
            print('\nGood initial fit: sum squares %lf' % (sum_squares))

        os.chdir('../')

print('\nFinished\n')
if len(was_fit)>0:
    print('Kept %d fits\n' % (len(was_fit)))
else:
    print('No torsion fits retained\n')

output_tor=[]
for val in was_fit:
    output_tor=output_tor+tor_lines('piece_'+str(val)+'/fitted_torsion.frcmod')

with open('fragment_fitted.frcmod','w') as f:
    f.write('fragment fits to XTB\n')
    f.write('DIHE\n')
    for val in output_tor:
        f.write(val)
    f.write('\n')

