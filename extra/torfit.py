#!/usr/bin/env python
import sys
import os
import copy
import numpy as np
import string
from operator import itemgetter
from rdkit import Chem
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

#
# Callum Dickson, NIBR CADD: callum.dickson@novartis.com
#

##############################################################################

def fit_criteria(qm_ene_profile,md_ene_profile):
    output=False

    sum_squares=np.sum(np.square(md_ene_profile-qm_ene_profile))
    # criteria 100 kcal/mol
    if sum_squares>100:
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

def clean_freeze(freeze_list,tor):

    out_list=[]
    for val in freeze_list:
        if (val[1]!=tor[1]) and (val[2]!=tor[2]):
            if (val[1]!=tor[2]) and (val[2]!=tor[1]):
                out_list.append(val)

    return out_list

##############################################################################
## MAIN ##
##############################################################################

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='TorFit: fit torsions to XTB scans')
    parser.add_argument("-i", type=str, help="Input SDF/MOL2 file",required=True)
    parser.add_argument("-p", type=str, help="Input PRMTOP of molecule",required=True)
    parser.add_argument("-substruct", type=str, help="Smiles / Smarts string of substructure to exclude from fitting",required=False)
    parser.add_argument("-sample", help="Flag to sample one extended conformer",required=False,action='store_true')
    parser.add_argument("-ensemble", help="Flag to sample multiple extended conformers",required=False,action='store_true')
    parser.add_argument("-abc", help="Write parameters in ABC format",required=False,action='store_true')
    parser.add_argument("-gbsa",dest="gbsa",help="Use implicit water model",action="store_true",required=False,default=False)
    args=vars(parser.parse_args())

    ##############################################################################
    ## MOL SETUP ##
    ##############################################################################

    #if not distutils.spawn.find_executable('parmed'):
    #    print('Error: parmed not loaded\n')
    #    sys.exit()

    # Load the SDF / MOL2 file
    if check_file(args['i']):

        input_type=args['i'].split('.')[-1]
        if input_type=='sdf':
            mol=Chem.SDMolSupplier(args['i'],removeHs=False)[0]
        elif input_type=='mol2':
            mol=Chem.MolFromMol2File(args['i'],removeHs=False)
        else:
            print('Error: cannot load %s\n' % (args['i']))
            sys.exit()

    else:
        print('Error: cannot find %s\n' % (args['i']))
        sys.exit()

    if not check_file(args['p']):
        print('Error: cannot find %s\n' % (args['p']))
        sys.exit()

    # Keep a copy
    rdmol=copy.deepcopy(mol)

    # Parse the substructure to exclude from fitting
    frag=None
    if args['substruct'] is not None:
        if Chem.MolFromSmiles(args['substruct']):
            frag=Chem.MolFromSmiles(args['substruct'])
        elif Chem.MolFromSmarts(args['substruct']):
            frag=Chem.MolFromSmarts(args['substruct'])
        else:
            print('Error: cannot convert %s to RDKit mol\n' % (args['substruct']))
            sys.exit()

    # Sample extended conformers
    if args['sample'] and args['ensemble']:
        print('Error: select either sample OR ensemble\n')
        sys.exit()
    elif args['sample'] and not args['ensemble']:
        print('\nSelecting one extended conformer for fitting\n')
        conf_list=generate_extended_conformer(mol,1)
    elif args['ensemble'] and not args['sample']:
        print('\nEnsemble of 3 extended conformers for fitting\n')
        conf_list=generate_extended_conformer(mol,3)
    else:
        print('\nUsing input conformer for fitting\n')
        conf_list=[mol]

    # Molecule_ff
    mol_ff=Molecule_ff(name='NUM')

    # output and save all data
    for info in ['Details','Bonds','Angles','Dihedrals']:
        get_info('details.in',info,args['p'])
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/parmed -i details.in > %s_gaff' % (info))
        os.system('rm details.in')

    # get all mol FF data
    mol_ff=save_atom(mol_ff,'Details_gaff',fields=['name','type','charge'])
    mol_ff=save_bond(mol_ff,'Bonds_gaff')
    mol_ff=save_angle(mol_ff,'Angles_gaff')
    mol_ff=save_dihed(mol_ff,'Dihedrals_gaff')
    os.system('rm Details_gaff Bonds_gaff Angles_gaff Dihedrals_gaff')

    # set negative dihedrals
    for dihe in mol_ff.dihedrals:
        dihe.set_negative_period()
        dihe.sort_dihed()

    # write num frcmod
    for at in mol_ff.atoms:
        at.atom_type=at.idx+1
    write_frcmod('num_param.frcmod',types=mol_ff.atoms,bonds=mol_ff.bonds,angles=mol_ff.angles,dihedrals=mol_ff.dihedrals)

    # leaprc file
    leaprc_out('leaprc.num',mol,mol_ff)

    # write a PDB
    write_rd_pdb(mol_ff,mol,mol_ff.name,'NUM.pdb')

    # write NUM.off
    make_off(mol_ff,'make_off.leap')
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_off.leap>out')
    os.system('rm make_off.leap out')

    # build prmtop for amber torsion scans
    build_parm(mol_ff.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
    os.system('rm build_parm.leap out inpcrd')

    ##############################################################################
    ## GET ROTATABLE TORSIONS ##
    ##############################################################################

    # list all torsions present in mol
    torsionList,frag_skip = enumerateSkipTorsions(mol,frag)

    # Get unique torsion list by atom type
    unique_mol_tor,unique_mol_id=unique_tor(torsionList,mol,mol_ff)

    # NB passing unique_mol_tor works since atom types are indices
    per_tor=fix_duplicate(unique_mol_tor,mol)

    # save torsions to fit
    torsions_to_fit=[]
    torsions_zero=[]
    for i in range(0,len(per_tor)):
        for j in range(0,len(per_tor[i])):
            if j==0:
                torsions_to_fit.append(sorted(per_tor[i],key=itemgetter(1),reverse=True)[j][0])
            else:
                torsions_zero.append(sorted(per_tor[i],key=itemgetter(1),reverse=True)[j][0])

    # repeat for fragment torsions not fitted
    frag_tor_skip=[]
    if frag is not None:
        unique_skip,unique_skip_id=unique_tor(frag_skip,mol,mol_ff)
        per_tor_skip=fix_duplicate(unique_skip,mol)
        for i in range(0,len(per_tor_skip)):
            for j in range(0,len(per_tor_skip[i])):
                frag_tor_skip.append(sorted(per_tor_skip[i],key=itemgetter(1),reverse=True)[j][0])

    ##############################################################################
    ## RUN XTB SCANS ## 
    ##############################################################################

    print('Found %d rotatable torsions\n' % (len(torsions_to_fit)))

    for i in range(0,len(conf_list)):
        conf=conf_list[i]

        os.mkdir('CONF_%d' % (i))
        os.chdir('CONF_%d' % (i))

        for tor in torsions_to_fit:
            print('\nXTB scan: %d %d %d %d\n' % (tor[0],tor[1],tor[2],tor[3]))
            freeze_list=list(torsions_to_fit)
            freeze_list.remove(tor)
            freeze_list=freeze_list+frag_tor_skip

            freeze_list=clean_freeze(freeze_list,tor)

            tor_dir=str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])
            os.mkdir(tor_dir)
            os.chdir(tor_dir)

            convert_sdf(conf)
            if args['gbsa']:
                xtb_torsion_scan(conf,'tmpAni.xyz',tor,freeze_list,36,10,'h2o',int(rdmolops.GetFormalCharge(conf)))
            else:
                xtb_torsion_scan(conf,'tmpAni.xyz',tor,freeze_list,36,10,None,int(rdmolops.GetFormalCharge(conf)))

            # convert sdf to mdcrd file
            create_mdcrd(conf,mol_ff,Chem.SDMolSupplier('torsion_mols_xtb.sdf',removeHs=False),'../../num.prmtop','amber_xtb.mdcrd')

            # run AMBER MM scan
            amber_scan('../../num.prmtop','amber_xtb.mdcrd',tor,'amber_xtb-scan.dat',gbsa_flag=args['gbsa'])

            os.chdir('../')
        os .chdir('../')

    ##############################################################################
    ## UPDATE ANGLE PARAMETERS ##
    ##############################################################################

    all_xtb_mols=[]
    for i in range(0,len(conf_list)):
        for tor in torsions_to_fit:
            tor_dir='CONF_'+str(i)+'/'+str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])
            suppl=Chem.SDMolSupplier(tor_dir+'/torsion_mols_xtb.sdf',removeHs=False)
            for mol in suppl:
                if float(mol.GetProp('xtb-kcal'))<2:
                    all_xtb_mols.append(mol)

    for angle in mol_ff.angles:
        a=int(angle.atom1.atom_type)-1
        b=int(angle.atom2.atom_type)-1
        c=int(angle.atom3.atom_type)-1

        angle_deg=[]
        for mol in all_xtb_mols:
            ang_val=rdMolTransforms.GetAngleDeg(mol.GetConformer(),a,b,c)
            angle_deg.append(ang_val)

        angle_avg_deg=round(float(np.average(angle_deg)),3)
        if abs(float(angle.ref_angle)-angle_avg_deg)>2.5:
            print('Setting angle %d %d %d from %lf to %lf\n' % (angle.atom1.atom_type,angle.atom2.atom_type,angle.atom3.atom_type,angle.ref_angle,angle_avg_deg))
            angle.ref_angle=angle_avg_deg

    write_frcmod('num_param.frcmod',types=mol_ff.atoms,bonds=mol_ff.bonds,angles=mol_ff.angles,dihedrals=mol_ff.dihedrals)

    build_parm(mol_ff.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
    os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
    os.system('rm build_parm.leap out inpcrd')

    ##############################################################################
    ## FIT ALL TORSIONS ##
    ##############################################################################

    copy_fitted_mol=[]
    for i in range(0,len(conf_list)):
        keep_fitted_terms=[]
        keep_fitted_torIdx=[]
        fit_run=[]

        print('Conformer %d\n' % (i))
        os.chdir('CONF_%d' % (i))
        os.system('cp ../leaprc.num .')
        os.system('cp ../NUM.off .')
        os.system('cp ../NUM.pdb .')

        for tor in torsions_to_fit:
            tor_dir=str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])
            qm_ene_profile=np.genfromtxt(tor_dir+'/torsion_profile_xtb.dat',skip_header=1,usecols=(1))
            md_ene_profile=np.genfromtxt(tor_dir+'/amber_xtb-scan.dat',skip_header=1,usecols=(1))

            if fit_criteria(qm_ene_profile,md_ene_profile):
                print('\nFitting torsion %d %d %d %d\n' % (tor[0],tor[1],tor[2],tor[3]))

                fit_run.append(tor)
                copy_mol=copy.deepcopy(mol_ff)
                torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(tor[0],tor[1],tor[2],tor[3]))
                copy_mol.dihedrals[torIdx].zero_dihed(torIdx)
                for extra in torsions_zero:
                    if (tor[1]==extra[1]) and (tor[2]==extra[2]):
                        torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(extra[0],extra[1],extra[2],extra[3]))
                        copy_mol.dihedrals[torIdx].zero_dihed(torIdx)
                    elif (tor[2]==extra[1]) and (tor[1]==extra[2]):
                        torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(extra[0],extra[1],extra[2],extra[3]))
                        copy_mol.dihedrals[torIdx].zero_dihed(torIdx)

                write_frcmod('num_param.frcmod',types=copy_mol.atoms,bonds=copy_mol.bonds,angles=copy_mol.angles,dihedrals=copy_mol.dihedrals)

                build_parm(copy_mol.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
                os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
                os.system('rm build_parm.leap out inpcrd')

                os.chdir(tor_dir)

                amber_scan('../num.prmtop','amber_xtb.mdcrd',tor,'amber_zero.dat',gbsa_flag=args['gbsa'])

                xtb_ene=np.loadtxt('torsion_profile_xtb.dat',usecols=(1,),skiprows=1)
                amber_zero=np.loadtxt('amber_zero.dat',usecols=(1,),skiprows=1)
                amber_angles=np.loadtxt('amber_zero.dat',usecols=(0,),skiprows=1)

                fitted=fit_GA_torsion(xtb_ene,amber_zero,amber_angles)
                frc_fitted=[round(float(fitted[8]),3),round(float(fitted[6]),3),round(float(fitted[4]),3),round(float(fitted[2]),3),round(float(fitted[0]),3)]
                phase_fitted=[round(float(fitted[9]),3),round(float(fitted[7]),3),round(float(fitted[5]),3),round(float(fitted[3]),3),round(float(fitted[1]),3)]

                torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(tor[0],tor[1],tor[2],tor[3]))
                copy_mol.dihedrals[torIdx].update_dihed_frc(frc_fitted)
                copy_mol.dihedrals[torIdx].update_dihed_phase(phase_fitted)

                keep_fitted_terms.append(copy_mol.dihedrals[torIdx])
                keep_fitted_torIdx.append(torIdx)

                os.chdir('../')

                # write out the fitted torsions, run a new torsion scan
                write_frcmod('num_param.frcmod',types=copy_mol.atoms,bonds=copy_mol.bonds,angles=copy_mol.angles,dihedrals=copy_mol.dihedrals)

                build_parm(copy_mol.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
                os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
                os.system('rm build_parm.leap out inpcrd')

                os.chdir(tor_dir)

                amber_scan('../num.prmtop','amber_xtb.mdcrd',tor,'amber_fitted.dat',gbsa_flag=args['gbsa'])

                os.chdir('../')

            else:
                print('\nSkipping torsion %d %d %d %d good initial fit\n' % (tor[0],tor[1],tor[2],tor[3])) 

        copy_mol=copy.deepcopy(mol_ff)
        for i in range(0,len(fit_run)):
            tor=fit_run[i]
            tor_dir=str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])
            qm_ene_profile=np.genfromtxt(tor_dir+'/torsion_profile_xtb.dat',skip_header=1,usecols=(1))
            md_orig=np.genfromtxt(tor_dir+'/amber_xtb-scan.dat',skip_header=1,usecols=(1))
            md_fitted=np.genfromtxt(tor_dir+'/amber_fitted.dat',skip_header=1,usecols=(1))
            sum_inital=np.sum(np.square(md_orig-qm_ene_profile))
            sum_fitted=np.sum(np.square(md_fitted-qm_ene_profile))

            if sum_fitted<sum_inital:
                print('\nKeeping fitted term %d %d %d %d: %lf -> %lf\n' % (tor[0],tor[1],tor[2],tor[3],sum_inital,sum_fitted))

                # update the dihedral term
                copy_mol.dihedrals[keep_fitted_torIdx[i]].zero_dihed(keep_fitted_torIdx[i])
                copy_mol.dihedrals[keep_fitted_torIdx[i]].update_dihed_frc(keep_fitted_terms[i].frc)
                copy_mol.dihedrals[keep_fitted_torIdx[i]].update_dihed_period(keep_fitted_terms[i].period)
                copy_mol.dihedrals[keep_fitted_torIdx[i]].update_dihed_phase(keep_fitted_terms[i].phase)

                for extra in torsions_zero:
                    if (tor[1]==extra[1]) and (tor[2]==extra[2]):
                        torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(extra[0],extra[1],extra[2],extra[3]))
                        copy_mol.dihedrals[torIdx].zero_dihed(torIdx)
                    elif (tor[2]==extra[1]) and (tor[1]==extra[2]):
                        torIdx=copy_mol.query_dihedral_idx(copy_mol.query_dihedral_type(extra[0],extra[1],extra[2],extra[3]))
                        copy_mol.dihedrals[torIdx].zero_dihed(torIdx)
            else:
                print('\nDiscarding %d %d %d %d: poor fit\n' % (tor[0],tor[1],tor[2],tor[3]))

        copy_fitted_mol.append(copy_mol)

        write_frcmod('num_param.frcmod',types=copy_mol.atoms,bonds=copy_mol.bonds,angles=copy_mol.angles,dihedrals=copy_mol.dihedrals)

        build_parm(copy_mol.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
        os.system('rm build_parm.leap out inpcrd')

        # CONF dir
        os.chdir('../')
       
    ##############################################################################
    ## KEEP BEST CONF FIT ##
    ##############################################################################

    if len(conf_list)==1:
        fit_mol=copy_fitted_mol[0]
        write_frcmod('num_param.frcmod',types=fit_mol.atoms,bonds=fit_mol.bonds,angles=fit_mol.angles,dihedrals=fit_mol.dihedrals)

        build_parm(copy_mol.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
        os.system('rm out')

        print('\nTorFit complete\n')

    else:
        print('\nEvaluating conformer fits across all scans\n')

        for i in range(0,len(conf_list)):
            os.chdir('CONF_%d' % (i))

            for tor in torsions_to_fit:
                tor_dir=str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])
                os.chdir(tor_dir)

                for j in range(0,len(conf_list)):
                    print('\nChecking conf_%d fit on conf_%d scan: torsion %d %d %d %d\n' % (j,i,tor[0],tor[1],tor[2],tor[3]))
                    amber_scan('../../CONF_'+str(j)+'/num.prmtop','amber_xtb.mdcrd',tor,'amber_final-conf'+str(j)+'.dat',gbsa_flag=args['gbsa'])

                os.chdir('../')

            os.chdir('../')

        conf_ene_vals=[]
        for i in range(0,len(conf_list)):
            run_ene=[]
            for j in range(0,len(conf_list)):

                for tor in torsions_to_fit:
                    tor_dir='CONF_'+str(j)+'/'+str(tor[0])+'-'+str(tor[1])+'-'+str(tor[2])+'-'+str(tor[3])

                    qm_ene_profile=np.genfromtxt(tor_dir+'/torsion_profile_xtb.dat',skip_header=1,usecols=(1))
                    amber_ene_profile=np.genfromtxt(tor_dir+'/amber_final-conf'+str(i)+'.dat',skip_header=1,usecols=(1))
                    sum_fitted=np.sum(np.square(amber_ene_profile-qm_ene_profile))
                    run_ene.append(sum_fitted)

            conf_ene_vals.append(np.average(run_ene))
            print('Conformer %d avg fit %lf\n' % (i,np.average(run_ene)))

        low_ene_fit=np.argmin(conf_ene_vals)

        print('\nLowest energy fit for conformer %d\n' % (low_ene_fit))

        fit_mol=copy_fitted_mol[low_ene_fit]

        write_frcmod('num_param.frcmod',types=fit_mol.atoms,bonds=fit_mol.bonds,angles=fit_mol.angles,dihedrals=fit_mol.dihedrals) 

        build_parm(fit_mol.name,'num','build_parm.leap',prmtop_name='num',frcmod_file='num_param.frcmod')
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f build_parm.leap>out')
        os.system('rm out')

    ##############################################################################
    ## WRITE ABC ##
    ##############################################################################

    if args['abc']:
        print('\nWriting ABC parameters\n')

        # make a dict - up to 676 atom types as ABC
        abc_list=[]
        for i in range(25,-1,-1):
            for j in range(0,26):
                abc_list.append('%s%s' % (string.ascii_lowercase[i],string.ascii_lowercase[j]))

        abc_dict={}
        for i in range(0,len(abc_list)):
            abc_dict.update({i:abc_list[i]})

        fit_mol.name='ABC'
    
        # modify atom type to abc
        for at in fit_mol.atoms:
            at.atom_type=abc_dict[at.idx]

        # leaprc file
        leaprc_out('leaprc.abc',mol,fit_mol)

        # write a PDB
        # NB I have used mol to iterate, so use original copy rdmol
        write_rd_pdb(fit_mol,rdmol,fit_mol.name,'ABC.pdb')

        # write the abc frcmod
        write_frcmod('abc_param.frcmod',types=fit_mol.atoms,bonds=fit_mol.bonds,angles=fit_mol.angles,dihedrals=fit_mol.dihedrals)

        # write ABC.off
        make_off(fit_mol,'make_off.leap')
        os.system('/usr/prog/cadd/amber_tools/alchemistry/bin/tleap -f make_off.leap>out')
        os.system('rm make_off.leap out')

        # leap file to build ligand
        #build_parm(fit_mol.name,'abc','build_abc.leap',prmtop_name='abc',frcmod_file='abc_param.frcmod')
        #os.system('tleap -f build_abc.leap>out')
        #os.system('rm out')

