# timber

import parmed
import numpy as np
from rdkit import Chem
from simtk import unit
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from .ligprep_tools import Info_OFF
from .molecule_ff import Molecule_ff,Atom_ff

##############################################################################
#
# tools to replace amber ligand FF with openFF, converting prmtop
# takes pre-calculated partial charges

# see:
# https://github.com/openforcefield/openff-toolkit/blob/master/examples/swap_amber_parameters/swap_existing_ligand_parameters.ipynb

# off_prmtop_converter('prot_UNL','UNL',xml='openff_unconstrained-2.0.0.offxml',gbsa=False)
#
# off_prmtop_converter('complex_ligands',['../core/LIG','../sec_lig/MOD'],xml='openff_unconstrained-2.0.0.offxml')
#
##############################################################################

def make_off_lig(piece,lig_file_path,force_field):

    rdmol=Chem.SDMolSupplier(lig_file_path+'.sdf',removeHs=False,sanitize=True)[0]

    # Make a Molecule_ff to get partial charges
    mol_ff=Molecule_ff(name='tmp')
    n_atoms=len(rdmol.GetAtoms())
    for at in rdmol.GetAtoms():
        x=rdmol.GetConformer().GetAtomPosition(at.GetIdx()).x
        y=rdmol.GetConformer().GetAtomPosition(at.GetIdx()).y
        z=rdmol.GetConformer().GetAtomPosition(at.GetIdx()).z

        mol_ff.add_atom(Atom_ff(idx=at.GetIdx(),atomic_num=at.GetAtomicNum(),atomic_weight=Chem.GetPeriodicTable().GetAtomicWeight(at.GetAtomicNum()),hybrid=at.GetHybridization(),bond_count=len(at.GetBonds()),x=x,y=y,z=z))

    mol_ff=Info_OFF(lig_file_path+'.off',mol_ff,n_atoms,fields=['name','type','charge'])

    # get the partial charges already assigned in .off amber file 
    partial_charges=np.array([float(mol_ff.atoms[i].atom_charge) for i in range(0,len(rdmol.GetAtoms()))])*unit.elementary_charge

    ligand_off_molecule = Molecule(lig_file_path+'.sdf')
    ligand_off_molecule.name = lig_file_path.split('/')[-1]
    ligand_off_molecule.partial_charges = partial_charges

    ligand_off_topology = ligand_off_molecule.to_topology()

    ligand_system = force_field.create_openmm_system(ligand_off_topology,charge_from_molecules=[ligand_off_molecule])
    new_ligand_structure = parmed.openmm.load_topology(ligand_off_topology.to_openmm(),ligand_system,xyz=piece.positions)

    # for some reason, atom names get appended with an 'x'
    # revert to original atom names
    for i in range(0,len(mol_ff.atoms)):
        new_ligand_structure.atoms[i].name=mol_ff.atoms[i].name

    return new_ligand_structure

def off_prmtop_converter(sys_name,lig_file,xml,gbsa=False):

    # lig_file: provide path to .sdf / .off files eg ["core/LIG","sec_lig/MOD"]
    # xml is open force field ligand xml name

    force_field = ForceField(xml)

    off_ligands=[]

    if type(lig_file)==str:
        lig_file=[lig_file]

    in_prmtop = sys_name+'.prmtop'
    in_crd = sys_name+'.inpcrd'
    orig_structure = parmed.amber.AmberParm(in_prmtop, in_crd)

    pieces = orig_structure.split()

    # assemble protein, ligand, spectators 
    complex_structure = parmed.Structure()

    for i in range(0,len(pieces)):
        resn=pieces[i][0].atoms[0].residue.name
        matches=[s for s in lig_file if resn in s]
        if len(matches)>0:
            struc=make_off_lig(pieces[i][0],matches[0],force_field)
            struc *= len(pieces[i][1])
            off_ligands.append(struc) 
        else:
            if len(pieces[i][1])==1:
                struc = pieces[i][0]
            else:
                struc = parmed.Structure()
                struc += pieces[i][0]
                struc *= len(pieces[i][1])

        # AmberParm desired
        # WARNING: with ff19SB, will write ChamberParm regardless (CMAPs)
        # TI is incompatible ChamberParm; use ff14SB
        complex_structure += parmed.amber.AmberParm.from_structure(struc)

    # Copy over the original coordinates and box vectors
    complex_structure.coordinates = orig_structure.coordinates
    if orig_structure.box_vectors:
        complex_structure.box_vectors = orig_structure.box_vectors
    complex_structure.update_dihedral_exclusions()
    for dihedral in complex_structure.dihedral_types:
        if 0.01<float(dihedral.scnb)<1.99:
            dihedral.scee = 1.2
            dihedral.scnb = 2.0

    # overwrite original files ...
    complex_structure.save(sys_name+'.prmtop', overwrite=True)
    complex_structure.save(sys_name+'.inpcrd', overwrite=True)

    # make the PBSA files; force ChamberParm type for ff19SB
    if len(lig_file)==1:
        if gbsa:
            com=parmed.Structure()
            rec=parmed.Structure()

            lig=parmed.amber.ChamberParm.from_structure(off_ligands[0])
        
            for i in range(0,len(pieces)):
                resn=pieces[i][0].atoms[0].residue.name
                matches=[s for s in lig_file if resn in s]
                if len(matches)>0:
                    com+=parmed.amber.ChamberParm.from_structure(off_ligands[0])
                
                elif resn not in ['WAT','Na+','K+','Cl-']:
                    spectator = parmed.Structure()
                    spectator += pieces[i][0]
                    spectator *= len(pieces[i][1])
                        
                    com += parmed.amber.ChamberParm.from_structure(spectator)
                    rec += parmed.amber.ChamberParm.from_structure(spectator)

            # overwrite original files ...
            # this was a pain to force correct PBRadii for MM-PBSA
            com.update_dihedral_exclusions()
            for dihedral in com.dihedral_types:
                if 0.01<float(dihedral.scnb)<1.99:
                    dihedral.scee = 1.2
                    dihedral.scnb = 2.0
            com=parmed.amber.ChamberParm.from_structure(com)
            action = parmed.tools.changeRadii(com,'mbondi2')
            action.execute()
            com.write_parm('com.prmtop')
            com.write_rst7('com.inpcrd')
            
            rec.update_dihedral_exclusions()
            for dihedral in rec.dihedral_types:
                if 0.01<float(dihedral.scnb)<1.99:
                    dihedral.scee = 1.2
                    dihedral.scnb = 2.0
            rec=parmed.amber.ChamberParm.from_structure(rec)
            action = parmed.tools.changeRadii(rec,'mbondi2')
            action.execute()
            rec.write_parm('rec.prmtop')
            rec.write_rst7('rec.inpcrd')

            for dihedral in lig.dihedral_types:
                if 0.01<float(dihedral.scnb)<1.99:
                    dihedral.scee = 1.2
                    dihedral.scnb = 2.0
            action = parmed.tools.changeRadii(lig,'mbondi2')
            action.execute()
            lig.write_parm('lig.prmtop')
            lig.write_rst7('lig.inpcrd')

