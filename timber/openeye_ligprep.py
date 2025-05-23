# timber

import os
from openeye.oechem import *
from openeye import oequacpac,oeomega

##############################################################################

def check_quacpac():
    return oequacpac.OEQuacPacIsLicensed('python')

def check_omega():
    return OEChemIsLicensed('python')

def openeye_charges(input_file,residue_name='UNL'):
    '''
    Generate clean SDF an PDB files
    Return partial charges
    '''

    if not check_quacpac():
        raise Exception('QuacPac not licensed\n')
    if not check_omega():
         raise Exception('OEChem  not licensed\n')

    # Read molecule
    ifs=oemolistream()

    oemol = OEMol()
    ifs = oemolistream(input_file)
    OEReadMolecule(ifs, oemol)
    oemol.SetTitle('%s' % (residue_name))

    # this re-orders atoms
    OEPDBOrderAtoms(oemol)
    OETriposAtomNames(oemol)

    # Try am1bcc elf10 charges
    charge_method=oequacpac.OEAM1BCCELF10Charges()

    # Charge fitting
    quacpac_status = oequacpac.OEAssignCharges(oemol, charge_method)

    # If this fails, fit regular am1bcc
    if not quacpac_status:
        charge_method=oequacpac.OEAM1BCCCharges()
        quacpac_status = oequacpac.OEAssignCharges(oemol, charge_method)
        if not quacpac_status:
            raise Exception('Error: charge fitting failed\n')

    output_charges=[]
    for oeatom in oemol.GetAtoms():
        charge = oeatom.GetPartialCharge()
        output_charges.append(float(charge))

    # Write SDF file for antechamber
    ofs=oemolostream(residue_name+'.sdf')
    OEWriteMolecule(ofs,oemol)

    return output_charges

