# timber

import math 
from .ligprep_tools import check_file

##############################################################################

def protein_charge(prot):

    # simple determination of protein net charge state
    # also return number of residues

    chg_dict={'ARG':1,'LYS':1,'HIP':1,'GLU':-1,'ASP':-1,'ASH':0,'GLH':0,'LYN':0}

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

