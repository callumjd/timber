# timber

import math 

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

