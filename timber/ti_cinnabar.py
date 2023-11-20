# timber

import os
import glob as glob
import copy
import collections
import networkx as nx
import pandas as pd
import numpy as np
from rdkit import Chem
from cinnabar import plotting, stats, femap, plotlying

##############################################################################

def run_cinnabar(df,data,col,units):

    if len(list(df.columns))==2:
        name1_col=list(df.columns)[0]
        name2_col=list(df.columns)[1]
        all_ligs=list(df[name1_col])+list(df[name2_col])
        rbfe=True
    elif len(list(df.columns))==1:
        name1_col=list(df.columns)[0]
        all_ligs=list(df[name1_col])
        rbfe=False
    else:
        raise Exception('Error: cannot parse mapping dataframe\n')

    if data:
        exp=pd.read_csv(data)
        name=list(exp.columns)[0]
        crossover=exp[exp[name].isin(all_ligs)]
        if len(crossover)==0:
            raise Exception('Error: data file supplied does not match any ligand transforms!\n')

        exp=exp.dropna(subset=[list(exp.columns)[col-1]])
        exp=exp.reset_index(drop=True)

    if rbfe:
        all_TI_data=gather_rbfe_data()
    else:
        all_TI_data=gather_abfe_data()

    if len(all_TI_data)==0:
        raise Exception('Error: cannot run cinnabar analysis\n')

    # prepare for cinnabar
    exp_mapping={}
    if data:
        for index,row in exp.iterrows():
            name_val=row[name]
            exp_val=float(row[list(exp.columns)[col-1]])

            dG=get_free_energy(exp_val,units)
            exp_mapping.update({name_val:dG})

    # run cinnabar RBFE
    if rbfe:

        make_arsenic_csv(exp_mapping,all_TI_data,all_ligs)

        # run cinnabar
        fe = femap.FEMap('data_arsenic.csv')
        fe.generate_absolute_values()

        print('Writing transformation pair data: pairs_amberTI_DDG.csv\n')
        print('Writing per-compound data: compounds_amberTI_DDG.csv\n')

        do_cycle_closure(all_TI_data,'cycle_closure_amberTI.csv')

        DG_plot_shift=output_rbfe_csv(exp_mapping,all_TI_data,fe,'pairs_amberTI_DDG.csv','compounds_amberTI_DG.csv')

        # remove fe nodes if exp data unknown
        to_remove=[]
        for x in fe.graph.nodes(data=True):
            exp_val=x[1]['exp_DG']
            if exp_val==0.0:
                to_remove.append(x[0])

        for val in to_remove:
            fe.graph.remove_node(val)

        if len(fe.graph.nodes)>1:
            print('Writing DDG and DG figures\n')

            plotting.plot_DDGs(fe.graph,filename='DDGs_fig.png',title='Amber TI',plotly=False)
            plotting.plot_DGs(fe.graph,filename='DGs_fig.png',shift=DG_plot_shift,title='Amber TI',guidelines=True,plotly=False)

    # ABFE 
    else:
        print('Writing per-compound data: compounds_amberTI_DDG.csv\n')
        output_abfe_csv(exp_mapping,all_TI_data,'compounds_amberTI_DG.csv')

    # output statistics
    compounds=pd.read_csv('compounds_amberTI_DG.csv')
    if check_file('pairs_amberTI_DDG.csv'):
        pairs=pd.read_csv('pairs_amberTI_DDG.csv')

    if rbfe:
        fe_statistics('stats_amberTI.csv',compounds,pairs,fe)
    else:
        fe_statistics('stats_amberTI.csv',compounds)

def fe_statistics(stats_csv,compounds,pairs=None,fe=None):

    df=pd.DataFrame(columns=['N_Lig','N_Pert','R2','KTAU','DG_MUE','DG_RMSE','DDG_MUE','DDG_RMSE'])

    df.at[0,'N_Lig']=len(compounds)
    exp_err=[0 for x in range(0,len(compounds))]

    for statistic in ['R2','KTAU','MUE','RMSE']:
        if fe:
            s=stats.bootstrap_statistic(list(compounds['exp_DG_kcal_mol']),list(compounds['RBFE_calc_DG_kcal_mol']),exp_error,list(compounds['RBFE_calc_dDG_kcal_mol']))
        else:
            s=stats.bootstrap_statistic(list(compounds['exp_DG_kcal_mol']),list(compounds['ABFE_calc_DG_kcal_mol']),exp_error,list(compounds['ABFE_calc_dDG_kcal_mol']))

        string=f"{s['mle']:.2f} [95% {s['low']:.2f},{s['high']:.2f}]"
        if statistic in ['R2','KTAU']:
            df.at[0,statistic]=string
        else:
            df.at[0,'DG_'+statistic]=string
    
    # rbfe
    if fe:
        df.at[0,'N_Pert']=len(pairs)
        exp_err=[0 for x in range(0,len(pairs))]

        for statistic in ['MUE','RMSE']:
             s=stats.bootstrap_statistic(list(pairs['exp_DDG_kcal_mol']),list(pairs['RBFE_calc_DDG_kcal_mol']),exp_error,list(compounds['RBFE_calc_dDDG_kcal_mol']))
             string=f"{s['mle']:.2f} [95% {s['low']:.2f},{s['high']:.2f}]"
             df.at[0,'DDG_'+statistic]=string

    df.to_csv(stats_csv,index=False)

def revert_free_energy(exp_val,units):

    temp=298.15
    R=0.001987204

    if units=='mM':
        unit_conv=1e-3
    elif units=='uM':
        unit_conv=1e-6
    elif units=='nM':
        unit_conv=1e-9
    elif units=='pM':
        unit_conv=1e-12
    elif units=='fM':
        unit_conv=1e-15

    output=np.exp(exp_val/(R*temp))*(1/unit_conv)

    return float(output)

def get_free_energy(exp_val,units):

    temp=298.15
    R=0.001987204

    if units=='mM':
        unit_conv=1e-3
    elif units=='uM':
        unit_conv=1e-6
    elif units=='nM':
        unit_conv=1e-9
    elif units=='pM':
        unit_conv=1e-12
    elif units=='fM':
        unit_conv=1e-15

    output=R*temp*np.log(exp_val*unit_conv)
    return float(output)

def get_diff(val1,val2):
    vals=sorted([val1,val2],reverse=True)

    return abs(float(val1-val2))

def star_count(all_pairs):
    name_count=[]
    for entry in all_pairs:
        name_count.append(entry.split('~')[0])
    return list(collections.Counter(name_count).keys())[0]

def gather_rbfe_data():
    all_pairs=[]
    all_dirs=glob.glob('*~*/')

    # sort by occurrence, in case of a star map
    star_name=star_count(all_dirs)
    all_dirs=sorted(all_dirs,key=lambda x: x.split('~')[0]!=star_name)

    for entry in all_dirs:
        entry=entry.replace('/','')
        lig1=entry.split('~')[0].strip()
        lig2=entry.split('~')[1].strip()

        if (lig1,lig2) not in all_pairs and (lig2,lig1) not in all_pairs:
            all_pairs.append((lig1,lig2))

    all_calc=collections.defaultdict(list)
    for pair in all_pairs:
        all_calc.update({pair:[]})

        lig1=pair[0]
        lig2=pair[1]

        forward=lig1+'~'+lig2
        reverse=lig2+'~'+lig1

        for direction in [forward,reverse]:

            if os.path.exists(direction):
                conv_files=glob.glob(direction+'/conv*dat')

                if len(conv_files)>0:

                    for myfile in conv_files:
                        with open(myfile,'r') as f:
                            data=f.readlines()
                            calc=float(data[-1].split()[1])

                            if direction==reverse:
                                calc=-1*calc

                            all_calc[pair].append(calc)
                else:
                    if pair in all_calc:
                        all_calc.pop(pair)

    return all_calc

def gather_abfe_data():

    all_lig=[]
    all_dirs=glob.glob('./*')
    
    for entry in all_dirs:
        entry=entry.replace('/','')
        lig1=entry.strip()

        if lig1 not in all_lig:
            all_lig.append(lig1)

    all_calc=collections.defaultdict(list)
    for lig in all_lig:
        conv_files=glob.glob(lig+'/conv*absolute*dat')

        if len(conv_files)>0:

            for myfile in conv_files:
                with open(myfile,'r') as f:
                    data=f.readlines()
                    calc=float(data[-1].split()[1])

                all_calc[lig].append(calc)

        else:
            if lig in all_calc:
                all_calc.pop(lig)

    return all_calc

def cycle_ddG(G,path,calc_DDG):

    running_calc=0
    for i in range(0,len(path)-1):
        ddG=float(G[path[i]][path[i+1]][calc_DDG])
        running_calc+=ddG

    return running_calc

def do_cycle_closure(all_calc,cycle_closure_csv):

    # Make a graph of the transforms
    G=nx.DiGraph()

    for key,value in all_calc.items():
        lig1=key[0]
        lig2=key[1]
        DDG=np.average(value)
        dDDG=np.std(value)

        if lig1 not in list(G.nodes()):
            G.add_node(lig1)
        if lig2 not in list(G.nodes()):
            G.add_node(lig2)

        G.add_edge(lig1,lig2,calc_DDG=DDG,calc_dDDG=dDDG)

        # add reverse if it wasn't run
        if (lig2,lig1) not in all_calc.keys():
            G.add_edge(lig2,lig1,calc_DDG=-1*DDG,calc_dDDG=dDDG)

    node_to_cycles = {}
    for source in G.nodes():
        paths = []
        for target in G.neighbors(source):
            paths += [l + [source] for l in list(nx.all_simple_paths(G, source=source, target=target)) if len(l) > 2]
        if len(paths)>0:
            node_to_cycles[source] = paths

    # the shortest cycles for each node
    all_cycle_closure={}
    for k,v in node_to_cycles.items():
        v.sort(key=lambda x: len(x))
        cyc=v[0]
        trim=sorted(cyc[:-1])

        trim_name=str(trim[0])
        for i in range(1,len(trim)):
            trim_name=trim_name+'~'+trim[i]

        if trim_name not in all_cycle_closure:
            all_cycle_closure.update({trim_name:[cyc,cycle_ddG(G,cyc,'calc_DDG')]})

    # final output
    od = collections.OrderedDict()

    for k in sorted(all_cycle_closure, key=len, reverse=False):
        trim_name=str(all_cycle_closure[k][0][0])
        for i in range(1,len(all_cycle_closure[k][0])-1):
            trim_name=trim_name+'~'+all_cycle_closure[k][0][i]

        od.update({trim_name:all_cycle_closure[k][1]})

    if len(od)>0:
        print('Writing cycle closure data: %s\n' % (cycle_closure_csv))

        with open(cycle_closure_csv,'w') as f:
            f.write('Ligand_cycle,Closure_term_kcal_mol\n')
            for k,v in od.items():
                f.write('%s,%lf\n' % (k,float(v)))

def make_arsenic_csv(exp_mapping,all_TI_data,all_ligs):

    with open('data_arsenic.csv','w') as f:
        f.write('# Experimental block\n')
        f.write('# Ligand, expt_DDG, expt_dDDG\n')

        for name in all_ligs:
            if name in exp_mapping:
                dG=float(exp_mapping[name])
                f.write('%s,%lf,0.0\n' % (name,dG))
            else:
                f.write('%s,0.0,0.0\n' % (name))
        f.write('\n')
        f.write('# Calculated block\n')
        f.write('# Ligand1,Ligand2, calc_DDG, calc_dDDG(MBAR),calc_dDDG(additional)\n')
        for key,value in all_TI_data.items():
            f.write('%s,%s,%lf,%lf,0.0\n' % (key[0],key[1],np.average(value),np.std(value)+0.00001)) # bug in arsenic if error is exactly zero ...
        f.write('\n')

def output_rbfe_csv(exp_mapping,all_TI_data,fe,pairs_csv,compounds_csv):

    # get fe with only exp values
    fe_copy=copy.deepcopy(fe)
    to_remove=[]
    for x in fe_copy.graph.nodes(data=True):
        exp_val=x[1]['exp_DG']
        if exp_val==0.0:
            to_remove.append(x[0])

    for val in to_remove:
        fe_copy.graph.remove_node(val)

    msd=[]
    for x in fe_copy.graph.nodes(data=True):
        DG=float(x[1]['calc_DG'])
        exp=float(x[1]['exp_DG'])
        msd.append(DG-exp)

    if len(msd)>0:
        DG_plot_shift=-1*(np.sum(msd)/len(msd))
    else:
        DG_plot_shift=0.0

    # could MolStandardize here if needed for smiles
    smi_dict={}
    for key,value in all_TI_data.items():
        lig1=key[0]
        lig2=key[1]

        pair_dir=lig1+'~'+lig2
        mol1=Chem.SDMolSupplier(pair_dir+'/core/LIG.sdf')[0]
        mol2=Chem.SDMolSupplier(pair_dir+'/sec_lig/MOD.sdf')[0]

        smi_dict.update({lig1:Chem.MolToSmiles(mol1)})
        smi_dict.update({lig2:Chem.MolToSmiles(mol2)})

    # pairs
    rbfe_pairs=pd.DataFrame(columns=['Ligand_1','Smiles_1','Ligand_2','Smiles_2','exp_DDG_kcal_mol','RBFE_calc_DDG_kcal_mol','RBFE_calc_dDDG_kcal_mol'])

    idx=0
    for key,value in all_TI_data.items():
        name1=key[0]
        name2=key[1]
        DDG=float(np.average(value))
        dDDG=float(np.std(value))

        smi1=smi_dict[name1]
        smi2=smi_dict[name2]

        rbfe_pairs.at[idx,'Ligand_1']=name1
        rbfe_pairs.at[idx,'Smiles_1']=smi1
        rbfe_pairs.at[idx,'Ligand_2']=name2
        rbfe_pairs.at[idx,'Smiles_2']=smi2
        rbfe_pairs.at[idx,'RBFE_calc_DDG_kcal_mol']=DDG
        rbfe_pairs.at[idx,'RBFE_calc_dDDG_kcal_mol']=dDDG

        if (name1 in exp_mapping) and (name2 in exp_mapping):
            fe1=exp_mapping[name1]
            fe2=exp_mapping[name2]
            rbfe_pairs.at[idx,'exp_DDG_kcal_mol']=fe2 - fe1
        else:
            rbfe_pairs.at[idx,'exp_DDG_kcal_mol']=None

        idx+=1

    for index,row in rbfe_pairs.iterrows():
        if row['exp_DDG_kcal_mol']:
            rbfe_pairs.at[index,'RBFE_diff_DDG']=get_diff(row['exp_DDG_kcal_mol'],row['RBFE_calc_DDG_kcal_mol'])

    rbfe_pairs=rbfe_pairs.sort_values(by=['RBFE_calc_DDG_kcal_mol'],ascending=[True])

    # write pairs CSV
    rbfe_pairs.to_csv(pairs_csv,index=False)

    # compounds
    rbfe_compounds=pd.DataFrame(columns=['Ligand','Smiles','exp_DG_kcal_mol','RBFE_calc_DG_kcal_mol','RBFE_calc_dDG_kcal_mol'])

    idx=0
    for x in fe.graph.nodes(data=True):
        name=str(x[1]['name'])
        my_DG=float(x[1]['calc_DG'])+DG_plot_shift
        my_dDG=float(x[1]['calc_dDG'])
        if my_dDG<0.00002:
            my_dDG=0.0

        smi=smi_dict[name]

        rbfe_compounds.at[idx,'Ligand']=name
        rbfe_compounds.at[idx,'Smiles']=smi
        if name in exp_mapping:
            rbfe_compounds.at[idx,'exp_DG_kcal_mol']=exp_mapping[name]
        else:
            rbfe_compounds.at[idx,'exp_DG_kcal_mol']=None
        rbfe_compounds.at[idx,'RBFE_calc_DG_kcal_mol']=my_DG
        rbfe_compounds.at[idx,'RBFE_calc_dDG_kcal_mol']=my_dDG
        idx+=1

    for index,row in rbfe_compounds.iterrows():
        if row['exp_DG_kcal_mol']:
            rbfe_compounds.at[index,'RBFE_diff_DG']=get_diff(row['exp_DG_kcal_mol'],row['RBFE_calc_DG_kcal_mol'])

    rbfe_compounds=rbfe_compounds.sort_values(by=['RBFE_calc_DG_kcal_mol'],ascending=[True])

    # write pairs CSV
    rbfe_compounds.to_csv(compounds_csv,index=False)

    return DG_plot_shift

def output_abfe_csv(exp_mapping,all_TI_data,compounds_csv):

    # could MolStandardize here if needed for smiles
    smi_dict={}
    for key,value in all_TI_data.items():
        lig1=key

        pair_dir=lig1
        mol1=Chem.SDMolSupplier(pair_dir+'/core/LIG.sdf')[0]

        smi_dict.update({lig1:Chem.MolToSmiles(mol1)})

    # compounds
    abfe_compounds=pd.DataFrame(columns=['Ligand','Smiles','exp_DG_kcal_mol','ABFE_calc_DG_kcal_mol','ABFE_calc_dDG_kcal_mol'])

    idx=0
    for key,value in all_TI_data.items():
        name=key
        my_DG=float(np.average(value))
        my_dDG=float(np.stdev(value))

        smi=smi_dict[name]

        abfe_compounds.at[idx,'Ligand']=name
        abfe_compounds.at[idx,'Smiles']=smi
        if name in exp_mapping:
            abfe_compounds.at[idx,'exp_DG_kcal_mol']=exp_mapping[name]
        else:
            abfe_compounds.at[idx,'exp_DG_kcal_mol']=None
        abfe_compounds.at[idx,'ABFE_calc_DG_kcal_mol']=my_DG
        abfe_compounds.at[idx,'ABFE_calc_dDG_kcal_mol']=my_dDG
        idx+=1

    for index,row in abfe_compounds.iterrows():
        if row['exp_DG_kcal_mol']:
            abfe_compounds.at[index,'ABFE_diff_DG']=get_diff(row['exp_DG_kcal_mol'],row['ABFE_calc_DG_kcal_mol'])

    abfe_compounds=abfe_compounds.sort_values(by=['ABFE_calc_DG_kcal_mol'],ascending=[True])

    # write pairs CSV
    abfe_compounds.to_csv(compounds_csv,index=False)

