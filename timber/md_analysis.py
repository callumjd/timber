# timber

import os
import sys
import copy
import numpy as np
import pandas as pd
import pytraj as pt
import glob as glob
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D

##############################################################################

def image_traj(pt_traj,ref_pdb=None):
    '''
    image the traj and rmsfit to first frame
    '''
    
    pt_traj=pt_traj.autoimage()
    if ref_pdb:
        # fit to the reference PDB if we can
        try:
            pt_traj=pt_traj.rmsfit(ref=ref_pdb,mask='@CA')
        except:
            pt_traj=pt_traj.rmsfit(ref=0,mask='@CA')
            
    else:
        pt_traj=pt_traj.rmsfit(ref=0,mask='@CA')
        
    return pt_traj

def get_rmsd(pt_traj,mask):
    '''
    compute RMSD, traj will NOT be fitted
    '''
 
    return pt.rmsd(pt_traj,mask=mask,ref=0,nofit=True,update_coordinate=False)

def average_frame(traj,avg_filename='avg.pdb',represent_filename='represent.pdb',lig_mask=':UNL'):
    '''
    make an average and represent PBD
    '''   
 
    # strip all not water
    strip_traj=pt.strip(traj,mask=':WAT,Na+,Cl-,K+,PA,OL,PC,CHL')
    strip_top=pt.strip(traj.top,mask=':WAT,Na+,Cl-,K+,PA,OL,PC,CHL')
 
    # single average frame of the full trajectory
    avg_frame=pt.mean_structure(strip_traj)
    
    # rmsd of the strip_traj to the average frame
    rms_to_avg=pt.rmsd(strip_traj,ref=avg_frame,mask='@CA,C,N,'+lig_mask)
    
    rms_min_index=int(np.argmin(rms_to_avg))
   
    with open('represent_idx.dat','w') as f:
        f.write('frame: %d\n' % (int(rms_min_index+1)))
 
    # represent frame
    if lig_mask!='':
        represent_frame=pt.closest(pt.strip(traj,mask=':Na+,Cl-,K+,PA,OL,PC,CHL'),mask=lig_mask,solvent_mask=':WAT',n_solvents=50,frame_indices=[rms_min_index],dtype='trajectory')
    else:
        represent_frame=pt.Trajectory(xyz=np.reshape(strip_traj[rms_min_index,:],(1,strip_traj.n_atoms,3)),top=strip_traj.top)

    # write the PDB files
    pt.write_traj(represent_filename,traj=represent_frame,frame_indices=[0],overwrite=True)
    
    # for the average structure, make a traj object
    avg_traj=pt.Trajectory(xyz=np.reshape(avg_frame.xyz,(1,strip_traj.n_atoms,3)),top=strip_traj.top)
 
    pt.write_traj(avg_filename,traj=avg_traj,overwrite=True)
   
def sse_to_num(sse):
    num = np.empty(sse.shape, dtype=int)
    num[sse == '0'] = 0
    num[sse == 'b'] = 1
    num[sse == 'B'] = 2
    num[sse == 'S'] = 3
    num[sse == 'T'] = 4
    num[sse == 'H'] = 5
    num[sse == 'G'] = 6
    num[sse == 'I'] = 7
    num[sse == 'C'] = 8
    return num

def dssp_plot(file_str,png_filename='dssp.png'):

    resn=pt.read_pickle(glob.glob(file_str+'resn_dssp.pkl')[0])

    dssp_data=[]
    for f in glob.glob(file_str+'dssp_raw_dssp.pkl'):
        dssp_data.append(pt.read_pickle(f))

    dssp_raw=np.concatenate(dssp_data)

    sse=sse_to_num(dssp_raw)

    color_assign = {
    r"none": "white",
    r"$\beta$-sheet": "red",
    r"$\beta$-bridge": "black",
    r"bend": "green",
    r"turn": "yellow",
    r"$\alpha$-helix": "blue",
    r"$3_{10}$-helix": "gray",
    r"$\pi$-helix": "purple",
    r"coil": "orange",
    }
    cmap = colors.ListedColormap(color_assign.values())

    # Plotting DSSP
    plt.figure(figsize=(8.0, 6.0))
    plt.imshow(sse.T, cmap=cmap, origin='lower')
    plt.xlabel("Time / ns")
    plt.ylabel("Residue")
    plt.xticks(ticks=np.arange(0,len(dssp_raw),50),labels=np.arange(0,len(dssp_raw)/100,0.5))
    plt.yticks(ticks=np.arange(0,len(resn),50),labels=[resn[x] for x in np.arange(0,len(resn),50)])

    # Custom legend below the DSSP plot
    custom_lines = [
    Line2D([0], [0], color=cmap(i), lw=4) for i in range(len(color_assign))
    ]
    plt.legend(
    custom_lines, color_assign.keys(), loc="upper center",
    bbox_to_anchor=(0.5, -0.15), ncol=len(color_assign), fontsize=8)
    plt.tight_layout()

    if png_filename:
        plt.savefig(png_filename,dpi=320,bbox_inches='tight')

    plt.clf()

def csv_to_png(filelist,title='plot'):
    
    df_list=[]
    for f in filelist:
        df=pd.read_csv(f)
        df_list.append(df)
        
    plt.figure()
    for df in df_list:
        plt.plot(list(df[list(df_list[0].columns)[0]]),list(df[list(df_list[0].columns)[1]]))
        
        plt.xlabel(list(df_list[0].columns)[0],fontsize=12)
        plt.ylabel(list(df_list[0].columns)[1],fontsize=12)
        
        if list(df_list[0].columns)[0]!='time_ns':
            plt.xticks(ticks=np.arange(0,len(df),50),labels=np.array(df[list(df_list[0].columns)[0]])[list(np.arange(0,len(df),50))])
        
    plt.title(title)
    plt.legend(['run_'+str(x+1) for x in range(0,len(df_list))])
        
    plt.savefig(title+'.png',dpi=320,bbox_inches='tight')

    plt.clf()
   
def sasa_figure(rdmol,file_str='run_*/sasa.csv',filename=None):

    atomic_sasa=[]

    frames=[]
    for f in glob.glob(file_str):
        local_df=pd.read_csv(f)
        frames.append(local_df)
   
    if len(frames)>0: 
        df=pd.concat(frames)
        for i in range(2,len(df.columns)):
            atomic_sasa.append(np.average(list(df[list(df.columns)[i]])))
    else:
        return None

    # set the matplotlib colormap
    norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(atomic_sasa))
    cmap = cm.cool

    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    # save rgb of the highlight colors
    atoms_to_highlight=[i for i in range(0,len(rdmol.GetAtoms()))]

    highlights={}
    for i in atoms_to_highlight:
        highlights.update({i:m.to_rgba(atomic_sasa[i])})

    rdDepictor.Compute2DCoords(rdmol)
    d2d = Draw.MolDraw2DCairo(450,450)
    d2d.DrawMolecule(rdmol,highlightAtoms=atoms_to_highlight,highlightAtomColors=highlights,highlightBonds=None)
    d2d.FinishDrawing()
    png_data = d2d.GetDrawingText()

    if filename:
        with open(filename, 'wb') as f:
            f.write(png_data)

import glob as glob

def torsion_figures(rdmol,file_str='run_*/torsion_*csv',hist_filename=None,highlight_filename=None):
    
    torsion_dict=dict()

    filelist=glob.glob(file_str)
   
    if len(filelist)==0:
        return None
 
    # get all torsion data into a dictionary
    for f in filelist:
        df=pd.read_csv(f)
        tor_str=list(df.columns)[1]
        
        a=int(tor_str.split('-')[0])
        b=int(tor_str.split('-')[1])
        c=int(tor_str.split('-')[2])
        d=int(tor_str.split('-')[3])
        
        
        if (a,b,c,d) in torsion_dict:
            data=torsion_dict[(a,b,c,d)]
            torsion_dict.update({(a,b,c,d):np.concatenate([data,np.array(list(df[tor_str]))])})
        else:
            torsion_dict.update({(a,b,c,d):np.array(list(df[tor_str]))})
            
    # histogram figure
    if len(torsion_dict)>1:
        fig, axs = plt.subplots(len(torsion_dict))
    else:
        fig, axs = plt.subplots(len(torsion_dict)+1)
    fig.subplots_adjust(hspace=1)
    fig.set_figheight(14)
    fig.set_figwidth(9)

    ctr=0
    for key,value in torsion_dict.items():
        
        hist,edges=np.histogram(value,range=(-180,180),bins=360)
        magnitude=(360-int(hist[np.where(hist==0)].size))/360

        axs[ctr].hist(value,range=(-180,180),color='red',bins=60)
        axs[ctr].set_xlabel('torsion angle (degrees)')
        axs[ctr].set_ylabel('count')
        axs[ctr].set_title('%d-%d-%d-%d' % (key[0],key[1],key[2],key[3]))
        axs[ctr].text(x=-190,y=2,s=str('{:.2f}'.format(magnitude)),fontsize=10)

        ctr+=1
        
    if hist_filename:
        plt.savefig(hist_filename)

    plt.clf()

    # glowing SASA figure
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    cmap = cm.bwr
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    # rdkit image
    hit_bonds=[]
    bond_cols={}

    for key,value in torsion_dict.items():
        
        hist,edges=np.histogram(value,range=(-180,180),bins=360)
        magnitude=(360-int(hist[np.where(hist==0)].size))/360
        b=key[1]
        c=key[2]

        mybond=rdmol.GetBondBetweenAtoms(b,c).GetIdx()
        hit_bonds.append(mybond)
        bond_cols.update({mybond:m.to_rgba(magnitude)})

    rdDepictor.Compute2DCoords(rdmol)

    for idx in range(0,len(rdmol.GetAtoms())):
        rdmol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',str(rdmol.GetAtomWithIdx(idx).GetIdx()))

    d2d = Draw.MolDraw2DCairo(450,450)
    d2d.DrawMolecule(rdmol,highlightAtoms=None,highlightBonds=hit_bonds,highlightBondColors=bond_cols)
    d2d.FinishDrawing()
    png_data = d2d.GetDrawingText()
    if highlight_filename:
        with open(highlight_filename, 'wb') as f:
            f.write(png_data)

def hbond_figure(file_str='run_*/prolif_all.csv',filename=None,cutoff=0.05):

    frames=[]
    for f in glob.glob(file_str):
        local_df=pd.read_csv(f,header=[0,1])
        frames.append(local_df)
   
    if len(frames)>0: 
        df=pd.concat(frames)
    else:
        return None    

    # make dict
    hb_tracker=OrderedDict()
    for res in list(df.columns.levels[0]):
        donor=np.average(df[res]['HBDonor'])
        acc=np.average(df[res]['HBAcceptor'])

        if donor>acc:
            if donor>cutoff:
                hb_tracker.update({res:list(df[res]['HBDonor'])})
        elif acc>donor:
            if acc>cutoff:
                hb_tracker.update({res:list(df[res]['HBAcceptor'])})
                
    if len(hb_tracker)==0:
        return None

    # sort by res
    hb_tracker_rename=OrderedDict()
    res_sort=sorted(hb_tracker, key=lambda x: int(x[3:]))
    for res in res_sort:
        hb_tracker_rename.update({res:hb_tracker[res]})

    plt.bar(range(len(hb_tracker_rename)), [np.average(x) for x in hb_tracker_rename.values()], yerr=[np.std(x) for x in hb_tracker_rename.values()], tick_label=list(hb_tracker_rename.keys()))
    plt.xlabel('residue name')
    plt.ylabel('hbond count')

    if filename:
        plt.savefig(filename)

    plt.clf()

    return hb_tracker_rename

