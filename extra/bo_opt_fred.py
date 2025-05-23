#!/usr/bin/env python
import pandas as pd
import os
import numpy as np
import time
import random
import itertools
import json

from molscoreforge.base import ScoringFunction
from molscoreforge.moo import EvaluationResults
from gats.evaluators import MultiObjectiveEvaluator,Evaluator
from gats.ga import GeneticAlgorithm
from gats.generators import *
from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import AllChem
from IPython.core.display import HTML
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole

from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.base import BaseEstimator
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import DataStructs

from openeye import oechem
from openeye import oeomega
from openeye import oequacpac
from openeye import oeshape
from openeye import oedocking

from joblib import Parallel,delayed

##############################################################################

def predict_posterior(model: RandomForestRegressor, X: np.ndarray):
    """
    Produces a predictive distribution from the trees in a random forest model.

    Arguments:
    - model: sklearn.ensemble.RandomForestRegressor
    - X: Input to the model. Shape: num_samples x num_features

    """
    # Because we want check_input=False below
    X = model._validate_X_predict(X)

    pred_dist = Parallel(n_jobs=-1, verbose=0, backend="threading")(
        delayed(tree.predict)(X, check_input=False) for tree in model.estimators_
    )

    # Each row corresponds to the tree predictions for a sample.
    # Shape: num_samples x num_trees
    pred_dist = np.transpose(np.vstack(pred_dist))

    # Compute the mean and average over the trees.
    mean = np.mean(pred_dist, axis=1)
    std = np.std(pred_dist, axis=1)
    return mean, std

def expected_improvement(
    model: BaseEstimator, X: np.ndarray,  model_name: str, max_observed: float, jitter: float = 0.01
):
    """Compute expected improvement.

    EI attempts to balance exploration and exploitation by accounting
    for the amount of improvement over the best observed value.

    Arguments:
    ----------
    - model: sklearn.ensemble.RandomForestRegressor
    - X: np.array input to the model
    - max_observed: Maximum experimental value observed.

    Returns:
    --------
    - ei: Acquisition values.
    - mean: Mean predictions.
    """
    if model_name == "RF":
        mean, std = predict_posterior(model, X)
        stdev = std + 1e-6
    if model_name == "Gaussian" or model_name == "Gaussian1":
        mean, stdev = model.predict(X, return_std = True)

    # EI parameter values
    z = (mean - max_observed - jitter) / stdev
    imp = mean - max_observed - jitter
    # exploitation + exploration
    acq_values = imp * stats.norm.cdf(z) + stdev * stats.norm.pdf(z)

    acq_values[stdev < jitter] = 0.0

    return acq_values, mean

def upper_confidence_bound(model, X, model_name, beta=1):
    """Compute upper confidence bound
    Arguments:
    ----------
    - model: sklearn.ensemble.RandomForestRegressor
    - X: np.array input to the model
    Returns:
    --------
    - Acquisition values
    - Mean predictions
    """

    # Mean and standard deviation
    if model_name == "RF":
        mean, std = predict_posterior(model, X)
        stdev = std + 1e-6
    if model_name == "Gaussian":
        mean, stdev = model.predict(X, return_std = True)
    acq_values = mean + beta * stdev

    return acq_values, mean

def greedy(model, X):
    """Compute greedy acquisition
    Arguments:
    ----------
    - model: sklearn.ensemble.RandomForestRegressor
    - X: np.array input to the model
    Returns:
    --------
    - Acquisition values
    - Mean predictions
    """

    # Mean and standard deviation
    mean = model.predict(X)
    acq_values = mean

    return acq_values, mean

class Generator(ABC):
    @abstractmethod
    def generate(self, *args, **kwargs) -> list[str]:
        pass

    def to_pickle(self, path: str):
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_pickle(cls, path: str):
        with open(path, "rb") as f:
            return pickle.load(f)

class StaticGenerator(Generator):
    """
    Static generator, which generates new compounds based on a static reaction and a list of reactants
    """

    def generate(self, indices: Sequence[Sequence[int]], *args, **kwargs) -> list[str]:
        """Generates products from a list of indices

        Args:
            indices (Sequence[Sequence[int]]): List of indices for each reactant

        Raises:
            ValueError: if the reaction is not one-to-one

        Returns:
            list[str]: List of SMILES strings for the products
        """
        rxn: ChemicalReaction = AllChem.ReactionFromSmarts(self.reaction)
        products = []
        all_reactants_smiles = []
        for combination in indices:
            reactants,reactants_smiles = self._get_reactants(combination)
            product = rxn.RunReactants(reactants)
            product_list = list(itertools.chain(*product))
            product_smis = [Chem.MolToSmiles(mol) for mol in product_list]
            unique_products = list(set(product_smis))
            if len(unique_products) == 0:
                unique_products = [None]
                
            products.append(unique_products[0])
            all_reactants_smiles.append(reactants_smiles)
        
        df_r=pd.DataFrame(all_reactants_smiles)
        df_p=pd.DataFrame(products,columns=['smiles'])
        
        return df_p.join(df_r)

    def _get_reactants(self, combination: Sequence[int]) -> list[Chem.Mol]:
        reactants = []
        reactants_smiles = []
        for idx, reactant_idx in enumerate(combination):
            reactant_smi = self.reactants[idx][reactant_idx]
            reactants_smiles.append(reactant_smi)
            reactants.append(Chem.MolFromSmiles(reactant_smi))
        return reactants,reactants_smiles

    @property
    @abstractmethod
    def reaction(self) -> str:
        "SMARTS reaction string"
        pass

    @property
    @abstractmethod
    def reactants(self) -> list[list[str]]:
        "List of building blocks for each reactant"
        pass

    @property
    def n_reagents(self) -> int:
        """Number of reactants in the reaction"""
        return len(self.reactants)

    @property
    def n_products(self) -> int:
        """Number of products available to the generation"""
        components = [len(reactant) for reactant in self.reactants]
        return np.prod(components)

class RDGenerator(StaticGenerator):
    def __init__(self, rxn: str, reactants: Sequence[Sequence[str]]):
        self._rxn = rxn
        self._reactans = reactants
        super().__init__(list(chain.from_iterable(reactants)))

    @property
    def reactants(self):
        return self._reactans

    @property
    def reaction(self):
        return self._rxn

def get_rf_models(score_array,fp_array):
    
    output_rf=[]
    
    for i in range(0,len(score_array)):
        
        rf = RandomForestRegressor(n_estimators=100)
        rf.fit(fp_array[i], score_array[i])
        
        output_rf.append(rf)
    
    return output_rf
 
def retain_scores(sample_idx,reagent_list,data_array,m_scores):
    
    # track the scores
    for idx in range(0,len(sample_idx)):
        for x in range(0,len(reagent_list)):
            data_array[x][sample_idx[idx][x]]=np.max([data_array[x][sample_idx[idx][x]],m_scores[idx]])
                
    return data_array

def return_reagent_space(reagent_list):
    
    pairs=[]
    if len(reagent_list)==2:
        for i,j in itertools.product(range(len(reagent_list[0])),range(len(reagent_list[1]))):
            pairs.append([i,j])

    elif len(reagent_list)==3:
        for i,j,k in itertools.product(range(len(reagent_list[0])),range(len(reagent_list[1])),range(len(reagent_list[2]))):
            pairs.append([i,j,k])
    
    else:
        raise Exception('Error: reagent size not supported')
    
    return pairs

def score_from_rank_df(rank_output,reagent_space,reagent_ei):
       
    score_list=[]
    for x in range(0,len(rank_output)):
        score_list.append(np.array(list(rank_output[x]['ei'])))
    
    if len(rank_output)==2:
        for index,(i,j) in enumerate(reagent_space):
            reagent_ei[index]=float(score_list[0][i])+float(score_list[1][j])
            
    elif len(rank_output)==3:
        for index,(i,j,k) in enumerate(reagent_space):
            reagent_ei[index]=float(score_list[0][i])+float(score_list[1][j])+float(score_list[2][k])
    
    return reagent_ei

class FredDock(ScoringFunction):
    def __init__(
        self,
        du_path: os.PathLike,
        desirability=None,
        weight: float = 1.0,
    ):
        name = "FredDock"
        self.du_path=du_path
        super().__init__(
            func=self._score,
            name=name,
            desirability=desirability,
            weight=weight,
            supports_batches=True,
        )

    def _score(self, mol):
        return fred_docking(self.du_path,mol)

    def _evaluate(self, mol):
        return fred_docking(self.du_path,mol)

    def evaluate(self, mols):
        score_list=[self._evaluate(mol) for mol in mols]
        # this is meant to get scaled
        return EvaluationResults(score_list,score_list,score_list,self.name)

def expand_confs(mol, max_confs=1000):
    """
    Set quacpac protonation state
    Expand tautomers
    Expand up to four stereocenters
    Enumerate max_confs omega conformers
    """
    output=[]

    rms = 0.5
    strict_stereo = False
    omega = oeomega.OEOmega()
    omega.SetRMSThreshold(rms)
    omega.SetStrictStereo(strict_stereo)
    omega.SetMaxConfs(max_confs)
    error_level = oechem.OEThrow.GetLevel()
    # Turn off OEChem warnings
    oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Error)

    tautomer_options=oequacpac.OETautomerOptions()
    tautomer_options.SetMaxTautomersGenerated(1000)
    tautomer_options.SetMaxTautomersToReturn(16)
    tautomer_options.SetCarbonHybridization(True)
    tautomer_options.SetMaxZoneSize(50)
    tautomer_options.SetApplyWarts(True)

    # max_centers=4, forceflip=True, enum_nitrogens=True
    for enantiomer in oeomega.OEFlipper(mol.GetActive(), 4, True, True):
        fmol = oechem.OEMol(enantiomer)
        # pKa_norm=True
        for tautomer in oequacpac.OEGetReasonableTautomers(fmol,tautomer_options,True):
            #oequacpac.OEGetReasonableProtomer(fmol)
            status = omega.Build(fmol)
            if status == oeomega.OEOmegaReturnCode_Success:
                output.append(fmol)

    oechem.OEThrow.SetLevel(error_level)
    return output

def read_design_unit(filename):
    """Read an OpenEye design unit
    :param filename: design unit filename (.oedu)
    :return: a docking grid
    """
    du = oechem.OEDesignUnit()
    rfs = oechem.oeifstream()
    if not rfs.open(filename):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % filename)

    du = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(rfs, du):
        oechem.OEThrow.Fatal("Failed to read design unit")
    if not du.HasReceptor():
        oechem.OEThrow.Fatal("Design unit %s does not contain a receptor" % du.GetTitle())
    dock_opts = oedocking.OEDockOptions()
    dock = oedocking.OEDock(dock_opts)
    dock.Initialize(du)
    return dock

def fred_docking(du,mol,max_confs=1000):

    dock = read_design_unit(du)

    smi = Chem.MolToSmiles(mol)
    mc_mol = oechem.OEMol()
    oechem.OEParseSmiles(mc_mol, smi)
    expand_list=expand_confs(mc_mol)
    score = 0.0

    output=[]
    if len(expand_list)>0:
        for m in expand_list:
            docked_mol = oechem.OEGraphMol()
            ret_code = dock.DockMultiConformerMolecule(docked_mol, m)

            if ret_code == oedocking.OEDockingReturnCode_Success:
                dock_opts = oedocking.OEDockOptions()
                sd_tag = oedocking.OEDockMethodGetName(dock_opts.GetScoreMethod())
                # this is a stupid hack, I need to figure out how to do this correctly
                oedocking.OESetSDScore(docked_mol, dock, sd_tag)
                #score = float(oechem.OEGetSDData(docked_mol, sd_tag))
                output.append(docked_mol)

    if len(output)>0:
        output=sorted(output, key=lambda x: float(oechem.OEGetSDData(x, sd_tag)))
        score=float(oechem.OEGetSDData(output[0], sd_tag)) # set to positive for optimizer
    return -1*score

def run_warmup(reagent_list,n_sample):
    
    # Prepare arrays of entire space
    reagent_space=return_reagent_space(reagent_list)
    reagent_ei=[0.0 for x in range(len(reagent_space))]
    
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    
    # Prepare reagent fps
    fp_array=[]
    for i in range(0,len(reagent_list)):
        fp_array.append([mfpgen.GetFingerprint((AllChem.MolFromSmiles(x))) for x in reagent_list[i]])
        
    # Prepare array that keeps the score for each individual reagent, every time it gets evaluated in a mol
    data_array=[]
    for i in range(0,len(reagent_list)):
        data_array.append(np.zeros(len(reagent_list[i])))
    
    # random sample of reagents to combine into n_sample
    sample_idx=random.sample(reagent_space,n_sample)
            
    # create the molecules
    m_df=graph.generate(sample_idx)
        
    # score
    m_scores=evaluator.evaluate([AllChem.MolFromSmiles(smi) for smi in list(m_df['smiles'])]).final_scores
    data_array=retain_scores(sample_idx,reagent_list,data_array,m_scores)
        
    # inspect how many ground truth we have identified
    m_df['score']=m_scores
    
    return data_array,fp_array,m_df,reagent_space,reagent_ei

def run_trial(reagent_list,n_sample,data_array,fp_array,reagent_space,reagent_ei,bo='ucb',n_random=0):

    rf_models=get_rf_models(data_array,fp_array)
        
    # run each RF model on the reagent space
    rank_output=[]
    for x in range(0,len(reagent_list)):

        if bo=='ei':
            ei,ei_m=expected_improvement(rf_models[x],fp_array[x],model_name='RF',max_observed=np.max(data_array[x]))
        elif bo=='ucb':
            ei,ei_m=upper_confidence_bound(rf_models[x],fp_array[x],model_name='RF')
        elif bo=='greedy':
            ei,ei_m=greedy(rf_models[x],fp_array[x])

        rank_df=pd.DataFrame(zip(ei,ei_m),columns=['ei','ei_m'])
        rank_output.append(rank_df)

    # update the reagent_ei selection scores
    reagent_ei=score_from_rank_df(rank_output,reagent_space,reagent_ei)

    # sample N reagents and score
    top_idx=list(np.argsort(reagent_ei)[::-1][:n_sample])
    sample_idx=[reagent_space[x] for x in top_idx]

    # add in random samples. As yes, these do not get popped out
    if n_random>0:
        sample_idx=sample_idx+random.sample(reagent_space,n_random)
    
    # remove from both the reagent space and score space
    for x in sorted(top_idx,reverse=True):
        reagent_space.pop(x)
        reagent_ei.pop(x)

    # create the molecules
    m_df=graph.generate(sample_idx)

    # score
    m_scores=evaluator.evaluate([AllChem.MolFromSmiles(smi) for smi in list(m_df['smiles'])]).final_scores
    data_array=retain_scores(sample_idx,reagent_list,data_array,m_scores)    
        
    # inspect how many ground truth we have identified
    m_df['score']=m_scores
    
    return m_df,top_idx,reagent_space,reagent_ei,data_array,fp_array

def get_HBD(smi):
    return int(rdMolDescriptors.CalcNumHBD(AllChem.MolFromSmiles(smi)))

##############################################################################

# Prepare the evaluator
# Evaluator
evaluator=FredDock(du_path='/mnt/tmplabdata/caddusers/dicksca3/MEMBRANE/TripleAgonist/virtual_screen/bo_run/glp1_local.oedu')

# Prepare the generator
# Reaction
rxn_smarts='([N;D1$(N-[#6]),D2$(N(-[#6])-[#6]);$(N-[#6])!$(N-C=[O,N,S]):1].[C&H3]C(-,:[C&H3])(-,:[C&H3])[#8]-[#6](-[#7:2])=O).[C;D1,$(C[#6]):3](=[OD1:4])[OD1,Cl].[C;D1,$(C[#6]):5](=[OD1:6])[OD1,Cl]>>([N:2][C:5](=[O:6]).[N:1][C:3](=[O:4]))'

# Reagents
reactant_1=pd.DataFrame(['CC(C)(C)OC(=O)N1CC[C@@H](N)C11CC1'],columns=['smiles'])

# load acids no amino
reactant_2=pd.read_csv('/usr/prog/cadd/amber_tools/thompson/RXN_DATA/275592/reagents_1.smi',delim_whitespace=True,names=['smiles','id'],header=None)
reactant_2=reactant_2.drop_duplicates(subset=['smiles'])

# acids with 1 or 2 hbd
reactant_3=pd.read_csv('/usr/prog/cadd/amber_tools/thompson/RXN_DATA/275592/reagents_1.smi',delim_whitespace=True,names=['smiles','id'],header=None)
reactant_3['HBD']=reactant_3['smiles'].apply(get_HBD)
# this count will include the acid -OH
reactant_3=reactant_3[(reactant_3['HBD']>1) & (reactant_3['HBD']<4)]
reactant_3=reactant_3.drop(columns=['HBD'])
reactant_3=reactant_3.drop_duplicates(subset=['smiles'])

# all reagents
reagent_list=[list(reactant_1['smiles']),list(reactant_2['smiles']),list(reactant_3['smiles'])]

# Initiate the rxn graph
graph = RDGenerator(rxn=rxn_smarts,reactants=reagent_list)

# Load the settings
with open('settings.json','r') as f:
    settings=json.load(f)

n_start=settings['n_start']
n_sample=settings['n_sample']
n_random=settings['n_random']
n_iters=settings['n_iters']
bo=settings['bo']

# Start
start=time.time()

# Results dataframe
output=pd.DataFrame(columns=['Iteration','Sampled','Max_score'])

# Collect all molecules sampled and score into list of running_df
running_df=[]

# Warm-up
ctr=0
data_array,fp_array,iter0_df,reagent_space,reagent_ei=run_warmup(reagent_list,n_start)
running_df.append(iter0_df)

output.at[ctr,'Iteration']=-1
output.at[ctr,'Sampled']=len(pd.concat(running_df))
output.at[ctr,'Max_score']=np.max(list(pd.concat(running_df)['score']))
ctr+=1

warmup_time=(time.time()-start)/60
print('Product space: %d' % (len(reagent_ei)))
print('Timing to run %d evaluations: %lf mins\n' % (n_start,warmup_time))
print('Estimated time for %d trials: %lf hrs\n' % (n_iters,(1/60)*(warmup_time/n_start)*(n_sample+n_random)*(n_iters)))

# Run the trials
for i in range(0,n_iters):
    
    iter_df,top_idx,reagent_space,reagent_ei,data_array,fp_array=run_trial(reagent_list,n_sample,data_array,fp_array,reagent_space,reagent_ei,bo,n_random)
        
    running_df.append(iter_df)
    pd.concat(running_df).to_csv('running_df.csv',index=False)
    
    output.at[ctr,'Iteration']=i
    output.at[ctr,'Sampled']=len(pd.concat(running_df))
    output.at[ctr,'Max_score']=np.max(list(pd.concat(running_df)['score']))
    print('Trial: %d Top score: %lf' % (i,np.max(list(pd.concat(running_df)['score']))))
    ctr+=1
 
end=time.time()

# Final results
print('%d cycles complete' % (n_iters))
print('Timing: %lf hrs' % ((end-start)/60.0/60))
print('Top score: %lf' % (output.loc[len(output)-1,'Max_score']))
output.to_csv('results.csv',index=False)

