# timber

from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

##############################################################################

def ClusterFps(fps,cutoff=0.2):

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

def get_cluster_center(mols):
    fps = [AllChem.GetMorganFingerprintAsBitVect(x,3,1024) for x in mols]
    
    cs=ClusterFps(fps,cutoff=0.5)
    all_names=[m.GetProp('_Name') for m in mols]
    return all_names[cs[0][0]]

