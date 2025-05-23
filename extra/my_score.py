#!/usr/bin/env python
import pandas as pd
import numpy as np
import json
import sys

# Read Smiles from stdin
#my_smiles=[smiles.strip() for smiles in sys.stdin]
my_smiles=['CC','CCC']

# Make a df
df=pd.DataFrame({"Structure":my_smiles})

# Assign scores by some scoring function
df['Score']=np.random.uniform(0,1.,len(df))

# Send back to reinvent stdout
data = {"version":1, "payload":{"predictions":list(df["Score"])}}
print(json.dumps(data))

