#!/usr/bin/env python
import math
import os
import torch
import torchani
import ase 
import numpy as np

def spe_to_int(ipt):
    opt = []
    #dic = 'HCNO'
    for s in ipt:
        for i, c in enumerate([1, 6, 7, 8, 16, 9, 17]):
            if s == c:
                opt.append(i)
                # print(dic[i], end=' ')
    return np.array(opt)

def pad_ase_molecule(mol_list, device):
    """
    input: ase molecule list
    return: [species, coordinates]
            padded species and coordinates input for torchani (in torch tensor)
    """
    species = []
    coordinates = []
    for mol in mol_list:
        species.append(torch.from_numpy(spe_to_int(mol.get_atomic_numbers())).to(device))
        coordinates.append(torch.from_numpy(mol.get_positions()).to(device))
    species = torch.nn.utils.rnn.pad_sequence(species,
                                              batch_first=True,
                                              padding_value=-1,)
    coordinates = torch.nn.utils.rnn.pad_sequence(coordinates,
                                              batch_first=True,
                                              padding_value=0.0, )

    return species, coordinates

# all of this is setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

Rcr = 5.1000e+00
Rca = 3.5000e+00
EtaR = torch.tensor([1.9700000e+01], device=device)
ShfR = torch.tensor([8.0000000e-01,1.0687500e+00,1.3375000e+00,1.6062500e+00,1.8750000e+00,2.1437500e+00,2.4125000e+00,2.6812500e+00,2.9500000e+00,3.2187500e+00,3.4875000e+00,3.7562500e+00,4.0250000e+00,4.2937500e+00,4.5625000e+00,4.8312500e+00], device=device)
Zeta = torch.tensor([1.4100000e+01], device=device)
ShfZ = torch.tensor([3.9269908e-01,1.1780972e+00,1.9634954e+00,2.7488936e+00], device=device)
EtaA = torch.tensor([1.2500000e+01], device=device)
ShfA = torch.tensor([8.0000000e-01,1.1375000e+00,1.4750000e+00,1.8125000e+00,2.1500000e+00,2.4875000e+00,2.8250000e+00,3.1625000e+00], device=device)
species_order = ['H','C','N','O','S','F','Cl']

num_species = len(species_order)

aev_computer = torchani.AEVComputer(Rcr, Rca, EtaR, ShfR, EtaA, Zeta, ShfA, ShfZ, num_species)

# create a methane molecule in ASE
coordinates = torch.tensor([[[0.03192167, 0.00638559, 0.01301679],
                             [-0.83140486, 0.39370209, -0.26395324],
                             [-0.66518241, -0.84461308, 0.20759389],
                             [0.45554739, 0.54289633, 0.81170881],
                             [0.66091919, -0.16799635, -0.91037834]]],
                           requires_grad=True)
methane = ase.Atoms(['C', 'H', 'H', 'H', 'H'], positions=coordinates.squeeze().detach().numpy())

# create a water molecule in ASE
d = 0.9575
t = math.pi / 180 * 104.51
coordinates = torch.tensor([[[d,0,0],
                             [d * math.cos(t), d * math.sin(t), 0],
                             [0,0,0]]],requires_grad=True)

water = ase.Atoms(['O','H','H'], positions=coordinates.squeeze().detach().numpy())

# first, pad the species and coords
species, coordinates=pad_ase_molecule([methane,water],device)

# get the AEV
aev_result=aev_computer((species, coordinates), cell=None, pbc=None)

#print(aev_computer.aev_length)

print('Species:',aev_result.species)
print('AEV:',aev_result.aevs)

# convert to numpy
np_aev=aev_result.aevs.detach().numpy()

print('Tensor shape:',np_aev.shape)

