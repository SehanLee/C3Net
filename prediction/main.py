import sys, os, time, glob
dir = os.getcwd()
dir = os.path.dirname(dir)
sys.path.append(dir)


from module.nn import Atomistic_Model, Atomwise
import torch
import torch.nn.functional as F
from torch.optim import Adam

import module
from module.nn import *
from module.data import *
import prediction.predictor as predictor
from module.input.input_generator import generator_predict

solvent = 'toluene'
solute = '../dataset/training_set/sdf/toluene_Solute_Prep_Conf5_RDKit.sdf'
input_path = '../dataset/training_set/npz'
generator_predict(solvent, solute, input_path)
#1=single conformer, 5=multiconformer
n_conformer = 1
checkpoint_path = './output/prediction'
data = MoleculeSet(input_path, checkpoint_path, n_conformer)
loader = MoleculeLoader(data, batch_size=2, num_workers=0)

# create model
reps = C3Net()
output = Atomwise(n_in=64)
model = Atomistic_Model(reps, output)

device = 'cpu'
predictor = predictor.Predictor(checkpoint_path, model, loader, device)

# start training
predictor.run(torch.device(device, 0))