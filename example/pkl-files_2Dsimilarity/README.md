download trainmols_fps.pkl from (https://drive.google.com/file/d/1TS7-jMJs3esuXmOPsXnyeprzl05u-oWt/view?usp=sharing)


import pickle
from glob import glob
import os
import os.path as osp
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def read_pkl(file_path):
    with open(file_path, 'rb') as file1:
        file = pickle.load(file1)
    return file

def write_pkl(file_path, data):
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)


def calculate_fps(mols):
    fps = []
    for mol in mols:
        try:
            fp = Chem.RDKFingerprint(mol)
            fps.append(fp)
        except Exception as e:
            print(e)
    return fps


def calculate_similarity_from_fps(fps_a, fps_b):
    sims = []
    for fp_a in fps_a:
        for fp_b in fps_b:
            sim = DataStructs.TanimotoSimilarity(fp_a, fp_b)
            sims.append(sim)
    return sims
