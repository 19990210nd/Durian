import pandas as pd
import os
import multiprocessing as mp
#from rdkit import Chem
import pickle
#from prody import *
from glob import glob
import os.path as osp
import numpy as np


# protein preparation
def prepare_protein(protein_file=None, protein_path=None):
    pred_protein_file = protein_file.replace('.pdb', '.maegz')
    # cmd
    cmdline = 'cd %s && module load schrodinger/2021-1 && prepwizard -rehtreat -watdist 0 -disulfides ' \
                           '-fillsidechains -ms 10 -s -propka_pH 7.0 -minimize_adj_h -epik_pH 7.0 -epik_pHt 0.0 ' \
                           '-delwater_hbond_cutoff 3 -r 0.3 -fix -f 2005 -NOJOBID -WAIT  %s %s' % (protein_path, protein_file, pred_protein_file)
    os.system(cmdline)
    pdb_name = protein_file[:4]
    cmdline = 'cd %s &&' % protein_path
    cmdline += 'rm -rf %s-epik* &&' % pdb_name
    cmdline += 'rm -rf %s-impref* &&' % pdb_name
    cmdline += 'rm -rf %s-protassign* &&' % pdb_name
    cmdline += 'rm -rf %s-missing-side-chains* &&' % pdb_name
    cmdline += 'rm -rf %s-missing-loops*' % pdb_name
    os.system(cmdline)
    # final_pred_pdb = pdb_name + '_pred.pdb'
    # cmdline = 'cd %s &&' % protein_path
    # cmdline += 'structconvert -imae %s -opdb %s' % (pred_protein_file, final_pred_pdb)
    #os.system(cmdline)

def struconvert(protein_file=None, protein_path=None):
    pdb_name = protein_file[:4]
    final_pred_pdb = pdb_name + '_pred.pdb'
    cmdline = 'cd %s &&' % protein_path
    cmdline += 'structconvert -imae %s -opdb %s' % (pred_protein_file, final_pred_pdb)
    os.system(cmdline)

def mol2centroid(mol2_file):
    #schrodinger的rdbase和rdkit似乎不兼容，因此要先prepare完蛋白再module load rdkit
    #cmdline = "module purge && module load rdkit"
    #os.system(cmdline)
    from rdkit import Chem
    mol = Chem.MolFromMol2File(mol2_file, sanitize=False)
    lig_xyz = mol.GetConformer().GetPositions()
    centroid_x, centroid_y, centroid_z = lig_xyz.mean(axis=0)
    return centroid_x, centroid_y, centroid_z

def write_glide_grid_in(mae_file, grid_center_x, grid_center_y, grid_center_z, output_file):
    content = f"""GRID_CENTER   {grid_center_x}, {grid_center_y}, {grid_center_z}
GRIDFILE   glide-grid.zip
INNERBOX   10, 10, 10
OUTERBOX   27.2567433467, 27.2567433467, 27.2567433467
RECEP_FILE   {mae_file}
""" 
    # 使用 with open 写入文件
    with open(output_file, 'w') as f:
        f.write(content)

def generate_grid(output_file):
    dir_path, file_name = os.path.split(output_file)
    cmdline = "module purge && module load schrodinger &&"
    cmdline += "cd %s &&" % dir_path
    cmdline += "\"${SCHRODINGER}/glide\" glide-gird.in -OVERWRITE -HOST \"localhost:48\" -TMPLAUNCHDIR -ATTACHED -WAIT"
    os.system(cmdline)



def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()
    if os.path.exists(output_file):
        print(f'{output_file}  Have been written!')
    else:
        print(f'{output_file}: Failed, not created!')



def glide_dock_score(grid_file, ligand_file, output_dir):
    #output_dir = os.path.dirname(com_prmtop)
    output_dir_final = output_dir + '/glide'
    cmdline = "mkdir %s" % output_dir_final
    os.system(cmdline)
    content1 = f"""GRIDFILE   {grid_file}
LIGANDFILE   {ligand_file}
POSE_OUTTYPE   ligandlib
POSTDOCK_NPOSE   1
PRECISION   SP
RINGCONFCUT   1.0
""" 
    content2 = f"""DOCKING_METHOD   inplace
GRIDFILE   {grid_file}
LIGANDFILE   {ligand_file}
POSE_OUTTYPE   ligandlib
POSTDOCK   False
PRECISION   SP
RINGCONFCUT   1.0
""" 
    if os.path.isdir(output_dir_final):
        glide_dock_file = os.path.join(output_dir_final,f'glide-dock_SP_flexible.in')
        glide_score_file = os.path.join(output_dir_final,f'glide-dock_SP_score.in')
        write_file(glide_dock_file, content1)
        write_file(glide_score_file, content2)
        cmdline = "module purge && module load schrodinger &&"
        cmdline += "cd %s &&" % output_dir_final
        cmdline +="\"${SCHRODINGER}/glide\" glide-dock_SP_flexible.in -OVERWRITE -NJOBS 1 -HOST \"localhost:48\" -TMPLAUNCHDIR -ATTACHED -WAIT &&"
        cmdline +="\"${SCHRODINGER}/glide\" glide-dock_SP_score.in -OVERWRITE -NJOBS 1 -HOST \"localhost:48\" -TMPLAUNCHDIR -ATTACHED -WAIT"
        os.system(cmdline)
        print("Finished: %s" % output_dir)
    else:
        print("ERROR: %s" % output_dir)


def addh_ligprep(target, ligand_file, out_file):
    #ligprep -isd crystal.sdf -osd crystal_addH.sdf -Rh -a -m 1 -NOJOBID -WAIT
    # ligprep -isd crystal.sdf -osd crystal_addH.sdf  -nt -ns -a -NOJOBID -WAIT
    cmdline = 'cd %s && module load schrodinger/2024-1 && ligprep -isd %s -osd %s -Rh -a -m 1 -NOJOBID -WAIT' % (target, ligand_file, out_file)
    os.system(cmdline)


# Function to read molecules from SDF
def read_sdf(sdf_file):
    from rdkit import Chem
    supp = Chem.SDMolSupplier(sdf_file, removeHs=False)  # 设置为False以保留氢原子
    mols_list = [mol for mol in supp if mol is not None]
    return mols_list

# Function to write unique molecules to SDF
def write_sdf(mols_list, output_file):
    from rdkit import Chem
    writer = Chem.SDWriter(output_file)
    for mol in mols_list:
        if mol is not None:
            writer.write(mol)
    writer.close()

# Function to remove duplicates based on SMILES
def remove_duplicates(input_sdf, output_sdf):
    from rdkit import Chem
    mols = read_sdf(input_sdf)
    unique_smiles = set()
    unique_mols = []
    for mol in mols:
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)  # Convert to SMILES without hydrogens
        if smiles not in unique_smiles:
            unique_smiles.add(smiles)
            unique_mols.append(mol)  # Keep original molecule with its hydrogens and conformation
    write_sdf(unique_mols, output_sdf)

source_dir3 = '/home/niedou/durain/real_world/durain_graphbp'
target_dir3 = '/home/niedou/durain/real_world/real-world_dataset'

#prepare protein
for subfolder in os.listdir(target_dir3):
    target_subfolder_path = os.path.join(target_dir3, subfolder)
    if os.path.isdir(target_subfolder_path):
            # prepare protein
             for file_name in os.listdir(target_subfolder_path):
                 if file_name.endswith('.pdb'):
                     dir_path, protein_name = os.path.split(file_name)
                     prepare_protein(file_name, target_subfolder_path)
                     print("%s has been prepared." % file_name)
                 #if file_name.endswith('.mol2'):
                 #    mol2_file = target_subfolder_path + '/' + file_name    
                 #    centroid_x, centroid_y, centroid_z = mol2centroid(mol2_file)
                 #generate_grid(mae_file, centroid_x, centroid_y, centroid_z, out_file)

#generate grid-infile
for subfolder in os.listdir(target_dir3):
    target_subfolder_path = os.path.join(target_dir3, subfolder)
    if os.path.isdir(target_subfolder_path):
        mae_file = None
        centroid_x = 0
        centroid_y = 0
        centroid_z = 0
        grid_infile = target_subfolder_path + '/glide-gird.in'
        for file_name in os.listdir(target_subfolder_path):
            #print(file_name)
            if file_name.endswith('.maegz'):
                   mae_file = file_name
            if file_name.endswith('.mol2'):
                    mol2_file = target_subfolder_path + '/' + file_name    
                    centroid_x, centroid_y, centroid_z = mol2centroid(mol2_file)
        write_glide_grid_in(mae_file, centroid_x, centroid_y, centroid_z, grid_infile)


#generate grid
for subfolder in os.listdir(target_dir3):
    target_subfolder_path = os.path.join(target_dir3, subfolder)
    if os.path.isdir(target_subfolder_path):
        grid_infile = target_subfolder_path + '/glide-gird.in'
        generate_grid(grid_infile)
        print("Finished: %s" % grid_infile)


#addH for ligands
methods = glob('/home/niedou/durain/real_world/durain_*')
#out_dir = '/home/niedou/durain/real_world/glide_result'
for method in methods:
    targets_path = method + '/*'
    targets = glob(targets_path)
    for target in targets:
        print(target)
        #ligand_path = target + '/SDF'
        gen_ligands = target + "/SDF/gen.sdf"
        gen_ligands_H = target + "/SDF/gen_addH.sdf"
        ori_ligand = target + "/crystal.sdf"
        ori_ligand_H = target + "/crystal_addH.sdf"
        addh_ligprep(target, ori_ligand, ori_ligand_H)
        addh_ligprep(target, gen_ligands, gen_ligands_H)


#Glide
#for realworld dataset
target_dir = '/home/niedou/durain/real_world/real-world_dataset'
methods = glob('/home/niedou/durain/real_world/durain_*')
for method in methods:
    targets_path = method + '/*'
    targets = glob(targets_path)
    for target in targets:
        #print(target)
        ligand_path = target + '/SDF'
        ligand_file = ligand_path + "/gen_addH.sdf"
        grid_file = target_dir + "/" + target.split('/')[-1] + "/glide-grid.zip"
        glide_dock_score(grid_file, ligand_file, target)
