#module load rdkit
#module load openbabel
import sys
sys.path.append('/home/niedou/durain/scripts')
import pickle
from glob import glob
import os
import os.path as osp
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import espsim
from mol_property import read_sdf
from mol_property import MoleculeProperties
import random
from sklearn.model_selection import train_test_split


def read_sdf(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file)
    mols_list = [i for i in supp]
    return mols_list

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

def read_mols(methods, train_targets, test_targets, k):
    path_marker = '/'
    print(k)
    vina_scores = {}
    mol_list = {}
    ori_mols = {}
    vina_scores2 = {}
    mol_list2 = {}
    ori_mols2 = {}
    mol_list_target = {}
    ori_mols_target = {}
    mol_list_target2 = {}
    ori_mols_target2 = {}
    for method in methods:
        mol_list[method] = []
        vina_scores[method] = []
        ori_mols[method] = []
        mol_list2[method] = []
        vina_scores2[method] = []
        ori_mols2[method] = []
        targets = glob(osp.join(method, '*'))
        try:
            #train_targets, test_targets = train_test_split(targets, test_size=0.5, random_state=random.randint(0, 10000))
            #print("OK! %s" % method)
            for target in tqdm(train_targets):
              base_path = target[:target.rfind('/')]
              if base_path == method:
                mol_list_target[target] = []
                ori_mols_target[target] = []         
                target_vina_mol = []
                try:
                    ori_mol2_file = glob(osp.join(target, 'qvina', 'ori', '*.mol2'))
                    ori_mol2 = Chem.MolFromMol2File(ori_mol2_file[0])
                    if ori_mol2 is not None:
                        ori_mols[method].append(ori_mol2)
                        ori_mols_target[target].append(ori_mol2)
                    else:
                        continue
                except Exception as e:
                    print(e)
                    continue
                for i in range(k):
                    try:
                        mol_vina = osp.join(target, 'qvina', f'{i}')
                        top1_score = glob(osp.join(mol_vina, 'top1_*'))
                        mol2_file = glob(osp.join(mol_vina, '*.mol2'))
                        if len(top1_score) > 0:
                            vina = float(np.loadtxt(top1_score[0]))
                            mol = Chem.MolFromMol2File(mol2_file[0])
                            target_vina_mol.append((vina, mol))
                    except Exception as e:
                        print(e)
                sorted_vina_mol = sorted(target_vina_mol, key=lambda x: x[0])
                sorted_mol = [i[1] for i in sorted_vina_mol]
                mol_list[method].append(sorted_mol)
                mol_list_target[target].append(sorted_mol)
            ####
            for target in tqdm(test_targets):
              base_path = target[:target.rfind('/')]
              if base_path == method:
                mol_list_target2[target] = []
                ori_mols_target2[target] = []               
                target_vina_mol = []
                try:
                    ori_mol2_file = glob(osp.join(target, 'qvina', 'ori', '*.mol2'))
                    ori_mol2 = Chem.MolFromMol2File(ori_mol2_file[0])
                    if ori_mol2 is not None:
                        ori_mols2[method].append(ori_mol2)
                        ori_mols_target2[target].append(ori_mol2)
                    else:
                        continue
                except Exception as e:
                    print(e)
                    continue
                for i in range(k):
                    try:
                        mol_vina = osp.join(target, 'qvina', f'{i}')
                        top1_score = glob(osp.join(mol_vina, 'top1_*'))
                        mol2_file = glob(osp.join(mol_vina, '*.mol2'))
                        if len(top1_score) > 0:
                            vina = float(np.loadtxt(top1_score[0]))
                            mol = Chem.MolFromMol2File(mol2_file[0])
                            target_vina_mol.append((vina, mol))
                    except Exception as e:
                        print(e)
                sorted_vina_mol = sorted(target_vina_mol, key=lambda x: x[0])
                sorted_mol = [i[1] for i in sorted_vina_mol]
                mol_list2[method].append(sorted_mol)
                mol_list_target2[target].append(sorted_mol)
        except:
            print("Error! %s" % method)
    # num
    all_mols = {}
    all_mols_wonan = {}
    all_mols2 = {}
    all_mols_wonan2 = {}
    #
    all_mols_target = {}
    all_mols_wonan_target = {}
    all_mols_target2 = {}
    all_mols_wonan_target2 = {}
    for method in methods:
            all_mols[method] = []
            all_mols_wonan[method] = []
            all_mols2[method] = []
            all_mols_wonan2[method] = []
            for target_mols in mol_list[method]:
                all_mols[method].extend(target_mols[:k])
            all_mols_wonan[method] = [i for i in all_mols[method] if i is not None]
            for target_mols in mol_list2[method]:
                all_mols2[method].extend(target_mols[:k])
            all_mols_wonan2[method] = [i for i in all_mols2[method] if i is not None]
            # print(method.split('/')[-1],'The number of Nan', len(all_mols[method]) - len(all_mols_wonan[method]))
            print(method.split('/')[-1], 'The number of train_targets Mol', len(all_mols_wonan[method]))
            print(method.split('/')[-1], 'The number of test_targets Mol', len(all_mols_wonan2[method]))
            for target in tqdm(train_targets):
                all_mols_target[target] = []
                all_mols_wonan_target[target] = []
                for target_mols in mol_list_target[target]:
                    all_mols_target[target].extend(target_mols[:k])
                all_mols_wonan_target[target] = [i for i in all_mols_target[target] if i is not None]
            for target in tqdm(test_targets):
                all_mols_target2[target] = []
                all_mols_wonan_target2[target] = []
                for target_mols in mol_list_target2[target]:
                    all_mols_target2[target].extend(target_mols[:k])
                all_mols_wonan_target2[target] = [i for i in all_mols_target2[target] if i is not None]
    sum1 = 0
    sum2 = 0
    for target in tqdm(train_targets):
        print(target.split('/')[-1], 'The number of train_targets Mol', len(all_mols_wonan_target[target]))
        if len(all_mols_wonan_target[target]) != 0:
         sum1 += len(all_mols_wonan_target[target])
    print(method.split('/')[-1], 'The number of train_targets Mol', sum1)
    for target in tqdm(test_targets):
        print(target.split('/')[-1], 'The number of test_targets Mol', len(all_mols_wonan_target2[target]))
        if len(all_mols_wonan_target2[target]) != 0:
            sum2 += len(all_mols_wonan_target2[target])
    print(method.split('/')[-1], 'The number of test_targets Mol', sum2)
    mol_file1 = out_dir + '/genmols_method1.pkl'
    mol_file2 = out_dir + '/genmols_method2.pkl'
    mol_file3 = out_dir + '/genmols_target1.pkl'
    mol_file4 = out_dir + '/genmols_target2.pkl'
    ori_file1 = out_dir + '/ori_mol_method1.pkl'
    ori_file2 = out_dir + '/ori_mol_method2.pkl'
    ori_file3 = out_dir + '/ori_mol_target1.pkl'
    ori_file4 = out_dir + '/ori_mol_target2.pkl'
    if not os.path.exists(mol_file1):
        write_pkl(mol_file1, all_mols_wonan)
    else:
        print("%s has been exsited." % mol_file1)
    if not os.path.exists(mol_file2):
        write_pkl(mol_file2, all_mols_wonan2)
    else:
        print("%s has been exsited." % mol_file2)
    if not os.path.exists(mol_file3):
        write_pkl(mol_file3, all_mols_wonan_target)
    else:
        print("%s has been exsited." % mol_file3)
    if not os.path.exists(mol_file4):
        write_pkl(mol_file4, all_mols_wonan_target2)
    else:
        print("%s has been exsited." % mol_file4)
    if not os.path.exists(ori_file1):
     write_pkl(ori_file1, ori_mols)
    else:
     print("%s has been exsited." % ori_file1)
    if not os.path.exists(ori_file2):
        write_pkl(ori_file2, ori_mols2)
    else:
        print("%s has been exsited." % ori_file2)
    if not os.path.exists(ori_file3):
     write_pkl(ori_file3, ori_mols_target)
    else:
     print("%s has been exsited." % ori_file3)
    if not os.path.exists(ori_file4):
     write_pkl(ori_file4, ori_mols_target2)
    else:
     print("%s has been exsited." % ori_file4)
    return ori_mols, ori_mols2, ori_mols_target, ori_mols_target2, all_mols_wonan, all_mols_wonan2, all_mols_wonan_target, all_mols_wonan_target2


def druglikeness(methods, all_mols_wonan, all_mols_wonan2, out_dir):
    calculator = MoleculeProperties()
    logps = {}
    qeds = {}
    sas = {}
    lipinskis = {}
    tpsas = {}
    mws = {}
    logps_mean = []
    qeds_mean = []
    sas_mean = []
    lipinskis_mean = []
    tpsas_mean = []
    mws_mean = []
    logps_mean_std = []
    qeds_mean_std = []
    sas_mean_std = []
    lipinskis_mean_std = []
    tpsas_mean_std = []
    mws_mean_std = []
    #
    logps2 = {}
    qeds2 = {}
    sas2 = {}
    lipinskis2 = {}
    tpsas2 = {}
    mws2 = {}
    logps_mean2 = []
    qeds_mean2 = []
    sas_mean2 = []
    lipinskis_mean2 = []
    tpsas_mean2 = []
    mws_mean2 = []
    logps_mean_std2 = []
    qeds_mean_std2 = []
    sas_mean_std2 = []
    lipinskis_mean_std2 = []
    tpsas_mean_std2 = []
    mws_mean_std2 = []
    for method in methods:
        logps[method] = [calculator.calculate_logp(i) for i in all_mols_wonan[method]]
        qeds[method] = [calculator.calculate_qed(i) for i in all_mols_wonan[method]]
        sas[method] = [calculator.calculate_sa(i) for i in all_mols_wonan[method]]
        lipinskis[method] = [calculator.calculate_lipinski(i) for i in all_mols_wonan[method]]
        tpsas[method] = [calculator.calculate_tpsa(i) for i in all_mols_wonan[method]]
        mws[method] = [calculator.calculate_mw(i) for i in all_mols_wonan[method]]
        #
        logps2[method] = [calculator.calculate_logp(i) for i in all_mols_wonan2[method]]
        qeds2[method] = [calculator.calculate_qed(i) for i in all_mols_wonan2[method]]
        sas2[method] = [calculator.calculate_sa(i) for i in all_mols_wonan2[method]]
        lipinskis2[method] = [calculator.calculate_lipinski(i) for i in all_mols_wonan2[method]]
        tpsas2[method] = [calculator.calculate_tpsa(i) for i in all_mols_wonan2[method]]
        mws2[method] = [calculator.calculate_mw(i) for i in all_mols_wonan2[method]]
    column_names = []
    for method in methods:
        method_name = method.split('/')[-1]
        column_names.append(method_name)
        logps_mean.append(np.mean(logps[method]))
        qeds_mean.append(np.mean(qeds[method]))
        sas_mean.append(np.mean(sas[method]))
        lipinskis_mean.append(np.mean(lipinskis[method]))
        tpsas_mean.append(np.mean(tpsas[method]))
        mws_mean.append(np.mean(mws[method]))
        #
        logps_mean_std.append(np.std(logps[method]))
        qeds_mean_std.append(np.std(qeds[method]))
        sas_mean_std.append(np.std(sas[method]))
        lipinskis_mean_std.append(np.std(lipinskis[method]))
        tpsas_mean_std.append(np.std(tpsas[method]))
        mws_mean_std.append(np.std(mws[method]))
        #
        logps_mean2.append(np.mean(logps2[method]))
        qeds_mean2.append(np.mean(qeds2[method]))
        sas_mean2.append(np.mean(sas2[method]))
        lipinskis_mean2.append(np.mean(lipinskis2[method]))
        tpsas_mean2.append(np.mean(tpsas2[method]))
        mws_mean2.append(np.mean(mws2[method]))
        #
        logps_mean_std2.append(np.std(logps2[method]))
        qeds_mean_std2.append(np.std(qeds2[method]))
        sas_mean_std2.append(np.std(sas2[method]))
        lipinskis_mean_std2.append(np.std(lipinskis2[method]))
        tpsas_mean_std2.append(np.std(tpsas2[method]))
        mws_mean_std2.append(np.std(mws2[method]))
    df = pd.DataFrame()
    df['method'] = column_names
    df['logps'] = logps_mean
    df['logps'] = logps_mean
    df['qeds'] = qeds_mean
    df['sas'] = sas_mean
    df['lipinskis'] = lipinskis_mean
    df['tpsas'] = tpsas_mean
    df['mws'] = mws_mean
    df['logps2'] = logps_mean2
    df['qeds2'] = qeds_mean2
    df['sas2'] = sas_mean2
    df['lipinskis2'] = lipinskis_mean2
    df['tpsas2'] = tpsas_mean2
    df['mws2'] = mws_mean2
    file_name = out_dir + '/druglikeness.csv'
    df.to_csv(file_name)
    df2 = pd.DataFrame()
    df2['method'] = column_names
    df2['logps'] = logps_mean
    df2['qeds'] = qeds_mean
    df2['sas'] = sas_mean
    df2['lipinskis'] = lipinskis_mean
    df2['tpsas'] = tpsas_mean
    df2['mws'] = mws_mean
    df2['logps_std'] = logps_mean_std
    df2['qeds_std'] = qeds_mean_std
    df2['sas_std'] = sas_mean_std
    df2['lipinskis_std'] = lipinskis_mean_std
    df2['tpsas_std'] = tpsas_mean_std
    df2['mws_std'] = mws_mean_std
    df2['logps2'] = logps_mean2
    df2['qeds2'] = qeds_mean2
    df2['sas2'] = sas_mean2
    df2['lipinskis2'] = lipinskis_mean2
    df2['tpsas2'] = tpsas_mean2
    df2['mws2'] = mws_mean2
    df2['logps_std2'] = logps_mean_std2
    df2['qeds_std2'] = qeds_mean_std2
    df2['sas_std2'] = sas_mean_std2
    df2['lipinskis_std2'] = lipinskis_mean_std2
    df2['tpsas_std2'] = tpsas_mean_std2
    df2['mws_std2'] = mws_mean_std
    file_name2 = out_dir + '/druglikeness_std.csv'
    df2.to_csv(file_name2)


def similarity_2d(methods, all_mols_wonan, all_mols_wonan2, out_dir):
    path_marker = '/'
    #2dsimilarity*change fps_files
    calculator = MoleculeProperties()
    train_path = '/home/niedou/durain/entire/analysis_20240922/fps/trainmols_fps.pkl'
    val_path = '/home/niedou/durain/entire/analysis_20240922/fps/valmols_fps.pkl'
    tra_fps = read_pkl(train_path)
    val_fps = read_pkl(val_path)
    cmdline = 'mkdir %s/fps_1 &&' % out_dir
    cmdline += 'mkdir %s/fps_2' % out_dir
    os.system(cmdline)
    for method in methods:
        mol_fps1 = calculate_fps(all_mols_wonan[method])
        mol_fps2 = calculate_fps(all_mols_wonan2[method])
        fps_file1 = out_dir + '/fps_1/%s_genmols_fps.pkl' % method.split('/')[-1]
        fps_file2 = out_dir + '/fps_2/%s_genmols_fps.pkl' % method.split('/')[-1]
        if not os.path.exists(fps_file1):
            write_pkl(fps_file1, mol_fps1)
        else:
            print("%s has been exsited." % fps_file1)
        if not os.path.exists(fps_file2):
            write_pkl(fps_file2, mol_fps2)
        else:
            print("%s has been exsited." % fps_file2)
    #for crossdock dataset
    fps_files1 = glob(out_dir +'/fps_1/*.pkl')
    fps_files2 = glob(out_dir +'/fps_2/*.pkl')
    fps = {}
    tra_smi = {}
    val_smi = {}
    fps_names = []
    tra_smi_mean = []
    val_smi_mean = []
    tra_smi_mean_std = []
    val_smi_mean_std = []
    tra_smi2 = {}
    val_smi2 = {}
    column_names = []
    for fps_file in fps_files1:
        fps = read_pkl(fps_file)
        fps_name = fps_file.split('/')[-1].split('.')[0]
        column_names.append(fps_name + '_1')
        fps_names.append(fps_file.split('/')[-1].split('.')[0])
        val_smi[fps_name] = calculate_similarity_from_fps(val_fps, fps)
        tra_smi[fps_name] = calculate_similarity_from_fps(tra_fps, fps)
        #k = len(tra_smi[fps_name])
        val_smi_mean.append(np.mean(val_smi[fps_name]))
        tra_smi_mean.append(np.mean(tra_smi[fps_name]))
        val_smi_mean_std.append(np.std(val_smi[fps_name]))
        tra_smi_mean_std.append(np.std(tra_smi[fps_name]))
    for fps_file in fps_files2:
        fps = read_pkl(fps_file)
        fps_name = fps_file.split('/')[-1].split('.')[0]
        column_names.append(fps_name+ '_2')
        #fps_names2.append(fps_file.split('/')[-1].split('.')[0])
        val_smi2[fps_name] = calculate_similarity_from_fps(val_fps, fps)
        tra_smi2[fps_name] = calculate_similarity_from_fps(tra_fps, fps)
        #k = len(tra_smi[fps_name])
        val_smi_mean.append(np.mean(val_smi2[fps_name]))
        tra_smi_mean.append(np.mean(tra_smi2[fps_name]))
        val_smi_mean_std.append(np.std(val_smi2[fps_name]))
        tra_smi_mean_std.append(np.std(tra_smi2[fps_name]))
    df= pd.DataFrame()
    df['Method'] = column_names
    df['Train. Sim'] = tra_smi_mean
    df['Train. Sim_std'] = tra_smi_mean_std
    df['Test. Sim'] = val_smi_mean
    df['Test. Sim_std'] = val_smi_mean_std
    df.to_csv(out_dir + '/smi2d.csv')


def similarity_3d(methods, train_targets, test_targets, ori_mols_target, ori_mols_target2, all_mols_wonan_target, all_mols_wonan_target2, k=100):
    path_marker = '/'
    calculator = MoleculeProperties()
    sha_smi = {}
    pha_smi = {}
    ele_smi = {}
    sha_smi_mean = {}
    pha_smi_mean = {}
    ele_smi_mean = {}
    sha_smi_mean_final = {}
    pha_smi_mean_final = {}
    ele_smi_mean_final = {}
    sha_smi_mean_out = []
    pha_smi_mean_out = []
    ele_smi_mean_out = []
    sha_smi_std = {}
    pha_smi_std = {}
    ele_smi_std = {}
    sha_smi_std_final = {}
    pha_smi_std_final = {}
    ele_smi_std_final = {}
    sha_smi_std_out = []
    pha_smi_std_out = []
    ele_smi_std_out = []
    column_names = []
    ori_mols_H = {}
    gen_mols_H = {}
    #
    sha_smi2 = {}
    pha_smi2 = {}
    ele_smi2 = {}
    sha_smi_mean2 = {}
    pha_smi_mean2 = {}
    ele_smi_mean2 = {}
    sha_smi_mean_final2 = {}
    pha_smi_mean_final2 = {}
    ele_smi_mean_final2 = {}
    sha_smi_mean_out2 = []
    pha_smi_mean_out2 = []
    ele_smi_mean_out2 = []
    sha_smi_std2 = {}
    pha_smi_std2 = {}
    ele_smi_std2 = {}
    sha_smi_std_final2 = {}
    pha_smi_std_final2 = {}
    ele_smi_std_final2 = {}
    sha_smi_std_out2 = []
    pha_smi_std_out2 = []
    ele_smi_std_out2 = []
    #column_names2 = []
    ori_mols_H2 = {}
    gen_mols_H2 = {}
    for method in methods:
        print(method)
        method_name = method.split('/')[-1]
        column_names.append(method_name)
        sha_smi_mean_final[method] = []
        pha_smi_mean_final[method] = []
        ele_smi_mean_final[method] = []
        sha_smi_std_final[method] = []
        pha_smi_std_final[method] = []
        ele_smi_std_final[method] = []
        sha_smi_mean_final2[method] = []
        pha_smi_mean_final2[method] = []
        ele_smi_mean_final2[method] = []
        sha_smi_std_final2[method] = []
        pha_smi_std_final2[method] = []
        ele_smi_std_final2[method] = []
        for target in train_targets:
          base_path = target[:target.rfind('/')]
          #print(target)
          #print(base_path)
          if base_path == method:
            print(target)
            sha_smi[target] = []
            pha_smi[target] = []
            ele_smi[target] = []
            sha_smi_mean[target] = []
            pha_smi_mean[target] = []
            ele_smi_mean[target] = []
            sha_smi_std[target] = []
            pha_smi_std[target] = []
            ele_smi_std[target] = []
            ori_mols_H[target] = []
            gen_mols_H[target] = []
            try:
                if ori_mols_target[target] is not None:
                    sdf_mol = ori_mols_target[target][0]
                    romol = Chem.Mol(sdf_mol)
                    romol_H = Chem.AddHs(romol)
                    ori_mols_H[target].append(romol_H)
                else:
                    continue
                for gen_mol_file in all_mols_wonan_target[target]:
                    try:
                    # gen_mol2 = Chem.MolFromMol2File(gen_mol_file)
                    # gen_mol2 = read_sdf(gen_mol_file)
                    #gen_mol_supplier = Chem.SDMolSupplier(gen_mol_file)
                    #gen_mol2 = [mol for mol in gen_mol_supplier][0]
                    # gen_mol2_H = gen_mol2
                        gen_mol2_H = Chem.AddHs(gen_mol_file)
                    # Chem.AllChem.ComputeGasteigerCharges(gen_mol2_H)
                        if gen_mol2_H is not None:
                            gen_mols_H[target].append(gen_mol2_H)
                    except Exception as e:
                        print(f"Error processing {gen_mol_file}: {str(e)}")
                        continue
            except Exception as e:
                print("error: %s" %target)
                print(e)
                continue
            # shape
            for i in gen_mols_H[target]:
                try:
                    sha_smi[target].append(calculator.calculate_shape(ori_mols_H[target][0], i))
                except:
                    sha_smi[target].append(0)
            for i in gen_mols_H[target]:
                try:
                    ele_tmp = calculator.calculate_esp(ori_mols_H[target][0], i)
                    import math
                    if math.isnan(ele_tmp):
                        ele_smi[target].append(0)
                    else:
                        ele_smi[target].append(ele_tmp)
                except:
                    ele_smi[target].append(0)
            for i in gen_mols_H[target]:
                try:
                    pha_smi[target].append(calculator.calculate_sc_score(ori_mols_H[target][0], i))
                except:
                    pha_smi[target].append(0)
            sha_smi_mean[target].append(np.mean(sha_smi[target]))
            ele_smi_mean[target].append(np.mean(ele_smi[target]))
            pha_smi_mean[target].append(np.mean(pha_smi[target]))
            sha_smi_mean_final[method].append(sha_smi_mean[target])
            ele_smi_mean_final[method].append(ele_smi_mean[target])
            pha_smi_mean_final[method].append(pha_smi_mean[target])
            sha_smi_std[target].append(np.std(sha_smi[target]))
            ele_smi_std[target].append(np.std(ele_smi[target]))
            pha_smi_std[target].append(np.std(pha_smi[target]))
            sha_smi_std_final[method].append(sha_smi_std[target])
            ele_smi_std_final[method].append(ele_smi_std[target])
            pha_smi_std_final[method].append(pha_smi_std[target])
        sha_smi_mean_out.append(np.mean(sha_smi_mean_final[method]))
        ele_smi_mean_out.append(np.mean(ele_smi_mean_final[method]))
        pha_smi_mean_out.append(np.mean(pha_smi_mean_final[method]))
        sha_smi_std_out.append(np.mean(sha_smi_std_final[method]))
        ele_smi_std_out.append(np.mean(ele_smi_std_final[method]))
        pha_smi_std_out.append(np.mean(pha_smi_std_final[method]))
        ####
        for target in test_targets:
          base_path = target[:target.rfind('/')]
          if base_path == method:
            sha_smi2[target] = []
            pha_smi2[target] = []
            ele_smi2[target] = []
            sha_smi_mean2[target] = []
            pha_smi_mean2[target] = []
            ele_smi_mean2[target] = []
            sha_smi_std2[target] = []
            pha_smi_std2[target] = []
            ele_smi_std2[target] = []
            ori_mols_H2[target] = []
            gen_mols_H2[target] = []
            try:
                if ori_mols_target2[target] is not None:
                    sdf_mol = ori_mols_target2[target][0]
                    romol = Chem.Mol(sdf_mol)
                    romol_H = Chem.AddHs(romol)
                    ori_mols_H2[target].append(romol_H)
                else:
                    continue
                for gen_mol_file in all_mols_wonan_target2[target]:
                    try:
                        gen_mol2_H = Chem.AddHs(gen_mol_file)
                        # Chem.AllChem.ComputeGasteigerCharges(gen_mol2_H)
                        if gen_mol2_H is not None:
                            gen_mols_H2[target].append(gen_mol2_H)
                    except Exception as e:
                        print(f"Error processing {gen_mol_file}: {str(e)}")
                        continue
            except Exception as e:
                print("error: %s" %target)
                print(e)
                continue
            # shape
            for i in gen_mols_H2[target]:
                try:
                    sha_smi2[target].append(calculator.calculate_shape(ori_mols_H2[target][0], i))
                except:
                    sha_smi2[target].append(0)
        # charge
            for i in gen_mols_H2[target]:
                try:
                    ele_tmp = calculator.calculate_esp(ori_mols_H2[target][0], i)
                    import math
                    if math.isnan(ele_tmp):
                        ele_smi2[target].append(0)
                    else:
                        ele_smi2[target].append(ele_tmp)
                except:
                    ele_smi2[target].append(0)
            # pharm
            for i in gen_mols_H2[target]:
                try:
                    pha_smi2[target].append(calculator.calculate_sc_score(ori_mols_H2[target][0], i))
                except:
                    pha_smi2[target].append(0)
            sha_smi_mean2[target].append(np.mean(sha_smi2[target]))
            ele_smi_mean2[target].append(np.mean(ele_smi2[target]))
            pha_smi_mean2[target].append(np.mean(pha_smi2[target]))
            sha_smi_mean_final2[method].append(sha_smi_mean2[target])
            ele_smi_mean_final2[method].append(ele_smi_mean2[target])
            pha_smi_mean_final2[method].append(pha_smi_mean2[target])
            sha_smi_std2[target].append(np.std(sha_smi2[target]))
            ele_smi_std2[target].append(np.std(ele_smi2[target]))
            pha_smi_std2[target].append(np.std(pha_smi2[target]))
            sha_smi_std_final2[method].append(sha_smi_std2[target])
            ele_smi_std_final2[method].append(ele_smi_std2[target])
            pha_smi_std_final2[method].append(pha_smi_std2[target])
        sha_smi_mean_out2.append(np.mean(sha_smi_mean_final2[method]))
        ele_smi_mean_out2.append(np.mean(ele_smi_mean_final2[method]))
        pha_smi_mean_out2.append(np.mean(pha_smi_mean_final2[method]))
        sha_smi_std_out2.append(np.mean(sha_smi_std_final2[method]))
        ele_smi_std_out2.append(np.mean(ele_smi_std_final2[method]))
        pha_smi_std_out2.append(np.mean(pha_smi_std_final2[method]))
    # outputs
    out_file = out_dir + '/smi3d_top100.csv'
    df = pd.DataFrame()
    df = df.append(pd.Series(sha_smi_mean_out, name='sha_smi'))
    df = df.append(pd.Series(sha_smi_std_out, name='sha_std'))
    df = df.append(pd.Series(ele_smi_mean_out, name='ele_smi'))
    df = df.append(pd.Series(ele_smi_std_out, name='ele_std'))
    df = df.append(pd.Series(pha_smi_mean_out, name='pha_smi'))
    df = df.append(pd.Series(pha_smi_std_out, name='pha_std'))
    df = df.append(pd.Series(sha_smi_mean_out2, name='sha_smi2'))
    df = df.append(pd.Series(sha_smi_std_out2, name='sha_std2'))
    df = df.append(pd.Series(ele_smi_mean_out2, name='ele_smi2'))
    df = df.append(pd.Series(ele_smi_std_out2, name='ele_std2'))
    df = df.append(pd.Series(pha_smi_mean_out2, name='pha_smi2'))
    df = df.append(pd.Series(pha_smi_std_out2, name='pha_std2'))
    df.rename(columns=dict(zip(df.columns, column_names)), inplace=True)
    df.to_csv(out_file)




#####read_mols
methods = glob('/home/niedou/durain/entire/crossdock/durian_*')
num_iterations = 10
k=100
path_marker = '/'
# 开始重复10次随机划分
for i in range(num_iterations):
    print(f"Iteration {i}:")
    out_dir = '/home/niedou/durain/entire/analysis_20241010/analysis_%s' % i
    #cmdline = 'mkdir %s' % out_dir
    #os.system(cmdline)
    train_targets_file = out_dir + '/train_targets.pkl'
    test_targets_file = out_dir + '/test_targets.pkl'
    train_targets = read_pkl(train_targets_file)
    test_targets = read_pkl(test_targets_file)
    ori_mols, ori_mols2, ori_mols_target, ori_mols_target2, all_mols_wonan, all_mols_wonan2, all_mols_wonan_target, all_mols_wonan_target2 = read_mols(methods, train_targets, test_targets, k)
    druglikeness(methods, all_mols_wonan, all_mols_wonan2, out_dir)
    similarity_2d(methods, ori_mols, ori_mols2, all_mols_wonan, all_mols_wonan2, out_dir)
    similarity_3d(methods, train_targets, test_targets, ori_mols_target, ori_mols_target2, all_mols_wonan_target, all_mols_wonan_target2, k=100)