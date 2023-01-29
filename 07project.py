import random
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import pickle
import pandas as pd
import torch

def load_pickle(path):
    with open(path,'rb') as f:
        return pickle.load(f)

scaled_data=load_pickle('scaled_data.pkl')
cubes = torch.zeros(3000,2,100,100,100)
def cubes_building(cubes, scaled_data, n):
    discard_num = [0, 0]
    overlap_num = [0, 0]
    # protein:
    for i in range(n):
        # protein:
        for j in range(len(scaled_data['proteins'][i][0][0])):
            x_p, y_p, z_p = scaled_data['proteins'][i][0][0][j], scaled_data['proteins'][i][0][1][j], \
                            scaled_data['proteins'][i][0][2][j]
            # Shift all the atoms by 25 such that the center could be 25,25,25
            # and convert the data type of (x,y,z) into int
            x_p, y_p, z_p = int(round(x_p)) + 50, int(round(y_p)) + 50, int(round(z_p)) + 50
            # crop the cube size with 100*100*100
            if x_p >= 100 or y_p >= 100 or z_p >= 100:
                discard_num[0] += 1
                continue
            # check the overlaping atoms
            if cubes[i, 0, x_p, y_p, z_p] != 0 or cubes[i, 1, x_p, y_p, z_p] != 0:
                overlap_num[0] += 1
                continue
            # first channel
            cubes[i, 0, x_p, y_p, z_p] = -1
            # second channel
            if scaled_data['proteins'][i][1][j] == 'h':
                cubes[i, 1, x_p, y_p, z_p] = 1
            else:
                cubes[i, 1, x_p, y_p, z_p] = -1

            # ligands:
        for k in range(len(scaled_data['ligands'][i][0][0])):
            x_l, y_l, z_l = scaled_data['ligands'][i][0][0][k], scaled_data['ligands'][i][0][1][k], \
                            scaled_data['ligands'][i][0][2][k]
            x_l, y_l, z_l = int(round(x_l)) + 50, int(round(y_l)) + 50, int(round(z_l)) + 50
            if x_l >= 100 or y_l >= 100 or z_l >= 100:
                discard_num[1] += 1
                continue
            if cubes[i, 0, x_l, y_l, z_l] != 0 or cubes[i, 1, x_l, y_l, z_l] != 0:
                overlap_num[1] += 1
                continue
            # first channel
            cubes[i, 0, x_l, y_l, z_l] = 1
            # second channel
            if scaled_data['ligands'][i][1][k] == 'h':
                cubes[i, 1, x_l, y_l, z_l] = 1
            else:
                cubes[i, 1, x_l, y_l, z_l] = -1
    return cubes, discard_num, overlap_num



