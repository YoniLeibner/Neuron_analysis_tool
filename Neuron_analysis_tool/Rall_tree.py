###########################################################
#
# author: Yoni Leibner
# description: script that creates swc Rall models with
#              depth as spasafized in hyper parameter depth
#              the swc will be saved to the data folder
#              with the name: Rall_tree<depth>.swc
# date of modification: 16.11.2022
#
###########################################################

import numpy as np
from copy import deepcopy

EPS = 1e-3
depth = 5
spread_factor = 30.0
EQ_DIAM=10.0


def section(f, start, end, diam, num, prev_num, type=3, n_samples=10):
    assert n_samples > 0
    # f.write(str(num)+' '+str(type)+' '+' '.join(np.array([round(start['x'], 4), round(start['y'], 4), round(start['z'], 4), round(diam/2.0, 4)]).astype(str))+' '+str(prev_num)+'\n')
    # prev_num=num
    # num += 1
    for x, y, z in zip(np.linspace(start['x'], end['x'], n_samples), np.linspace(start['y'], end['y'], n_samples),
                       np.linspace(start['z'], end['z'], n_samples)):
        f.write(str(num) + ' ' + str(type) + ' ' + ' '.join(np.array([round(x, 4), round(y, 4), round(z, 4), round(diam / 2.0, 4)]).astype(str)) + ' ' + str(prev_num) + '\n')
        prev_num = num
        num += 1
    point = {'x': x, 'y': y, 'z': z}
    return point, num


def create_tree(max_depth=3, Rm=10000.0, Ra=100.0, E_length=0.25, n_samples=10):
    f = open('data/Rall_tree' + str(max_depth) + '.swc', 'w')
    f.write('# created by Yoni Leibner, Rall tree depth=' + str(max_depth) + ' , E_length=' + str(E_length) + '\n')
    f.write('# this need to have: Rm=' + str(Rm) + ', Ra=' + str(Ra) + '\n')
    f.write('# id, type, x, y, z, diam, prev\n')
    helping(f, 0, max_depth=max_depth, start=dict(x=0, y=0, z=0), num_to_create=1, diam=EQ_DIAM, Rm=Rm, Ra=Ra,
            n_samples=n_samples, E_length=E_length)
    f.close()


def get_length(Rm, Ra, diam, E_length=0.25):
    lamda = ((Rm / Ra) * (diam / 10000.0 / 4.0)) ** 0.5
    return lamda * E_length * 10000.0


def helping(f, depth, max_depth=3, start={'x': 0, 'y': 0, 'z': 0}, num_to_create=1, diam=1.0, Rm=10000.0, Ra=100.0,
            n_samples=10, E_length=0.25, num=0, prev_num=-1, types=[3, 4]):
    if depth >= max_depth: return num
    l = get_length(Rm, Ra, diam, E_length=0.25)
    if num_to_create == 1:
        end = {'x': start['x'] + l, 'y': start['y'], 'z': start['z']}
        start['z'] += EPS
        point, num = section(f, start, end, diam, num + 1, prev_num, type=1, n_samples=n_samples)
        num = helping(f, depth + 1, max_depth=max_depth, start=point, num_to_create=2,
                      diam=((diam ** 1.5) / 2.0) ** (2.0 / 3.0), Rm=Rm, Ra=Ra, n_samples=n_samples,
                      E_length=E_length, num=num - 1, prev_num=num - 1, types=[4, 4])
    else:
        elevation = 2**(max_depth - depth-1) / 2.0 * spread_factor
        l2 = (l ** 2 - elevation ** 2) ** 0.5
        end1 = {'x': start['x'] + l2, 'y': start['y'] + elevation, 'z': start['z']}
        start['z'] += EPS
        point, num = section(f, deepcopy(start), deepcopy(end1), diam, num + 1, prev_num, type=types[0], n_samples=n_samples)
        num = helping(f, depth + 1, max_depth=max_depth, start=point, num_to_create=2,
                      diam=((diam ** 1.5) / 2.0) ** (2.0 / 3.0), Rm=Rm, Ra=Ra, n_samples=n_samples,
                      E_length=E_length, num=num - 1, prev_num=num - 1, types=[types[0], types[0]])
        end2 = {'x': start['x'] + l2, 'y': start['y'] - elevation, 'z': start['z']}
        start['z'] -= 2 * EPS
        point, num = section(f, deepcopy(start), deepcopy(end2), diam, num + 1, prev_num, type=types[1], n_samples=n_samples)
        num = helping(f, depth + 1, max_depth=max_depth, start=point, num_to_create=2,
                      diam=((diam ** 1.5) / 2.0) ** (2.0 / 3.0), Rm=Rm, Ra=Ra, n_samples=n_samples,
                      E_length=E_length, num=num - 1, prev_num=num - 1, types=[types[1], types[1]])
    return num


create_tree(max_depth=depth, Rm=10000.0, Ra=100.0, E_length=0.25, n_samples=10)
