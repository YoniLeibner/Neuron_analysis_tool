###########################################################
#
# author: Yoni Leibner
# description: plot the morphology in 2d based on the
#              samplaed points of the neuron and color each
#              segment with the color func class
# date of modification: 16.11.2022
#
###########################################################

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as path_effects
import matplotlib as mpl
from tqdm import tqdm


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    else:
        return np.eye(3) #cross of all zeros only occurs on identical directions

def get_rotation_matrix_2d(theta):
    theta = np.radians(theta)
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def split_point(prev_point, next_point, split_pos):
    direction_vec = prev_point - next_point
    return prev_point + direction_vec * split_pos

def get_parts(sec, color_func, soma_loc, rotation_matrix = np.eye(3), rotation_matrix_2d = np.eye(2)):
    pos = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])-soma_loc
    pos = rotation_matrix @ pos
    pos_2d = rotation_matrix_2d @ pos[:2]
    res = []
    accumulative_len = 0
    seg_len = sec.L/sec.nseg
    segs = list(sec)
    for i in range(1, sec.n3d()):
        next_pos = rotation_matrix @ (np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])-soma_loc)
        next_pos_2d = rotation_matrix_2d @ next_pos[:2]
        curent_len=np.linalg.norm(pos-next_pos)

        if accumulative_len//seg_len == (accumulative_len+curent_len)//seg_len: # this part is only on one segment
            c, d = color_func(segs[int(accumulative_len//seg_len)])
            res.append(dict(x=[pos_2d[0], next_pos_2d[0]], y=[pos_2d[1], next_pos_2d[1]], seg=segs[int(accumulative_len//seg_len)], color=c, diam=d))
        else: # we need to split the points linearly

            missing_len = seg_len - accumulative_len % seg_len
            split_pos = missing_len / curent_len
            mid_points = split_point(pos, next_pos, split_pos)
            mid_points_2d = rotation_matrix_2d @ mid_points[:2]
            c1, d1 = color_func(segs[int(accumulative_len // seg_len)])
            c2, d2 = color_func(segs[min(int((accumulative_len+curent_len) // seg_len), len(segs)-1)])
            res.append(dict(x=[pos_2d[0], mid_points_2d[0]], y=[pos_2d[1], mid_points_2d[1]], seg=segs[int(accumulative_len // seg_len)], color=c1, diam=d1))
            res.append(dict(x=[mid_points_2d[0], next_pos_2d[0]], y=[mid_points_2d[1], next_pos_2d[1]], seg=segs[min(int((accumulative_len+curent_len)//seg_len), len(segs)-1)], color=c2, diam=d2))

        accumulative_len+=curent_len
        pos = next_pos
        pos_2d = next_pos_2d
    return res

def get_point_segs(sec, soma_loc, color_func, rotation_matrix=np.eye(3), rotation_matrix_2d=np.eye(2)):
    try:
        parts = get_parts(sec, color_func, soma_loc, rotation_matrix=rotation_matrix, rotation_matrix_2d=rotation_matrix_2d)
    except:
        return []
    return parts

def cell_to_points(cell, color_func, rotation_matrix=np.eye(3), rotation_matrix_2d=np.eye(2), ignore_sections=[]):
    number_of_soma_points = cell.soma[0].n3d()
    soma_loc = [cell.soma[0].x3d(number_of_soma_points // 2),
                cell.soma[0].y3d(number_of_soma_points // 2),
                cell.soma[0].z3d(number_of_soma_points // 2)]
    all_points = dict()
    for sec in cell.all:
        if sec in ignore_sections: continue
        sec_points = get_point_segs(sec, soma_loc, color_func=color_func, rotation_matrix=rotation_matrix, rotation_matrix_2d=rotation_matrix_2d)
        all_points[sec]=sec_points
    return all_points


def plot(ax, all_points, add_nums=False, seg_to_indicate={}, counter=None, diam_factor=None, ignore_soma=False):
    lines=[]
    segs = []
    for sec in all_points:
        if ignore_soma and sec in sec.cell().soma:
            continue
        for i in all_points[sec]:
            x=i['x']
            y=i['y']
            c = i['color']
            # d = i['diam']
            seg = i['seg']
            diam = seg.diam
            if diam_factor is None:
                diam=1
            else:
                diam*=diam_factor
            lines.append(ax.plot(x, y, color=c, linewidth=diam, zorder=1)[0])
            segs.append(seg)
            # prev_point = cur_point

            if seg in seg_to_indicate.keys():
                ax.scatter(np.mean(x), np.mean(y), s=seg_to_indicate[seg]['size'], color=seg_to_indicate[seg]['color'], zorder=2, alpha=seg_to_indicate[seg]['alpha'])

        if counter and counter.do_count(sec) and add_nums:
            num, color = counter.get_num_and_color(seg)
            txt = ax.text(x=x[1] + 1, y=y[1] + 1, z=z[1] + 1, s=str(num), color=color, fontsize=5)
            txt.set_path_effects([path_effects.withStroke(linewidth=1, foreground='k')])

    return ax, lines, segs

def get_norm(all_vals):
    norm = mpl.colors.Normalize(vmin=np.min(all_vals), vmax=np.max(all_vals))
    return norm

def plot_morph(cell, color_func, scatter=False, add_nums=False, seg_to_indicate={},
               counter=None, fig=None,  ax=None,
               sec_to_change =None, diam_factor=None, plot_color_bar=True, theta=0, ignore_sections=[], ignore_soma=False, color_bar_idx = [0.9, 0.2, 0.02, 0.6]): #cmap = plt.cm.coolwarm, norm_colors=True,


    all_points_arr = []
    for sec in cell.all:
        for i in range(sec.n3d()):
            all_points_arr.append([sec.x3d(i), sec.y3d(i), sec.z3d(i)])

    from sklearn.decomposition import PCA
    X = np.array(all_points_arr)
    pca = PCA()
    pca.fit(X)
    V0 = pca.components_[0]
    V1 = pca.components_[1]
    rotation_matrix = rotation_matrix_from_vectors(V0, V1)
    rotation_matrix_2d = get_rotation_matrix_2d(theta)

    all_points = cell_to_points(cell, color_func=color_func, rotation_matrix=rotation_matrix, rotation_matrix_2d=rotation_matrix_2d, ignore_sections=ignore_sections)
    if sec_to_change is not None:
        for sec in all_points.keys():
            if sec in sec_to_change.keys():
                from_seg = sec_to_change[sec]['from_']
                from_dict = all_points[from_seg.sec][int((len(all_points[from_seg.sec])-1)*from_seg.x)]
                from_vals = dict(x=np.array(from_dict['x']).mean(),
                                 y=np.array(from_dict['y']).mean())
                to_seg = sec_to_change[sec]['to']
                to_dict = all_points[to_seg.sec][int((len(all_points[to_seg.sec]) - 1) * to_seg.x)]
                to_vals = dict(x=np.array(to_dict['x']).mean(),
                               y=np.array(to_dict['y']).mean())
                for i in range(len(all_points[sec])):
                    all_points[sec][i]['x'] = all_points[sec][i]['x'] - from_vals['x'] + to_vals['x']
                    all_points[sec][i]['y'] = all_points[sec][i]['y'] - from_vals['y'] + to_vals['y']


    if ax is None:
        fig = plt.figure(figsize=(20, 20))
        ax = plt.axes()
    all_points2 = dict()
    for sec in tqdm(all_points, desc='optimizing lines'):
        if sec in cell.soma:
            pass
        all_points2[sec] = list()
        if len(all_points[sec])>0:
            x=all_points[sec][0]['x']
            y=all_points[sec][0]['y']
            color = all_points[sec][0]['color']
            current_seg = all_points[sec][0]['seg']
            for point in all_points[sec]:
                if point['seg'] == current_seg:
                    assert point['color'] == color
                    x.append(point['x'][1])
                    y.append(point['y'][1])
                else:
                    all_points2[sec].append(dict(x=x, y=y, seg=current_seg, color=color))
                    x = point['x']
                    y = point['y']
                    color = point['color']
                    current_seg = point['seg']
            all_points2[sec].append(dict(x=x, y=y, seg=current_seg, color=color))
    if ignore_soma:
        sec=cell.soma[0]
        soma_diam = cell.soma[0].diam * (diam_factor if diam_factor else 1)
        soma_length = cell.soma[0].L * (diam_factor if diam_factor else 1)
        soma_x = abs(all_points[sec][0]['x'][0] - all_points[sec][-1]['x'][-1])
        soma_y = abs(all_points[sec][0]['y'][0] - all_points[sec][-1]['y'][-1])
        soma_angle = np.rad2deg(np.arctan(soma_y/soma_x))
    all_points=all_points2
    ax, lines, segs = plot(ax, all_points, add_nums=add_nums, seg_to_indicate=seg_to_indicate,counter=counter, diam_factor=diam_factor, ignore_soma=ignore_soma)
    if ignore_soma:
        segs.append(list(cell.soma[0])[len(list(cell.soma[0]))//2])
        c, _ = color_func(segs[-1])
        lines.append(ax.add_patch(mpl.patches.Ellipse((0, 0), soma_length*2, soma_diam*2, soma_angle, color=c, clip_on=False, zorder=2)))
        for seg in seg_to_indicate.keys():
            if seg.sec == cell.soma[0]:
                ax.scatter(0, 0, s=seg_to_indicate[seg]['size'], color=seg_to_indicate[seg]['color'], zorder=3, alpha=seg_to_indicate[seg]['alpha'])
                break
        # ax.add_patch(lines[-1])
    # todo change x, y values to electrical units

    return ax, lines, segs