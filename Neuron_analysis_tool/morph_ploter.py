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
from Neuron_analysis_tool.distance import Distance


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
    """
    get a rotation matrix in 2d space
    :param theta:
    :return:
    """
    theta = np.radians(theta)
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def split_point(prev_point, next_point, split_pos):
    """
    split the space between two points to match the segment length
    the dots are from reconstruction and don't have to alien with the segments location, here we take care of that
    by spliting the distance between to dots into the length of the segment
    :param prev_point:
    :param next_point:
    :param split_pos:
    :return:
    """
    direction_vec = next_point-prev_point
    return prev_point + direction_vec * split_pos

def get_parts_3d(sec, soma_loc, rotation_matrix = np.eye(3)):
    """
    take a section and split its x, y, z cordinated into list of cordinates per segment.
    :param sec: the section to split into segment paths
    :param soma_loc: the location of the soma (to move to (0,0,0)
    :param rotation_matrix: the 3d rotationmatrix to aply on the dots
    :return: list of dictinary of points there color and there attached segment
    """
    pos = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])-soma_loc
    pos = rotation_matrix @ pos
    res = []
    seg_len = sec.L/sec.nseg
    segs = list(sec)
    seg_lengths = [0] * len(segs)
    current_seg_idx = 0
    x = [pos[0]]
    y = [pos[1]]
    z = [pos[2]]
    l=0
    for i in range(1, sec.n3d()):
        next_pos = rotation_matrix @ (np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])-soma_loc)
        curent_len=np.linalg.norm(pos-next_pos)
        l+=curent_len
        if curent_len+seg_lengths[current_seg_idx]<=seg_len:
            x.append(next_pos[0])
            y.append(next_pos[1])
            z.append(next_pos[2])
            seg_lengths[current_seg_idx]+=curent_len
            if seg_lengths[current_seg_idx]>=seg_len:
                res.append(dict(x=np.array(x), y=np.array(y), z=np.array(z),seg=segs[current_seg_idx]))
                x = x[-1]
                y = y[-1]
                z = z[-1]
                current_seg_idx+=1
        else:
            # create_new pos in between
            while current_seg_idx < len(seg_lengths) and curent_len+seg_lengths[current_seg_idx]>seg_len:
                mid_points = split_point(pos, next_pos, (seg_len - seg_lengths[current_seg_idx]) / curent_len)
                seg_lengths[current_seg_idx] = seg_len

                x.append(mid_points[0])
                y.append(mid_points[1])
                z.append(mid_points[2])
                res.append(dict(x=np.array(x), y=np.array(y), z=np.array(z), seg=segs[current_seg_idx]))
                curent_len = np.linalg.norm(mid_points-next_pos)
                current_seg_idx+=1
                pos=mid_points

            if current_seg_idx < len(seg_lengths):
                seg_lengths[current_seg_idx] = curent_len
            elif curent_len>10:
                print('not using length of:', curent_len)
            x = [x[-1]]
            y = [y[-1]]
            z = [z[-1]]
            x.append(next_pos[0])
            y.append(next_pos[1])
            z.append(next_pos[2])
        pos = next_pos
    last_seg_in = False
    for data in res:
        if data['seg'] == segs[-1]:
            last_seg_in=True
    if not last_seg_in:
        res.append(dict(x=np.array(x), y=np.array(y), z=np.array(z), seg=segs[-1]))
    return res

def cell_to_points_3d(cell, rotation_matrix=np.eye(3), ignore_sections=[]):
    """
    take a cell and split its x, y, z cordinated into list of cordinates per segment.
    :param cell:
    :param rotation_matrix:
    :param rotation_matrix_2d:
    :param ignore_sections: sectins to exude from the plot
    :return:
    """
    number_of_soma_points = cell.soma[0].n3d()
    soma_loc = [cell.soma[0].x3d(number_of_soma_points // 2),
                cell.soma[0].y3d(number_of_soma_points // 2),
                cell.soma[0].z3d(number_of_soma_points // 2)]
    all_points = dict()
    for sec in cell.all:
        if sec in ignore_sections: continue
        try:
            sec_points = get_parts_3d(sec, soma_loc, rotation_matrix=rotation_matrix)
            all_points[sec]=sec_points
        except:
            all_points[sec] = []
            print('skiping section:', sec)
    return all_points


def get_parts(sec, color_func, soma_loc, rotation_matrix = np.eye(3), rotation_matrix_2d = np.eye(2)):
    """
    take a section and split its x, y, z cordinated into list of cordinates per segment.
    :param sec: the section to split into segment paths
    :param color_func: color func that gives color to each segment
    :param soma_loc: the location of the soma (to move to (0,0,0)
    :param rotation_matrix: the 3d rotationmatrix to aply on the dots
    :param rotation_matrix_2d: the 2d rotation matrix to aply on the dots
    :return: list of dictinary of points there color and there attached segment
    """
    pos = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])-soma_loc
    pos = rotation_matrix @ pos
    pos_2d = rotation_matrix_2d @ pos[:2]
    res = []
    seg_len = sec.L/sec.nseg
    segs = list(sec)
    seg_lengths = [0] * len(segs)
    current_seg_idx = 0
    x = [pos_2d[0]]
    y = [pos_2d[1]]
    l=0
    for i in range(1, sec.n3d()):
        next_pos = rotation_matrix @ (np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])-soma_loc)
        next_pos_2d = rotation_matrix_2d @ next_pos[:2]
        curent_len=np.linalg.norm(pos-next_pos)
        l+=curent_len
        if curent_len+seg_lengths[current_seg_idx]<=seg_len:
            x.append(next_pos_2d[0])
            y.append(next_pos_2d[1])
            seg_lengths[current_seg_idx]+=curent_len
            if seg_lengths[current_seg_idx]>=seg_len:
                c, d = color_func(segs[current_seg_idx])
                res.append(dict(x=np.array(x), y=np.array(y),seg=segs[current_seg_idx], color=c, diam=d))
                x = x[-1]
                y = y[-1]
                current_seg_idx+=1
        else:
            # create_new pos in between
            while current_seg_idx < len(seg_lengths) and curent_len+seg_lengths[current_seg_idx]>seg_len:
                mid_points = split_point(pos, next_pos, (seg_len - seg_lengths[current_seg_idx]) / curent_len)
                mid_points_2d = rotation_matrix_2d @ mid_points[:2]
                seg_lengths[current_seg_idx] = seg_len

                x.append(mid_points_2d[0])
                y.append(mid_points_2d[1])
                c, d = color_func(segs[current_seg_idx])
                res.append(dict(x=np.array(x), y=np.array(y), seg=segs[current_seg_idx], color=c, diam=d))
                curent_len = np.linalg.norm(mid_points-next_pos)
                current_seg_idx+=1
                pos=mid_points

            if current_seg_idx < len(seg_lengths):
                seg_lengths[current_seg_idx] = curent_len
            elif curent_len>10:
                print('not using length of:', curent_len)
            x = [x[-1]]
            y = [y[-1]]
            x.append(next_pos_2d[0])
            y.append(next_pos_2d[1])

            # current_seg_idx+=1
        pos = next_pos
        pos_2d = next_pos_2d
    last_seg_in = False
    for data in res:
        if data['seg'] == segs[-1]:
            last_seg_in=True
    if not last_seg_in :
        c, d = color_func(segs[- 1])
        res.append(dict(x=np.array(x), y=np.array(y), seg=segs[-1], color=c, diam=d))
    return res


def cell_to_points(cell, color_func, rotation_matrix=np.eye(3), rotation_matrix_2d=np.eye(2), ignore_sections=[]):
    """
    take a cell and split its x, y, z cordinated into list of cordinates per segment.
    :param cell:
    :param color_func:
    :param rotation_matrix:
    :param rotation_matrix_2d:
    :param ignore_sections: sectins to exude from the plot
    :return:
    """
    number_of_soma_points = cell.soma[0].n3d()
    soma_loc = [cell.soma[0].x3d(number_of_soma_points // 2),
                cell.soma[0].y3d(number_of_soma_points // 2),
                cell.soma[0].z3d(number_of_soma_points // 2)]
    all_points = dict()
    for sec in cell.all:
        if sec in ignore_sections: continue
        try:
            sec_points = get_parts(sec, color_func, soma_loc, rotation_matrix=rotation_matrix, rotation_matrix_2d=rotation_matrix_2d)
            all_points[sec]=sec_points
        except:
            all_points[sec] = []
            print('skiping section:', sec)
    return all_points


def plot(ax, all_points, add_nums=False, seg_to_indicate={}, counter=None, diam_factor=None, ignore_soma=False, diam_const=1, diam_min=0):
    """
    plot the morphology and return list of ploted lines and corisponding list of segment so line[i] corespond to seg[i], note that one segment can have multipl lines!!!
    :param ax:
    :param all_points:
    :param add_nums:
    :param seg_to_indicate:
    :param counter:
    :param diam_factor:
    :param ignore_soma:
    :return:
    """
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
                diam=diam_const
            else:
                diam*=diam_factor
            if diam<diam_min:
                diam=diam_min
            lines.append(ax.plot(x, y, color=c, linewidth=diam, zorder=1)[0])
            segs.append(seg)
            # prev_point = cur_point

            if seg in seg_to_indicate.keys():
                ax.scatter(np.mean(x), np.mean(y), s=seg_to_indicate[seg]['size'], color=seg_to_indicate[seg]['color'], zorder=2, alpha=seg_to_indicate[seg]['alpha'])

        if counter and counter.do_count(sec) and add_nums:
            num, color = counter.get_num_and_color(seg)
            txt = ax.text(x=x[1] + 1, y=y[1] + 1, s=str(num), color=color, fontsize=5)
            txt.set_path_effects([path_effects.withStroke(linewidth=1, foreground='k')])

    return ax, lines, segs

def electrical_move_helper(all_points, sec, start_point, distance):
    """
    helper function for electrical_move
    :param all_points:
    :param sec:
    :param start_point:
    :param distance:
    :return:
    """
    for seg in sec:
        for idx, sec_data in enumerate(all_points[sec]):
            if sec_data['seg'] == seg:
                dx = sec_data['x'][:-1] - sec_data['x'][1:]
                dy = sec_data['y'][:-1] - sec_data['y'][1:]
                current_lens = (dx**2 + dy**2)**0.5
                wanted_len = distance.get_length(seg, electrical=True)
                f = wanted_len/sum(current_lens)
                # fixing the shift
                all_points[sec][idx]['x'] += start_point['x']-sec_data['x'][0]
                all_points[sec][idx]['y'] += start_point['y']-sec_data['y'][0]

                # fixing the length
                for point_idx in range(1, len(sec_data['x']), 1):
                    new_x = sec_data['x'][point_idx-1] + (sec_data['x'][point_idx]-sec_data['x'][point_idx-1])*f
                    new_y = sec_data['y'][point_idx-1] + (sec_data['y'][point_idx]-sec_data['y'][point_idx-1])*f
                    all_points[sec][idx]['x'][point_idx:] += (new_x-sec_data['x'][point_idx])
                    all_points[sec][idx]['y'][point_idx:] += (new_y-sec_data['y'][point_idx])
                start_point = dict(x=sec_data['x'][-1], y=sec_data['y'][-1])
    for son in sec.children():
        all_points = electrical_move_helper(all_points, son, start_point, distance)
    return all_points


def electrical_move(all_points, cell, distance):
    """
    move the morphology dots acording to there electrical property's
        (we move them so there distance in um in the plot will be there electrical distance in lamda)
    :param all_points:
    :param cell:
    :param distance:
    :return:
    """
    all_points = electrical_move_helper(all_points, cell.soma[0], dict(x=0, y=0), distance)
    return all_points

def electrical_move_helper_3d(all_points, sec, start_point, distance):
    """
    helper function for electrical_move
    :param all_points:
    :param sec:
    :param start_point:
    :param distance:
    :return:
    """
    for seg in sec:
        if sec not in all_points: continue
        for idx, sec_data in enumerate(all_points[sec]):
            if sec_data['seg'] == seg:
                dx = sec_data['x'][:-1] - sec_data['x'][1:]
                dy = sec_data['y'][:-1] - sec_data['y'][1:]
                dz = sec_data['z'][:-1] - sec_data['z'][1:]
                current_lens = (dx**2 + dy**2 + dz**2)**0.5
                wanted_len = distance.get_length(seg, electrical=True)
                f = wanted_len/sum(current_lens)
                # fixing the shift
                all_points[sec][idx]['x'] += start_point['x']-sec_data['x'][0]
                all_points[sec][idx]['y'] += start_point['y']-sec_data['y'][0]
                all_points[sec][idx]['z'] += start_point['z']-sec_data['z'][0]

                # fixing the length
                for point_idx in range(1, len(sec_data['x']), 1):
                    new_x = sec_data['x'][point_idx-1] + (sec_data['x'][point_idx]-sec_data['x'][point_idx-1])*f
                    new_y = sec_data['y'][point_idx-1] + (sec_data['y'][point_idx]-sec_data['y'][point_idx-1])*f
                    new_z = sec_data['z'][point_idx-1] + (sec_data['z'][point_idx]-sec_data['z'][point_idx-1])*f
                    all_points[sec][idx]['x'][point_idx:] += (new_x-sec_data['x'][point_idx])
                    all_points[sec][idx]['y'][point_idx:] += (new_y-sec_data['y'][point_idx])
                    all_points[sec][idx]['z'][point_idx:] += (new_y-sec_data['z'][point_idx])
                start_point = dict(x=sec_data['x'][-1], y=sec_data['y'][-1], z=sec_data['z'][-1])
    for son in sec.children():
        all_points = electrical_move_helper_3d(all_points, son, start_point, distance)
    return all_points


def electrical_move_3d(all_points, cell, distance):
    """
    move the morphology dots acording to there electrical property's
        (we move them so there distance in um in the plot will be there electrical distance in lamda)
    :param all_points:
    :param cell:
    :param distance:
    :return:
    """
    all_points = electrical_move_helper_3d(all_points, cell.soma[0], dict(x=0, y=0, z=0), distance)
    return all_points

def get_dots(cell, electrical=False, distance=None, more_conductances=None, time=None, dt=1, dt_func= lambda x: np.mean(x)):

    if electrical and distance is None:
        distance = Distance(cell, more_conductances, dt_func=dt_func)
        distance.compute(time=time, dt=dt)
    if electrical:
        assert (more_conductances is not None) or (distance is not None)
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

    all_points = cell_to_points_3d(cell, rotation_matrix=rotation_matrix, ignore_sections=[])

    if electrical:
        all_points = electrical_move_3d(all_points, cell, distance)
    return all_points

def plot_morph(cell, color_func, add_nums=False, seg_to_indicate={},
               counter=None, ax=None,
               sec_to_change =None, diam_factor=None, theta=0,
               ignore_sections=[], ignore_soma=False,
               electrical=False, distance=None, more_conductances=None, time=None, dt=1,
               dt_func= lambda x: np.mean(x), diam_const=1, diam_min=0):
    """
    plot the morphology
    :param cell: the cell to plot
    :param color_func: function that get a segment and return its color and part name
    :param seg_to_indicate: dictinary of {seg {color, size, alpha} to scatter on the segment location
    :param add_nums: if to add numbers to the acording to the counter
    :param counter: class that have a function get_num_and_color that get a section and return num, color for the text to add,
    and aa function do_count that get a section and retun if he is part off the counting (if to add text next to it or not)
    :param ax: axis to plot on
    :param sec_to_change: section to change there location from one parent seg ot ather arent seg {sec: {from_: seg, to: seg}}
    :param diam_factor: factor for the diameter, if None that consatnt diameter is used for all the tree
    :param theta: rotation angle in 2d
    :param ignore_sections: section not to plot
    :param ignore_soma: if True the soma is ploted as a circle
    :param electrical: if True, the morphology will be change to show distance in electical units
    :param distance: distance of the cell (from a givin POV)
    :param more_conductances: the conductances to use if changing into electrical units
    :param time: time to take the more_conductances from
    :param dt: dt for the more_conductances
    :param dt_func: dt_fuction for the more_conductances
    :return:
    """
    if electrical and distance is None:
        distance = Distance(cell, more_conductances, dt_func=dt_func)
        distance.compute(time=time, dt=dt)
    if electrical:
        assert (more_conductances is not None) or (distance is not None)
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
        ax = plt.gca()
    if ignore_soma:
        sec=cell.soma[0]
        soma_diam = cell.soma[0].diam * (diam_factor if diam_factor else diam_const)
        soma_length = cell.soma[0].L * (diam_factor if diam_factor else diam_const)
        if soma_diam< diam_min:
            soma_diam=diam_min
        if soma_length< diam_min:
            soma_length=diam_min
        soma_x = abs(all_points[sec][0]['x'][0] - all_points[sec][-1]['x'][-1])
        soma_y = abs(all_points[sec][0]['y'][0] - all_points[sec][-1]['y'][-1])
        soma_angle = np.rad2deg(np.arctan(soma_y/soma_x))

    if electrical:
        soma_diam = soma_length = sum([distance.get_length(seg, electrical=True) for seg in cell.soma[0]])*2*(diam_factor if diam_factor else diam_const)
        if soma_diam< diam_min:
            soma_diam=diam_min
        all_points = electrical_move(all_points, cell, distance)
    ax, lines, segs = plot(ax, all_points, add_nums=add_nums, seg_to_indicate=seg_to_indicate,counter=counter, diam_factor=diam_factor, ignore_soma=ignore_soma, diam_const=diam_const, diam_min=diam_min)
    if ignore_soma:
        segs.append(list(cell.soma[0])[len(list(cell.soma[0]))//2])
        c, _ = color_func(segs[-1])
        lines.append(ax.add_patch(mpl.patches.Ellipse((0, 0), soma_length*2, soma_diam*2, soma_angle, color=c, clip_on=False, zorder=2)))
        for seg in seg_to_indicate.keys():
            if seg.sec == cell.soma[0]:
                ax.scatter(0, 0, s=seg_to_indicate[seg]['size'], color=seg_to_indicate[seg]['color'], zorder=3, alpha=seg_to_indicate[seg]['alpha'])
                break
    return ax, lines, segs