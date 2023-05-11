#############################################################
#
# author: Yoni Leibner
# description: calculate and plot cables from a
#              givin start seg
#              cable options are:
#                       d3_2 - equivalent cable
#                       electric - dA/dX cable
#                       dist - dA/dx cable
# date of modification: 11.05.2023
#
#############################################################

import numpy as np
from Neuron_analysis_tool.distance import Distance
import matplotlib.pyplot as plt
import numbers


def get_part(part_dict, seg):
    """
    get the direction of the segment (in order to get all the cables based on there parts, basal, oblique, tuft, etc.)
    :param part_dict: the part dict as {part name: list of segments}
    :param seg: the current seg
    :return: the
    """
    for part in part_dict.keys():
        if seg in part_dict[part]:
            return part
    raise Exception('seg have no part, '+str(seg)+', '+str(seg.sec))#this shold not bereached


def get_cable(cell, start_seg, factor_e_space=100, factor_m_space=10, more_conductances={}, seg_dist_dict={},
              part_dict=dict(), ignore_sections = [], distance=None, dt_func=lambda x: np.mean(x)):
    """
    calculate the cables.
    for each direction in 'sons', 'parent' and for each part from part dict and all
    :param cell: the cell neuron model
    :param start_seg: the seg to start the cable from
    :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
    :param factor_m_space: the bin size in physical units, every factor_m_space um is a bin
    :param more_conductances: more_conductances to initiate the distance if distance is None
    :param seg_dist_dict: dict in the form {direction: {seg: []}} (segs to get there electric distane
    :param part_dict: the part dict as {part name: list of segments}
    :param ignore_sections: section to skip in the cable
    :param distance: a distance for the segments
                    (if not givin its uses the default distance from the init more_conductances function)
    :param dt_func: function for dt in the more_conductances
    :return:
    """
    if (distance is None) or (not distance.start_seg == start_seg):
        distance = Distance(cell, more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)

    cross_dist_dict = dict(sons=[], parent=[])
    total_res = dict()
    for direction in seg_dist_dict:
        for seg in seg_dist_dict[direction]:
            if distance.get_direction(seg) == direction:
                distance.get_mid_point(seg, electrical=True)
                seg_dist_dict[direction][seg].append(dict(dist_m=distance.get_mid_point(seg, electrical=False),
                                                          dist_e=distance.get_mid_point(seg, electrical=True)))

    for direction in ['sons', 'parent']:
        total_res[direction] = dict()
        total_res[direction]['all'] = dict(dist=np.zeros(1000), electric= np.zeros(1000), d3_2 = np.zeros(1000))
        for part in part_dict.keys():
            total_res[direction][part] = dict(dist=np.zeros(1000), electric= np.zeros(1000), d3_2 = np.zeros(1000))
            total_res[direction][part]['dist'][0] = total_res[direction][part]['electric'][0] = start_seg.area()
    e_threshs = np.arange(0, 1000, 1)/factor_e_space
    m_threshs = np.arange(0, 1000, 1)*factor_m_space

    for sec in distance.distance_dict:
        if sec in ignore_sections: continue
        for seg in distance.distance_dict[sec]['segs']:
            e_dist = distance.get_mid_point(seg, electrical=True)
            m_dist = distance.get_mid_point(seg, electrical=False)
            m_idx = np.where(m_threshs<=m_dist)[0][-1]
            e_idx = np.where(e_threshs<=e_dist)[0][-1]
            direction = distance.get_direction(seg)
            start_end = distance.get_start_end(seg, electrical=True)

            part = get_part(part_dict, seg)
            total_res[direction]['all']['dist'][m_idx] += seg.area()
            total_res[direction][part]['dist'][m_idx] += seg.area()
            total_res[direction]['all']['electric'][e_idx]+=seg.area()
            total_res[direction][part]['electric'][e_idx]+=seg.area()
            total_res[direction]['all']['d3_2'][np.logical_and(e_threshs >= start_end['start'], e_threshs < start_end['end'])] += seg.diam ** 1.5
            total_res[direction][part]['d3_2'][np.logical_and(e_threshs >= start_end['start'], e_threshs < start_end['end'])] += seg.diam ** 1.5


    for direction in total_res.keys():
        for part in total_res[direction].keys():
            total_res[direction][part]['d3_2'] = np.power(total_res[direction][part]['d3_2'], 2.0/3.0)
    return total_res, seg_dist_dict, dict()


def plot_cable(cell, start_seg, ax=None, cable_type='d3_2', factor_e_space=25, factor_m_space=10,
               more_conductances=dict(), part_dict=dict(), colors_dict=dict(), ignore_sections = [], distance=None,
               segs_to_indecate=dict(), start_loc=0, vertical=True, dots_size=10, start_color='k', shift=None,
               plot_legend=True, cable_factor=1, labal_start=None, return_shift=False, dt_func=lambda x: np.mean(x)):
    """
    calculate and plot a cable.
    :param cell: the cell neuron model
    :param start_seg: the seg to start the cable from
    :param ax: the axes to plot on
    :param cable_type: the cable type to plot from:
                            d3_2 - equivalent cable
                            electric - dA/dX cable
                            dist - dA/dx cable
    :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
    :param factor_m_space: the bin size in physical units, every factor_m_space um is a bin
    :param more_conductances: more_conductances to initiate the distance if distance is None
    :param part_dict: the part dict as {part name: list of segments}
    :param colors_dict: the color dict as {part name: color}
    :param ignore_sections: section to skip in the cable
    :param distance: a distance for the segments
                    (if not givin its uses the default distance from the init more_conductances function)
    :param segs_to_indecate: dict in the form {direction: {seg: []}} (segs to get there electric distane
    :param start_loc: the start location on the axes
    :param vertical: if to plot the cable vertically or horizontally
    :param dots_size: the size of the dot at the origin
    :param start_color: the color of the dot in the origin
    :param shift: a constant shift or a chosen from the data
    :param plot_legend:  if True it will dd a legend with the parts_dict keys and colors
    :param cable_factor: factor to add to the cable to match the axes
    :param labal_start: the label for the legend of the start seg
    :param return_shift: if return the result shift
    :param dt_func: function for dt in the more_conductances
    :return: the axes ploted on and shift if return_shift is True
    """
    if ax is None:
        ax = plt.gca()
    if start_seg is None:
        start_seg = cell.soma[0](0.5)

    seg_dist_dict = dict()
    for direction in ['sons', 'parent']:
        seg_dist_dict[direction] = dict()
        for seg_ in segs_to_indecate.keys():
            seg_dist_dict[direction][seg_] = []
    results, seg_dist, cross_dist_dict = get_cable(cell,
                                                   start_seg=start_seg,
                                                   factor_e_space=factor_e_space,
                                                   factor_m_space=factor_m_space,
                                                   more_conductances=more_conductances,
                                                   seg_dist_dict=seg_dist_dict,
                                                   part_dict=part_dict,
                                                   ignore_sections=ignore_sections,
                                                   distance=distance,
                                                   dt_func=dt_func)
    max_cable = 0
    for direction in results:
        for part in results[direction]:
            for cable_type in results[direction][part]:
                results[direction][part][cable_type] *= cable_factor

    if not isinstance(shift, numbers.Number):
        shift = 0

        for part, direction in zip(results.keys(), [1, -1]):
            cable = results[part]['all'][cable_type].flatten()
            if max_cable < cable.max() / 2.0:
                max_cable = cable.max() / 2.0
                shift = start_loc + cable.max() / 2.0

    for part, direction in zip(results, [1, -1]):
        cable = results[part]['all'][cable_type].flatten()
        if cable_type == 'd3_2':
            befor_d_3_2 = np.power(cable, 3.0 / 2.0)
        y = np.arange(0, len(cable), 1) / factor_e_space

        start_pos = -cable / 2.0 + shift
        for morph_part in part_dict.keys():
            remove_start_diam = False
            part_cable = results[part][morph_part][cable_type].flatten()
            if cable_type == 'd3_2':
                part_cable = results[part][morph_part][cable_type].flatten()
                part_cable_befor_d_3_2 = np.power(part_cable, 3.0 / 2.0)
                part_cable = cable * (part_cable_befor_d_3_2 / befor_d_3_2)
            if part_cable[1] > 0 and part_cable[0] == 0:
                part_cable[0] = (results[part]['all'][cable_type].flatten())[0]
                remove_start_diam = True
                temp_ = start_pos[0]
                start_pos[0] = -cable[0] / 2.0 + shift
            plot_cable = part_cable[part_cable > 0]
            start_pos_temp = start_pos[part_cable > 0]
            if vertical:
                ax.fill_betweenx(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable,
                                 label=morph_part, color=colors_dict[morph_part])
            else:
                ax.fill_between(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable,
                                label=morph_part, color=colors_dict[morph_part])
            if remove_start_diam:
                part_cable[0] = 0
                start_pos[0] = temp_
            start_pos += part_cable
    if vertical:
        for seg_ in segs_to_indecate:
            if len(seg_dist_dict['parent'][seg_]) > 0:
                dist = -seg_dist_dict['parent'][seg_][0]['dist_e'] if cable_type in ['electric', 'd3_2'] else \
                seg_dist_dict['parent'][seg_][0]['dist_m']
                ax.scatter(shift, dist, s=segs_to_indecate[seg_]['size'], color=segs_to_indecate[seg_]['color'])
            else:
                dist = seg_dist_dict['sons'][seg_][0]['dist_e'] if cable_type in ['electric', 'd3_2'] else \
                seg_dist_dict['sons'][seg_][0]['dist_m']
                ax.scatter(shift, dist, s=segs_to_indecate[seg_]['size'], color=segs_to_indecate[seg_]['color'])
        if labal_start is None:
            ax.scatter(shift, 0, s=dots_size, color=start_color, edgecolor='k', linewidths=2)
        else:
            ax.scatter(shift, 0, s=dots_size, color=start_color, label=labal_start, edgecolor='k', linewidths=2)
    else:
        for seg_ in segs_to_indecate:
            if len(seg_dist_dict['parent'][seg_]) > 0:
                dist = -seg_dist_dict['parent'][seg_][0]['dist_e'] if cable_type in ['electric', 'd3_2'] else \
                seg_dist_dict['parent'][seg_][0]['dist_m']
                ax.scatter(dist, shift, s=segs_to_indecate[seg_]['size'], color=segs_to_indecate[seg_]['color'])
            else:
                dist = seg_dist_dict['sons'][seg_][0]['dist_e'] if cable_type in ['electric', 'd3_2'] else \
                seg_dist_dict['sons'][seg_][0]['dist_m']
                ax.scatter(dist, shift, s=segs_to_indecate[seg_]['size'], color=segs_to_indecate[seg_]['color'])
        if labal_start is None:
            ax.scatter(0, shift, s=dots_size, color=start_color, edgecolor='k', linewidths=2)
        else:
            ax.scatter(0, shift, s=dots_size, color=start_color, label=labal_start, edgecolor='k', linewidths=2)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if plot_legend:
        ax.legend(by_label.values(), by_label.keys(), loc='best')
    if return_shift:
        return ax, shift
    return ax
