#############################################################
#
# author: Yoni Leibner
# description: plot the dendogram of the cell in micro-meters
#              or electricl units
# date of modification: 16.11.2022
#
#############################################################

import numpy as np
from Neuron_analysis_tool.distance import Distance
from neuron import h
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Neuron_analysis_tool.utils import LAMDA, MICRO

FIX_DIAM=1
BRANCH_SPACE=2


def plot_recursive(sec, distance, ax, color_func, x_pos, lines, segs, ignore_sections=[], segs_to_indecate=dict(), electrical=True, diam_factor=None, BRANCH_SPACE_=BRANCH_SPACE):
    # start_end = distance.get_sec_start_end(sec, electrical=electrical)
    [color, part_name]=color_func.get_seg_color(list(sec)[0])
    # colors, lengths = color_func(sec, distance.get_sec_parent(sec).sec) # change to seg colors!!!!!!!
    mid_points = [x_pos]
    mean_points=[]
    sons = distance.get_sons(sec)
    if not distance.is_terminal(sec):
        for son in sons:
            if son not in ignore_sections:
                pos, mean_point, lines, segs = plot_recursive(son, distance, ax, color_func, mid_points[-1],lines, segs,
                                                              ignore_sections=ignore_sections, segs_to_indecate=segs_to_indecate,
                                                              electrical=electrical, diam_factor=diam_factor, BRANCH_SPACE_=BRANCH_SPACE_)
                mid_points.append(pos)
                mean_points.append(mean_point)
        if len(sons)>1:
            m_point = (mid_points[0]+ mid_points[-1]-BRANCH_SPACE_)/2.0
        else:
            m_point = mid_points[-1] - BRANCH_SPACE_ # single cheldren
    else:
        m_point = x_pos # terminal

    mul = 1 if distance.get_part(list(sec)[0]) == 'sons' else -1
    for seg in distance.get_segs(sec):
        start_end = distance.get_start_end(seg, electrical=electrical)
        [color, part_name]=color_func.get_seg_color(seg)
        lines.append(ax.plot([m_point] * 2,  mul * np.array([start_end['start'], start_end['end']]), color=color, linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)[0])
        segs.append(seg)
        if seg in segs_to_indecate:
            ax.scatter(m_point, mul * (start_end['start']+(start_end['end']-start_end['start'])/2.0), color=segs_to_indecate[seg]['color'], s=segs_to_indecate[seg]['size'], alpha=segs_to_indecate[seg]['alpha'], zorder=3)
    if len(sons) > 1:
        sec_start_end = distance.get_sec_start_end(sec, electrical=electrical)
          # plot vertical at end
        segs.append(seg)
        lines.append(ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']* mul] * 2 , color=color, linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)[0])

    if not distance.is_terminal(sec):
        return mid_points[-1], m_point, lines, segs
    return mid_points[-1] + BRANCH_SPACE_, m_point, lines, segs

def plot_dendogram(cell, start_seg, more_conductances, color_func, ax=None, plot_legend=False, ignore_sections=[],
                   segs_to_indecate=dict(), electrical=True, diam_factor=None, distance=None, BRANCH_SPACE_=None, dt_func= lambda x: np.mean(x)):
    if BRANCH_SPACE_ is None:
        BRANCH_SPACE_=BRANCH_SPACE
    if ax is None:
        ax = plt.gca()
    if (distance is None) or (not distance.start_seg == start_seg):
        distance = Distance(cell, more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)
    x_pos = 0.0
    [color, _] = color_func.get_seg_color(start_seg)
    mid_points = [x_pos]
    total_means = []
    lines = []
    segs = []
    for part in ['sons', 'parent']:
        mean_points=[]
        mul = 1 if part == 'sons' else -1
        mid_points = [mid_points[-1]]
        for son in distance.get_sons(start_seg.sec):
            if (distance.get_part(list(son)[0])==part) and (son not in ignore_sections):
                pos, mean_point, lines, segs = plot_recursive(son, distance, ax, color_func, mid_points[-1], lines, segs,
                                                              ignore_sections=ignore_sections, segs_to_indecate=segs_to_indecate,
                                                              electrical=electrical, diam_factor=diam_factor, BRANCH_SPACE_=BRANCH_SPACE_)
                mid_points.append(pos)
                mean_points.append(mean_point)
        sec_start_end = distance.get_sec_start_end_part(start_seg.sec, part=part, electrical=electrical)
        if len(mean_points) > 1:
            if part == 'sons':
                segs.append(list(start_seg.sec)[-1])
            else:
                segs.append(list(start_seg.sec)[0])
            lines.append(ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']*mul] * 2, color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)[0]) # plot horizental at end
        # change to segs and add scatter
        for seg in start_seg.sec:
            if distance.get_part(seg) == part or seg == start_seg:
                start_end = distance.get_start_end(seg, electrical=electrical)
                segs.append(seg)
                lines.append(ax.plot([np.mean(mean_points)]*2, [start_end['start']*mul, start_end['end']*mul], color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)[0])  # plot vertical at end
                if seg in segs_to_indecate:
                    ax.scatter(np.mean(mean_points), mul * (start_end['start'] + (start_end['end'] - start_end['start']) / 2.0), color=segs_to_indecate[seg]['color'], s=segs_to_indecate[seg]['size'], alpha=segs_to_indecate[seg]['alpha'], zorder=3)
        total_means.append(np.mean(mean_points))
    if len(total_means) > 1:
        segs.append(start_seg)
        lines.append(ax.plot([np.min(total_means), np.max(total_means)], [0] * 2, color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)[0]) # plot vertical at end

    if plot_legend:
        legend_elements = [Line2D([0], [0], color=color_func.color_dict[label], lw=2, label=label) for label in color_func.color_dict]
        ax.legend(handles=legend_elements, loc="best")
    max_y = ax.get_ylim()[1]
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if electrical:
        ax.set_ylabel('distance from origin ('+LAMDA+')')
    else:
        ax.set_ylabel('distance from origin ('+MICRO+'m)')

    return max_y, mid_points[-1], lines, segs



