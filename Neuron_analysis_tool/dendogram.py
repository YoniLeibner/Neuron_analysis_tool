import numpy as np
from Neuron_analysis_tool.distance import Distance
from neuron import h
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

FIX_DIAM=1
BRANCH_SPACE=2


def plot_recursive(sec, distance, ax, color_func, x_pos, ignore_sections=[], segs_to_indecate=dict(), electrical=True, diam_factor=None):
    # start_end = distance.get_sec_start_end(sec, electrical=electrical)
    [color, part_name]=color_func.get_seg_color(list())
    # colors, lengths = color_func(sec, distance.get_sec_parent(sec).sec) # change to seg colors!!!!!!!
    mid_points = [x_pos]
    mean_points=[]
    sons = distance.get_sons(sec)
    if not distance.is_terminal(sec):
        for son in sons:
            if son not in ignore_sections:
                pos, mean_point = plot_recursive(son, distance, ax, color_func, mid_points[-1], ignore_sections=ignore_sections, segs_to_indecate=segs_to_indecate, electrical=electrical, diam_factor=diam_factor)
                mid_points.append(pos)
                mean_points.append(mean_point)
        if len(sons)>1:
            m_point = (mid_points[0]+ mid_points[-1]-BRANCH_SPACE)/2.0
        else:
            m_point = mid_points[-1] - BRANCH_SPACE # single cheldren
    else:
        m_point = x_pos # terminal

    mul = 1 if distance.get_part(list(sec)[0]) == 'sons' else -1
    for seg in distance.get_segs(sec):
        start_end = distance.get_start_end(seg, electrical=electrical)
        [color, part_name]=color_func.get_seg_color(seg)
        ax.plot([m_point] * 2,  mul * np.array([start_end['start'], start_end['end']]), color=color, linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)
        if seg in segs_to_indecate:
            ax.scatter(m_point, mul * (start_end['start']+(start_end['end']-start_end['start'])/2.0), color=segs_to_indecate[seg]['color'], s=segs_to_indecate[seg]['size'], alpha=segs_to_indecate[seg]['alpha'], zorder=2)
    if len(sons) > 1:
        sec_start_end = distance.get_sec_start_end(sec, electrical=electrical)
        ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']* mul] * 2 , color=color, linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)  # plot vertical at end

    if not distance.is_terminal(sec):
        return mid_points[-1], m_point
    return mid_points[-1] + BRANCH_SPACE, m_point

def plot_dendogram(cell, start_seg, more_conductances, color_func, ax=None, plot_legend=False, ignore_sections=[], segs_to_indecate=dict(), electrical=True, diam_factor=None, distance=None):
    if ax is None:
        ax = plt.gca()
    if (distance is None) or (not distance.start_seg==start_seg):
        distance = Distance(cell, more_conductances)
        distance.compute(start_seg=start_seg)
    x_pos = 0.0
    [color, part_name] = color_func.get_seg_color(start_seg)
    mid_points = [x_pos]
    total_means = []
    for part in ['sons', 'parent']:
        mean_points=[]
        mul = 1 if part == 'sons' else -1
        mid_points = [mid_points[-1]]
        for son in distance.get_sons(start_seg.sec):
            if (distance.get_part(list(son)[0])==part) and (son not in ignore_sections):
                pos, mean_point = plot_recursive(son, distance, ax, color_func, mid_points[-1], ignore_sections=ignore_sections, segs_to_indecate=segs_to_indecate, electrical=electrical, diam_factor=diam_factor)
                mid_points.append(pos)
                mean_points.append(mean_point)
        sec_start_end = distance.get_sec_start_end_part(start_seg.sec, part=part, electrical=electrical)
        if len(mean_points) > 1:
            ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']*mul] * 2, color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)  # plot horizental at end
        # change to segs and add scatter
        for seg in start_seg.sec:
            if distance.get_part(seg) == part:
                start_end = distance.get_start_end(seg, electrical=electrical)
                ax.plot([np.mean(mean_points)]*2, [start_end['start']*mul, start_end['end']*mul], color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)  # plot vertical at end
                if seg in segs_to_indecate:
                    ax.scatter(np.mean(mean_points), mul * (start_end['start'] + (start_end['end'] - start_end['start']) / 2.0), color=segs_to_indecate[seg]['color'], s=segs_to_indecate[seg]['size'], alpha=segs_to_indecate[seg]['alpha'], zorder=2)
        # ax.plot([np.mean(mean_points)]*2, [sec_start_end['start']*mul, sec_start_end['end']*mul], color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)  # plot vertical at end
        total_means.append(np.mean(mean_points))
    if len(total_means) > 1:
        ax.plot([np.min(total_means), np.max(total_means)], [0] * 2, color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)  # plot vertical at end
    if plot_legend:
        legend_elements = [Line2D([0], [0], color=color_func.color_dict[label], lw=2, label=label) for label in color_func.color_dict]
        ax.legend(handles=legend_elements, loc="best")
    max_y = ax.get_ylim()[1]
    return max_y, mid_points[-1]


