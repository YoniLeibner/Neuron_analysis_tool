#############################################################
#
# author: Yoni Leibner
# description: plot the dendogram of the cell in micro-meters
#              or electricl units
# date of modification: 11.05.2023
#
#############################################################

# from neuron import h
import numpy as np
from Neuron_analysis_tool.distance import Distance
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Neuron_analysis_tool.utils import LAMDA, MICRO

FIX_DIAM=1
BRANCH_SPACE=2


def plot_recursive(sec, distance, ax, color_func, x_pos, lines, segs, ignore_sections=[], segs_to_indecate=dict(),
                   electrical=True, diam_factor=None, BRANCH_SPACE_=BRANCH_SPACE):
    """
    plot the dendogram recursivly
    :param sec: the current section to plot
    :param distance: a pre calculated distance class
    :param ax: the axes to plot on
    :param color_func: color class givin a color and part name for each segment.
    :param x_pos: the current x_pos for plot
    :param lines: list of ploted lines
    :param segs: list of ploted line segment ref
    :param ignore_sections: section not to plot
    :param segs_to_indecate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
    :param electrical: if True it will plot the dendogram in electrical units
    :param diam_factor:  factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
    :param BRANCH_SPACE_: the dx to modve between terminal branches
    :return: mid point of the current sub-tree, the last x_pos, list of ploted lines, list of ploted line segment ref
    """
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

    mul = 1 if distance.get_direction(list(sec)[0]) == 'sons' else -1
    for seg in distance.get_segs(sec):
        start_end = distance.get_start_end(seg, electrical=electrical)
        [color, part_name]=color_func.get_seg_color(seg)
        lines.append(ax.plot([m_point] * 2,  mul * np.array([start_end['start'], start_end['end']]), color=color,
                             linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)[0])
        segs.append(seg)
        if seg in segs_to_indecate:
            ax.scatter(m_point, mul * (start_end['start']+(start_end['end']-start_end['start'])/2.0),
                       color=segs_to_indecate[seg]['color'], s=segs_to_indecate[seg]['size'],
                       alpha=segs_to_indecate[seg]['alpha'], zorder=3)
    if len(sons) > 1:
        sec_start_end = distance.get_sec_start_end(sec, electrical=electrical)
          # plot vertical at end
        segs.append(seg)
        lines.append(ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']* mul] * 2 , color=color,
                             linewidth= FIX_DIAM if diam_factor is None else seg.diam * diam_factor, zorder=2)[0])

    if not distance.is_terminal(sec):
        return mid_points[-1], m_point, lines, segs
    return mid_points[-1] + BRANCH_SPACE_, m_point, lines, segs


def plot_dendogram(cell, start_seg, more_conductances, color_func, ax=None, plot_legend=False, ignore_sections=[],
                   segs_to_indecate=dict(), electrical=True, diam_factor=None, distance=None, BRANCH_SPACE_=None,
                   dt_func= lambda x: np.mean(x)):
    """
    plot a dendogram of a cell from a givin start segment
    :param cell: the cell to plot (Neuron model)
    :param start_seg: the seg to start the dendogram from
    :param more_conductances: more_conductances to initiate the distance if distance is None
    :param color_func: a givin color func that gives color to each segment via get_seg_color func
    :param ax: the axes to plot on
    :param plot_legend: if True it will dd a legend with the parts_dict keys and colors
    :param ignore_sections: section not to plot
    :param segs_to_indecate:  dictinary {seg: dict(size=100, color='b', alpha=0.75)}
    :param electrical: if True it will plot the dendogram in electrical units
    :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
    :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
    :param BRANCH_SPACE_: spacing between the terminal branches (makes the dendogram more/less sparse)
    :param dt_func: function for dt in the more_conductances
    :return: max ylim of the axes, the larget x_pos, list of ploted lines, list of ploted line segment ref
    """
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
    for direction in ['sons', 'parent']:
        mean_points=[]
        mul = 1 if direction == 'sons' else -1
        mid_points = [mid_points[-1]]
        for son in distance.get_sons(start_seg.sec):
            if (distance.get_direction(list(son)[0])==direction) and (son not in ignore_sections):
                pos, mean_point, lines, segs = plot_recursive(son, distance, ax, color_func, mid_points[-1], lines, segs,
                                                              ignore_sections=ignore_sections, segs_to_indecate=segs_to_indecate,
                                                              electrical=electrical, diam_factor=diam_factor, BRANCH_SPACE_=BRANCH_SPACE_)
                mid_points.append(pos)
                mean_points.append(mean_point)
        sec_start_end = distance.get_sec_start_end_direction(start_seg.sec, direction=direction, electrical=electrical)
        if len(mean_points) > 1:
            if direction == 'sons':
                segs.append(list(start_seg.sec)[-1])
            else:
                segs.append(list(start_seg.sec)[0])
            lines.append(ax.plot([mean_points[0], mean_points[-1]], [sec_start_end['end']*mul] * 2, color=color, linewidth= FIX_DIAM if diam_factor is None else start_seg.diam * diam_factor, zorder=2)[0]) # plot horizental at end
        # change to segs and add scatter
        for seg in start_seg.sec:
            if distance.get_direction(seg) == direction or seg == start_seg:
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



