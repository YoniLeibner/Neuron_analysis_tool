#############################################################
#
# author: Yoni Leibner
# description: calculate and plot the attanuation in a cell
#              following a givin protocol x axis can be shown
#              in micro-meters or electricl units
# date of modification:  11.05.2023
#
#############################################################

import numpy as np
from Neuron_analysis_tool.record import record_all, sec_name, seg_name
from Neuron_analysis_tool.distance import Distance
from neuron import h
import matplotlib.pyplot as plt


def record_to_value(rec):
    """
    defult function to change from record to the plotted value (takes the max diffrence from baseline)
    :param rec: np array of a recorded value
    :return:
    """
    return rec.max() - rec[0]


class color_func:
    """
    default color function color every thing in black
    """
    def get_seg_color(seg):
        return 'k', 'all'


def plot_attenuation(cell, start_seg, protocol, more_conductances, color_func=color_func, ax=None, record_name='v',
                     cut_start_ms=None, record_to_value=record_to_value, norm_by = None, norm=True, electrical=True,
                     seg_to_indicate=dict(), distance=None, records=None, ignore_sections=[], start_time=0.0,
                     end_time=None, dt_func=lambda x: np.mean(x), direction_dist_factors=dict(sons=1, parent=1), **kwargs):
    """
    ploting the attentuation along distance
    :param cell: the cell model
    :param start_seg: the egment to set the distance from (and also to norm if norm by is None and norm is True)
    :param protocol: the protocol to run in case records are not provided
    :param more_conductances: more_conductances to initiate the distance if distance is None
    :param color_func: color class givin a color and part name for each segment. default is color_func (all black)
    :param ax: the axes to plot on
    :param record_name: the name of the value o record
    :param cut_start_ms: the start time to cut from the recorded protocol (if records is None)
    :param record_to_value: a function taking a record into a single value per seg. default is record_to_value above
    :param norm_by: a numarical value to norm the records by.
    :param norm: bollean if to norm or display the original value
    :param electrical: if True it will plot the dendogram in electrical units
    :param seg_to_indicate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
    :param distance: distance class
    :param records: the records from a pre-calc protocol, if this is None a protocol will run.
    :param ignore_sections: section not to plot
    :param start_time: the start time to take the record from in ms
    :param end_time: the end time to take the record from in ms
    :param dt_func: function for dt in the more_conductances
    :param kwargs: extra kwargs for the each line in matplotlib such lw, ls etc.
    :return: the axes ploted on, the norm_by value, the lines ploted, the segments for each line, the records used in the plot
    """
    if ax is None:
        ax = plt.gca()
    if start_seg is None:
        segs = list(cell.soma[0])
        start_seg = segs [len(segs)//2]
    if (distance is None) or (not distance.start_seg==start_seg):
        distance = Distance(cell, more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)
    if records is None:
        records = record_all(cell, record_name=record_name)
        delay, extra = protocol(cell, start_seg)
        if cut_start_ms is None:
            cut_start_ms = delay - 50
        records.extract(lambda x:np.array(x)[int(cut_start_ms / h.dt):])

    if end_time is None:
        end_time = records.time[-1]
    if norm:
        if norm_by is None:
            norm_by = records.get_record_at_dt(start_seg, start_time, end_time, dt_func = record_to_value)
    else:
        norm_by=1.0

    lines = []
    segs=[]
    for sec in cell.all:
        if sec in cell.axonal: continue
        if sec in ignore_sections: continue
        for seg in sec:
            if seg == start_seg: continue
            parent_seg = distance.get_seg_parent(seg)
            start_end = distance.get_start_end(seg, electrical=electrical)
            factor = direction_dist_factors[distance.get_direction(seg)]
            x = [factor*start_end['start'], factor*start_end['end']]

            y=[records.get_record_at_dt(parent_seg, start_time, end_time, dt_func = record_to_value)/norm_by,
               records.get_record_at_dt(seg, start_time, end_time, dt_func = record_to_value)/norm_by]
            color, part_name = color_func.get_seg_color(seg)
            segs.append(seg)
            lines.append(ax.plot(x[-2:], y[-2:], color=color, zorder=1, **kwargs)[0])
            if seg in seg_to_indicate.keys():
                ax.scatter(x[-1], y[-1], color=seg_to_indicate[seg]['color'], s=seg_to_indicate[seg]['size'],
                           alpha=seg_to_indicate[seg]['alpha'], zorder=3)
    if start_seg in seg_to_indicate.keys():
        ax.scatter(0, records.get_record_at_dt(start_seg, start_time, end_time, dt_func = record_to_value)/norm_by,
                   color=seg_to_indicate[start_seg]['color'],
                   s=seg_to_indicate[start_seg]['size'],
                   alpha=seg_to_indicate[start_seg]['alpha'],
                   zorder=3)
    return ax, norm_by, lines, segs, records
