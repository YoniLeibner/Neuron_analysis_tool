#############################################################
#
# author: Yoni Leibner
# description: calculate and plot the attanuation in a cell
#              following a givin protocol x axis can be shown
#              in micro-meters or electricl units
# date of modification: 16.11.2022
#
#############################################################

import numpy as np
from Neuron_analysis_tool.record import record_all
from Neuron_analysis_tool.distance import Distance
from neuron import h
import matplotlib.pyplot as plt

def record_to_value(rec):
    return rec.max() - rec[0]

class color_func:
    def get_seg_color(seg):
        return 'k', 'all'

def plot_attenuation(cell, start_seg, protocol, more_conductances, color_func=color_func, ax=None, record_name='v',
                     cut_start_ms=None, record_to_value=record_to_value, norm_by = None, norm=True, electrical=True,
                     seg_to_indicate=dict(), distance=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    if start_seg is None:
        segs = list(cell.soma[0])
        start_seg = segs [len(segs)//2]
    if (distance is None) or (not distance.start_seg==start_seg):
        distance = Distance(cell, more_conductances)
        distance.compute(start_seg=start_seg)
    records = record_all(cell, record_name=record_name)
    tstop, delay, dur, amp, extra = protocol(cell, start_seg)
    if cut_start_ms is None:
        cut_start_ms = delay - 50
    records.extract(lambda x:np.array(x)[int(cut_start_ms / h.dt):])

    attanuation_vals = records.get_vals(func = record_to_value)
    if norm:
        if norm_by is None:
            norm_by = attanuation_vals[start_seg.sec][start_seg]
    else:
        norm_by=1.0


    for sec in attanuation_vals:
        if sec in cell.axonal: continue
        for seg in attanuation_vals[sec]:
            if seg == start_seg: continue
            parent_seg = distance.get_seg_parent(seg)
            start_end = distance.get_start_end(seg, electrical=electrical)
            x = [start_end['start'], start_end['end']]
            y=[attanuation_vals[parent_seg.sec][parent_seg]/norm_by, attanuation_vals[seg.sec][seg]/norm_by]
            color, part_name = color_func.get_seg_color(seg)
            ax.plot(x[-2:], y[-2:], color=color, zorder=1, **kwargs)
            if seg in seg_to_indicate.keys():
                ax.scatter(x[-1], y[-1], color=seg_to_indicate[seg]['color'], s=seg_to_indicate[seg]['size'], alpha=seg_to_indicate[seg]['alpha'], zorder=3)
    if start_seg in seg_to_indicate.keys():
        ax.scatter(0, attanuation_vals[start_seg.sec][start_seg]/norm_by, color=seg_to_indicate[start_seg]['color'], s=seg_to_indicate[start_seg]['size'], alpha=seg_to_indicate[start_seg]['alpha'], zorder=3)
    return ax, norm_by