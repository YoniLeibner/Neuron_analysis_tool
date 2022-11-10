from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from copy import deepcopy


def Rin_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)


# analyser = Analyzer(type='L5PC')

# dir_path = os.path.abspath('')
# morph_path=os.path.join(dir_path,'data/morph.ASC')
# analyser = Analyzer(type='ASC', morph_path=morph_path)

#
# record_dict, time = analyser.record_protocol(cut_start_ms=1000.0, record_name='v')
# animation = analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=10, clip_name='test1', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-90)
# animation.ipython_display(fps=10, loop=True, autoplay=True)
# # analyser.plot_morph(scale=500, diam_factor=1, theta=-90, ignore_soma=True)
# # plt.show()

analyser = Analyzer(type='Rall_tree')
colors_dict  = analyser.colors_dict
colors_dict['soma']='r'
colors_dict['basal']='pink'

from Neuron_analysis_tool.attenuation import plot_attenuation, record_to_value
s=list(analyser.cell.apic[20])[-1]
s2=list(analyser.cell.apic[5])[-1]
s3=list(analyser.cell.soma[0])[0]
s4=list(analyser.cell.apic[10])[0]
segs_to_indicate = {s:dict(color='green', alpha=0.5, size=50), s2:dict(color='b', alpha=0.75, size=20), s3:dict(color='r', alpha=0.3, size=20), s4:dict(color='k', alpha=0.25, size=20)}
start_seg=list(analyser.cell.soma[0])[0]
start_seg=list(analyser.cell.apic[29])[-1]
ax, norm_by = plot_attenuation(analyser.cell, start_seg, protocol=long_pulse_protocol, more_conductances=analyser.more_conductances, color_func=analyser.colors, ax=None, record_name='v',
                     cut_start_ms=None, record_to_value=record_to_value, norm_by = None, norm=True, electrical=True,
                     seg_to_indicate=segs_to_indicate, ls='-')
y_ticks = np.array([10**0, 10**-1])
ax.set_yticks(y_ticks)
# fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.apic[29])[-1], diam_factor=1, factor_e_space=100)
# fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1, factor_e_space=100)
plt.show()
#
#
# show_records_from = dict()
# nexus_sec = list(analyser.cell.apic[29])
# show_records_from[list(analyser.cell.soma[0])[0]] = dict(label='soma', alpha=0.75, color='lime', size=50)
# show_records_from[list(nexus_sec)[len(nexus_sec)//2]] = dict(label='nexus', alpha=0.75, color='grey', size=50)
#
# plt.close()
# records = analyser.record_protocol(cut_start_ms=1000.0, record_name='v')
# animation = analyser.create_movie_from_rec(records=records, slow_down_factor=50,
#                                            func_for_missing_frames=np.max, theta=-90, diam_factor=0.5,
#                                            show_records_from=show_records_from,
#                                            seg_to_indicate_dict=dict())
# animation.ipython_display(fps=10, loop=True, autoplay=True)