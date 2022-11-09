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

dir_path = os.path.abspath('')
morph_path=os.path.join(dir_path,'data/morph.ASC')
analyser = Analyzer(type='ASC', morph_path=morph_path)

#
# record_dict, time = analyser.record_protocol(cut_start_ms=1000.0, record_name='v')
# animation = analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=10, clip_name='test1', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-90)
# animation.ipython_display(fps=10, loop=True, autoplay=True)
# # analyser.plot_morph(scale=500, diam_factor=1, theta=-90, ignore_soma=True)
# # plt.show()

# analyser = Analyzer(type='Rall_tree')
# colors_dict  = analyser.colors_dict
# colors_dict['soma']='r'
# colors_dict['basal']='pink'

show_records_from = [list(analyser.cell.soma[0])[0], list(analyser.cell.apic[29])[-1]]
show_records_from = []
# h.dt=0.5
# h.steps_per_ms = 1.0/h.dt
record_dict, time = analyser.record_protocol(cut_start_ms=1000.0, record_name='v')
animation = analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=10, clip_name='bAP_rall_tree', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-90, diam_factor=0.5, show_records_from=show_records_from)
animation.ipython_display(fps=10, loop=True, autoplay=True)


# analyser.change_color_dict(colors_dict)
fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1)
plt.show()
analyser.create_card(scale=100, start_seg=list(analyser.cell.apic[29])[-1], diam_factor=1)
plt.show()