from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt


analyser = Analyzer(type='Rall_tree')
# analyser.plot_morph()

terninals = [sec for sec in analyser.cell.all if len(sec.children())==0]

def test1_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

analyser.plot_morph_with_value_func(func = test1_func, run_time=1000)
analyser.plot_dendogram()
plt.show()

#
# analyser.plot_dendogram(start_seg=terninals[0](1))
# plt.show()
exit(0)
print('run')
# analyser.plot_cable(start_seg=None, ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='d3_2', start_loc=0, x_axis=True, plot_legend=False)
# plt.show()
#
# analyser.plot_cable(start_seg=cell.apic[60](0.5), ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='dist', start_loc=0, x_axis=True)
# plt.show()


# analyser.plot_attanuation(protocol=long_pulse_protocol, ax=None, seg_to_indicate=[bif_seg],
#                           start_seg =None, record_to_value_func=None, norm=True)
# plt.show()
#
# analyser.plot_attanuation(protocol=long_pulse_protocol, ax=None, seg_to_indicate=[bif_seg],
#                           start_seg =cell.apic[60](0.5), record_to_value_func=None, norm=True)
# plt.show()

# analyser.create_morph_movie(cut_start_ms=1998.0, fps=1, clip_name='clip_3')

record_dict, time = analyser.record_protocol(cut_start_ms=1000.0)

import timeit
print('create_movie_from_rec, in seconds:',timeit.timeit(lambda:analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=1000, clip_name='spikes_land_mark_optim', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-75), number=1))

# print('create_morph_movie, in seconds:',timeit.timeit(lambda:analyser.create_morph_movie(cut_start_ms=1000.0, fps=10, clip_name='spikes_new_9', threads=4, slow_down_factor=100, func_for_missing_frames=np.max, theta=-75), number=1))
# print('create_morph_movie2, in seconds:', timeit.timeit(lambda:analyser.create_morph_movie2(cut_start_ms=1000.0, fps=10, clip_name='spikes_new_5', threads=4, slow_down_factor=100, func_for_missing_frames=np.max, theta=-75), number=1))


# analyser.create_morph_movie(cut_start_ms=1000.0, fps=100, clip_name='spikes_new_4', threads=4, slow_down_factor=100, func_for_missing_frames=np.mean)
# analyser.create_morph_movie2(cut_start_ms=1000.0, fps=100, clip_name='spikes_new_4', threads=4, slow_down_factor=100, func_for_missing_frames=np.mean)

