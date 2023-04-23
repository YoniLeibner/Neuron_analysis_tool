from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt




analyser = Analyzer(type='Rall_tree')
colors_dict  = analyser.colors_dict
colors_dict['soma']='r'
colors_dict['basal']='pink'

analyser.change_color_dict(colors_dict)

record_dict, extra = analyser.record_protocol(cut_start_ms=1000.0, record_name='v')

animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=50,
                                           func_for_missing_frames=np.max, theta=0, diam_factor=0.5,
                                           show_records_from=dict(), electrical=False,
                                           plot_all_records=True, distance_factor=10, plot_every=0.5)
animation.ipython_display(fps=10, loop=True, autoplay=True)

a=1
def Ca_spike_protocol(cell, start_seg=None):
    delay = 200.0
    stim = h.IClamp(0.5, sec=cell.soma[0])
    stim.dur = 5
    stim.delay = delay

    syn = h.epsp(cell.apic[36](0.9))
    syn.tau0 = 0.5
    syn.tau1 = 5
    syn.onset = stim.delay + 5
    syn.imax = 0.5
    stim.amp = 1.9

    h.tstop = 400
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 200, {}

def resting_protocol2(cell, start_seg=None):
    h.tstop = 50.0
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}

analyser = Analyzer(type='L5PC')

# animation=analyser.dancing_morph(protocol=Ca_spike_protocol, seg_to_indicate_dict=dict(), diam_factor=1,
#                             sec_to_change=None, ignore_sections=[], theta=-90, scale=0.25,
#                             slow_down_factor=1, figsize=(5,5))
# animation.ipython_display(fps=10, loop=True, autoplay=True)
record_dict, extra = analyser.record_protocol(protocol=Ca_spike_protocol, cut_start_ms=None,
                                              record_name='v', compute_more_condunctances=True)

seg_to_indicate_dict={list(analyser.cell.dend[30])[-2]:dict(color='cyan', alpha=0.5, size=20),
                     list(analyser.cell.apic[70])[-2]:dict(color='green', alpha=0.5, size=20)}

animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=20, scale=0.25,
                                           func_for_missing_frames=np.max, theta=-90, diam_factor=None,
                                           show_records_from=dict(), draw_funcs=extra['draw_funcs'],
                                           electrical=True, dancing=True,
                                           more_conductances_=extra['more_conductances'],
                                          base_plot_type='dendogram',figsize=(7,7),
                                           seg_to_indicate_dict=seg_to_indicate_dict)
animation.ipython_display(fps=10, loop=True, autoplay=True)



record_dict, extra = analyser.record_protocol(cut_start_ms=1000.0, record_name='v', compute_more_condunctances=True)

animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=1, scale=0.25,
                                           func_for_missing_frames=np.max, theta=-90, diam_factor=0.5,
                                           show_records_from=dict(), draw_funcs=extra['draw_funcs'],
                                          electrical=True, dancing=True, more_conductances_=extra['more_conductances'],
                                          base_plot_type='dendogram')
animation.ipython_display(fps=10, loop=True, autoplay=True)
# animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=1, scale=0.25,
#                                            func_for_missing_frames=np.max, theta=-90, diam_factor=0.5,
#                                            show_records_from=dict(), draw_funcs=extra['draw_funcs'],
#                                           electrical=True, dancing=True, more_conductances_=extra['more_conductances'])
# animation.ipython_display(fps=10, loop=True, autoplay=True)



#####################################################################################################
_, _, _ = analyser.plot_morph(theta=-90, seg_to_indicate_dict=dict(),
                          scale=0.25, diam_factor=1,
                          ignore_soma=True, distance=None, electrical=True)
# plt.show()
plt.figure()
#####################################################################################################
_,_,_ = analyser.plot_morph(scale=500, diam_factor=0.5, theta=-90, ignore_soma=True)
plt.show()
#####################################################################################################

show_records_from = dict()
show_records_from[list(analyser.cell.soma[0])[0]] = dict(label='soma', alpha=0.75, color='lime', size=50)
show_records_from[list(analyser.cell.apic[36])[0]] = dict(label='dend1', alpha=0.75, color='grey', size=40)

record_dict, _ = analyser.record_protocol(protocol=Ca_spike_protocol, cut_start_ms=None, record_name='v')
animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=50,
                                           func_for_missing_frames=np.max, theta=-90, diam_factor=0.5,
                                           show_records_from=show_records_from, seg_to_indicate_dict=dict(),
                                           base_plot_type='dendogram')
animation.ipython_display(fps=10, loop=True, autoplay=True)
