from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from Neuron_analysis_tool.color_func import color_func

analyser = Analyzer(type='Rall_tree')
# analyser.plot_morph()

def test1_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

# analyser.plot_morph_with_value_func(func = test1_func, run_time=1000)
# plt.show()
#
analyser.plot_dendogram(electrical=True, initial_seg=analyser.cell.soma[0](0))
s_pos = 0
plt.axhline(s_pos, color='r', ls='--')
for i in range(1, 6):
    plt.axhline(s_pos + i*0.25, color='r', ls='--')

plt.show()

analyser.plot_dendogram(initial_seg=analyser.cell.apic[29](1))
plt.show()

analyser.plot_dendogram(initial_seg=analyser.cell.apic[29](1), electrical=True)
plt.show()
exit(0)
# a=1

# print('run')
analyser.plot_cable(initial_seg=analyser.cell.soma[0](0), ax=None,
            factor_e_space=150,#25, #this is a problem do to the spacing!!! w ecan get the same branch twise!!!
            factor_m_space=10,
            dots_loc_seg=[], ignore_sections=[],
            cable_type='d3_2', start_loc=0, x_axis=True, plot_legend=False)
plt.show()

# analyser.plot_cable(initial_seg=cell.apic[60](0.5), ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='dist', start_loc=0, x_axis=True)
# plt.show()


fig, ax = plt.subplots(1, 1, figsize=(10, 10))
colors_dict2 = {'soma':'k', 'basal':'b', 'apical':'b', 'axon':'b', 'else':'b'}

analyser.change_color_dict(colors_dict2)

ax, norm_by = analyser.plot_attanuation(protocol=long_pulse_protocol, ax=ax, seg_to_indicate=[bif_seg],
                          initial_seg =list(analyser.cell.apic[29])[-1], record_to_value_func=None, norm=True)

colors_dict3 = {'soma':'pink', 'basal':'r', 'apical':'r', 'axon':'r', 'else':'r'}
analyser.change_color_dict(colors_dict3)
analyser.plot_attanuation(protocol=long_pulse_protocol, ax=ax, seg_to_indicate=[bif_seg],
                          initial_seg =list(analyser.cell.soma[0])[0], record_to_value_func=None, norm=False, norm_by=norm_by, ls='--', dashes=(1, 200))
plt.show()
# analyser.change_color_dict(analyser.colors_dict)
exit(0)

# analyser.create_morph_movie(cut_start_ms=1998.0, fps=1, clip_name='clip_3')

record_dict, time = analyser.record_protocol(cut_start_ms=1000.0)
