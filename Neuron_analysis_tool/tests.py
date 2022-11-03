from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt


def test1_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)


dir_path = os.path.dirname(os.path.realpath(__file__))
morph_path=os.path.join(dir_path,'data/morph.ASC')

# analyser = Analyzer(type='ASC', morph_path=morph_path)
analyser = Analyzer(type='L5PC')#, morph_path=morph_path)
cell = analyser.cell
apic = list(cell.apical)
basal = list(cell.basal)

parts_dict = {'soma':[], 'basal':[], 'apical':[], 'axon':[], 'else':[]}
colors_dict = {'soma':'k', 'basal':'r', 'apical':'b', 'axon':'green', 'else':'cyan'}
for sec in cell.all:
    for seg in sec:
        if sec in cell.soma:
            parts_dict['soma'].append(seg)
        elif sec in basal:
            parts_dict['basal'].append(seg)
        elif sec in apic:
            parts_dict['apical'].append(seg)
        elif sec in cell.axon:
            parts_dict['axon'].append(seg)
        else:
            parts_dict['else'].append(seg)

analyser.change_parts_dict(parts_dict=parts_dict, colors_dict=colors_dict)

analyser.create_card(theta=-75)
plt.show()
analyser.create_card(initial_seg=list(cell.apic[101])[-1], theta=-75)
plt.show()
analyser.plot_morph_with_value_func(func = test1_func, run_time=1000, theta=-75)
# analyser.plot_morph_with_value_func(run_time=1000, theta=-75)
# plt.show()

# analyser.plot_dendogram(electrical=True)
# s_pos = 0.25/2
# plt.axhline(s_pos, color='r', ls='--')
# for i in range(1, 5):
#     plt.axhline(s_pos + i*0.25, color='r', ls='--')
#
# plt.show()

# analyser.plot_dendogram(initial_seg=cell.apic[60](0.5))
# plt.show()
# a=1

# print('run')
# analyser.plot_cable(initial_seg=soma(0), ax=None,
#             factor_e_space=50,#25, #this is a problem do to the spacing!!! w ecan get the same branch twise!!!
#             factor_m_space=10,
#             dots_loc_seg=[], ignore_sections=[],
#             cable_type='d3_2', start_loc=0, x_axis=True, plot_legend=False)
# plt.show()

# analyser.plot_cable(initial_seg=cell.apic[60](0.5), ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='dist', start_loc=0, x_axis=True)
# plt.show()


# fig, ax = plt.subplots(1, 1, figsize=(10, 10))
# colors_dict2 = {'soma':'k', 'basal':'b', 'apical':'b', 'axon':'b', 'else':'b'}
#
# analyser.change_color_dict(colors_dict2)
#
# ax, norm_by = analyser.plot_attanuation(protocol=long_pulse_protocol, ax=ax, seg_to_indicate=[bif_seg],
#                           initial_seg =list(apic[29])[-1], record_to_value_func=None, norm=True)
#
# colors_dict3 = {'soma':'pink', 'basal':'r', 'apical':'r', 'axon':'r', 'else':'r'}
# analyser.change_color_dict(colors_dict3)
# analyser.plot_attanuation(protocol=long_pulse_protocol, ax=ax, seg_to_indicate=[bif_seg],
#                           initial_seg =list(soma)[0], record_to_value_func=None, norm=False, norm_by=norm_by, ls='--', dashes=(1, 200))
# plt.show()
# analyser.change_color_dict(colors_dict)
# exit(0)

# analyser.create_morph_movie(cut_start_ms=1998.0, fps=1, clip_name='clip_3')


# record_dict, time = analyser.record_protocol(cut_start_ms=1000.0)
# #
# import timeit
# print('create_movie_from_rec, in seconds:',timeit.timeit(lambda:analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=1000, clip_name='spikes_land_mark_optim2', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-75), number=1))
