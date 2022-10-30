from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from Neuron_analysis_tool.color_func import color_func

dir_path = os.path.dirname(os.path.realpath(__file__))
morph_path=os.path.join(dir_path,'data/morph.ASC')
# morph_path=os.path.join(dir_path,'data/Rall_tree5.swc')
# morph_path=os.path.join(dir_path,'data\\morph.ASC')
hoc_file_name = 'allen_model.hoc'

h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h("objref cell, tobj") #neuron object
h.init()
h.load_file(os.path.join(dir_path,'allen_model.hoc').replace('\\', '/'))


h.execute("cell = new " + hoc_file_name[:-4] + "()")  # replace?
if morph_path.endswith('ASC'):
    nl = h.Import3d_Neurolucida3()
    nl.quiet = 1
elif morph_path.endswith('swc'):
    nl = h.Import3d_SWC_read()
else:
    print ('problem 2')
    raise BaseException("the morphology path isn't correct,\n check the path")
nl.input(morph_path)
i3d = h.Import3d_GUI(nl, 0)
i3d.instantiate(h.cell)
cell=h.cell

#this is proc of the hoc template do I need them?
cell.geom_nseg()
cell.delete_axon()
cell.biophys()

soma = cell.soma[0]
apic = list(cell.apical)
basal = list(cell.basal)
all = cell.all
bif_seg = cell.soma[0](0)
bif_seg = cell.apic[101](0.99)
#segment Rc cercate every 5 um of morphology
for sec in cell.all:
    sec.nseg = int(sec.L/10) + 1
    sec.e_pas=-70
    sec.Ra=100
    sec.g_pas = 1.0/10000.0
    sec.cm=1

# soma.insert('hh')
# f=10
# soma.gnabar_hh = 0.12 * f
# soma.gkbar_hh = 0.036 * f
# print(soma.gnabar_hh)
# print(soma.gkbar_hh)
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

analyser = Analyzer(cell, parts_dict, colors_dict)
# analyser.plot_morph(theta=-75)
# plt.show()
def test1_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)
    # return imp.transfer(bif_seg.x, sec=bif_seg.sec)

# analyser.plot_morph_with_value_func(func = test1_func, run_time=1000, theta=-75)
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


record_dict, time = analyser.record_protocol(cut_start_ms=1000.0)

import timeit
print('create_movie_from_rec, in seconds:',timeit.timeit(lambda:analyser.create_movie_from_rec2(record_dict=record_dict, time=time, fps=1000, clip_name='spikes_land_mark_optim2', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-75), number=1))
