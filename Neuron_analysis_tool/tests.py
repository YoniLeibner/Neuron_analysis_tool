from Neuron_analysis_tool.load import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt

dir_path = os.path.dirname(os.path.realpath(__file__))
morph_path=os.path.join(dir_path,'data/morph.ASC')
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
apic = cell.apic
basal = cell.dend
all = cell.all
bif_seg = cell.apic[27](0.99)
#segment Rc cercate every 5 um of morphology
for sec in cell.all:
    sec.nseg = int(sec.L/40) + 1
    sec.e_pas=-70

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
        elif sec in cell.dend:
            parts_dict['basal'].append(seg)
        elif sec in cell.apic:
            parts_dict['apical'].append(seg)
        elif sec in cell.axon:
            parts_dict['axon'].append(seg)
        else:
            parts_dict['else'].append(seg)

analyser = Analyzer(cell, parts_dict, colors_dict)
# analyser.plot_morph()

def test1_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.transfer(bif_seg.x, sec=bif_seg.sec)

# analyser.plot_morph_with_value_func(func = test1_func, run_time=1000)
# analyser.plot_dendogram()
# plt.show()
#
# analyser.plot_dendogram(initial_seg=cell.apic[60](0.5))
# plt.show()
# a=1
print('run')
# analyser.plot_cable(initial_seg=None, ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='d3_2', start_loc=0, x_axis=True, plot_legend=False)
# plt.show()
#
# analyser.plot_cable(initial_seg=cell.apic[60](0.5), ax=None,
#            factor_e_space=25, factor_m_space=10,
#            dots_loc_seg=[], ignore_sections=[],
#            cable_type='dist', start_loc=0, x_axis=True)
# plt.show()


# analyser.plot_attanuation(protocol=long_pulse_protocol, ax=None, seg_to_indicate=[bif_seg],
#                           initial_seg =None, record_to_value_func=None, norm=True)
# plt.show()
#
# analyser.plot_attanuation(protocol=long_pulse_protocol, ax=None, seg_to_indicate=[bif_seg],
#                           initial_seg =cell.apic[60](0.5), record_to_value_func=None, norm=True)
# plt.show()

# analyser.create_morph_movie(cut_start_ms=1998.0, fps=1, clip_name='clip_3')

record_dict, time = analyser.record_protocol(cut_start_ms=1000.0)

import timeit
print('create_movie_from_rec, in seconds:',timeit.timeit(lambda:analyser.create_movie_from_rec(record_dict=record_dict, time=time, fps=1000, clip_name='spikes_land_mark_cp', threads=4, slow_down_factor=50, func_for_missing_frames=np.max, theta=-75), number=1))

# print('create_morph_movie, in seconds:',timeit.timeit(lambda:analyser.create_morph_movie(cut_start_ms=1000.0, fps=10, clip_name='spikes_new_9', threads=4, slow_down_factor=100, func_for_missing_frames=np.max, theta=-75), number=1))
# print('create_morph_movie2, in seconds:', timeit.timeit(lambda:analyser.create_morph_movie2(cut_start_ms=1000.0, fps=10, clip_name='spikes_new_5', threads=4, slow_down_factor=100, func_for_missing_frames=np.max, theta=-75), number=1))


# analyser.create_morph_movie(cut_start_ms=1000.0, fps=100, clip_name='spikes_new_4', threads=4, slow_down_factor=100, func_for_missing_frames=np.mean)
# analyser.create_morph_movie2(cut_start_ms=1000.0, fps=100, clip_name='spikes_new_4', threads=4, slow_down_factor=100, func_for_missing_frames=np.mean)

