#testing script dont use this!!!!
# for examples look on jupyter notebook

from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from copy import deepcopy



analyser = Analyzer(type='Rall_tree')
colors_dict  = analyser.colors_dict
colors_dict['soma']='r'
colors_dict['basal']='pink'

analyser.change_color_dict(colors_dict)

for sec in analyser.cell.all:
    lamda = ((1.0/sec.g_pas/sec.Ra) * (sec.diam/10000.0/4.0))**0.5 * 10000.0
    print('name:',sec, ', e_length:',round(sec.L/lamda, 3), ', diam:', sec.diam)
#########################################################################################################################################################################

fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1)
plt.show()

#########################################################################################################################################################################

fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.apic[29])[-1], diam_factor=1, factor_e_space=100)
plt.show()

#########################################################################################################################################################################

def Rin_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

plt.title('Rin with color code')
ax, color_bar, colors = analyser.plot_morph_with_value_func(func = Rin_func, run_time=1000, theta=0, diam_factor=1, scale=500)
color_bar.set_ylabel('Rin (M ohm)')
plt.show()
#########################################################################################################################################################################
analyser.plot_dendogram_with_value_func(func = None, diam_factor=1, colors=colors)
plt.show()
#########################################################################################################################################################################
def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

colors_dict = deepcopy(analyser.colors_dict)
plt.figure()
colors_dict1 = {key:'k' for key in colors_dict.keys()}
analyser.change_color_dict(colors_dict1)

seg_to_indicate_dict =dict()
seg_to_indicate_dict[list(analyser.cell.soma[0])[0]]=dict(size=100, color='r', alpha=0.75)
seg_to_indicate_dict[list(analyser.cell.apic[29])[-1]]=dict(size=100, color='b', alpha=0.75)
analyser.plot_morph(ax=plt.gca(), theta=0, seg_to_indicate_dict=seg_to_indicate_dict, scale=200, diam_factor=1)


fig, ax = plt.subplots(1, 1, figsize=(10, 10))
colors_dict2 = {key:'b' for key in colors_dict.keys()}
analyser.change_color_dict(colors_dict2)

ax, norm_by, _, _ = analyser.plot_attenuation(protocol=long_pulse_protocol, ax=ax, start_seg=list(analyser.cell.apic[29])[-1], seg_to_indicate_dict=seg_to_indicate_dict, label='forward')

colors_dict3 =  {key:'r' for key in colors_dict.keys()}
analyser.change_color_dict(colors_dict3)
ax, norm_by, _, _ = analyser.plot_attenuation(protocol=long_pulse_protocol, ax=ax, start_seg=list(analyser.cell.soma[0])[0], norm=True, norm_by=norm_by, seg_to_indicate_dict=seg_to_indicate_dict, ls='--', dashes=(1, 100), label='backward')
analyser.change_color_dict(colors_dict)

ax.legend()
legend_without_duplicate_labels(ax)
plt.show()
#########################################################################################################################################################################