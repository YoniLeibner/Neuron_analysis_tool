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

# analyser.create_card(scale=100, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1)
# plt.show()
analyser.create_card(scale=100, start_seg=list(analyser.cell.apic[29])[-1], diam_factor=1)
plt.show()
