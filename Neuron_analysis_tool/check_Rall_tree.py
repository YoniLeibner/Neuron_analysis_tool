#testing script dont use this!!!!
# for examples look on jupyter notebook

from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol, short_pulse_protocol, spike_protocol2, spike_protocol3, spike_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from copy import deepcopy


def add_syn(seg, g_AMPA=0.0004, g_NMDA=0.0004):
    netstim = h.NetStim()
    netstim.interval = 400  # fr of 1
    netstim.start = 400
    netstim.noise = 1.0
    netstim.number = 0
    # AMPA part
    AMPA = h.Exp2Syn(seg.x, sec=seg.sec)
    AMPA_con = h.NetCon(netstim, AMPA)
    AMPA.e = 0
    AMPA.tau1 = 0.3
    AMPA.tau2 = 1.5
    AMPA_con.weight[0] = g_AMPA
    AMPA_con.delay = 0

    # NMDA part
    NMDA = h.NMDA(seg.x, sec=seg.sec)
    NMDA_con = h.NetCon(netstim, NMDA)
    NMDA.e = 0
    NMDA.tau_r_NMDA = 8
    NMDA.tau_d_NMDA = 35
    NMDA.n_NMDA = 0.27
    NMDA.gama_NMDA = 0.076
    NMDA_con.weight[0] = g_NMDA
    NMDA_con.delay = 0
    return [AMPA, AMPA_con], [NMDA, NMDA_con], netstim


def add_isyn(seg, g_GABA=0.0004):
    netstim = h.NetStim()
    netstim.interval = 500  # fr of 1
    netstim.start = 400
    netstim.noise = 1.0
    netstim.number = 0
    # AMPA part
    AMPA = h.Exp2Syn(seg.x, sec=seg.sec)
    AMPA_con = h.NetCon(netstim, AMPA)
    AMPA.e = -100
    AMPA.tau1 = 0.2
    AMPA.tau2 = 4
    AMPA_con.weight[0] = g_GABA
    AMPA_con.delay = 0

    return [AMPA, AMPA_con], netstim


def random_syn_protocol(cell, start_seg):
    syns = []
    isyns = []
    segs_e = []
    syn_e_times = []

    segs_i = []
    syn_i_times = []

    total_time = 1400
    start_time = 400
    number_of_e_syns = 5000
    number_of_i_syns = 2000

    amp = 0
    e_FR = 2  # Hz
    i_FR = 2  # Hz
    # difining the synapses times
    sim_time = total_time - start_time
    for sec in np.random.choice(list(cell.all), number_of_e_syns):
        seg_num = np.random.randint(0, len(list(sec)))
        segs_e.append(list(sec)[seg_num])
        syns.append(add_syn(list(sec)[seg_num]))
        syn_time = np.random.rand(int(np.ceil(sim_time / 1000.0 * e_FR))) * sim_time
        syn_e_times.append(syn_time)

    for sec in np.random.choice(list(cell.all), number_of_i_syns):
        seg_num = np.random.randint(0, len(list(sec)))
        segs_i.append(list(sec)[seg_num])
        isyns.append(add_isyn(list(sec)[seg_num]))
        syn_time = np.random.rand(int(np.ceil(sim_time / 1000.0 * i_FR))) * sim_time
        syn_i_times.append(syn_time)

    # function to insert synapse time to neuron
    def event_setter():
        for syn, times in zip(syns, syn_e_times):
            for t in times:
                syn[0][1].event(t)  # AMPA netcon
                syn[1][1].event(t)  # NMDA netcon
        for syn, times in zip(isyns, syn_i_times):
            for t in times:
                syn[0][1].event(t)  # GABA netcon

    # this function is return to the movie maker and is called every frame generation, than one can add stuff to the frame.
    # you must return all the new created items, so the function can call .remove() to remove this items in the next frame

    def draw_func1(start_time, end_time, segs, lines, ax, records):
        elements = []
        for seg, syn, times in zip(segs_e, syns, syn_e_times):
            for t in times:
                if t > start_time and t < end_time:
                    if seg in segs:
                        l = lines[segs == seg][0]
                        try:
                            elements.append(
                                ax.scatter(l.get_xdata().mean(), l.get_ydata().mean(), color='r', alpha=0.75, s=10))
                        except:
                            pass
        for seg, syn, times in zip(segs_i, isyns, syn_i_times):
            for t in times:
                if t > start_time and t < end_time:
                    if seg in segs:
                        l = lines[segs == seg][0]
                        try:
                            elements.append(
                                ax.scatter(l.get_xdata().mean(), l.get_ydata().mean(), color='b', alpha=0.75, s=10))
                        except:
                            pass
        return elements

    h.tstop = total_time
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    fih = h.FInitializeHandler(0, event_setter)
    h.run()
    return start_time, dict(draw_funcs=[draw_func1])

def Ca_spike_protocol(cell, start_seg=None):
    delay = 1000.0
    stim = h.IClamp(0.5, sec=cell.soma[0])
    stim.dur = 5
    stim.delay = delay

    syn = h.epsp(cell.apic[36](0.9))
    syn.tau0 = 0.5
    syn.tau1 = 5
    syn.onset = stim.delay + 5
    syn.imax = 0.5
    stim.amp = 1.9

    stim2 = h.IClamp(0.5, sec=cell.soma[0])
    stim2.dur = 5
    stim2.delay = delay + 200

    syn2 = h.epsp(cell.apic[36](0.9))
    syn2.tau0 = 0.5
    syn2.tau1 = 5
    syn2.onset = stim.delay + 200 + 5
    syn2.imax = 0.5
    stim2.amp = 1.9

    h.tstop = 1300
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 400, {}


analyser = Analyzer(type='L5PC')

# analyser = Analyzer(type='Rall_tree')
# colors_dict  = analyser.colors_dict
# colors_dict['soma']='r'
# colors_dict['basal']='pink'

# analyser.change_color_dict(colors_dict)

# for sec in analyser.cell.all:
#     lamda = ((1.0/sec.g_pas/sec.Ra) * (sec.diam/10000.0/4.0))**0.5 * 10000.0
#     print('name:',sec, ', e_length:',round(sec.L/lamda, 3), ', diam:', sec.diam)
# #########################################################################################################################################################################
# _, _, _ = analyser.plot_morph(theta=0, seg_to_indicate_dict=dict(),
#                           scale=0.25, diam_factor=1,
#                           ignore_soma=not analyser.type.startswith('Rall_tree'),
#                               distance=None, electrical=True)
# plt.show()
#
# #########################################################################################################################################################################
#
#
# fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1)
# plt.show()
#
# #########################################################################################################################################################################
#
# fig, ax = analyser.create_card(scale=500, start_seg=list(analyser.cell.apic[29])[-1], diam_factor=1, factor_e_space=100)
# plt.show()
#
# #########################################################################################################################################################################
#
# def Rin_func(seg):
#     imp = h.Impedance(seg.x, sec=seg.sec)
#     imp.loc(seg.x, sec=seg.sec)
#     imp.compute(0, 1)
#     return imp.input(seg.x, sec=seg.sec)
#
# plt.title('Rin with color code')
# ax, color_bar, colors, _, _ = analyser.plot_morph_with_value_func(func = Rin_func, run_time=1000, theta=0, diam_factor=1, scale=500)
# color_bar.set_ylabel('Rin (M ohm)')
# plt.show()
# #########################################################################################################################################################################
#
# analyser.plot_dendogram_with_value_func(func = None, diam_factor=1, colors=colors)
# plt.show()
# #########################################################################################################################################################################

# start_seg=None
# record_dict, extra = analyser.record_protocol(protocol=short_pulse_protocol,cut_start_ms=1990.0, record_name='v', start_seg=start_seg)
# animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=500,
#                                            func_for_missing_frames=np.max, theta=0, diam_factor=0.5,
#                                            show_records_from=dict(), electrical=True,
#                                            base_plot_type='attenuation', start_seg=start_seg)
# # animation.ipython_display(fps=10, loop=True, autoplay=True)
# animation.write_videofile('./check_Rall.mp4',fps=10, threads=4,audio=False)

# start_seg=list(analyser.cell.apic[29])[-1]
# record_dict, extra = analyser.record_protocol(protocol=short_pulse_protocol,cut_start_ms=1990.0, record_name='v', start_seg=start_seg)
# animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=500,
#                                            func_for_missing_frames=np.max, theta=0, diam_factor=0.5,
#                                            show_records_from=dict(), electrical=True,
#                                            base_plot_type='attenuation', start_seg=start_seg)
# # animation.ipython_display(fps=10, loop=True, autoplay=True)
# animation.write_videofile('./check_Rall2.mp4',fps=10, threads=4,audio=False)


record_dict, extra = analyser.record_protocol(protocol=Ca_spike_protocol,cut_start_ms=970.0, record_name='v')
animation = analyser.create_movie_from_rec(records=record_dict, slow_down_factor=500,
                                           func_for_missing_frames=np.max, theta=0, diam_factor=0.5,
                                           show_records_from=dict(), electrical=True,
                                           base_plot_type='attenuation')
# animation.ipython_display(fps=10, loop=True, autoplay=True)
# animation.write_videofile('./check_Rall_spikes4.mp4',fps=10, threads=4,audio=False)
animation.write_videofile('./Ca_spike_protocol.mp4',fps=10, threads=4,audio=False)
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