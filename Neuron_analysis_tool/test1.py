from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Neuron_analysis_tool.utils import LAMDA, MICRO, plot_shring_axes

from glob import glob
fig, ax = plt.subplots(1, 2)
f='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/DATA/P9/Hipp5/H22.29.209.21.03.04.ASC'
analyser2 = Analyzer(type='L5PC')

# f='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/DATA/P9/Hipp6/H22.29.209.21.03.04_obliquesremoved_metbasals.ASC'
folder_ = '/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/DATA/P9/Hipp1/'
print(glob(folder_+'*.asc'))
print(glob(folder_+'*.ASC'))
# f=glob(folder_+'*.asc')[0]
# f=glob(folder_+'*.ASC')[0]
# f='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/DATA/P9/Hipp4/H21.29.198.21.02.09_eline_edit.asc'
analyser = Analyzer(type='ASC', morph_path=f)
_=analyser.plot_morph(theta=50, scale=500, ax=ax[0])
_=analyser2.plot_morph(theta=50, scale=500, ax=ax[1])

plt.show()
exit(0)
# analyser = Analyzer(type='Rall_tree')
# colors_dict  = analyser.colors_dict
# colors_dict['soma']='r'
# colors_dict['basal']='pink'
#
# analyser.change_color_dict(colors_dict)
# fig, ax = analyser.create_morpho_card(scale=500, start_seg=list(analyser.cell.soma[0])[0], diam_factor=1)
# plt.show()
# exit(0)

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

    start_time_ = 100
    total_time = 300 + start_time_  # 1400

    number_of_e_syns = 1000
    number_of_i_syns = 500

    amp = 0
    e_FR = 4  # Hz
    i_FR = 4  # Hz
    # difining the synapses times
    sim_time = total_time - start_time_
    for sec in np.random.choice(list(cell.all), number_of_e_syns):
        seg_num = np.random.randint(0, len(list(sec)))
        segs_e.append(list(sec)[seg_num])
        syns.append(add_syn(list(sec)[seg_num]))
        syn_time = np.random.rand(int(np.ceil(sim_time / 1000.0 * e_FR))) * sim_time + start_time_
        syn_e_times.append(syn_time)

    for sec in np.random.choice(list(cell.all), number_of_i_syns):
        seg_num = np.random.randint(0, len(list(sec)))
        segs_i.append(list(sec)[seg_num])
        isyns.append(add_isyn(list(sec)[seg_num]))
        syn_time = np.random.rand(int(np.ceil(sim_time / 1000.0 * i_FR))) * sim_time + start_time_
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
        start_time += start_time_
        end_time += start_time_

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
    return start_time_, dict(draw_funcs=[draw_func1])


# analyser = Analyzer(type='Rall_tree', seg_every=10)
# colors_dict  = analyser.colors_dict
# colors_dict['soma']='r'
# colors_dict['basal']='pink'
# analyser.change_color_dict(colors_dict)
analyser = Analyzer(type='L5PC', seg_every=10)
# dir_path = os.path.abspath('')
# morph_path=os.path.join(dir_path,'data/morph.ASC')
# analyser = Analyzer(type='ASC', morph_path=morph_path)

print('start protocol')

records, extra = analyser.record_protocol(protocol=random_syn_protocol, record_names=['g_total_g_total','g_syn_g_total', 'v'], compute_more_condunctances=True)

# fig, ax = plt.subplots(1, 5, figsize=(25,5))
fig, ax = plt.subplots(1, 3, figsize=(25,5))
soma_seg = list(analyser.cell.soma[0])[0]

record_dict = records.all_records['g_syn_g_total']
record_dict2 = records.all_records['v']
print('min, max:', record_dict.get_min(), record_dict.get_max())
print('min, max:', record_dict2.get_min(), record_dict2.get_max())
draw_funcs = extra['draw_funcs']
more_conductances_ = extra['more_conductances']

plot_kwargs = [
    # dict(ax=ax[0], seg = soma_seg, records=record_dict, electrical=True, plot_type='morph',
    #      plot_color_bar=True, theta=0, draw_funcs=draw_funcs, dancing=True, more_conductances_=more_conductances_),
    dict(ax=ax[0], seg = soma_seg, records=record_dict2, electrical=True, plot_type='dendogram',
         plot_color_bar=False, draw_funcs=draw_funcs, dancing=True, more_conductances_=more_conductances_),
    dict(ax=ax[1], seg = soma_seg, records=record_dict, distance_factor=0, plot_every=0.25, plot_type='all_records'),
    dict(ax=ax[2], seg = soma_seg, records=record_dict2, distance_factor=0, plot_every=0.25, plot_type='all_records'),
    # dict(ax=ax[1], seg = soma_seg, records=record_dict, electrical=True, plot_type='attenuation'),
    # dict(ax=ax[2], seg = soma_seg, records=record_dict2, electrical=True, plot_type='attenuation'),
    # dict(ax=ax[3], seg = soma_seg, records=record_dict, distance_factor=0, plot_every=0.25,
    #      plot_type='all_records'),
    # dict(ax=ax[4], seg = soma_seg, records=record_dict, plot_type='single_record'),
              ]
f = lambda: plt.tight_layout()
slow_down_factor=100
videos_folder = 'videos/synapses/'
video_name = 'dancing_synapses_rall_tree_many_synapses3.mp4'
# video_name = 'dancing_synapses.mp4'
os.makedirs(videos_folder, exist_ok=True)
analyser.save_movie_from_rec(fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs, func_before_run=[f],
                             save_to=videos_folder, clip_name=video_name, fps=10,
                             threads=16, preset='ultrafast')
# video_player(Path.cwd(), videos_folder+video_name)