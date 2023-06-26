from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Neuron_analysis_tool.utils import LAMDA, MICRO, plot_shring_axes



analyser = Analyzer(type='L5PC')

def single_spike_protocol(cell, start_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt = spike_data.T[0][1] - spike_data.T[0][0]
    v = spike_data.T[1][800:2400]# the recorded spike voltage is for 40 ms
    v = v[:800] #crope only 20 ms
    start_time = 400
    extra_end = 0

    V = np.concatenate([np.zeros(int(start_time / dt)) + v[0]] + [v] + [np.zeros(int(extra_end / dt)) + v[-1]])
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)
    h.dt=dt
    h.steps_per_ms = 1.0 / h.dt
    clamp = h.SEClamp(start_seg.x, sec=start_seg.sec)
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    spike_vec.play(clamp._ref_amp1, spike_data.T[0][1]-spike_data.T[0][0])
    h.tstop = T[-1]
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return start_time-1, {}

records2, extra2 = analyser.record_protocol(protocol=single_spike_protocol, record_names=['v', 'gCa_LVAst_Ca_LVAst'], compute_more_condunctances=True)

more_conductances_ = extra2['more_conductances']

fig, ax = plt.subplots(2, 5, figsize=(16, 8), gridspec_kw={'width_ratios': [1, 0.5, 0.5, 0.5, 0.5]})
f = lambda: plt.tight_layout()
ax_morph = ax[0,0]
ax_dendogram = ax[1,0]
ax_v = ax[0,1]
ax_soma_cable = ax[0,2]
ax_oblique1_cable = ax[0,3]
ax_oblique2_cable = ax[0,4]
ax_basal1_cable = ax[1,1]
ax_basal2_cable = ax[1,2]
ax_tuft1_cable = ax[1,3]
ax_tuft2_cable = ax[1,4]

# fig = plt.figure(constrained_layout=True, figsize=(16, 8))
# gs = GridSpec(2, 5, figure=fig, width_ratios=[1, 0.5, 0.5, 0.5, 0.5])
# ax_morph = fig.add_subplot(gs[0, 0])
# ax_dendogram = fig.add_subplot(gs[1, 0])
#
# ax_v = fig.add_subplot(gs[0, 1])
#
# ax_soma_cable = fig.add_subplot(gs[0, 2])
# ax_oblique1_cable = fig.add_subplot(gs[0, 3])
# ax_oblique2_cable = fig.add_subplot(gs[0, 4])
# ax_basal1_cable = fig.add_subplot(gs[1, 1])
# ax_basal2_cable = fig.add_subplot(gs[1, 2])
# ax_tuft1_cable = fig.add_subplot(gs[1, 3])
# ax_tuft2_cable = fig.add_subplot(gs[1, 4])
# f = lambda: plt.subplots_adjust(wspace=0.6)
# f()

soma_seg = list(analyser.cell.soma[0])[0]
basal_tip1 = list(analyser.cell.dend[28])[-1]
basal_tip2 = list(analyser.cell.dend[62])[-1]

oblique_tip1 = list(analyser.cell.apic[102])[-1]
oblique_tip2 = list(analyser.cell.apic[19])[-1]
tuft_tip1 = list(analyser.cell.apic[67])[-1]
tuft_tip2 = list(analyser.cell.apic[52])[-1]

seg_to_indicate_dict = dict()
dot_size = 20
seg_to_indicate_dict[soma_seg] = dict(label='soma', alpha=0.75, color='k', size=dot_size)
seg_to_indicate_dict[basal_tip1] = dict(label='basal', alpha=0.75, color='r', size=dot_size)
seg_to_indicate_dict[basal_tip2] = dict(label='basal', alpha=0.75, color='r', size=dot_size)

seg_to_indicate_dict[oblique_tip1] = dict(label='oblique', alpha=0.75, color='cyan', size=dot_size)
seg_to_indicate_dict[oblique_tip2] = dict(label='oblique', alpha=0.75, color='cyan', size=dot_size)
seg_to_indicate_dict[tuft_tip1] = dict(label='tuft', alpha=0.75, color='b', size=dot_size)
seg_to_indicate_dict[tuft_tip2] = dict(label='tuft', alpha=0.75, color='b', size=dot_size)

plot_kwargs = [
    dict(ax=ax_morph, seg=soma_seg, records=records2.all_records['v'], electrical=True, plot_type='morph',
         plot_color_bar=True, theta=270, dancing=True, more_conductances_=more_conductances_,
         seg_to_indicate_dict=seg_to_indicate_dict, ),
    dict(ax=ax_dendogram, seg=soma_seg, records=records2.all_records['v'], electrical=True, plot_type='dendogram',
         plot_color_bar=False, dancing=True, more_conductances_=more_conductances_,
         seg_to_indicate_dict=seg_to_indicate_dict, ),
    dict(ax=ax_v, seg=soma_seg, records=records2.all_records['v'], plot_type='single_record', color='grey',
         title='soma voltage', ylabel='V (mV)'),

    dict(ax=ax_soma_cable, seg=soma_seg, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='soma cable', start_color='k', dots_size=dot_size),

    dict(ax=ax_basal1_cable, seg=basal_tip1, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='basal cable 1', start_color='r', dots_size=dot_size),

    dict(ax=ax_basal2_cable, seg=basal_tip2, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='basal cable 2', start_color='r', dots_size=dot_size),

    dict(ax=ax_oblique1_cable, seg=oblique_tip1, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='oblique cable 1', start_color='cyan', dots_size=dot_size),

    dict(ax=ax_oblique2_cable, seg=oblique_tip2, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='oblique cable 2', start_color='cyan', dots_size=dot_size),

    dict(ax=ax_tuft1_cable, seg=tuft_tip1, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='tuft cable 1', start_color='b', dots_size=dot_size),

    dict(ax=ax_tuft2_cable, seg=tuft_tip2, plot_type='cable', cable_type='d3_2', shift=0,
         factor_e_space=25, more_conductances_=more_conductances_,
         scales=None, plot_legend=False, title='tuft cable 2', start_color='b', dots_size=dot_size),

]
slow_down_factor = 1000
videos_folder = 'videos/L5PC/'
video_name = 'bAP_cables_new_more_condactances3.mp4'
os.makedirs(videos_folder, exist_ok=True)

cable_limets_func = plot_shring_axes(plot_kwargs, [3, 4, 5, 6, 7, 8, 9])
analyser.save_movie_from_rec(fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs, func_before_run=[f], func_during_run=[cable_limets_func],
                             save_to=videos_folder, clip_name=video_name, fps=10,
                             threads=16, preset='ultrafast')
