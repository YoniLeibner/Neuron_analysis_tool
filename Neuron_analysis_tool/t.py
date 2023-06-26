from neuron import h, gui
import matplotlib.pyplot as plt
import pickle
import os
import sys
import numpy as np
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new')
from utiles_func import get_cell2
from utils import fig_names as names
from Neuron_analysis_tool.Analyzer import Analyzer
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from Neuron_analysis_tool.record import record_all
from Neuron_analysis_tool.attenuation import record_to_value
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib as mpl


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

a='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/'

def long_pulse_protocol(cell, start_seg):
    delay = 2000.0
    dur = 1000.0
    amp = .1
    h.tstop = delay + dur + 500.0
    clamp = h.IClamp(start_seg.x, sec=start_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return delay, {}


def short_pulse_protocol(cell, start_seg):
    delay = 2000.0
    dur = 10.0
    amp = .1
    h.tstop = delay + dur + 500.0
    clamp = h.IClamp(start_seg.x, sec=start_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return delay, {}

def very_short_pulse_protocol(cell, start_seg):
    delay = 2000.0
    dur = 2.0
    amp = .1
    h.tstop = delay + dur + 500.0
    clamp = h.IClamp(start_seg.x, sec=start_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return delay, {}

def resting_protocol2(cell, start_seg=None):
    h.tstop = 1500.0
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}


def run_cell(number, ax, load_num=0, dur=1000.0, remove_obliques=False, remove_Ih=False, plot_morph=False, save_folder=None):
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)
    cell = bbp_loader.get_model()
    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]
    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    tuft_sec = np.unique([seg.sec for seg in tuft])
    basal_terminals = [sec for sec in cell.dend if len(sec.children()) == 0]
    obliqes_terminals = [sec for sec in obliqes_sec if len(sec.children()) == 0]
    tuft_terminals = [sec for sec in tuft_sec if len(sec.children()) == 0]

    if remove_obliques:

        obliques_start = []
        for sec in obliqes_sec:
            if sec.parentseg().sec not in obliqes_sec:
                obliques_start.append(sec)
        for sec in obliques_start:
            h.disconnect(sec=sec)
        ignore_sec = obliqes_sec
        # obliqes = []
    else:
        ignore_sec=[]

    if remove_Ih:
        for sec in cell.all:
            try:
                sec.gIhbar_Ih_human_linear = 0
                sec.gIhbar_exp_Ih_human_linear = 0
                sec.dist_Ih_human_linear = 0
            except:
                print(sec, 'no Ih here')


    parts_dict = dict(tuft=tuft, trunk=trunk, oblique=obliqes, basal=basal, axon=axon, soma=[seg for seg in cell.soma[0]])

    colors_dict = {"soma": "silver",
                   "tuft": "#35BD99",
                   "oblique": "#4A9ABE",
                   "trunk": "#1F1A54",
                   "basal": "#C12328",
                   "axon": "green"}


    analyser = Analyzer(cell=cell, parts_dict=parts_dict, colors_dict=colors_dict,
                        more_conductances_protocol=resting_protocol2)

    soma_seg = list(cell.soma[0])
    soma_seg = soma_seg[len(soma_seg) // 2]
    dots_size = 80
    seg_to_indicate_ = {
        hot_spot[0]: dict(size=dots_size, color='#B37214', alpha=1, lw=2),
        hot_spot[1]: dict(size=dots_size, color='#FCB653', alpha=1, lw=2),
        soma_seg: dict(size=dots_size, color='k', alpha=1, lw=2),
    }
    if plot_morph:
        ax_morph = inset_axes(ax, width="25%", height="80%", loc=1, borderpad=0)
        analyser.plot_morph(ax=ax_morph, seg_to_indicate_dict = seg_to_indicate_, diam_factor=None, sec_to_change=None,
                   ignore_sections=ignore_sec,
                   theta=more_azim, scale=0, scale_text=False, ignore_soma=True, distance=None, electrical=False)

    assert dur in [1000.0, 10.0, 2]
    if dur==1000.0:
        protocol = long_pulse_protocol
    elif dur==10.0:
        protocol = short_pulse_protocol
    else:
        protocol = very_short_pulse_protocol

    save_name = str(dur)
    if remove_obliques:
        save_name += '_no_obliques'
    if remove_Ih:
        save_name += '_no_Ih'

    load_flag=os.path.isfile(save_folder + save_name+'.p')
    records=None
    if load_flag:
        records = record_all(cell, record_name='v')
        records.load(save_folder + save_name+ '.p')
    ax, norm_by, lines, segs, records = analyser.plot_attenuation(protocol=protocol, ax=ax, seg_to_indicate_dict=seg_to_indicate_, start_seg = None, record_to_value_func=record_to_value,
                                                                  norm=True, record_name='v', norm_by=None, electrical=True, records=records, ignore_sections=ignore_sec)
    if not load_flag and save_folder is not None:
        records.save(save_folder + save_name+'.p')


    all_tip_attenuation = dict()
    for terminals, tip_name in zip([basal_terminals, obliqes_terminals, tuft_terminals], ['basal', 'oblique', 'tuft']):
        tip_attenuation = []
        for terminal in terminals:
            tip_attenuation.append(record_to_value(records.get_record(list(terminal)[-1])) / norm_by)
        all_tip_attenuation[tip_name]=tip_attenuation
    bbp_loader.model.destroy(bbp_loader.sim)
    return all_tip_attenuation


if __name__ == "__main__":
    if len(sys.argv)>1:
        cell_number=int(sys.argv[1])
    else:
        cell_number=17
    print('runing cell num:', cell_number)
    try:
        os.mkdir('sup_fig1/'+str(cell_number))
    except:
        pass
    # short_pulse_dur=10.0
    short_pulse_dur=2.0

    save_folder = 'sup_fig1/'+str(cell_number)+'/'
    fig = plt.figure(figsize=(5*4, 5*2))
    ax = fig.add_gridspec(2, 4, width_ratios=[1,1,1,1])
    # morphs
    ax_ss = fig.add_subplot(ax[0, 0])
    ax_ss_no_obliques = fig.add_subplot(ax[0, 1], sharex=ax_ss, sharey=ax_ss)
    ax_ss_no_Ih = fig.add_subplot(ax[0, 2], sharex=ax_ss, sharey=ax_ss)

    ax_short = fig.add_subplot(ax[1, 0], sharex=ax_ss, sharey=ax_ss)
    ax_short_no_obliques = fig.add_subplot(ax[1, 1], sharex=ax_ss, sharey=ax_ss)
    ax_short_no_Ih = fig.add_subplot(ax[1, 2], sharex=ax_ss, sharey=ax_ss)

    ax_stats_ss = fig.add_subplot(ax[0, 3])
    ax_stats_short = fig.add_subplot(ax[1, 3])#, sharex=ax_stats_ss, sharey=ax_stats_ss)

    all_tip_attenuation1 = run_cell(cell_number, ax_ss, load_num=0, dur=1000.0, remove_obliques=False, remove_Ih=False, plot_morph=False, save_folder=save_folder)
    all_tip_attenuation2 = run_cell(cell_number, ax_ss_no_obliques, load_num=0, dur=1000.0, remove_obliques=True, remove_Ih=False, plot_morph=False, save_folder=save_folder)
    all_tip_attenuation3 = run_cell(cell_number, ax_ss_no_Ih, load_num=0, dur=1000.0, remove_obliques=False, remove_Ih=True, plot_morph=False, save_folder=save_folder)

    all_tip_attenuation4 = run_cell(cell_number, ax_short, load_num=0, dur=short_pulse_dur, remove_obliques=False, remove_Ih=False, plot_morph=True, save_folder=save_folder)
    all_tip_attenuation5 = run_cell(cell_number, ax_short_no_obliques, load_num=0, dur=short_pulse_dur, remove_obliques=True, remove_Ih=False, plot_morph=True, save_folder=save_folder)
    all_tip_attenuation6 = run_cell(cell_number, ax_short_no_Ih, load_num=0, dur=short_pulse_dur, remove_obliques=False, remove_Ih=True, plot_morph=False, save_folder=save_folder)

    for a in [ax_ss_no_obliques, ax_ss_no_Ih, ax_short_no_obliques, ax_short_no_Ih]:
        a.set_ylabel('')
    ax_ss.set_ylabel('attenuation (%)')
    ax_short.set_ylabel('attenuation (%)')
    for a in [ax_ss, ax_ss_no_obliques, ax_ss_no_Ih, ax_short, ax_short_no_Ih]:
        a.set_xlabel('')
    for a in [ax_ss, ax_ss_no_obliques, ax_ss_no_Ih, ax_short, ax_short_no_obliques, ax_short_no_Ih, ax_stats_ss, ax_stats_short]:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)

    colors_dict = {"soma": "k",
                   "tuft": "#35BD99",
                   "oblique": "#4A9ABE",
                   "trunk": "#1F1A54",
                   "basal": "#C12328",
                   "hot-spot":'#B37214',
                   "hot-spot2":'#FCB653'
                   }
    legend_elements = [Line2D([0], [0], color=colors_dict[label], lw=2, label=label) for label in colors_dict]
    ax_short_no_Ih.legend(handles=legend_elements, loc="best")

    pickle.dump(dict(
        short_pulse=dict(
            control=all_tip_attenuation4,
            no_obliques=all_tip_attenuation5,
            no_ih=all_tip_attenuation6,
        ),
        long_pulse=dict(
            control=all_tip_attenuation1,
            no_obliques=all_tip_attenuation2,
            no_ih=all_tip_attenuation3,
        ),
    ), open(save_folder+'tip_data2.p', 'wb'))
    plt.tight_layout()
    plt.savefig(save_folder+'temp2.png')
    plt.savefig(save_folder+'temp.pdf')
