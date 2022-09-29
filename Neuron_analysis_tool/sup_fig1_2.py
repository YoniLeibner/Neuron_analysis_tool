from neuron import h, gui
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from Load_BBP_Model import Load_BBP_Model
from open_morph_with_obliques import morph_reader, get_section_nums
from cables import get_cable
from cables_bi_direction import get_cable as bi_direction_get_cable
from dendogram_class import Dendogram, more_conductances, color_func, get_segment_length_lamda,  get_segment_length_um, more_conductances_fake
from morph_plotter_2d import plot_morph
from dendogram_class_bi_direction import Dendogram as bi_direction_Dendogram
from extraClasses import eiline_resonance
from attenuation_plotter import run_attenuation_ploter
import json
from utils import colors_dict
from utils import fig_names as names
from glob import glob
import seaborn as sns
import pandas as pd

a='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/'
# a='C:/Users/USER1/Documents/Hippocampus_new/'
def get_cell2(person, cell_name, seg_every_um=10, load_num=0):
    morph_opener = morph_reader(df_loc = a+'new_data/table.xlsx',
                                base_folder=a+'new_data/morphs/')
    morphology_dict = morph_opener.get_cell(person=person, cell_name=cell_name)
    morphology_dirr = morphology_dict['morph_path']

    parameter_folder = 'sag_models/'
    sag_model_add = '_linear'
    bbp_loader = Load_BBP_Model(morphology_dirr,
                                parameter_folder + 'mechanism' + sag_model_add + '.txt',
                                parameter_folder + 'parameters' + sag_model_add + '.txt',
                                chnnels_loc=None)
    bbp_loader.load_hall_of_fame(load_num, base_folder + person + '_' + cell_name + '/hall_of_fame.p')
    cell = bbp_loader.get_model()
    for sec in cell.all:
        sec.nseg = max(min(int(sec.L / seg_every_um), 32767), 1)

    obliques = [seg for sec in get_section_nums(cell, morphology_dict['obliques']) for seg in sec]

    tuft = []
    for sec in get_section_nums(cell, morphology_dict['tuft']):
        if int(sec.name().split('[')[-1][:-1]) in morphology_dict['tuft-trunk']: continue
        for seg in sec:
            tuft.append(seg)

    trunk = []
    for sec in get_section_nums(cell, morphology_dict['trunk']):
        if int(sec.name().split('[')[-1][:-1]) in morphology_dict['tuft-trunk']: continue
        for seg in sec:
            trunk.append(seg)
    hot_spot = []
    for sec_num in morphology_dict['tuft-trunk']:
        sec = cell.apic[sec_num]
        hot_spot.append(sec(morphology_dict['tuft-trunk'][sec_num]['x']))
        for seg in sec:
            if seg.x < morphology_dict['tuft-trunk'][sec_num]['x']:
                trunk.append(seg)
            else:
                tuft.append(seg)
    return bbp_loader, obliques, trunk, tuft, hot_spot


def plot_morph_(cell, colors, seg_to_indicate, colors_dict, num=0, more_elev=0, more_azim = 0, ax = None, fig=None,
                diam_factor = None, sec_to_change =None, bounds=None, cmap=plt.cm.turbo, norm_colors=False, plot_color_bar=True, ignore_sections=[]):
    if ax is None:
        fig = plt.figure(figsize=(20, 20))
        ax = plt.axes(projection='3d')

    seg_to_indicate_dict = {seg: dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in seg_to_indicate}
    seg_to_indicate_dict[list(cell.soma[0])[0]] = dict(size=30, color=colors_dict['soma'], alpha=1)

    fig, ax, color_bar, points_dict = plot_morph(cell, color_func=colors.get_seg_color, scatter=False, add_nums=False,
                                                 seg_to_indicate=seg_to_indicate_dict,
                                                 norm_colors=norm_colors, fig=fig, ax=ax, diam_factor = diam_factor,
                                                 sec_to_change =sec_to_change, bounds=bounds, cmap=cmap, plot_color_bar=plot_color_bar,
                                                 theta=more_azim, ignore_sections=ignore_sections)
    # remove grid and axis
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()
    # ax.plot([addition[0], addition[0]+direction_vector[0]], [addition[1], addition[1]+direction_vector[1]], [addition[2], addition[2]+direction_vector[2]], color='k')
    ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0]+500], color='k')

    # ax.view_init(elev=90+more_elev, azim=90+more_azim)
    # plt.draw()
    return fig, ax


def get_branch_init(seg, seg_list, prev):
    if seg.sec not in seg_list:
        return seg, prev
    return get_branch_init(seg.sec.parentseg(), seg_list, prev=seg.sec)


def run_cell(person, cell_name, fig, ax_arr, load_num=0, dur=1000.0):

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)
    cell = bbp_loader.get_model()
    json.dump(bbp_loader.all_res_hall[load_num], open(os.path.join(save_folder, person+'_'+cell_name+'_params.json'), 'w'), indent=4, sort_keys=True)
    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]


    parts_dict = dict(tuft=tuft, trunk=trunk, oblique=obliqes, basal=basal, axon=axon, soma=[seg for seg in cell.soma[0]])
    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    basal_sec = np.unique([seg.sec for seg in basal])
    tuft_sec = np.unique([seg.sec for seg in tuft])

    obliqes_terminals = []
    for sec in obliqes_sec:
        if len(sec.children()) == 0:
            obliqes_terminals.append(sec)

    obliques_start = []
    for sec in obliqes_sec:
        if sec.parentseg().sec not in obliqes_sec:
            obliques_start.append(sec)

    for sec in obliques_start:
        h.disconnect(sec=sec)
    # colors_dict['oblique'] = 'w'
    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)
    fig, ax_arr[0] = plot_morph_(cell, colors, seg_to_indicate=hot_spot, colors_dict=colors_dict, num=0, more_elev=more_elev, more_azim = more_azim, ax = ax_arr[0], diam_factor=None, ignore_sections=obliqes_sec)
    START_TIME = 1500
    more_conductances_ = more_conductances(cell, run_time=3000)

    seg = cell.soma[0](0.5)
    ax_arr[1] = run_attenuation_ploter(cell=cell,
                                       seg_start=list(cell.soma[0])[0],
                                       color_func=colors,
                                       seg_length_function = get_segment_length_lamda,
                                       more_conductances_=more_conductances_,
                                       param_to_record='v',
                                       record_to_value_func=None,
                                       norm=True,
                                       delay=2000.0,
                                       dur=dur,
                                       amp=0.1,
                                       seg_to_indicate={seg:dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot},
                                       ax = ax_arr[1])

    ax_arr[1].spines['top'].set_visible(False)
    ax_arr[1].spines['right'].set_visible(False)
    bbp_loader.model.destroy(bbp_loader.sim)

    return fig, ax_arr

def run_cell_no_Ih(person, cell_name, fig, ax_arr, load_num=0, dur=1000.0):

    colors_dict['oblique'] = '#4A9ABE'

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)
    cell = bbp_loader.get_model()
    json.dump(bbp_loader.all_res_hall[load_num], open(os.path.join(save_folder, person+'_'+cell_name+'_params.json'), 'w'), indent=4, sort_keys=True)
    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]


    parts_dict = dict(tuft=tuft, trunk=trunk, oblique=obliqes, basal=basal, axon=axon, soma=[seg for seg in cell.soma[0]])
    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    basal_sec = np.unique([seg.sec for seg in basal])
    tuft_sec = np.unique([seg.sec for seg in tuft])

    obliqes_terminals = []
    for sec in obliqes_sec:
        if len(sec.children()) == 0:
            obliqes_terminals.append(sec)

    obliques_start = []
    for sec in obliqes_sec:
        if sec.parentseg().sec not in obliqes_sec:
            obliques_start.append(sec)

    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)
    fig, ax_arr[0] = plot_morph_(cell, colors, seg_to_indicate=hot_spot, colors_dict=colors_dict, num=0, more_elev=more_elev, more_azim = more_azim, ax = ax_arr[0], diam_factor=None)
    START_TIME = 1500
    more_conductances_ = more_conductances(cell, run_time=3000)

    seg = cell.soma[0](0.5)
    ax_arr[1] = run_attenuation_ploter(cell=cell,
                                       seg_start=list(cell.soma[0])[0],
                                       color_func=colors,
                                       seg_length_function = get_segment_length_lamda,
                                       more_conductances_=more_conductances_,
                                       param_to_record='v',
                                       record_to_value_func=None,
                                       norm=True,
                                       delay=2000.0,
                                       dur=dur,
                                       amp=0.1,
                                       seg_to_indicate={seg:dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot},
                                       ax = ax_arr[1])

    ax_arr[1].spines['top'].set_visible(False)
    ax_arr[1].spines['right'].set_visible(False)


    for sec in cell.all:
        try:
            sec.gIhbar_Ih_human_linear = 0
        except:
            print('no Ih at:', sec.name())
    more_conductances_ = more_conductances_fake(cell)
    ax_arr[2] = run_attenuation_ploter(cell=cell,
                                       seg_start=list(cell.soma[0])[0],
                                       color_func=colors,
                                       seg_length_function = get_segment_length_lamda,
                                       more_conductances_=more_conductances_,
                                       param_to_record='v',
                                       record_to_value_func=None,
                                       norm=True,
                                       delay=2000.0,
                                       dur=dur,
                                       amp=0.1,
                                       seg_to_indicate={seg:dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot},
                                       ax = ax_arr[2])

    ax_arr[2].spines['top'].set_visible(False)
    ax_arr[2].spines['right'].set_visible(False)

    bbp_loader.model.destroy(bbp_loader.sim)

    return fig, ax_arr

def plot_lamda_parts(ax):
    df = pd.read_csv('lamda_length_terminals_Ih_f0.csv')
    df['name'] = df.person + df.cell_name
    df = df[~(df.part=='trunk')]
    sns.swarmplot(data=df, x="part", y="L", hue="name", ax=ax[0],
                  palette=sns.color_palette("hls", len(pd.unique(df.name))))
    ax[0].legend([], [], frameon=False)
    ax[0].set_title('control')
    print('control')
    for part in pd.unique(df.part):
        print(part, df.L[df.part == part].mean(), df.L[df.part == part].std())

    df = pd.read_csv('lamda_length_terminals_passive_f0.csv')
    df['name'] = df.person + df.cell_name
    df = df[~(df.part == 'trunk')]
    sns.swarmplot(data=df, x="part", y="L", hue="name", ax=ax[1],
                  palette=sns.color_palette("hls", len(pd.unique(df.name))))
    ax[1].legend([], [], frameon=False)
    ax[1].set_ylabel('')
    ax[1].set_title('ZD')
    print('ZD')
    for part in pd.unique(df.part):
        print(part, df.L[df.part == part].mean(), df.L[df.part == part].std())

    for ax_ in ax:
        ax_.spines['top'].set_visible(False)
        ax_.spines['right'].set_visible(False)

if __name__ == "__main__":

    base_folder = 'optim_ZAP_power_new6_Ih_linear/'
    parameter_folder = 'sag_models/'
    load_num=0
    dur = 1000.0
    # dur = 20.0
    # dur = 10.0
    # dur = 2.0
    try: os.mkdir('graphs_for_fig1')
    except: pass
    print('dur=', dur)
    save_folder = 'graphs_for_fig1/'
    # morph morph_no obliques, attenuation normal, attenuation no Ih, attenuation no obliques, distance lamdas
    fig, ax_arr = plt.subplots(2, 6, figsize=(5*6, 5*2), gridspec_kw={'width_ratios': [2, 2, 2, 2, 2, 1]}, sharex='col', sharey='col')
    if dur == 1000.0:
        plot_lamda_parts(ax_arr[:, -1])
    number = 1
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]
    run_cell(person , cell_name, fig=fig, ax_arr=ax_arr[0, [1, 4]], load_num=0, dur=dur)
    run_cell_no_Ih(person , cell_name, fig=fig, ax_arr=ax_arr[0, [0, 2, 3]], load_num=0, dur=dur)

    number = 5
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]
    run_cell(person, cell_name, fig=fig, ax_arr=ax_arr[1, [1, 4]], load_num=0, dur=dur)
    run_cell_no_Ih(person, cell_name, fig=fig, ax_arr=ax_arr[1, [0, 2, 3]], load_num=0, dur=dur)


    for ax in ax_arr[0, [2, 3, 4]]:
        ax.set_ylim(10**-2, 10**0)
        if not dur == 1000.0:
            ax.set_xlim(-0.1, 2)
    plt.savefig(os.path.join(save_folder, 'sup_fig1_'+str(dur)+'.png'))
    plt.savefig(os.path.join(save_folder, 'sup_fig1_'+str(dur)+'.pdf'))
    plt.close()
