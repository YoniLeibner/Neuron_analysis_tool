from neuron import h, gui
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from Load_BBP_Model import Load_BBP_Model
from open_morph_with_obliques import morph_reader, get_section_nums
from cables import get_cable
from cables_bi_direction import get_cable as bi_direction_get_cable
from dendogram_class import Dendogram, more_conductances, color_func, get_segment_length_lamda,  get_segment_length_um
from morph_plotter_2d import plot_morph
from dendogram_class_bi_direction import Dendogram as bi_direction_Dendogram
from extraClasses import eiline_resonance
from attenuation_plotter import run_attenuation_ploter
import json
from utils import colors_dict
from utils import fig_names as names
from glob import glob

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


def cable_option1_parts(ax, results, factor, factor_e_space, x_pos, cable_type='electric', s=50):
    # max_x = 0
    max_cable = 0
    shift = 0
    for part, direction in zip(results.keys(), [1, -1]):
        cable = results[part]['all'][cable_type].flatten()
        if max_cable < cable.max() / 2:
            max_cable = cable.max() / 2
            # max_x = x_pos + cable.max() / 2
            shift = x_pos + cable.max() / factor / 2 + 5
    for part, direction in zip(results, [1, -1]):
        cable = results[part]['all'][cable_type].flatten()
        if cable_type == 'd3_2':
            befor_d_3_2 = np.power(cable, 3.0/2.0)
        cable /= factor
        y = np.arange(0, len(cable), 1) / factor_e_space

        start_pos = -cable / 2.0 + shift

        for morph_part in ['basal', 'tuft', 'trunk', 'oblique']:
            part_cable = results[part][morph_part][cable_type].flatten() / factor
            if cable_type == 'd3_2':
                part_cable = results[part][morph_part][cable_type].flatten()
                part_cable_befor_d_3_2=np.power(part_cable, 3.0/2.0)
                part_cable = cable*(part_cable_befor_d_3_2/befor_d_3_2)
            plot_cable = part_cable[part_cable > 0]
            start_pos_temp = start_pos[part_cable > 0]
            ax.fill_betweenx(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable,
                             label=morph_part, color=colors_dict[morph_part])
            start_pos += part_cable

    ax.scatter(shift, 0, s=s, color='k', label='start')
    # add_scale(ax, x_=int(50 * 100 * factor / (factor_e_space)), y_=1)  # for 100 x and 1 y

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='best')
    # ax.legend(loc='best')
    ax.set_ylabel('distanse')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()


def cable_option1(ax, results, factor, factor_e_space, x_pos, cable_type = 'electric'):
    cable = results['all'][cable_type].flatten()

    cable /= factor
    y = np.arange(0, len(cable), 1) / factor_e_space
    shift = x_pos + cable.max() / 2 + 5
    max_x = x_pos + cable.max() / 2
    start_pos = -cable / 2.0 + shift

    for part in ['tuft', 'trunk', 'oblique', 'basal']:
        part_cable = results[part][cable_type].flatten() / factor
        plot_cable = part_cable[part_cable > 0]
        start_pos_temp = start_pos[part_cable > 0]
        ax.fill_betweenx(y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable, label=part,
                         color=colors_dict[part])
        start_pos += part_cable

    ax.scatter(shift, 0, s=150, color='k', label='soma')
    # add_scale(ax, x_=int(50 * 100 * factor / (factor_e_space)), y_=1)  # for 100 x and 1 y

    # ax.legend(loc='best')
    ax.set_ylabel('distanse')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()

def plot_morph_(cell, colors, seg_to_indicate, colors_dict, num=0, more_elev=0, more_azim = 0, ax = None, fig=None,
                diam_factor = None, sec_to_change =None, bounds=None, cmap=plt.cm.turbo, norm_colors=False, plot_color_bar=True):
    if ax is None:
        fig = plt.figure(figsize=(20, 20))
        ax = plt.axes(projection='3d')
    if type(seg_to_indicate)==list:
        seg_to_indicate_dict = {seg: dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in seg_to_indicate}
        seg_to_indicate_dict[list(cell.soma[0])[0]] = dict(size=30, color=colors_dict['soma'], alpha=1)
    elif type(seg_to_indicate)==dict:
        print(seg_to_indicate)
        seg_to_indicate_dict = seg_to_indicate
    else:
        raise BaseException('seg_to_indicate is in the wrong type')
    fig, ax, color_bar, points_dict = plot_morph(cell, color_func=colors.get_seg_color, scatter=False, add_nums=False,
                                                 seg_to_indicate=seg_to_indicate_dict,
                                                 norm_colors=norm_colors, fig=fig, ax=ax, diam_factor = diam_factor,
                                                 sec_to_change =sec_to_change, bounds=bounds, cmap=cmap, plot_color_bar=plot_color_bar,
                                                 theta=more_azim)
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


def run_cell(person, cell_name, fig, ax_arr=None, load_num=0):
    if ax_arr is None:
        fig, ax_arr = plt.subplots(1, 4, figsize=(20, 5), gridspec_kw={'width_ratios': [1, 1, 1, 1]})
        ax_arr[0].remove()
        ax_arr[0] = fig.add_subplot(1, 4, 1, projection='3d')

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


    fig, ax_arr[0] = plot_morph_(cell, colors, seg_to_indicate=hot_spot, colors_dict=colors_dict, num=0, more_elev=more_elev, more_azim = more_azim, ax = ax_arr[0], diam_factor=None)
    # START_TIME = 10
    # more_conductances_ = more_conductances(cell, run_time=10)
    START_TIME = 1500
    more_conductances_ = more_conductances(cell, run_time=1500)

    seg = cell.soma[0](0.5)
    d = bi_direction_Dendogram(cell,
                               seg_length_function=get_segment_length_lamda,
                               color_func=colors,
                               dots_loc=[[seg_.sec.name().split('.')[-1], seg_.x] for seg_ in hot_spot],
                               more_conductances=more_conductances_,
                               diam_factor=None, s=10, fix_diam=1.)

    d.cumpute_distances(seg)
    max_y, x_pos = d.plot(ax=ax_arr[2], plot_legend=False, ignore_sections=[])
    seg_dist_dict = dict()
    for part in ['sons', 'parent']:
        seg_dist_dict[part] = dict()
        seg_dist_dict[part][cell.soma[0](0.5)] = []
        # seg_dist = {cell.soma[0](0.5): []}
        for hh in hot_spot:
            seg_dist_dict[part][hh] = []
    cross_dist_dict = []
    factor_e_space=25
    results, seg_dist, cross_dist_dict = bi_direction_get_cable(cell,
                                                                factor_e_space=factor_e_space,
                                                                factor_m_space=10,
                                                                start_section=seg.sec,
                                                                x_start=seg.x,
                                                                more_conductions=more_conductances_,
                                                                seg_dist_dict=seg_dist_dict,
                                                                cross_dist_dict=cross_dist_dict,
                                                                part_dict=parts_dict,
                                                                ignore_sections=[])

    # the bin size is 1/25=0.04
    # i reduce the diameter by a factor of 100 for visualizations
    # so dA/dX is 25*100 = 2500
    # I use scale bar of 10, so the its 10*2500=25000 = 0.25mm^2 (factor of 1e-6)
    cable_option1_parts(ax_arr[2], results, 0.04*25*100, factor_e_space, x_pos, cable_type='electric')
    # cable_option1_parts(plt.gca(), results,1, 25, x_pos, cable_type='d3_2', s=25)

    ax_arr[2].plot([-10, -10], [0, 1], color='k')
    # ax_arr[2].text(-15, 0.25, '1 lamda', rotation=90)

    ax_arr[2].plot([x_pos+40, x_pos+50], [-0.2, -0.2], color='k')
    # ax_arr[2].text(x_pos+35, -0.8, 'dA/dX\n0.025 mm^2')


    clamp = h.IClamp(0.5, sec=cell.soma[0])
    #set records
    soma_V = h.Vector()
    soma_V.record(cell.soma[0](0.5)._ref_v)
    h_T = h.Vector()
    h_T.record(h._ref_t)

    IV_folder = 'IV_traces/'
    IV_data = pickle.load(open(IV_folder + person + '_' + cell_name + '/data.p', 'rb'))
    IV_T = IV_data['T']
    IV_I = IV_data['I']
    IV_V = IV_data['V']
    T_filter = []
    I_filter = []
    V_filter = []
    for i, I in enumerate(IV_I):  # filter negative currents
        if I.max() < 10:
            T_filter.append(IV_T[i].copy())
            I_filter.append(IV_I[i].copy())
            t_v = IV_V[i].copy()
            t_v[:100] = t_v[101]
            V_filter.append(t_v)

    all_E_PAS = []
    for T, I, V in zip(T_filter, I_filter, V_filter):
        E_PAS = np.median(V[:100])
        all_E_PAS.append(E_PAS)
    E_PAS = np.mean(all_E_PAS)
    # run long pulses
    amps = []
    max_x = 2000
    for T_, I_, V_ in zip(T_filter, I_filter, V_filter):
        h.tstop = max_x + START_TIME
        I_abs = np.abs(I_)
        if I_abs.max() == 0: continue
        T_bound = T_[I_abs > (I_abs.max() / 3.0)]
        amp = np.min(I_)
        if I_abs.max() == I_.max():  # positive curent inj
            amp = np.max(I_)
        amps.append(float(amp))
        clamp.dur = T_bound[-1] - T_bound[0]
        clamp.delay = START_TIME + T_bound[0]
        clamp.amp = amp/1000.0
        # print(amp, clamp.dur, clamp.amp)
        h.v_init = E_PAS
        # h.cvode_active(1)
        h.run()

        T = np.array(h_T)
        V = np.array(soma_V)[T >= START_TIME]
        T = T[T >= START_TIME]
        T -= T[0]
        ax_arr[1].plot(T_, V_-np.median(V_[:500])+V[0], color='k', zorder=1)
        ax_arr[1].plot(T, V, color='r', ls='--', zorder=2)
        # break
    print("done_amps:", amps)
    json.dump(dict(
                amps=amps,
                area = dict(
                    all = sum([seg.area() for sec in cell.all for seg in sec]),
                    basal = sum([seg.area() for seg in basal]),
                    oblique = sum([seg.area() for seg in obliqes]),
                    tuft = sum([seg.area() for seg in tuft]),
                    trunk = sum([seg.area() for seg in trunk]),
                )
            ), open(os.path.join(save_folder, person+'_'+cell_name+'_extra_data.json'), 'w'), indent=4, sort_keys=True)

    ax_arr[1].set_xlim(0, max_x)
    ax_arr[1].plot([max_x-550, max_x-50], [E_PAS-10, E_PAS-10], color='k')
    ax_arr[1].plot([max_x-550, max_x-550], [E_PAS-5, E_PAS-10], color='k')
    ax_arr[1].set_xticks([])
    ax_arr[1].set_yticks([])
    ax_arr[1].set_axis_off()

    ax_arr[3] = run_attenuation_ploter(cell=cell,
                           seg_start=list(cell.soma[0])[0],
                           color_func=colors,
                           seg_length_function = get_segment_length_lamda,
                           more_conductances_=more_conductances_,
                           param_to_record='v',
                           record_to_value_func=None,
                           norm=True,
                           delay=2000.0,
                           dur=1000.0,
                           amp=0.1,
                           seg_to_indicate={seg:dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot},
                           ax = ax_arr[3])

    ax_arr[3].spines['top'].set_visible(False)
    ax_arr[3].spines['right'].set_visible(False)
    bbp_loader.model.destroy(bbp_loader.sim)
    return max(np.where(results['sons']['all']['d3_2'].flatten()>0)[0][-1], np.where(results['parent']['all']['d3_2'].flatten()>0)[0][-1])/factor_e_space

if __name__ == "__main__":

    base_folder = 'optim_ZAP_power_new6_Ih_linear/'
    parameter_folder = 'sag_models/'
    load_num=0

    try: os.mkdir('graphs_for_fig1')
    except: pass

    save_folder = 'graphs_for_fig1/'

    fig, ax_arr = plt.subplots(2, 4, figsize=(5*4, 5*2), gridspec_kw={'width_ratios': [3, 1, 3, 2]}, sharex='col', sharey='col')

    number = 1
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]
    max_cable1 = run_cell(person , cell_name, fig=fig, ax_arr=ax_arr[0], load_num=0)

    number = 5
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]
    max_cable2 = run_cell(person , cell_name, fig=fig, ax_arr=ax_arr[1], load_num=0)

    ax_arr[0, 3].set_xlim(-0.1, max(max_cable1, max_cable2)+0.1)
    plt.savefig(os.path.join(save_folder, 'temp2'))
    plt.savefig(os.path.join(save_folder, 'temp2.pdf'))
    plt.close()
