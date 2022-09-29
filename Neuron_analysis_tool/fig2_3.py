from neuron import h, gui
import matplotlib.pyplot as plt
import pickle
import os
import numpy as np
from Load_BBP_Model import Load_BBP_Model
from open_morph_with_obliques import morph_reader, get_section_nums
from dendogram_class import color_func
from utils import colors_dict
from utils import fig_names as names
from fig1_2 import plot_morph_


a='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/'
# a='C:/Users/USER1/Documents/Hippocampus_new/'
def get_cell2(person, cell_name, seg_every_um=10, load_num=0):
    morph_opener = morph_reader(df_loc = a+'new_data/table.xlsx',
                                base_folder=a+'new_data/morphs/')
    morphology_dict = morph_opener.get_cell(person=person, cell_name=cell_name)
    morphology_dirr = morphology_dict['morph_path']
    parameter_folder = 'sag_models/'

    linear_Ih = 'linear' in base_folder
    sag_model_add = '' if not linear_Ih else '_linear'
    bbp_loader = Load_BBP_Model(morphology_dirr,
                                parameter_folder + 'mechanism'+sag_model_add+'.txt',
                                parameter_folder + 'parameters'+sag_model_add+'.txt',
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

    ax.scatter(shift, 0, s=s, color='forestgreen', label='start')

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='best')
    ax.set_ylabel('distanse')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()

freqs = np.concatenate([np.arange(0, 1, 0.01), np.arange(1, 10, 0.05), np.arange(10, 1000,1)])
def transfer_freq(seg1, seg2, freqs = freqs):
    imp = h.Impedance(seg1.x, sec=seg1.sec)
    imp.loc(seg1.x, sec=seg1.sec)
    transfer = []
    for freq in freqs:
        imp.compute(freq)
        # imp.compute(freq, 1)
        transfer.append(imp.transfer(seg2.x, sec=seg2.sec))
    return freqs, np.array(transfer)

def get_branch_init(seg, seg_list, prev):
    if seg.sec not in seg_list:
        return seg, prev
    return get_branch_init(seg.sec.parentseg(), seg_list, prev=seg.sec)

parameter_folder = 'sag_models/'
base_folder = 'optim_ZAP_power_new6_Ih_linear/'


def get_trunk_init(cell, trunk_sec):
    res = []
    for sec in cell.soma[0].children():
        if sec in trunk_sec:
            res.append(sec)
    return res

def sons_in_list(sec, sec_list):
    sons = []
    for son in sec.children():
        if son in sec_list:
            sons.append(son)
    return sons

def split_trunk_bif(trunk_init, trunk_sec, sec_cross, crossed = False, befor=[], after = []):
    if not crossed:
        befor.append(trunk_init)
    else:
        after.append(trunk_init)
    trunk_sons = sons_in_list(trunk_init, trunk_sec)
    if trunk_init == sec_cross:
        for son in trunk_sons:
            befor, after=split_trunk_bif(son, trunk_sec, sec_cross, crossed=True, befor=befor, after=after)
    else:
        for son in trunk_sons:
            befor, after=split_trunk_bif(son, trunk_sec, sec_cross, crossed=crossed, befor=befor, after=after)
    return befor, after

def all_sons_in_list(sec, sec_list, res=[]):
    if sec in sec_list:
        res.append(sec)
    for son in sec.children():
        res= all_sons_in_list(son, sec_list, res)
    return res

def sec_list_to_segs(sec_list):
    res=[]
    for sec in sec_list:
        res.extend(list(sec))
    return res


def add_syn(seg, g_AMPA, g_NMDA, netstim):
    # AMPA part
    AMPA = h.Exp2Syn(seg.x, sec=seg.sec)
    AMPA_con= h.NetCon(netstim, AMPA)
    AMPA.e = 0
    AMPA.tau1 = 0.3
    AMPA.tau2 = 1.5
    AMPA_con.weight[0] = g_AMPA
    AMPA_con.delay = 0

    # NMDA part
    NMDA=h.NMDA(seg.x, sec=seg.sec)
    NMDA_con = h.NetCon(netstim, NMDA)
    NMDA.e = 0
    NMDA.tau_r_NMDA = 8
    NMDA.tau_d_NMDA = 35
    NMDA.n_NMDA = 0.27
    NMDA.gama_NMDA = 0.076
    NMDA_con.weight[0] = g_NMDA
    NMDA_con.delay = 0
    return [AMPA, AMPA_con], [NMDA, NMDA_con]

def run_cell(number, ax_arr, fig=None, norm = False, syn = True):

    cell_name = names[number][1]
    person = names[number][0]
    load_num = 0
    obliqes_nums = names[number][2]
    trunk_initial_num = names[number][3]
    more_elev = names[number][4]
    more_azim = names[number][5]

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)

    cell = bbp_loader.get_model()

    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]

    parts_dict = dict(tuft=tuft, trunk=trunk, oblique=obliqes, basal=basal, axon=axon,
                      soma=[seg for seg in cell.soma[0]])
    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    basal_sec = np.unique([seg.sec for seg in basal])
    tuft_sec = np.unique([seg.sec for seg in tuft])
    trunk_sec = np.unique([seg.sec for seg in trunk])

    obliqes_terminals = []
    for sec in obliqes_sec:
        if len(sec.children()) == 0:
            obliqes_terminals.append(sec)

    obliques_start = []
    for sec in obliqes_sec:
        if sec.parentseg().sec not in obliqes_sec:
            obliques_start.append(sec)


    seg_to_indicate_ = {seg: dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot}
    seg_to_indicate_[list(cell.soma[0])[0]] = dict(size=30, color=colors_dict['soma'], alpha=1)

    trunk_init = [cell.apic[trunk_initial_num]]
    assert len(trunk_init) == 1
    befor, after = split_trunk_bif(cell.apic[0], trunk_sec, trunk_init[0])
    befor_all, after_all = split_trunk_bif(cell.apic[0], cell.apic, trunk_init[0])
    obliques_to_change = []
    for sec in obliques_start:
        if sec.parentseg().sec in after:
            obliques_to_change.append(sec)
    connect_to = trunk_init[0](1)
    bif_sons = trunk_init[0].children()
    assert len(bif_sons) == 2
    first_side_all = all_sons_in_list(bif_sons[0], cell.apic, res=[])
    first_side = all_sons_in_list(bif_sons[0], obliqes_terminals, res=[])
    second_side_all = all_sons_in_list(bif_sons[1], cell.apic, res=[])
    second_side = all_sons_in_list(bif_sons[1], obliqes_terminals, res=[])
    first_side_sec = all_sons_in_list(bif_sons[0], obliqes_sec, res=[])
    second_side_sec = all_sons_in_list(bif_sons[1], obliqes_sec, res=[])
    oblique_before = []
    for sec in obliqes_sec:
        if sec not in first_side_sec and sec not in second_side_sec:
            oblique_before.append(sec)

    first_side_seg = sec_list_to_segs(first_side_sec)
    second_side_seg = sec_list_to_segs(second_side_sec)
    oblique_before_seg = sec_list_to_segs(oblique_before)
    hot_spot_dict = dict(first_side=[], second_side=[])
    for seg in hot_spot:
        if seg.sec in first_side_all:
            hot_spot_dict['first_side'].append(seg)
        elif seg.sec in second_side_all:
            hot_spot_dict['second_side'].append(seg)
    print(hot_spot)
    hot_spot = [hot_spot_dict['first_side'][0], hot_spot_dict['second_side'][0]]
    print(hot_spot)
    assert len(hot_spot) == 2, 'the figure cn take only 2 hot-spots'
    parts_dict = dict(soma=trunk+tuft+basal+axon+[seg for seg in cell.soma[0]],
                      oblique_before=oblique_before_seg,
                      oblique_right=first_side_seg,
                      oblique_left=second_side_seg, )

    colors_dict['oblique_before'] = 'b'
    colors_dict['oblique_right'] = 'r'
    colors_dict['oblique_left'] = '#38E31E'
    colors_dict['oblique_moved'] = '#BD1EE3'

    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)
    fig_, ax_arr[0] = plot_morph_(cell, colors, seg_to_indicate=hot_spot,
                          colors_dict=colors_dict, num=0, more_elev=more_elev,
                          more_azim=more_azim, ax=ax_arr[0], diam_factor=None)

    sec_to_change = dict()
    for sec in obliques_start:
        if sec in after_all:
            from_seg = sec.parentseg()
            sec_to_change[sec] = dict(from_=from_seg, to=connect_to)
            for son in all_sons_in_list(sec, obliqes_sec, res=[]):
                sec_to_change[son] = dict(from_=from_seg, to=connect_to)

    parts_dict = dict(soma=tuft+trunk+basal+axon+[seg for seg in cell.soma[0]],
                      oblique_before=oblique_before_seg,
                      oblique_moved=first_side_seg+second_side_seg,)

    colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)
    fig_, ax_arr[1] = plot_morph_(cell, colors, seg_to_indicate=hot_spot,
                                 colors_dict=colors_dict, num=0, more_elev=more_elev,
                                 more_azim=more_azim, ax=ax_arr[1], diam_factor=None, sec_to_change =sec_to_change )

    # return
    obliqes_terminals_segs = [list(sec)[-1] for sec in obliqes_terminals]

    ###### NMDA spikes ######

    record_dict = dict()
    color_dict_voltage = dict()
    parts_dict = dict()
    START_TIME = 3000
    h.tstop = START_TIME + 1500

    # START_TIME = 10
    # h.tstop = START_TIME + 10

    for sec in cell.all:
        for seg in sec:
            key = sec.name() + '_' + str(seg.x)
            parts_dict[key] = [seg]
            record_dict[key] = h.Vector()
            record_dict[key].record(seg._ref_v)
    if syn:
        netstim = h.NetStim()
        netstim.interval = 1
        netstim.start = START_TIME
        netstim.noise = 0
        netstim.number = 1
        syns = []
        for seg in obliqes_terminals_segs:
            number_of_syns = 15
            syns.append(add_syn(seg, g_AMPA=0.0004 * number_of_syns, g_NMDA=0.0012 * number_of_syns, netstim=netstim))

    else:
        print('simulated_NMDA')

        NMDA = np.loadtxt(
            '/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/graphs_for_fig2/voltages/NMDA.txt').T
        NMDA_T = NMDA[0]
        NMDA_V = NMDA[1]
        dt = NMDA_T[1] - NMDA_T[0]
        start = np.zeros(int(START_TIME / dt)) + NMDA[0]
        V_vec = h.Vector(np.concatenate([start, NMDA]))
        Vclamps = []
        for seg in obliqes_terminals_segs:
            Vclamps.append(h.SEClamp(seg.x, sec=seg.sec))
            Vclamps[-1].rs = 1e-3
            Vclamps[-1].dur1 = 1e9
            V_vec.play(Vclamps[-1]._ref_amp1, dt)
    h.run()

    for key in record_dict:
        V = np.array(record_dict[key])[int((START_TIME - 10) / h.dt):]
        color_dict_voltage[key] = V.max() - V[0]
    colors = color_func(parts_dict=parts_dict, color_dict=color_dict_voltage)
    fig_, ax_arr[-1] = plot_morph_(cell, colors, seg_to_indicate=hot_spot,
                                   colors_dict=colors_dict, num=0, more_elev=more_elev,
                                   more_azim=more_azim, fig=fig, ax=ax_arr[-1], diam_factor=None, bounds=[0, 75],
                                   norm_colors=True, cmap=plt.cm.turbo, plot_color_bar=fig is not None)
    if syn:
        syns = None
    soma_seg = list(cell.soma[0])[0]
    initial_res = dict()
    for seg in hot_spot + [soma_seg]:
        for seg2 in obliqes_terminals_segs:

            f, R_tr = transfer_freq(seg, seg2)
            initial_res[seg.sec.name() + '_' + seg2.sec.name()] = dict(f=f, R_tr=R_tr)

    # trunk_init = get_trunk_init(cell, trunk_sec)
    trunk_init = [cell.apic[trunk_initial_num]]
    assert len(trunk_init) == 1
    befor, after = split_trunk_bif(cell.apic[0], trunk_sec, trunk_init[0])
    obliques_to_change = []
    for sec in obliques_start:
        if sec.parentseg().sec in after:
            obliques_to_change.append(sec)

    connect_to = trunk_init[0](1)

    print(connect_to.sec.children())
    for sec in obliques_to_change:
        h.disconnect(sec=sec)
        sec.connect(connect_to, 1.0, 0.0)
    print(connect_to.sec.children())

    for i, seg in enumerate(hot_spot+[soma_seg]):
        for seg2 in obliqes_terminals_segs:
            name = seg.sec.name() + '_' + seg2.sec.name()

            r_tr = initial_res[name]['R_tr']
            if norm:
                r_tr /= r_tr[0]
            color = colors_dict['oblique_before']
            if seg2.sec in first_side:
                color = colors_dict['oblique_right']
            elif seg2.sec in second_side:
                color = colors_dict['oblique_left']
            initial_res[name]['color'] = color
            ax_arr[2+i].plot(initial_res[name]['f'], r_tr, label='initial', color=color, alpha=0.5)

            if (seg2 in first_side_seg) or (seg2 in second_side_seg):
                f, R_tr = transfer_freq(seg, seg2)
                initial_res[name]['R_tr_changed'] = R_tr
                initial_res[name]['color_changed'] = colors_dict['oblique_moved']
                if norm:
                    R_tr /= R_tr[0]
                ax_arr[2 + i].plot(f, R_tr, label='changed', color=colors_dict['oblique_moved'], alpha=0.5)
            else:
                initial_res[name]['R_tr_changed'] = None #np.zeros(r_tr.shape)
        ax_arr[2 + i].set_xscale('log')
        ax_arr[2 + i].spines['top'].set_visible(False)
        ax_arr[2 + i].spines['right'].set_visible(False)
        ax_arr[2 + i].get_xaxis().tick_bottom()
        ax_arr[2 + i].get_yaxis().tick_left()
    return initial_res



def run_cell_NMDA(number, ax, syn = True):

    cell_name = names[number][1]
    person = names[number][0]
    load_num = 0
    obliqes_nums = names[number][2]
    trunk_initial_num = names[number][3]
    more_elev = names[number][4]
    more_azim = names[number][5]

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)

    cell = bbp_loader.get_model()

    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]

    parts_dict = dict(tuft=tuft, trunk=trunk, oblique=obliqes, basal=basal, axon=axon,
                      soma=[seg for seg in cell.soma[0]])

    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    basal_sec = np.unique([seg.sec for seg in basal])
    tuft_sec = np.unique([seg.sec for seg in tuft])
    trunk_sec = np.unique([seg.sec for seg in trunk])

    obliqes_terminals = []
    for sec in obliqes_sec:
        if len(sec.children()) == 0:
            obliqes_terminals.append(sec)

    obliques_start = []
    for sec in obliqes_sec:
        if sec.parentseg().sec not in obliqes_sec:
            obliques_start.append(sec)


    seg_to_indicate_ = {seg: dict(size=30, color=colors_dict['synapse'], alpha=1) for seg in hot_spot}
    seg_to_indicate_[list(cell.soma[0])[0]] = dict(size=30, color=colors_dict['soma'], alpha=1)

    trunk_init = [cell.apic[trunk_initial_num]]
    assert len(trunk_init) == 1
    befor, after = split_trunk_bif(cell.apic[0], trunk_sec, trunk_init[0])
    befor_all, after_all = split_trunk_bif(cell.apic[0], cell.apic, trunk_init[0])


    colors_dict['oblique_before'] = 'b'
    colors_dict['oblique_right'] = 'r'
    colors_dict['oblique_left'] = '#38E31E'
    colors_dict['oblique_moved'] = '#BD1EE3'

    # return
    obliqes_terminals_segs = [list(sec)[-1] for sec in obliqes_terminals]

    ###### NMDA spikes ######

    record_dict = dict()
    color_dict_voltage = dict()
    parts_dict = dict()
    START_TIME = 3000
    h.tstop = START_TIME + 1500

    # START_TIME = 10
    # h.tstop = START_TIME + 10

    for sec in cell.all:
        for seg in sec:
            key = sec.name() + '_' + str(seg.x)
            parts_dict[key] = [seg]
            record_dict[key] = h.Vector()
            record_dict[key].record(seg._ref_v)
    if syn:
        netstim = h.NetStim()
        netstim.interval = 1
        netstim.start = START_TIME
        netstim.noise = 0
        netstim.number = 1
        syns = []
        for seg in obliqes_terminals_segs:
            number_of_syns = 15
            syns.append(add_syn(seg, g_AMPA=0.0004 * number_of_syns, g_NMDA=0.0012 * number_of_syns, netstim=netstim))

    else:
        print('simulated_NMDA')

        NMDA = np.loadtxt(
            '/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/graphs_for_fig2/voltages/NMDA.txt').T
        NMDA_T = NMDA[0]
        NMDA_V = NMDA[1]
        dt = NMDA_T[1] - NMDA_T[0]
        start = np.zeros(int(START_TIME / dt)) + NMDA[0]
        V_vec = h.Vector(np.concatenate([start, NMDA]))
        Vclamps = []
        for seg in obliqes_terminals_segs:
            Vclamps.append(h.SEClamp(seg.x, sec=seg.sec))
            Vclamps[-1].rs = 1e-3
            Vclamps[-1].dur1 = 1e9
            V_vec.play(Vclamps[-1]._ref_amp1, dt)
    h.run()

    for key in record_dict:
        V = np.array(record_dict[key])[int((START_TIME - 10) / h.dt):]
        color_dict_voltage[key] = V.max() - V[0]
    colors = color_func(parts_dict=parts_dict, color_dict=color_dict_voltage)
    fig_, ax = plot_morph_(cell, colors, seg_to_indicate=hot_spot,
                                   colors_dict=colors_dict, num=0, more_elev=more_elev,
                                   more_azim=more_azim, fig=fig, ax=ax, diam_factor=None, bounds=[0, 75],
                                   norm_colors=True, cmap=plt.cm.turbo, plot_color_bar=fig is not None)
    if syn:
        syns = None


if __name__ == "__main__":
    row_font_pad = 50
    fig = plt.figure(figsize=(6*4, 12))
    ax = fig.add_gridspec(6, 4)
    ax1 = fig.add_subplot(ax[0:3, 0]) #, projection='3d'
    ax2 = fig.add_subplot(ax[0:3, 1], sharex=ax1)
    ax3 = fig.add_subplot(ax[0, 2])
    ax4 = fig.add_subplot(ax[1, 2], sharex=ax3)
    ax5 = fig.add_subplot(ax[2, 2], sharex=ax3)
    # ax2 = fig.add_subplot(ax[0:3, 3], projection='3d')

    ax11 = fig.add_subplot(ax[3:6, 0], sharex=ax1)
    ax12 = fig.add_subplot(ax[3:6, 1], sharex=ax1)
    ax13 = fig.add_subplot(ax[3:6, 2], sharex=ax1)

    ax11.set_title(' ', fontsize=row_font_pad )
    ax12.set_title(' ', fontsize=row_font_pad )
    ax13.set_title(' ', fontsize=row_font_pad )
    try:
        os.mkdir('gra phs_for_fig3')
    except:
        pass
    save_folder = 'graphs_for_fig3/'

    try:
        os.mkdir(save_folder)
    except: pass
    Z_rt_dict = run_cell(1, ax_arr = [ax1, ax2, ax3, ax4, ax5, ax11], norm = False, fig=fig)
    pickle.dump(Z_rt_dict , open(save_folder+'Z_rt_dict.p', 'wb'))
    run_cell_NMDA(5, ax = ax12, syn=True)
    run_cell_NMDA(4, ax = ax13, syn=True)

    plt.savefig(os.path.join(save_folder, 'temp4'))
    plt.savefig(os.path.join(save_folder, 'temp4.pdf'))