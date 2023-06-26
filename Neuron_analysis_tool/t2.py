import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new')
from neuron import h, gui
import os
import numpy as np
# from utils import fig_names as names
from utiles_func import get_cell2
from Neuron_analysis_tool.Analyzer import Analyzer
import pandas as pd
import matplotlib as mpl
from utils import name_to_cluster_dict
from Neuron_analysis_tool.distance import Distance


names = [
    ['P1', 'Hipp1', [30,61,77], 49, 0, 225, 'cluster2'], #0
    ['P1', 'Hipp2', [23, 61, 90], 53, 0, 10, 'cluster2'], #1
    ['P1', 'Hipp3', [25, 50, 75], 38, 0, 200, 'cluster2'], #2
    ['P1', 'Hipp4', [101, 90, 16], 0, 0, 110, 'cluster1'], #3
    ['P1', 'Hipp8', [17,58,79], 36, 0, 60, 'cluster2'], #4
    ['P1', 'Hipp13', [18, 42, 104], 27, 0, -75, 'cluster2'], #5
    ['P6', 'Hipp1', [104, 79, 36], 5, 0, 155, 'cluster2'], #6
    ['P6', 'Hipp3', [118, 39, 107], 1, 0, 180, 'cluster1'], #7
    ['P6', 'Hipp5', [66, 28, 60], 13, 60, 180, 'cluster1'], #8 #change elev to 150
    ['P6', 'Hipp7', [115, 43, 86], 27, 0, 0, 'cluster1'], #9

    ['P1', 'Hipp7', [115, 43, 86], 36, 0, -100, 'cluster1'], #10
    ['P1', 'Hipp9', [115, 43, 86], 10, 0, -60, 'cluster1'], #11
    ['P1', 'Hipp10', [115, 43, 86], 40, 0, -90, 'cluster2'], #12
    ['P1', 'Hipp11', [115, 43, 86], 4, 0, -90, 'cluster1'], #13
    ['P1', 'Hipp12', [115, 43, 86], 0, 0, 180, 'cluster1'], #14 # there is a morphology error here
    ['P6', 'Hipp8', [115, 43, 86], 19, 0, 20, 'cluster1'], #15
    ['P6', 'Hipp6', [115, 43, 86], 30, 0, 20, 'cluster3'], #16

    ['P9', 'Hipp1', [115, 43, 86], 13, 0, 20, 'cluster3'], #17
    ['P9', 'Hipp2', [115, 43, 86], 30, 0, 20, 'cluster3'], #18 # pass no tuft
    ['P9', 'Hipp3', [115, 43, 86], 30, 0, 20, 'cluster3'], #19 # pass no tuft
    ['P9', 'Hipp4', [115, 43, 86], 41, 0, 20, 'cluster3'], #20
    ['P9', 'Hipp5', [115, 43, 86], 21, 0, 20, 'cluster3'], #21
]

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

a='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Hippocampus_new/'

def get_segment_lamda(seg, more_conductances, time=None, dt=1):
    """
    return the segment  e_length
    :param seg_len:
    :param RM:
    :param RA:
    :return:
    """
    sec = seg.sec
    seg_len = sec.L/sec.nseg
    d = seg.diam
    R_total = more_conductances.cumpute(seg, time=time, dt=dt)
    lamda_0 = np.sqrt((R_total / sec.Ra) * (d / 10000.0) / 4.0)
    # return lamda_0
    # for 500 Hz
    tau = seg.cm * R_total
    f=1
    lamda_500 = lamda_0/(np.real(np.sqrt(1+2j*np.pi*f*tau)))
    return lamda_500

def get_segment_length(seg):
    """
    return the segment  e_length
    :param seg_len:
    :param RM:
    :param RA:
    :return:
    """
    sec = seg.sec
    return sec.L/sec.nseg

def resting_protocol2(cell, start_seg=None):
    h.tstop = 1500.0
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}


def run_cell(number, load_num=0, remove_obliques=False, remove_Ih=False, save_folder=None):
    cell_name = names[number][1]
    person = names[number][0]
    more_elev = names[number][4]
    more_azim = names[number][5]

    bbp_loader, obliqes, trunk, tuft, hot_spot = get_cell2(person, cell_name, load_num=load_num)
    cell = bbp_loader.get_model()
    # print(number, 'TDL=', sum([sec.L for sec in cell.all]), 'mean_diam=',np.mean([sec.diam for sec in cell.all]))
    # return
    basal = [seg for sec in cell.basal for seg in sec]
    axon = [seg for sec in cell.axon for seg in sec]
    obliqes_sec = np.unique([seg.sec for seg in obliqes])
    tuft_sec = np.unique([seg.sec for seg in tuft])
    basal_terminals = [sec for sec in cell.basal if len(sec.children())==0]
    tuft_terminals = [sec for sec in tuft_sec if len(sec.children())==0]
    obliqe_terminals = [sec for sec in obliqes_sec if len(sec.children())==0]
    if remove_obliques:
        obliques_start = []
        for sec in obliqes_sec:
            if sec.parentseg().sec not in obliqes_sec:
                obliques_start.append(sec)
        for sec in obliques_start:
            h.disconnect(sec=sec)
        ignore_sec = obliqes_sec
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

    analyser = Analyzer(cell=cell, more_conductances_protocol=resting_protocol2) #, parts_dict=parts_dict, colors_dict=colors_dict

    soma_seg = list(analyser.cell.soma[0])
    soma_seg = soma_seg[len(soma_seg) // 2]

    distance = Distance(cell=analyser.cell, more_conductances=analyser.more_conductances)
    # from Neuron_analysis_tool.more_conductances import more_conductances_fake
    # distance = Distance(cell=analyser.cell, more_conductances=more_conductances_fake(analyser.cell))
    distance.compute(start_seg=soma_seg)

    df = pd.DataFrame()
    # for part_name, part_seg_list in parts_dict.items():
    for terminals_sec_list , part_name in zip([basal_terminals, obliqe_terminals, tuft_terminals], ['basal_terminals', 'obliqe_terminals', 'tuft_terminals']):
        for sec in terminals_sec_list:
            df = df.append(dict(
                cell_num=number,
                cell_name=cell_name,
                person=person,
                cluster=name_to_cluster_dict[person+'_'+cell_name],
                part_name=part_name,
                electrical_length_lamda=distance.get_start_end(list(sec)[-1], electrical=True)['end'],
                length_um=distance.get_start_end(list(sec)[-1], electrical=False)['end'],
            ), ignore_index=True)
    bbp_loader.model.destroy(bbp_loader.sim)
    return df


if __name__ == "__main__":
    try:
        os.mkdir('sup_fig1/')
    except:
        pass
    save_folder = 'sup_fig1/'

    dfs = []
    for cell_number in range(len(names)):
    # for cell_number in [4, 5,6,17, 20, 21]:
        if cell_number < 20: continue
        if cell_number in [12, 13, 14, 18, 19]: continue
        print('cell_number:',cell_number)
        dfs.append(run_cell(cell_number, load_num=0, remove_obliques=False, remove_Ih=False, save_folder=save_folder))
        df = pd.concat(dfs)
        df.to_csv(save_folder+'tips_length_3.csv')
        # a=1
    # pickle.dump(dict(
    #     short_pulse=dict(
    #         control=all_tip_attenuation4,
    #         no_obliques=all_tip_attenuation5,
    #         no_ih=all_tip_attenuation6,
    #     ),
    #     long_pulse=dict(
    #         control=all_tip_attenuation1,
    #         no_obliques=all_tip_attenuation2,
    #         no_ih=all_tip_attenuation3,
    #     ),
    # ), open(save_folder+'tip_data2.p', 'wb'))
    #