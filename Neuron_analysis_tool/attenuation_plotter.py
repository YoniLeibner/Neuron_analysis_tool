from neuron import h
import matplotlib.pyplot as plt
import numpy as np
from Neuron_analysis_tool.more_conductances import more_conductances
from Neuron_analysis_tool.utils import *
# from dendogram_class import more_conductances, get_segment_length_lamda, get_segment_length_um
# import os
# from neuron import gui
# from scipy.signal import find_peaks
# from glob import glob
# from tqdm import tqdm
# import efel
# import pickle
# from scipy.stats import sem
# from scipy import stats
# import pandas as pd
# import seaborn as sns
# import math
# from matplotlib.lines import Line2D
# import matplotlib.colors as colors
# import matplotlib.cm as cmx

class attenuation:

    def __init__(self, cell, color_func, seg_length_function = get_segment_length_lamda, more_conductances = more_conductances, param_to_record = None, record_to_value_func=None):
        self.cell=cell
        self.record_dict = dict()
        self.seg_length_function = seg_length_function
        self.color_func = color_func
        self.more_conductances = more_conductances
        if not param_to_record  is None:
            self.set_recordings(param_to_record)
        self.start_seg=None
        self.distance_dict=dict()
        if record_to_value_func is None:
            self.record_to_value_func = self.record_to_value
        else:
            self.record_to_value_func=record_to_value_func

    def set_recordings(self, param_to_record):
        self.param_to_record = param_to_record
        for sec in self.cell.all:
            self.record_dict[sec]=dict()
            for seg in sec:
                self.record_dict[sec][seg] = h.Vector()
                self.record_dict[sec][seg].record(getattr(seg, '_ref_%s' % self.param_to_record))

    def record_to_value(self, rec):
        return rec.max() - rec[0]

    def compute_distances_helper(self, sec, parent_seg, done, reverse=False):
        if sec in done:
            return done
        done.add(sec)
        segs = list(sec)
        if reverse:
            segs = segs[::-1]
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=parent_seg, done=done)
        self.distance_dict[sec] = dict()
        for seg in segs:
            self.distance_dict[sec][seg] = dict(l=self.seg_length_function(seg, self.more_conductances) + self.distance_dict[parent_seg.sec][parent_seg]['l'], parent=parent_seg)
            parent_seg = seg

        if reverse:
            if not sec.parentseg() is None:
                done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=seg, done=done, reverse=True)
        else:
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=seg, done=done)
        return done

    def compute_distances(self, start_seg):
        if start_seg == self.start_seg:
            return
        self.start_seg = start_seg
        self.distance_dict = dict()
        sec = start_seg.sec
        self.distance_dict[sec] = dict()
        h.distance(0, start_seg.x, sec=sec)
        segs = [seg for seg in sec]
        done = set()
        done.add(sec)
        parent_seg = start_seg
        self.distance_dict[sec][start_seg] = dict(l=0, parent=None)
        for seg in segs:
            if seg.x > start_seg.x:
                self.distance_dict[sec][seg] = dict(l=self.seg_length_function(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['l'], parent=parent_seg)
                parent_seg = seg
        #now we go to the sones and
        for son in sec.children():
            done = self.compute_distances_helper(son, parent_seg=seg, done=done)

        parent_seg = start_seg
        for seg in segs[::-1]:
            if seg.x < start_seg.x:
                self.distance_dict[sec][seg] =dict(l= self.seg_length_function(seg, self.more_conductances)+ self.distance_dict[sec][parent_seg]['l'], parent=parent_seg)
                parent_seg = seg
        if not sec.parentseg() is None:
            done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=seg, done=done, reverse=True)

    def plot_helper(self, sec, parent_seg, done, ax, reverse=False, cut_start=20000, norm=1.0, seg_to_indicate=dict()):
        if sec in done:
            return done
        segs = list(sec)
        done.add(sec)

        if reverse and (sec is not self.cell.soma[0]):
            segs = segs[::-1]
            for son in sec.children():
                done = self.plot_helper(son, parent_seg=parent_seg, done=done, ax=ax, cut_start=cut_start, norm=norm, seg_to_indicate=seg_to_indicate)

        x, y = [self.distance_dict[parent_seg.sec][parent_seg]['l']], [self.record_to_value_func(np.array(self.record_dict[parent_seg.sec][parent_seg])[cut_start:])/norm]
        for seg in segs:
            c, _ = self.color_func.get_seg_color(seg)
            x.append(self.distance_dict[sec][seg]['l'])
            y.append(self.record_to_value_func(np.array(self.record_dict[sec][seg])[cut_start:])/norm)
            if seg in seg_to_indicate.keys():
                ax.scatter(x[-1], y[-1], color=seg_to_indicate[seg]['color'], s=seg_to_indicate[seg]['size'], zorder=10)
            ax.plot(x[-2:], y[-2:], color=c)
        if reverse:
            if not sec.parentseg() is None:
                done = self.plot_helper(sec.parentseg().sec, parent_seg=seg, done=done, reverse=True, ax=ax, cut_start=cut_start, norm=norm, seg_to_indicate=seg_to_indicate)
        else:
            for son in sec.children():
                done = self.plot_helper(son, parent_seg=seg, done=done, ax=ax, cut_start=cut_start, norm=norm, seg_to_indicate=seg_to_indicate)
        return done


    def plot(self, start_seg=None, ax=None, cut_start=20000, norm=False, seg_to_indicate=dict()):
        if start_seg is None:
            segs = list(self.cell.soma[0])
            start_seg = segs [len(segs)//2]
        self.compute_distances(start_seg)
        if ax is None:
            ax = plt.gca()

        sec = start_seg.sec
        segs = list(sec)
        done=set()
        done.add(sec)
        norm_by=1.0
        if norm:
            norm_by = self.record_to_value_func(np.array(self.record_dict[sec][start_seg])[cut_start:])
            print('norm by :', norm_by)
        x, y = [self.distance_dict[sec][start_seg]['l']], [self.record_to_value_func(np.array(self.record_dict[sec][start_seg])[cut_start:])/norm_by]
        for seg in segs:
            if seg.x > start_seg.x:
                c, _ = self.color_func.get_seg_color(seg)
                x.append(self.distance_dict[sec][seg]['l'])
                y.append(self.record_to_value_func(np.array(self.record_dict[sec][seg])[cut_start:])/norm_by)
                if seg in seg_to_indicate.keys():
                    ax.scatter(x[-1], y[-1], color=seg_to_indicate[seg]['color'], s=seg_to_indicate[seg]['size'], zorder=10)
                ax.plot(x[-2:], y[-2:], color=c)
        #now we go to the sones and
        for son in sec.children():
            done = self.plot_helper(son, parent_seg=seg, done=done, ax=ax, cut_start=cut_start, norm=norm_by, seg_to_indicate=seg_to_indicate)

        x, y = [self.distance_dict[sec][start_seg]['l']], [self.record_to_value_func(np.array(self.record_dict[sec][start_seg])[cut_start:])/norm_by]
        for seg in segs[::-1]:
            if seg.x < start_seg.x:
                c, _ = self.color_func.get_seg_color(seg)
                x.append(self.distance_dict[sec][seg]['l'])
                y.append(self.record_to_value_func(np.array(self.record_dict[sec][seg])[cut_start:])/norm_by)
                if seg in seg_to_indicate.keys():
                    ax.scatter(x[-1], y[-1], color=seg_to_indicate[seg]['color'], s=seg_to_indicate[seg]['size'], zorder=10)
                ax.plot(x[-2:], y[-2:], color=c)
        if not sec.parentseg() is None:
            done = self.plot_helper(sec.parentseg().sec, parent_seg=seg, done=done, reverse=True, ax=ax, cut_start=cut_start, norm=norm_by, seg_to_indicate=seg_to_indicate)
        return ax




def run_attenuation(cell, seg_start, color_func, seg_length_function, more_conductances_, param_to_record='v', record_to_value_func=None, delay=2000.0, dur=1000.0, amp=0.1):
    att = attenuation(cell, color_func=color_func, seg_length_function=seg_length_function, more_conductances=more_conductances_, param_to_record=param_to_record, record_to_value_func=record_to_value_func)
    h.tstop = delay+dur+250
    clamp = h.IClamp(seg_start.x, sec = seg_start.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return att


def run_attenuation_to_segs(cell, seg_start, param_to_record='v', record_to_value_func=None, norm=False, delay=2000.0, dur=1000.0, amp=0.1, seg_to_indicate=dict()):
    att = run_attenuation(cell, seg_start, None, None, None, param_to_record=param_to_record, record_to_value_func=record_to_value_func, delay=delay, dur=dur, amp=amp)
    cut_start = int((delay - 1) / h.dt)
    norm_by=att.record_to_value_func(np.array(att.record_dict[seg_start.sec][seg_start])[cut_start:])
    res = dict()
    for seg in seg_to_indicate:
        if norm:
            res[seg] = att.record_to_value_func(np.array(att.record_dict[seg.sec][seg])[cut_start:])/norm_by
        else:
            res[seg] = att.record_to_value_func(np.array(att.record_dict[seg.sec][seg])[cut_start:])
    return res



def run_attenuation_ploter(cell, seg_start, color_func, seg_length_function, more_conductances_, param_to_record='v', record_to_value_func=None, norm=False, delay=2000.0, dur=1000.0, amp=0.1, seg_to_indicate=dict(), ax=None):
    att = run_attenuation(cell, seg_start, color_func, seg_length_function, more_conductances_, param_to_record=param_to_record, record_to_value_func=record_to_value_func, delay=delay, dur=dur, amp=amp)
    # att = attenuation(cell, color_func=color_func, seg_length_function=seg_length_function, more_conductances=more_conductances_, param_to_record=param_to_record, record_to_value_func=record_to_value_func)
    # h.tstop = delay+dur
    # clamp = h.IClamp(seg_start.x, sec = seg_start.sec)
    # clamp.delay = delay
    # clamp.dur = dur
    # clamp.amp = amp
    # h.v_init = cell.soma[0].e_pas
    # h.celsius = 37
    # h.run()
    ax = att.plot(start_seg = seg_start, norm=norm, cut_start=int((delay-1)/h.dt), seg_to_indicate=seg_to_indicate, ax=ax)
    ax.set_yscale('log')
    return ax


def add_syn(seg, g_AMPA, g_NMDA, delay=2000.0):
    netstim = h.NetStim()
    netstim.interval = 1
    netstim.start = delay
    netstim.noise = 0
    netstim.number = 1

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

    return [AMPA, AMPA_con], [NMDA, NMDA_con], netstim

def get_voltages():
    import os
    current_path = os.path.dirname(os.path.realpath(__file__))
    NMDA = np.loadtxt(os.path.join(current_path, 'data/NMDA.txt')).T
    AMPA = np.loadtxt(os.path.join(current_path, 'data/AMPA.txt')).T
    return NMDA[1], AMPA[1]

def run_attenuation_ploter_syn(cell, seg_start, color_func, seg_length_function, more_conductances_, param_to_record='v', record_to_value_func=None, norm=False, delay=2000.0, do_NMDA=True, seg_to_indicate=dict()):
    att = attenuation(cell, color_func=color_func, seg_length_function=seg_length_function, more_conductances=more_conductances_, param_to_record=param_to_record, record_to_value_func=record_to_value_func)
    NMDA, AMPA = get_voltages()

    h.tstop = delay+300# on the safe side
    if do_NMDA:
        start = np.zeros(int(delay/h.dt))+NMDA[0]
        V_vec=h.Vector(np.concatenate([start, NMDA]))
    else:
        start = np.zeros(int(delay / h.dt)) + AMPA[0]
        V_vec=h.Vector(np.concatenate([start, AMPA]))
    # AMPA, NMDA, netstim = add_syn(seg_start, g_AMPA, g_NMDA, delay=delay)
    Vclamp = h.SEClamp(seg_start.x, sec=seg_start.sec)
    Vclamp.rs = 1e-3
    Vclamp.dur1 = 1e9
    V_vec.play(Vclamp._ref_amp1, h.dt)
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    ax = att.plot(start_seg = seg_start, norm=norm, cut_start=int(delay/h.dt), seg_to_indicate=seg_to_indicate)
    plt.yscale('log')
    return ax, np.array(att.record_dict[seg_start.sec][seg_start])[int((delay-10)/h.dt):]

