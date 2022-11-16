###########################################################
#
# author: Yoni Leibner
# description: class to calc R_total currently seort only
#              R_total at rest (after some run time)
#              can be changed using the record class to
#              record all conductances but its costy in
#              memory to record all conductances, check to
#              change the dt
# date of modification: 16.11.2022
#
###########################################################

from neuron import h
import numpy as np

#todo we need to generalize to all channels and allow to add protocol run and record all channels

def get_condactance(mechanisms):
    try:
        return getattr(mechanisms, 'g'+str(mechanisms))
    except:
        return 0


class more_conductances():
    def __init__(self, cell, run_time=3000, record_names=['gIhbar_Ih_human_linear'], is_resting=True):
        self.name = 'Ih_check'
        self.cell=cell
        self.run_time = run_time
        self.is_resting=is_resting
        self.record_names = []
        if not self.is_resting:
            self.record_names = ['_ref_'+name for name in record_names]
        self.run_resting()

    def record_condactances(self):
        record_dict = dict()
        for sec in self.cell.all:
            record_dict[sec] = dict()
            for i, seg in enumerate(sec):
                record_dict[sec][seg] = dict()
                for record_name in self.record_names:
                    try:
                        record_dict[sec][seg][record_name]=h.Vector()
                        record_dict[sec][seg][record_name].record(getattr(sec(seg.x), record_name))
                    except:
                        record_dict[sec][seg]=0 # no Ih hare
        return record_dict

    def run_resting(self):
        if not self.is_resting:
            self.record_dict = self.record_condactances()
        h.tstop = self.run_time
        h.run()
        if not self.is_resting:
            for sec in self.record_dict.keys():
                for seg in self.record_dict[sec].keys():
                    for record_name in self.record_dict[sec][seg].keys():
                        if self.record_dict[sec][seg][record_name] == 0: continue
                        self.record_dict[sec][seg][record_name] = np.array(self.record_dict[sec][seg][record_name])[-1] # stady state opening
        else:
            self.record_dict = dict()
            for sec in self.cell.all:
                self.record_dict[sec] = dict()
                for i, seg in enumerate(sec):
                    self.record_dict[sec][seg] = dict()
                    for mechanisms in seg:
                        self.record_dict[sec][seg]['g'+str(mechanisms)] = get_condactance(mechanisms)

    def cumpute(self, seg):
        sec= seg.sec
        g_total = seg.g_pas + sum([self.record_dict[sec][seg][record_name] for record_name in self.record_dict[sec][seg]])
        return 1.0/g_total


class more_conductances_fake():
    def __init__(self, cell):
        self.name='fake'

    def cumpute(self, seg):
        return 1.0/seg.g_pas
