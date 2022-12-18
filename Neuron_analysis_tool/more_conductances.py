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
from Neuron_analysis_tool.record import record_all
from Neuron_analysis_tool.protocols import resting_protocol

#todo we need to generalize to all channels and allow to add protocol run and record all channels

def get_condactance(mechanisms):
    try:
        return getattr(mechanisms, 'g'+str(mechanisms))
    except:
        return 0

class more_conductances():
    def __init__(self, cell, is_resting=True, extraction_func=None, protocol = resting_protocol):
        self.name = 'Ih_check'
        self.cell=cell
        self.is_resting=is_resting
        if extraction_func is None:
            self.extraction_func = lambda x: np.array(x)
        else:
            self.extraction_func=extraction_func
        self.protocol = protocol
        self.record_names=[]
        for sec in self.cell.all:
            for i, seg in enumerate(sec):
                for mechanisms in seg:
                    if hasattr(seg, 'g'+str(mechanisms)):
                        self.record_names.append('g' + str(mechanisms))# = get_condactance(mechanisms)
                    if hasattr(seg, 'g_'+str(mechanisms)):
                        self.record_names.append('g_' + str(mechanisms))# = get_condactance(mechanisms)

                    if hasattr(seg, 'g'+str(mechanisms)+'_'+str(mechanisms)):
                        self.record_names.append('g' + str(mechanisms)+'_'+str(mechanisms))# = get_condactance(mechanisms)
                    if hasattr(seg, 'g_'+str(mechanisms)+'_'+str(mechanisms)):
                        self.record_names.append('g_' + str(mechanisms)+'_'+str(mechanisms))# = get_condactance(mechanisms)
        self.record_names = list(set(self.record_names))
        self.set()
        if protocol is not None:
            delay, _ = self.protocol(self.cell, None)
            if extraction_func is None:
                self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]
            self.extract()

    def set(self):
        self.record_dict = dict()
        if not self.is_resting:
            for record_name in self.record_names:
                self.record_dict[record_name] = record_all(self.cell, record_name=record_name)

        # delay, _ = self.protocol(self.cell, None)
        # if self.extraction_func is None:
        #     self.extraction_func = lambda x: np.array(x)[int(delay/h.dt):]
        # self.time = np.arange(0, h.tstop, h.dt)
        # if not self.is_resting:
        #     for record_name in self.record_dict.keys():
        #         self.record_dict[record_name].extract(self.extraction_func)
        #         self.time = self.record_dict[record_name].time
        # else:
        #     self.record_dict = dict()
        #     for sec in self.cell.all:
        #         self.record_dict[sec] = dict()
        #         for i, seg in enumerate(sec):
        #             self.record_dict[sec][seg] = dict()
        #             for mechanisms in seg:
        #                 self.record_dict[sec][seg]['g'+str(mechanisms)] = get_condactance(mechanisms)

    def extract(self, extraction_func=None):
        if extraction_func is None:
            extraction_func =  self.extraction_func
        self.time = np.arange(0, h.tstop, h.dt)
        if not self.is_resting:
            for record_name in self.record_dict.keys():
                self.record_dict[record_name].extract(extraction_func)
                self.time = self.record_dict[record_name].time
        else:
            self.record_dict = dict()
            for sec in self.cell.all:
                self.record_dict[sec] = dict()
                for i, seg in enumerate(sec):
                    self.record_dict[sec][seg] = dict()
                    for mechanisms in seg:
                        self.record_dict[sec][seg]['g'+str(mechanisms)] = get_condactance(mechanisms)


    def cumpute(self, seg, time=None, dt=1):
        sec= seg.sec
        if self.is_resting:
            g_total = seg.g_pas + sum([self.record_dict[sec][seg][record_name] for record_name in self.record_dict[sec][seg]])
        else:
            if time is None:
                g_total = seg.g_pas + sum(
                    [record.get_record_at_dt(seg, t1=record.time[-2], t2=record.time[-1], dt_func = lambda x: x[-1]) for
                     record in self.record_dict.values()])
            else:
                g_total = seg.g_pas + sum(
                    [record.get_record_at_dt(seg, t1=time-dt/2, t2=time+dt/2, dt_func=lambda x: x[-1]) for
                     record in self.record_dict.values()])
        return 1.0/g_total


class more_conductances_fake():
    def __init__(self, cell):
        self.name='fake'

    def cumpute(self, seg):
        return 1.0/seg.g_pas
