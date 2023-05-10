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
from Neuron_analysis_tool.utils import sec_name, seg_name
from Neuron_analysis_tool.protocols import resting_protocol
import pickle

#todo we need to generalize to all channels and allow to add protocol run and record all channels


def get_condactance_name(mechanisms):
    try:
        valids = []
        for param in list(mechanisms):
            if param.name()[0] == 'g' and 'bar' not in param.name():
                valids.append(param.name())
        assert len(valids)==1, 'check the more_conductances, and change the get_condactance func to seport your case'
        return valids[0]
        # return getattr(mechanisms, 'g_'+str(mechanisms))
    except:
        return ''

def get_condactance(mechanisms):
    try:
        return getattr(mechanisms, 'g_'+str(mechanisms))
    except:
        return 0

class more_conductances():
    def __init__(self, cell, is_resting=True, extraction_func=None, protocol = resting_protocol):
        # assert protocol is not None, 'check your protocol in more_conductances'
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
                    if mechanisms.name() in ['ca_ion', 'k_ion', 'na_ion','CaDynamics_E2']: continue
                    mechanisms_g_name = get_condactance_name(mechanisms)
                    if mechanisms_g_name =='':
                        raise Exception('un recognized mechanisms:'+mechanisms.name()+', check get_condactance_name func in more_conductances')
                    self.record_names.append(mechanisms_g_name)
                    # if hasattr(seg, 'g'+str(mechanisms)):
                    #     self.record_names.append('g' + str(mechanisms))# = get_condactance(mechanisms)
                    # if hasattr(seg, 'g'+str(mechanisms)+'_'+str(mechanisms)):
                    #     self.record_names.append('g' + str(mechanisms)+'_'+str(mechanisms))# = get_condactance(mechanisms)

        self.record_names = list(set(self.record_names))
        self.set(extraction_func)

        # delay, _ = self.protocol(self.cell, None)
        # if extraction_func is None:
        #     self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]
        # self.extract()

    def set(self, extraction_func):
        self.record_dict = dict()
        if not self.is_resting:
            for record_name in self.record_names:
                self.record_dict[record_name] = record_all(self.cell, record_name=record_name)
        if self.protocol is not None:
            delay, _ = self.protocol(self.cell, None)
            if extraction_func is None:
                self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]

            self.extract()



    def extract(self, extraction_func=None):
        if extraction_func is None:
            extraction_func =  self.extraction_func
        # self.time = np.arange(0, h.tstop, h.dt)
        if not self.is_resting:
            for record_name in self.record_dict.keys():
                self.record_dict[record_name].extract(extraction_func)
                # self.time = self.record_dict[record_name].time
        else:
            for record_name in self.record_names:
                temp_record_dict = dict()
                for sec in self.cell.all:
                    temp_record_dict[sec] = dict()
                    for i, seg in enumerate(sec):
                        try:
                            temp_record_dict[sec][seg] = np.array([getattr(seg, record_name)])
                        except:
                            # print('record name:'+record_name+', is not in:'+str(sec.name())+' '+str(seg))
                            temp_record_dict[sec][seg] = np.array([0])
                self.record_dict[record_name] = record_all(self.cell, record_name=record_name)
                self.record_dict[record_name].push_records(temp_record_dict, np.array([0]))

            # self.record_dict = dict()
            # for sec in self.cell.all:
            #     self.record_dict[sec] = dict()
            #     for i, seg in enumerate(sec):
            #         self.record_dict[sec][seg] = dict()
            #         for mechanisms in seg:
            #             self.record_dict[sec][seg]['g'+str(mechanisms)] = get_condactance(mechanisms)

    def cumpute(self, seg, time=None, dt=1, dt_func = lambda x: np.mean(x)): #  lambda x: x[-1]
        sec= seg.sec
        if self.is_resting:
            # g_total = seg.g_pas + sum([self.record_dict[sec_name(sec)][seg_name(seg)][record_name] for record_name in self.record_dict[sec_name(sec)][seg_name(seg)]])

            g_total = seg.g_pas + sum([record.get_record(seg)[0] for record in self.record_dict.values()])
        else:
            if time is None:
                g_total = seg.g_pas + sum(
                    [record.get_record_at_dt(seg, t1=record.time[-2], t2=record.time[-1], dt_func = dt_func) for
                     record in self.record_dict.values()])
            else:
                g_total = seg.g_pas + sum(
                    [record.get_record_at_dt(seg, t1=time-dt/2, t2=time+dt/2, dt_func=dt_func) for
                     record in self.record_dict.values()])
        return 1.0/g_total


    def save(self, save_dir):
        pickle.dump(dict(
            is_resting = self.is_resting,
            record_names=self.record_names,
            record_dict=self.record_dict,
        ), open(save_dir, 'wb'))

    def load(self, save_dir):
        pickle.dump(dict(
            is_resting=self.is_resting,
            record_names=self.record_names,
            record_dict=self.record_dict,
        ), open(save_dir, 'wb'))


class more_conductances_fake():
    def __init__(self, cell):
        pass

    def cumpute(self, seg, time=None, dt=None):
        return 1.0/seg.g_pas
