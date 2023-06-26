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

from neuron import h, nrn
import numpy as np
from Neuron_analysis_tool.record import record_all
from Neuron_analysis_tool.utils import sec_name, seg_name
from Neuron_analysis_tool.protocols import resting_protocol
import pickle

#todo we need to generalize to all channels and allow to add protocol run and record all channels



class more_conductances_fake():
    def __init__(self, cell):
        pass

    def cumpute(self, seg, time=None, dt=None, dt_func = lambda x: np.mean(x)):
        return 1.0/seg.g_pas

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


def get_condactance_point_prosses_name(point_prosses):
    if hasattr(point_prosses, 'g'):
        return True, 'g'
    name = point_prosses.hname().split('[')[0]
    if hasattr(point_prosses, 'g_'+name):
        return True, 'g_'+name
    return False, 0

def get_seg_condactances(cell):
    res = dict()
    for sec in cell.all:
        for seg in sec:
            res[seg] = []
            for mechanisms in seg:
                if mechanisms.name() == 'g_total': continue
                if mechanisms.is_ion(): continue
                if mechanisms.name() in ['CaDynamics_E2']: continue
                record_name = get_condactance_name(mechanisms)
                res[seg].append([seg, record_name])
            for point_prosses in seg.point_processes():
                check, g_name = get_condactance_point_prosses_name(point_prosses)
                if check:
                    res[seg].append([point_prosses, g_name])
    return res
#
def callback(analyzer):
    if h.t==0: #start of new simulation
        set_refs(analyzer.cell)

def set_refs(cell):
    seg_conductances = get_seg_condactances(cell)

    for seg in seg_conductances :
        g_ref_count = 0
        g_syn_count = 0
        for [p, g_name] in seg_conductances[seg]:
            if type(p) == nrn.Segment:
                if g_ref_count > 19:
                    print('to many conductances to record, please update the g_total mod file.\n this conductance is not recorded:'+g_name, p.sec.name, p.x)
                else:
                    h.setpointer(getattr(p, '_ref_'+g_name), 'g_ref'+str(g_ref_count), p.g_total)
                    g_ref_count+=1
            else:
                if g_syn_count > 39:
                    print(
                    'to many synapses iin this segment, please update the g_total mod file.\n this synaptic conductance is not recorded:' + g_name, p.sec.name, p.x)
                else:
                    h.setpointer(getattr(p, '_ref_'+g_name), 'g_syn'+str(g_syn_count), seg.g_total)
                    g_syn_count+=1
        for i in range(g_ref_count, 20):
            h.setpointer(seg._ref_zero_val_g_total, 'g_ref' + str(i), seg.g_total)
        for i in range(g_syn_count, 40):
            h.setpointer(seg._ref_zero_val_g_total, 'g_syn' + str(i), seg.g_total)

class more_conductances():
    def __init__(self, cell, is_resting=True, extraction_func=None, protocol = resting_protocol):

        self.cell=cell
        self.is_resting=is_resting
        if extraction_func is None:
            self.extraction_func = lambda x: np.array(x)
        else:
            self.extraction_func=extraction_func
        self.protocol = protocol
        self.set(extraction_func)

    def set(self, extraction_func):
        self.g_total_rec = record_all(self.cell, record_name='g_total_g_total')
        self.g_syn_rec = record_all(self.cell, record_name='g_syn_g_total')
        if self.protocol is not None:
            delay, _ = self.protocol(self.cell, None)
            if extraction_func is None:
                self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]
            self.extract()

    def extract(self, extraction_func=None):
        if extraction_func is None:
            extraction_func =  self.extraction_func
        self.g_total_rec.extract(extraction_func)
        self.g_syn_rec.extract(extraction_func)

    def cumpute(self, seg, time=None, dt=1, dt_func = lambda x: np.mean(x)):
        # g_total is in S/cm^2
        seg_L = seg.sec.L/seg.sec.nseg
        if time is None:
            g_total = self.g_total_rec.get_record_at_dt(seg, t1=self.g_total_rec.time[-2],
                                                        t2=self.g_total_rec.time[-1], dt_func = dt_func)*(np.pi * seg.diam)*10**-2 # this is in uS/um first change to uS/um^2 and then divide by the Circumference
            g_total += self.g_syn_rec.get_record_at_dt(seg, t1=self.g_total_rec.time[-2],
                                                        t2=self.g_total_rec.time[-1], dt_func=dt_func)/seg_L #here its in uS and we divid by the length of the segment so uS/um
        else:
            g_total = self.g_total_rec.get_record_at_dt(seg, t1=time-dt/2, t2=time+dt/2, dt_func=dt_func)*(np.pi * seg.diam)*10**-2  # this is in uS/um first change to uS/um^2 and then divide by the Circumference
            g_total += self.g_syn_rec.get_record_at_dt(seg, t1=time-dt/2, t2=time+dt/2, dt_func=dt_func)/seg_L  #here its in uS and we divid by the length of the segment so uS/um
        return 1.0/g_total

    def save(self, save_dir):
        pickle.dump(dict(
            is_resting = self.is_resting,
            record_dict=self.g_total_rec,
        ), open(save_dir, 'wb'))

    def load(self, save_dir):
        pickle.dump(dict(
            is_resting=self.is_resting,
            record_dict=self.g_total_rec,
        ), open(save_dir, 'wb'))

#
# class more_conductances():
#     def __init__(self, cell, is_resting=True, extraction_func=None, protocol = resting_protocol):
#         # assert protocol is not None, 'check your protocol in more_conductances'
#         self.cell=cell
#         self.is_resting=is_resting
#         if extraction_func is None:
#             self.extraction_func = lambda x: np.array(x)
#         else:
#             self.extraction_func=extraction_func
#         self.protocol = protocol
#         self.record_names=[]
#         for sec in self.cell.all:
#             for i, seg in enumerate(sec):
#                 for mechanisms in seg:
#                     if mechanisms.is_ion(): continue
#                     if mechanisms.name() in ['CaDynamics_E2', 'g_total']: continue # ion pumps not listed
#                     mechanisms_g_name = get_condactance_name(mechanisms)
#                     if mechanisms_g_name == '':
#                         raise Exception('un recognized mechanisms:'+mechanisms.name()+', check get_condactance_name func in more_conductances')
#                     self.record_names.append(mechanisms_g_name)
#
#
#         self.record_names = list(set(self.record_names))
#         self.set(extraction_func)
#
#         # delay, _ = self.protocol(self.cell, None)
#         # if extraction_func is None:
#         #     self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]
#         # self.extract()
#
#     def set(self, extraction_func):
#         self.record_dict = dict()
#         if not self.is_resting:
#             for record_name in self.record_names:
#                 self.record_dict[record_name] = record_all(self.cell, record_name=record_name)
#         if self.protocol is not None:
#             delay, _ = self.protocol(self.cell, None)
#             if extraction_func is None:
#                 self.extraction_func = lambda x: np.array(x)[int(delay / h.dt):]
#
#             self.extract()
#
#
#
#     def extract(self, extraction_func=None):
#         if extraction_func is None:
#             extraction_func =  self.extraction_func
#         # self.time = np.arange(0, h.tstop, h.dt)
#         if not self.is_resting:
#             for record_name in self.record_dict.keys():
#                 self.record_dict[record_name].extract(extraction_func)
#                 # self.time = self.record_dict[record_name].time
#         else:
#             for record_name in self.record_names:
#                 temp_record_dict = dict()
#                 for sec in self.cell.all:
#                     temp_record_dict[sec] = dict()
#                     for i, seg in enumerate(sec):
#                         try:
#                             temp_record_dict[sec][seg] = np.array([getattr(seg, record_name)])
#                         except:
#                             # print('record name:'+record_name+', is not in:'+str(sec.name())+' '+str(seg))
#                             temp_record_dict[sec][seg] = np.array([0])
#                 self.record_dict[record_name] = record_all(self.cell, record_name=record_name)
#                 self.record_dict[record_name].push_records(temp_record_dict, np.array([0]))
#
#             # self.record_dict = dict()
#             # for sec in self.cell.all:
#             #     self.record_dict[sec] = dict()
#             #     for i, seg in enumerate(sec):
#             #         self.record_dict[sec][seg] = dict()
#             #         for mechanisms in seg:
#             #             self.record_dict[sec][seg]['g'+str(mechanisms)] = get_condactance(mechanisms)
#
#     def cumpute(self, seg, time=None, dt=1, dt_func = lambda x: np.mean(x)): #  lambda x: x[-1]
#         sec= seg.sec
#         if self.is_resting:
#             # g_total = seg.g_pas + sum([self.record_dict[sec_name(sec)][seg_name(seg)][record_name] for record_name in self.record_dict[sec_name(sec)][seg_name(seg)]])
#
#             g_total = seg.g_pas + sum([record.get_record(seg)[0] for record in self.record_dict.values()])
#         else:
#             if time is None:
#                 g_total = seg.g_pas + sum(
#                     [record.get_record_at_dt(seg, t1=record.time[-2], t2=record.time[-1], dt_func = dt_func) for
#                      record in self.record_dict.values()])
#             else:
#                 g_total = seg.g_pas + sum(
#                     [record.get_record_at_dt(seg, t1=time-dt/2, t2=time+dt/2, dt_func=dt_func) for
#                      record in self.record_dict.values()])
#         return 1.0/g_total
#
#
#     def save(self, save_dir):
#         pickle.dump(dict(
#             is_resting = self.is_resting,
#             record_names=self.record_names,
#             record_dict=self.record_dict,
#         ), open(save_dir, 'wb'))
#
#     def load(self, save_dir):
#         pickle.dump(dict(
#             is_resting=self.is_resting,
#             record_names=self.record_names,
#             record_dict=self.record_dict,
#         ), open(save_dir, 'wb'))
# #
