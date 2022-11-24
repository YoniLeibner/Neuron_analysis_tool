#########################################################
#
# author: Yoni Leibner
# description: record class to record from a segment and
#                           from all the neuron segments
# date of modification: 16.11.2022
#
#########################################################


import neuron
from neuron import h
import numpy as np

class record:

    def __init__(self, seg, record_name):
        assert hasattr(seg, '_ref_' + record_name), 'wrong recording name'
        self.seg = seg
        self.record_name = record_name
        self._record = h.Vector()
        self._record.record(getattr(seg, '_ref_' + record_name))


    def extract(self, extraction_func):
        self._record = extraction_func(self._record)

    def get_val(self, func=lambda x: np.mean(x)):
        if type(self._record) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return func(self._record)

    def set_record(self, record_to_set):
        self._record = record_to_set


class record_all:
    def __init__(self, cell, record_name='v'):
        self.cell=cell
        self.record_name = record_name
        self.restart()

    def restart(self):
        self.record_dict = dict()
        for sec in self.cell.all:
            self.record_dict[sec] = dict()
            for seg in sec:
                if hasattr(seg, '_ref_' + self.record_name):
                    self.record_dict[sec][seg] = record(seg, self.record_name)
                else:
                    self.record_dict[sec][seg] = 'non_exsisting'
        self.time = h.Vector()
        self.time.record(h._ref_t)

    def push_records(self, record_dict, time):
        for sec in self.cell.all:
            for seg in sec:
                self.record_dict[sec][seg].set_record(record_dict[sec][seg])
        self.time = time

    def push_records_seg(self, record_dict, time):
        for sec in self.cell.all:
            for seg in sec:
                self.record_dict[sec][seg].set_record(record_dict[seg])
        self.time = time

    def extract(self, extraction_func):
        for sec in self.record_dict:
            for seg in self.record_dict[sec]:
                if not self.record_dict[sec][seg] == 'non_exsisting':
                    self.record_dict[sec][seg].extract(extraction_func)
        self.time = extraction_func(self.time)
        self.time -= self.time[0]

    def get_vals(self, func=lambda x: np.mean(x), default_res=0):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        res=dict()
        for sec in self.record_dict:
            res[sec] = dict()
            for seg in self.record_dict[sec]:
                if not self.record_dict[sec][seg] == 'non_exsisting':
                    res[sec][seg] = self.record_dict[sec][seg].get_val(func)
                else:
                    res[sec][seg] = default_res
        return res

    def get_vals_at_t(self, t, default_res=0):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        res=dict()
        indexs = np.where(self.time>=t)[0]
        assert len(indexs)>0, 'the time bin dont exsists, make sure you got the correct time between 0 and '+str(self.time[-1])

        func = lambda x: x[indexs[0]]
        return self.get_vals(func=func, default_res=default_res)

    def get_vals_at_dt(self, t1, t2, default_res=0, dt_func = lambda x: np.max(x)):
        assert t1<t2, 'the time bins have no are not corect, t1='+str(t1)+', t2='+str(t2)
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        indexs1 = np.where(self.time >= t1)[0]
        indexs2 = np.where(self.time >= t2)[0]
        assert len(indexs1) > 0, 'the time bin (t1'+str(t1)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        assert len(indexs2) > 0, 'the time bin (t2'+str(t2)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        func = lambda x: dt_func(x[indexs1[0]:indexs2[0]])
        return self.get_vals(func=func, default_res=default_res)

    def get_max(self):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return np.nanmax([np.nan if self.record_dict[sec][seg]=='non_exsisting' else self.record_dict[sec][seg].get_val(lambda x:np.max(x)) for sec in self.record_dict for seg in self.record_dict[sec]])

    def get_min(self):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return np.nanmin([np.nan if self.record_dict[sec][seg]=='non_exsisting' else self.record_dict[sec][seg].get_val(lambda x:np.min(x)) for sec in self.record_dict for seg in self.record_dict[sec]])

    def get_bounds(self):
        return [self.get_min(), self.get_max()]

    def get_record(self, seg):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        if self.record_dict[seg.sec][seg]=='non_exsisting':
            return np.zeros(self.time.shape)
        return self.record_dict[seg.sec][seg]._record.copy()

    def get_record_at_dt(self, seg, t1, t2, dt_func = lambda x: np.max(x)):
        assert seg.sec in self.record_dict, 'seg not valid'
        assert seg in self.record_dict[seg.sec], 'seg not valid'
        assert t1<t2, 'the time bins have no are not corect, t1='+str(t1)+', t2='+str(t2)
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        indexs1 = np.where(self.time >= t1)[0]
        indexs2 = np.where(self.time >= t2)[0]

        assert len(indexs1) > 0, 'the time bin (t1'+str(t1)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        assert len(indexs2) > 0, 'the time bin (t2'+str(t2)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        func = lambda x: dt_func(x[indexs1[0]:indexs2[0]])
        return self.record_dict[seg.sec][seg].get_val(func)

    def is_existing(self, seg):
        if seg.sec in self.record_dict:
            if seg in self.record_dict[seg.sec]:
                if not self.record_dict[seg.sec][seg] == 'non_exsisting':
                    return True
        return False


