import neuron
from neuron import h
import numpy as np

class record:

    def __init__(self, seg, record_name):
        assert getattr(seg, '_ref_' + record_name), 'wrong recording name'
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


class record_all:
    def __init__(self, cell, record_name):
        self.cell=cell
        self.record_dict = dict()
        for sec in cell.all:
            self.record_dict[sec]=dict()
            for seg in sec:
                if getattr(seg, '_ref_' + record_name):
                    self.record_dict[sec][seg.x] = record(seg, record_name)
                else:
                    self.record_dict[sec][seg.x] = 'non_exsisting'
        self.record_name = record_name
        self.time = h.Vector()
        self.time.record(h._ref_t)

    def extract(self, extraction_func):
        for sec in self.record_dict:
            for seg in self.record_dict[sec]:
                if not self.record_dict[sec][seg.x] == 'non_exsisting':
                    self.record_dict[sec][seg.x] = self.record_dict[sec][seg.x].extract(extraction_func)
        self.time = extraction_func(self.time)
        self.time -= self.time[0]

    def get_vals(self, func=lambda x: np.mean(x), default_res=0):
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        res=dict()
        for sec in self.record_dict:
            res[sec] = dict()
            for seg in self.record_dict[sec]:
                if not self.record_dict[sec][seg.x] == 'non_exsisting':
                    res[sec][seg] = self.record_dict[sec][seg.x].get_val(func)
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
        assert t1<=t2, 'the time bins have no are not corect, t1='+str(t1)+', t2='+str(t2)
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        res = dict()
        indexs1 = np.where(self.time >= t1)[0]
        indexs2 = np.where(self.time >= t2)[0]
        assert len(indexs1) > 0, 'the time bin (t1'+str(t1)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        assert len(indexs2) > 0, 'the time bin (t2'+str(t2)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        func = lambda x: dt_func(x[indexs1[0]:indexs2[0]])
        return self.get_vals(func=func, default_res=default_res)




