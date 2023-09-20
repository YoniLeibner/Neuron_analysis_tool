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
import pickle
import os
from .utils import seg_name, sec_name
import matplotlib.pyplot as plt

class record:
    """
    record a single value from a seg
    """

    def __init__(self, seg, record_name):
        assert hasattr(seg, '_ref_' + record_name), 'wrong recording name'
        self.seg_name = seg_name(seg)
        self.sec = sec_name(seg.sec)
        self.record_name = record_name
        self._record = h.Vector()
        self._record.record(getattr(seg, '_ref_' + record_name))


    def extract(self, extraction_func):
        """
        extract the recorded value from h.Vector using extraction_func
        :param extraction_func: function that get vector and return numpy array (can cut the start)
        :return:
        """
        self._record = extraction_func(self._record)

    def get_val(self, func=lambda x: np.mean(x)):
        """
        return the record after func
        :param func:
        :return:
        """
        if type(self._record) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return func(self._record)

    def set_record(self, record_to_set):
        """
        set the record (for loading data)
        :param record_to_set:
        :return:
        """
        self._record = np.array(record_to_set)

    def plot(self, ax, time=None, elev=0, x_shift=0, color='k', **kwargs):
        """
        plot the record into the ax
        :param ax:
        :param time:
        :param elev:
        :param x_shift:
        :param color:
        :param kwargs:
        :return:
        """
        if type(self._record) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        if time is None:
            time = np.arange(0, len(self._record), 1)*h.dt
        ax.plot(time+x_shift, self._record+elev, color=color, **kwargs)
        return ax


class record_all:
    """
    record the same value for all the segment in a givin neuron
    """
    def __init__(self, cell, record_name='v'):
        self.cell = cell
        self.record_name = record_name
        self.restart()
        self.extraction_func = lambda x: x

    def restart(self):
        """
        seting a dictinary of empty vectors befor runing a protocol
        :return:
        """
        self.record_dict = dict()
        for sec in self.cell.all:
            self.record_dict[sec_name(sec)] = dict()
            for seg in sec:
                if hasattr(seg, '_ref_' + self.record_name):
                    self.record_dict[sec_name(sec)][seg_name(seg)] = record(seg, self.record_name)
                else:
                    self.record_dict[sec_name(sec)][seg_name(seg)] = 'non_exsisting'
        self.time = h.Vector()
        self.time.record(h._ref_t)

    def push_records(self, record_dict, time):
        """
        push prerecorded valued into the records
        :param record_dict: dictinary of {sec: seg: record as np array}}
        :param time: time in the simulation (np array)
        :return:
        """
        for sec in self.cell.all:
            for seg in sec:
                if self.record_dict[sec_name(sec)][seg_name(seg)] == 'non_exsisting':
                    continue
                self.record_dict[sec_name(sec)][seg_name(seg)].set_record(record_dict[sec][seg])
        self.time = time

    def push_records_seg(self, record_dict, time):
        """
        push prerecorded valued into the records
        :param record_dict: dictinary of {seg: record as np array}}
        :param time: time in the simulation (np array)
        :return:
        """
        for sec in self.cell.all:
            for seg in sec:
                self.record_dict[sec_name(sec)][seg_name(seg)].set_record(record_dict[seg])
        self.time = time

    def extract(self, extraction_func):
        """
        extract the recorded values from h.Vector using extraction_func
        :param extraction_func: function that get vector and return numpy array (can cut the start)
        :return:
        :param extraction_func:
        :return:
        """
        if type(self.time) == np.ndarray: return
        for sec in self.record_dict:
            for seg in self.record_dict[sec]:
                if not self.record_dict[sec][seg] == 'non_exsisting':
                    self.record_dict[sec][seg].extract(extraction_func)
        self.time = extraction_func(self.time)
        self.extraction_func=extraction_func
        self.time -= self.time[0]

    def plot_seg(self, seg, ax, elev, x_shift, color='k', **kwargs):
        """
        plot the record of a single segment into the ax
        :param seg:
        :param ax:
        :param elev:
        :param x_shift:
        :param color:
        :param kwargs:
        :return:
        """
        return self.plot_seg_str(sec_name(seg.sec), seg_name(seg), ax, elev, x_shift, color=color, **kwargs)

    def plot_seg_str(self, sec_name_, seg_name_, ax, elev, x_shift, color='k', **kwargs):
        """
        plot the record of a single segment into the ax
        :param sec_name_:
        :param seg_name_:
        :param ax:
        :param elev:
        :param x_shift:
        :param color:
        :param kwargs:
        :return:
        """
        self.record_dict[sec_name_][seg_name_].plot(ax, self.time, elev=elev, x_shift=x_shift, color=color, **kwargs)

    def plot_all(self, analyzer, ax, distance=None, distance_factor=1, plot_every=0.25, electrical=True,
                              color_distance=True, cmap=plt.cm.turbo, bounds=None,color_bar_kwarts=dict(shrink=0.6),
                              dt_func=np.mean, on_title=True, on_ylabel=False, **kwargs):
        """
        plot all the records of all the segments into the ax, you can use distance_factor to seperate the records in the y axis
        or color_distance to show the distance in diffrent colors
        :param analyzer:
        :param ax:
        :param distance:
        :param distance_factor:
        :param plot_every:
        :param electrical:
        :param color_distance:
        :param cmap:
        :param bounds:
        :param color_bar_kwarts:
        :param dt_func:
        :param kwargs:
        :return:
        """
        if distance is None:
            from Neuron_analysis_tool.distance import Distance
            distance = Distance(self.cell, analyzer.more_conductances, dt_func=dt_func)
            soma_sec = list(self.cell.soma[0])
            soma_seg = soma_sec[len(soma_sec) // 2]
            distance.compute(start_seg=soma_seg)
        t1 = self.time.copy()
        if color_distance:
            from Neuron_analysis_tool.color_func import color_func_norm
            colors = color_func_norm(distance.get_all_mid_distances(), bounds=bounds, cmap=cmap)
        else:
            colors = analyzer.colors
        for sec_name_ in self.record_dict.keys():
            for seg_name_ in self.record_dict[sec_name_].keys():
                try:
                    start_end = distance.get_start_end_str(sec_name_, seg_name_)
                    if (start_end['start'] // plot_every == start_end['end'] // plot_every) and (not sec_name_ == 'soma[0]') and (not (sec_name_==sec_name(distance.start_seg.sec) and seg_name_==seg_name(distance.start_seg))):  # show the soma and start seg in all cases!!
                        continue
                    d = distance.get_mid_point_str(sec_name_, seg_name_, electrical=electrical) * distance_factor
                    if distance.get_direction_str(sec_name_, seg_name_) == 'parent':
                        d = -d
                    color, _ = colors.get_seg_color_str(sec_name_, seg_name_)
                    self.record_dict[sec_name_][seg_name_].plot(ax, self.time, elev=d, color=color, zorder = 2 if sec_name_==sec_name(distance.start_seg.sec) and seg_name_==seg_name(distance.start_seg) else 1, **kwargs)
                except:
                    pass
        if color_distance:
            im = plt.cm.ScalarMappable(norm=colors.norm, cmap=colors.cmap)
            color_bar = plt.colorbar(im, ax=ax, **color_bar_kwarts)
            cax = color_bar.ax
            fig = ax.figure
            fig.add_axes(cax)
            if on_title:
                cax.set_title('distance')
            elif on_ylabel:
                cax.set_ylabel('distance')
            return ax, cax
        return ax

    def get_vals(self, func=lambda x: np.mean(x), default_res=0):
        """
        get the recorded values for all the sections givin a func to change from a list of record into a single vlaue
        :param func:
        :param default_res:
        :return:
        """
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
        """
        get the recorded values for all the sections in a givin time givin a func to change from a list of record into a single vlaue
        :param t:
        :param default_res:
        :return:
        """
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        res=dict()
        indexs = np.where(self.time>=t)[0]
        assert len(indexs)>0, 'the time bin dont exsists, make sure you got the correct time between 0 and '+str(self.time[-1])

        func = lambda x: x[indexs[0]]
        return self.get_vals(func=func, default_res=default_res)

    def get_vals_at_dt(self, t1, t2, default_res=0, dt_func = lambda x: np.max(x)):
        """
        get the recorded values for all the sections in a givin time window givin a func to change from a list of record into a single vlaue
        :param t1:
        :param t2:
        :param default_res:
        :param dt_func:
        :return:
        """
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
        """
        get the maximal value of all records
        :return:
        """
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return np.nanmax([np.nan if self.record_dict[sec][seg]=='non_exsisting' else self.record_dict[sec][seg].get_val(lambda x:np.max(x)) for sec in self.record_dict for seg in self.record_dict[sec]])

    def get_min(self):
        """
        get the minimal value of all the recors
        :return:
        """
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        return np.nanmin([np.nan if self.record_dict[sec][seg]=='non_exsisting' else self.record_dict[sec][seg].get_val(lambda x:np.min(x)) for sec in self.record_dict for seg in self.record_dict[sec]])

    def get_bounds(self):
        """
        get the num, max value of all the rcors
        :return:
        """
        return [self.get_min(), self.get_max()]

    def get_record_str(self, sec_name_, seg_name_):
        """
        get a specific record sec name and seg name
        :param sec_name_:
        :param seg_name_:
        :return:
        """
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        if self.record_dict[sec_name_][seg_name_]=='non_exsisting':
            return np.zeros(self.time.shape)
        return self.record_dict[sec_name_][seg_name_]._record.copy()

    def get_record(self, seg):
        """
        get a specific record seg
        :param seg:
        :return:
        """
        return self.get_record_str(sec_name(seg.sec), seg_name(seg))

    def get_record_at_dt_str(self, sec_name_, seg_name_, t1, t2, dt_func = lambda x: np.max(x)):
        """
        get a specific record sec name and seg name in a time window
        :param sec_name_:
        :param seg_name_:
        :param t1:
        :param t2:
        :param dt_func:
        :return:
        """
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        assert sec_name_ in self.record_dict, 'seg not valid '+sec_name_
        assert seg_name_ in self.record_dict[sec_name_], 'seg not valid'+seg_name_
        if t1>=self.time[-1] or t2>self.time[-1]:
            t1 = self.time[-2]
            t2 = self.time[-1]
        assert t1<t2, 'the time bins have no are not corect, t1='+str(t1)+', t2='+str(t2)

        indexs1 = np.where(self.time >= t1)[0]
        indexs2 = np.where(self.time >= t2)[0]

        assert len(indexs1) > 0, 'the time bin (t1'+str(t1)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        assert len(indexs2) > 0, 'the time bin (t2'+str(t2)+') dont exsists, make sure you got the correct time between 0 and ' + str(self.time[-1])
        func = lambda x: dt_func(x[indexs1[0]:indexs2[0]])
        if self.record_dict[sec_name_][seg_name_] == 'non_exsisting':
            return 0
        return self.record_dict[sec_name_][seg_name_].get_val(func)

    def get_record_at_dt(self, seg, t1, t2, dt_func = lambda x: np.max(x)):
        """
        get a specific record seg in a time window
        :param seg:
        :param t1:
        :param t2:
        :param dt_func:
        :return:
        """
        return self.get_record_at_dt_str(sec_name(seg.sec), seg_name(seg), t1, t2, dt_func = dt_func)

    def is_existing_str(self, sec_name_, seg_name_):
        """
        check if a section name and segment name exsit in the record
        :param sec_name_:
        :param seg_name_:
        :return:
        """
        if sec_name_ in self.record_dict:
            if seg_name_ in self.record_dict[sec_name_]:
                if not self.record_dict[sec_name_][seg_name_] == 'non_exsisting':
                    return True
        return False

    def is_existing(self, seg):
        """
        check if a segment exsit in the record
        :param seg:
        :return:
        """
        return self.is_existing_str(sec_name(seg.sec), seg_name(seg))

    def save(self, save_dir='records'):
        """
        save all the records into pickle
        :param save_dir: save path
        :return:
        """
        os.makedirs(os.path.basename(save_dir), exist_ok=True)
        if type(self.time) == neuron.hoc.HocObject:
            self.extract(lambda x: np.array(x))
        pickle.dump(dict(
            time = self.time,
            record_name=self.record_name,
            records=self.record_dict,
        ), open(save_dir, 'wb'))

    def load(self, save_dir='records'):
        """
        load all the recors from presaved recors
        :param save_dir: path to load from
        :return:
        """
        data = pickle.load(open(save_dir, 'rb'))
        self.time = data['time']
        self.record_name = data['record_name']
        self.record_dict = data['records']

class multi_record_all:
    """
    same as record_all but for multipul recording (recording all the record names
    """
    def __init__(self, cell, record_names=['v']):
        self.cell=cell
        self.record_names = record_names
        self.all_records = dict()
        for record_name in record_names:
            self.all_records[record_name] = record_all(cell, record_name=record_name)
        self.extraction_func = lambda x: x

        self.restart()

    def restart(self):
        for record_name in self.record_names:
            self.all_records[record_name].restart()

    def push_records(self, record_dict, time, record_name):
        assert record_name in self.all_records
        self.all_records[record_name].push_records(record_dict, time)

    def push_records_seg(self, record_dict, time, record_name):
        assert record_name in self.all_records
        self.all_records[record_name].push_records_seg(record_dict, time)

    def extract(self, extraction_func):
        for record_name in self.record_names:
            self.all_records[record_name].extract(extraction_func)

    def get_vals(self, func=lambda x: np.mean(x), default_res=0):
        res = dict()
        if not type(func) == list:
            func = [func]*len(self.record_names)
        for i, record_name in enumerate(self.record_names):
            res[record_name]=self.all_records[record_name].get_vals(func=func[i], default_res=default_res)
        return res

    def get_vals_at_t(self, t, default_res=0):
        res = dict()
        for record_name in self.record_names:
            res[record_name] = self.all_records[record_name].get_vals_at_t(t=t, default_res=default_res)
        return res

    def get_vals_at_dt(self, t1, t2, default_res=0, dt_func = lambda x: np.max(x)):
        res = dict()
        if not type(dt_func) == list:
            dt_func = [dt_func]*len(self.record_names)
        for i, record_name in enumerate(self.record_names):
            res[record_name] = self.all_records[record_name].get_vals_at_dt(t1=t1, t2=t2, default_res=default_res, dt_func = dt_func[i])
        return res

    def get_max(self):
        res = dict()
        for record_name in self.record_names:
            res[record_name] = self.all_records[record_name].get_max()
        return res

    def get_min(self):
        res = dict()
        for record_name in self.record_names:
            res[record_name] = self.all_records[record_name].get_min()
        return res

    def get_bounds(self):
        res = dict()
        for record_name in self.record_names:
            res[record_name] = self.all_records[record_name].get_bounds()
        return res

    def get_records(self, seg):
        res = dict()
        for record_name in self.record_names:
            res[record_name] = self.all_records[record_name].get_record(seg=seg)
        return res

    def get_record(self, seg, record_name):
        assert record_name in self.record_names
        return self.all_records[record_name].get_record(seg=seg)

    def get_record_at_dt(self, seg, t1, t2, dt_func = lambda x: np.max(x)):
        res = dict()
        if not type(dt_func) == list:
            dt_func = [dt_func] * len(self.record_names)
        for i, record_name in enumerate(self.record_names):
            res[record_name] = self.all_records[record_name].get_record_at_dt(seg=seg, t1=t1, t2=t2, dt_func=dt_func[i])
        return res

    def is_existing(self, seg):
        for record_name in self.record_names:
            if not self.all_records[record_name].is_existing(seg=seg):
                return False
        return True

    def save(self, save_dir='records'):
        os.makedirs(os.path.basename(save_dir), exist_ok=True)
        for record_name in self.record_names:
            self.all_records[record_name].save(save_dir=os.path.join(save_dir, record_name, 'records.p'))

    def load(self, save_dir='records', record_names=[]):
        self.record_names=record_names
        for record_name in self.record_names:
            self.all_records[record_name].load(save_dir=os.path.join(save_dir, record_name, 'records.p'))
