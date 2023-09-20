#########################################################
#
# author: Yoni Leibner
# description: class that calculate the distance of each
#              segment from the start seg both in
#              micro-meters and in electrical units
# date of modification: 11.05.2023
#
#########################################################

from neuron import h
import numpy as np
from Neuron_analysis_tool.utils import get_segment_length_lamda, get_segment_length_um
from .utils import seg_name, sec_name


class Distance:
    """
    class that calculate the physical and electrical distance from a givin segment
    """
    def __init__(self, cell, more_conductances, dt_func= lambda x: np.mean(x)):
        """

        :param cell: the cell to plot (Neuron model)
        :param more_conductances: more_conductances to initiate the distance if distance is None
        :param dt_func: function for dt in the more_conductances
        """
        self.cell=cell
        self.distance_dict=dict()
        self.start_seg = None
        self.more_conductances=more_conductances
        self.dt_func=dt_func

    def compute(self, start_seg=None, time=None, dt=1, dt_func= None):
        """
        compute the distance for each point on the tree
        :param start_seg: the segment to set as the origin
        :param time: time to get the more_conductances from
        :param dt: the dt around this time
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        if dt_func is None:
            dt_func=self.dt_func
        else:
            self.dt_func=dt_func
        if start_seg == None:
            start_seg = list(self.cell.soma[0])
            start_seg = start_seg[len(start_seg)//2]
        if start_seg == self.start_seg:
            return
        self.start_seg=start_seg
        self.distance_dict = dict()

        sec = start_seg.sec
        cell=sec.cell()
        if sec==cell.soma[0]:
            sons = list(cell.apical)
            parent = list(cell.basal)
        else:
            sons = sec.children()
            if sec.parentseg() is None:
                parent = []
            else:
                parent = [sec.parentseg().sec]
        h.distance(0, start_seg.x, sec=sec)
        segs = [seg for seg in sec]
        done = set()
        done.add(sec)
        parent_seg = start_seg
        self.distance_dict[sec_name(sec)] = dict(segs=dict(), parent_seg=None, sec_sons=[])
        self.distance_dict[sec_name(sec)]['segs'][seg_name(start_seg)] = dict(
                                        length=dict(start=0, end=get_segment_length_um(start_seg)),
                                        electrical_length=dict(start=0,
                                                               end=get_segment_length_lamda(start_seg,
                                                                                            self.more_conductances,
                                                                                            time=time, dt=dt,
                                                                                            dt_func=dt_func)),
                                        parent=parent_seg, direction='sons')
        for seg in segs:
            if seg.x > start_seg.x:
                self.distance_dict[sec_name(sec)]['segs'][seg_name(seg)] = dict(
                    length=dict(start = self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['length']['end'],
                                end = get_segment_length_um(seg) +
                                      self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['length']['end']),
                    electrical_length=dict(start = self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances, time=time, dt=dt,
                                                                        dt_func=dt_func) +
                                               self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end']),
                    parent=parent_seg,direction='sons')
                parent_seg = seg
        # now we go to the sones and
        for son in sons:# sec.children():
            done = self.compute_distances_helper(son, parent_seg=seg, done=done, direction='sons', time=time, dt=dt)

        parent_seg = start_seg
        for seg in segs[::-1]:
            if seg.x < start_seg.x:
                self.distance_dict[sec_name(sec)]['segs'][seg_name(seg)] = dict(
                    length=dict(start=self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['length']['end'],
                                end=get_segment_length_um(seg) +
                                    self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['length']['end']),
                    electrical_length=dict(start=self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances, time=time, dt=dt,
                                                                        dt_func=dt_func) +
                                               self.distance_dict[sec_name(sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end']),
                    parent=parent_seg, direction='parent')
                parent_seg = seg
        if len(parent)>0:
        # if not sec.parentseg() is None:
            for son in parent:
                done = self.compute_distances_helper(son, parent_seg=seg, done=done, reverse=True, direction='parent', time=time, dt=dt)


    def compute_distances_helper(self, sec, parent_seg, done, reverse=False, direction='sons', time=None, dt=1):
        """
        helper function that goes on all the tree to get the distances
        :param sec: current section
        :param parent_seg: the parent segmant (depand on the segment origin)
        :param done: section already calculated
        :param reverse: if this direction is in reverse (we got here from one of our son sections)
        :param part:
        :param time:
        :param dt:
        :return:
        """
        if sec in done:
            return done
        self.distance_dict[sec_name(parent_seg.sec)]['sec_sons'].append(sec)
        done.add(sec)
        segs = list(sec)
        if reverse:
            segs = segs[::-1]
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=parent_seg, done=done, direction=direction,
                                                     time=time, dt=dt)
        self.distance_dict[sec_name(sec)] = dict(segs=dict(), parent_seg=parent_seg, sec_sons=[])
        for seg in segs:
            self.distance_dict[sec_name(sec)]['segs'][seg_name(seg)] = dict(
                length=dict(start=self.distance_dict[sec_name(parent_seg.sec)]['segs'][seg_name(parent_seg)]['length']['end'],
                            end=get_segment_length_um(seg) +
                                self.distance_dict[sec_name(parent_seg.sec)]['segs'][seg_name(parent_seg)]['length']['end']),
                electrical_length=dict(start=self.distance_dict[sec_name(parent_seg.sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end'],
                                       end=get_segment_length_lamda(seg, self.more_conductances, time=time, dt=dt,
                                        dt_func=self.dt_func) +
                                           self.distance_dict[sec_name(parent_seg.sec)]['segs'][seg_name(parent_seg)]['electrical_length']['end']),
                parent=parent_seg, direction=direction)
            parent_seg = seg

        if reverse:
            if not sec.parentseg() is None:
                done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=segs[-1], done=done, reverse=True,
                                                     direction=direction, time=time, dt=dt)
        else:
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=segs[-1], done=done, direction=direction,
                                                     time=time, dt=dt)
        return done

    def get_all_mid_distances(self):
        results = dict()
        for sec_name_ in self.distance_dict:
            results[sec_name_]=dict()
            for seg_name_ in self.distance_dict[sec_name_]['segs']:
                results[sec_name_][seg_name_] = self.get_mid_point_str(sec_name_, seg_name_)
        return results

    def get_start_end(self, seg, electrical=True):
        """
        get the start and end distance of a givin segment from the start segment
        :param seg: the givin segment
        :param electrical: if True it will plot the dendogram in electrical units
        :return: dictionary of {start: value, end: value}
        """
        return self.get_start_end_str(sec_name(seg.sec), seg_name(seg), electrical=electrical)
        # if self.start_seg is None:
        #     print('you forgot to compute')
        #     self.compute()
        # if (seg.sec not in self.distance_dict) or (seg not in self.distance_dict[sec_name(seg.sec)]['segs']):
        #     return dict(start=0, end=0)
        # if electrical:
        #     return self.distance_dict[sec_name(seg.sec)]['segs'][seg_name(seg)]['electrical_length']
        # return self.distance_dict[sec_name(seg.sec)]['segs'][seg_name(seg)]['length']

    def get_start_end_str(self, sec_name_, seg_name_, electrical=True):
        """
        get the start and end distance of a givin segment from the start segment
        :param seg: the givin segment
        :param electrical: if True it will plot the dendogram in electrical units
        :return: dictionary of {start: value, end: value}
        """
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        if (sec_name_ not in self.distance_dict) or (seg_name_ not in self.distance_dict[sec_name_]['segs']):
            return dict(start=0, end=0)
        if electrical:
            return self.distance_dict[sec_name_]['segs'][seg_name_]['electrical_length']
        return self.distance_dict[sec_name_]['segs'][seg_name_]['length']

    def get_sec_start_end(self, sec, electrical=True):
        """
        get the start and end distance of a givin section from the start segment
        :param sec: the givin section
        :param electrical: if True it will plot the dendogram in electrical units
        :return: dictionary of {start: value, end: value}
        """
        starts=[]
        ends=[]
        for seg in sec:
            start_end = self.get_start_end(seg, electrical=electrical)
            starts.append(start_end['start'])
            ends.append(start_end['end'])

        return dict(start=np.min(starts), end=np.max(ends))

    def get_sec_start_end_direction(self, sec, direction, electrical=True):
        """
        get the start and end of a section based on its direction
        :param sec: the givin section
        :param direction: the deriction to look in
        :param electrical:if True it will plot the dendogram in electrical units
        :return: dictionary of {start: value, end: value}
        """
        starts = []
        ends = []
        for seg in sec:
            if self.get_direction(seg)==direction:
                start_end = self.get_start_end(seg, electrical=electrical)
                starts.append(start_end['start'])
                ends.append(start_end['end'])
        try:
            return dict(start=np.min(starts), end=np.max(ends))
        except:
            return self.get_start_end(self.start_seg, electrical=electrical)

    def get_mid_point(self, seg, electrical=True):
        """
        geting the distance for middle point of a givin segment
        :param seg: the givin segment
        :param electrical: if True it will plot the dendogram in electrical units
        :return: distance for middle point of a givin segment
        """
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        start_end = self.get_start_end(seg, electrical=electrical)
        return (start_end['end']+start_end['start'])/2.0

    def get_mid_point_str(self, sec_name_, seg_name_, electrical=True):
        """
        geting the distance for middle point of a givin segment
        :param seg: the givin segment
        :param electrical: if True it will plot the dendogram in electrical units
        :return: distance for middle point of a givin segment
        """
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        start_end = self.get_start_end_str(sec_name_, seg_name_, electrical=electrical)
        return (start_end['end']+start_end['start'])/2.0

    def get_length(self, seg, electrical=True):
        """
        get the length of a givin segment
        :param seg: the givin segment
        :param electrical: if True it will plot the dendogram in electrical units
        :return: the segment length
        """
        start_end = self.get_start_end(seg, electrical=electrical)
        return start_end['end']-start_end['start']

    def get_direction(self, seg):
        """
        get the direction of a givin segment
        :param seg: the givin segment
        :return: direction of a givin segment
        """
        return self.get_direction_str(sec_name(seg.sec), seg_name(seg))
        # if self.start_seg is None:
        #     print('you forgot to compute')
        #     self.compute()
        # return self.distance_dict[sec_name(seg.sec)]['segs'][seg_name(seg)]['direction']

    def get_direction_str(self, sec_name_, seg_name_):
        """
        get the direction of a givin segment
        :param seg: the givin segment
        :return: direction of a givin segment
        """
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        return self.distance_dict[sec_name_]['segs'][seg_name_]['direction']

    def get_max(self, electrical=True):
        """
        get the maximal distance by direction
        :param electrical: if True it will plot the dendogram in electrical units
        :return: dictionary of {parent: value, sons: value}
        """
        max_dict = dict(parent=0, sons=0)
        for sec in self.distance_dict:
            for seg in self.distance_dict[sec]['segs']:
                if max_dict[self.get_direction_str(sec, seg)] < self.get_start_end_str(sec, seg, electrical=electrical)['end']:
                    max_dict[self.get_direction_str(sec, seg)] = self.get_start_end_str(sec, seg, electrical=electrical)['end']
        return max_dict

    def is_terminal(self, sec):
        """
        check if a givin section is a terminal section
        :param sec: thegivin section
        :return: boolean, if the section is a terminal section
        """
        if sec_name(sec) in self.distance_dict:
            return len(self.distance_dict[sec_name(sec)]['sec_sons'])==0
        return True

    def get_sons(self, sec):
        """
        get all the section sons
        :param sec:
        :return:
        """
        return self.distance_dict[sec_name(sec)]['sec_sons']

    def get_sec_parent(self, sec):
        """
        get a section parent
        :param sec:
        :return: the section parent
        """
        return self.distance_dict[sec_name(sec)]['parent_seg']

    def get_seg_parent(self, seg):
        """
        get segment parent
        :param seg:
        :return:
        """
        return self.distance_dict[sec_name(seg.sec)]['segs'][seg_name(seg)]['parent']

    def get_segs(self, sec):
        """
        :param sec:
        :return: list of dictionary of {start: value, end: value} for all the segment in a givin section
        """
        if sec_name(sec) in self.distance_dict:
            return [sec(float(x)) for x in self.distance_dict[sec_name(sec)]['segs']]
        return []
