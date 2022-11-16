#########################################################
#
# author: Yoni Leibner
# description: class that calculate the distance of each
#              segment from the start seg both in
#              micro-meters and in electrical units
# date of modification: 16.11.2022
#
#########################################################

from neuron import h
import numpy as np
from Neuron_analysis_tool.utils import get_segment_length_lamda, get_segment_length_um


class Distance:
    def __init__(self, cell, more_conductances):
        self.cell=cell
        self.distance_dict=dict()
        self.start_seg = None
        self.more_conductances=more_conductances

    def compute(self, start_seg=None):
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
        self.distance_dict[sec] = dict(segs=dict(), parent_seg=None, sec_sons=[])
        self.distance_dict[sec]['segs'][start_seg] = dict(
                                        length=dict(start=0, end=get_segment_length_um(start_seg, self.more_conductances)),
                                        electrical_length=dict(start=0, end=get_segment_length_lamda(start_seg, self.more_conductances)),
                                        parent=parent_seg, part='sons')
        for seg in segs:
            if seg.x > start_seg.x:
                self.distance_dict[sec]['segs'][seg] = dict(
                    length=dict(start = self.distance_dict[sec]['segs'][parent_seg]['length']['end'],
                                end = get_segment_length_um(seg, self.more_conductances) + self.distance_dict[sec]['segs'][parent_seg]['length']['end']),
                    electrical_length=dict(start = self.distance_dict[sec]['segs'][parent_seg]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[sec]['segs'][parent_seg]['electrical_length']['end']),
                    parent=parent_seg,part='sons')
                parent_seg = seg
        # now we go to the sones and
        for son in sons:# sec.children():
            done = self.compute_distances_helper(son, parent_seg=seg, done=done, part='sons')

        parent_seg = start_seg
        for seg in segs[::-1]:
            if seg.x < start_seg.x:
                self.distance_dict[sec]['segs'][seg] = dict(
                    length=dict(start=self.distance_dict[sec]['segs'][parent_seg]['length']['end'],
                                end=get_segment_length_um(seg, self.more_conductances) + self.distance_dict[sec]['segs'][parent_seg]['length']['end']),
                    electrical_length=dict(start=self.distance_dict[sec]['segs'][parent_seg]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[sec]['segs'][parent_seg]['electrical_length']['end']),
                    parent=parent_seg, part='parent')
                parent_seg = seg
        if len(parent)>0:
        # if not sec.parentseg() is None:
            for son in parent:
                done = self.compute_distances_helper(son, parent_seg=seg, done=done, reverse=True, part='parent')


    def compute_distances_helper(self, sec, parent_seg, done, reverse=False, part='sons'):
        if sec in done:
            return done
        self.distance_dict[parent_seg.sec]['sec_sons'].append(sec)
        done.add(sec)
        segs = list(sec)
        if reverse:
            segs = segs[::-1]
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=parent_seg, done=done, part=part)
        self.distance_dict[sec] = dict(segs=dict(), parent_seg=parent_seg, sec_sons=[])
        for seg in segs:
            self.distance_dict[sec]['segs'][seg] = dict(
                length=dict(start=self.distance_dict[parent_seg.sec]['segs'][parent_seg]['length']['end'],
                            end=get_segment_length_um(seg, self.more_conductances) + self.distance_dict[parent_seg.sec]['segs'][parent_seg]['length']['end']),
                electrical_length=dict(start=self.distance_dict[parent_seg.sec]['segs'][parent_seg]['electrical_length']['end'],
                                       end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[parent_seg.sec]['segs'][parent_seg]['electrical_length']['end']),
                parent=parent_seg, part=part)
            parent_seg = seg

        if reverse:
            if not sec.parentseg() is None:
                done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=segs[-1], done=done, reverse=True, part=part)
        else:
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=segs[-1], done=done, part=part)
        return done

    def get_start_end(self, seg, electrical=True):
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        if electrical:
            return self.distance_dict[seg.sec]['segs'][seg]['electrical_length']
        return self.distance_dict[seg.sec]['segs'][seg]['length']

    def get_sec_start_end(self, sec, electrical=True):
        starts=[]
        ends=[]
        for seg in sec:
            start_end = self.get_start_end(seg, electrical=electrical)
            starts.append(start_end['start'])
            ends.append(start_end['end'])

        return dict(start=np.min(starts), end=np.max(ends))

    def get_sec_start_end_part(self, sec, part, electrical=True):
        starts = []
        ends = []
        for seg in sec:
            if self.get_part(seg)==part:
                start_end = self.get_start_end(seg, electrical=electrical)
                starts.append(start_end['start'])
                ends.append(start_end['end'])
        try:
            return dict(start=np.min(starts), end=np.max(ends))
        except:
            return self.get_start_end(self.start_seg, electrical=electrical)

    def get_mid_point(self, seg, electrical=True):
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        start_end = self.get_start_end(seg, electrical=electrical)
        return (start_end['end']+start_end['start'])/2.0

    def get_length(self, seg, electrical=True):
        start_end = self.get_start_end(seg, electrical=electrical)
        return start_end['end']-start_end['start']

    def get_part(self, seg):
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        return self.distance_dict[seg.sec]['segs'][seg]['part']

    def is_terminal(self, sec):
        return len(self.distance_dict[sec]['sec_sons'])==0

    def get_sons(self, sec):
        return self.distance_dict[sec]['sec_sons']

    def get_sec_parent(self, sec):
        return self.distance_dict[sec]['parent_seg']

    def get_seg_parent(self, seg):
        return self.distance_dict[seg.sec]['segs'][seg]['parent']

    def get_segs(self, sec):
        return self.distance_dict[sec]['segs']

