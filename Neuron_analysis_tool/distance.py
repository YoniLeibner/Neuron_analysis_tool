import neuron
from neuron import h
import numpy as np
from Neuron_analysis_tool.utils import get_segment_length_lamda, get_segment_length_um


class distance:
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

        h.distance(0, start_seg.x, sec=sec)
        segs = [seg for seg in sec]
        done = set()
        done.add(sec)
        parent_seg = start_seg
        self.distance_dict[sec] = dict(segs=dict(), parent_seg=None, sec_sons=[])
        self.distance_dict[sec]['segs'][start_seg] = dict(
                                        length=dict(start=0, end=get_segment_length_um(start_seg, self.more_conductances)),
                                        electrical_length=dict(start=0, end=get_segment_length_lamda(start_seg, self.more_conductances)),
                                        parent=parent_seg)
        for seg in segs:
            if seg.x > start_seg.x:
                self.distance_dict[sec]['segs'][seg] = dict(
                    length=dict(start = self.distance_dict[sec][parent_seg]['length']['end'],
                                end = get_segment_length_um(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['length']['end']),
                    electrical_length=dict(start = self.distance_dict[sec][parent_seg]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['electrical_length']['end']),
                    parent=parent_seg)
                parent_seg = seg
        # now we go to the sones and
        for son in sec.children():
            done = self.compute_distances_helper(son, parent_seg=seg, done=done)

        parent_seg = start_seg
        for seg in segs[::-1]:
            if seg.x < start_seg.x:
                self.distance_dict[sec]['segs'][start_seg] = dict(
                    length=dict(start=self.distance_dict[sec][parent_seg]['length']['end'],
                                end=get_segment_length_um(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['length']['end']),
                    electrical_length=dict(start=self.distance_dict[sec][parent_seg]['electrical_length']['end'],
                                           end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['electrical_length']['end']),
                    parent=parent_seg)
                parent_seg = seg
        if not sec.parentseg() is None:
            done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=seg, done=done, reverse=True)


    def compute_distances_helper(self, sec, parent_seg, done, reverse=False):
        if sec in done:
            return done
        self.distance_dict[parent_seg.sec]['sons'].append(sec)
        done.add(sec)
        segs = list(sec)
        if reverse:
            segs = segs[::-1]
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=parent_seg, done=done)
        self.distance_dict[sec] = dict(segs=dict(), parent_seg=None, sec_sons=[])
        for seg in segs:
            self.distance_dict[sec]['segs'][seg] = dict(
                length=dict(start=self.distance_dict[sec][parent_seg]['length']['end'],
                            end=get_segment_length_um(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['length']['end']),
                electrical_length=dict(start=self.distance_dict[sec][parent_seg]['electrical_length']['end'],
                                       end=get_segment_length_lamda(seg, self.more_conductances) + self.distance_dict[sec][parent_seg]['electrical_length']['end']),
                parent=parent_seg
            )
            parent_seg = seg

        if reverse:
            if not sec.parentseg() is None:
                done = self.compute_distances_helper(sec.parentseg().sec, parent_seg=segs[-1], done=done, reverse=True)
        else:
            for son in sec.children():
                done = self.compute_distances_helper(son, parent_seg=segs[-1], done=done)
        return done

    def get_start_end(self, seg, electrical=True):
        if self.start_seg is None:
            print('you forgot to compute')
            self.compute()
        if electrical:
            return self.distance_dict[seg.sec]['segs'][seg]['electrical_length']
        return self.distance_dict[seg.sec]['segs'][seg]['length']

    def get_length(self, seg, electrical=True):
        start_end = self.get_start_end(seg, electrical=electrical)
        return start_end['end']-start_end['start']
