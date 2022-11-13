import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class color_func:
    def __init__(self, parts_dict, color_dict, defult_name='else'):
        self.parts_dict=parts_dict
        self.color_dict=color_dict
        self.defult_name=defult_name

    def get_seg_color(self, seg):
        for part in self.parts_dict:
            if seg in self.parts_dict[part]:
                return [self.color_dict[part], part]
        return [self.color_dict[self.defult_name], self.defult_name]

    def __call__(self, sec, parent=None, *args, **kwargs):
        if type(sec) == str:
            if sec in self.color_dict:
                return [self.color_dict[sec], sec], [1]
            return [self.color_dict[self.defult_name], self.defult_name], [1]

        colors=[]
        lengths = {self.defult_name:0}
        segs = list(sec)
        if (not sec.parentseg() is None) and (not parent == sec.parentseg().sec):
            segs = segs[::-1]
        seg_L = sec.L/sec.nseg
        for seg in segs:
            colors.append(self.get_seg_color(seg))
            part = colors[-1][1]
            if not part in lengths:
                lengths[part]=0
            lengths[part]+=seg_L

        colors_only = np.array(colors)[:,0]
        indexes = np.unique(colors_only, return_index=True)[1]
        colors = [colors[index] for index in sorted(indexes)]
        assert sum(lengths.values()) - sec.L < 1e-3, str(sum(lengths.values()))+'!='+str(sec.L)
        lengths = [lengths[part]/sec.L for [color, part] in colors]
        lengths = np.array(lengths)
        lengths /= lengths.sum()
        return colors, lengths



class color_func_norm:
    def __init__(self, value_dict, bounds=None, cmap=plt.cm.coolwarm):
        self.value_dict=value_dict
        if bounds is None:
            bounds = []
            for key in value_dict:
                bounds+=list(value_dict[key].values())
        self.cmap=cmap
        self.norm = mpl.colors.Normalize(vmin=np.min(bounds), vmax=np.max(bounds))
        self.color_dict = dict()
        for key in value_dict:
            self.color_dict[key] = dict()
            for key2 in value_dict[key]:
                self.color_dict[key][key2]= cmap(self.norm(value_dict[key][key2]))

    def get_seg_color(self, seg):
        return [self.color_dict[seg.sec][seg], '']

    def __call__(self, sec, parent=None, *args, **kwargs):
        return ['error'], ['error']


class color_func_by_func:
    def __init__(self, cell, func, bounds=None, cmap=plt.cm.coolwarm):
        self.value_dict=dict()
        for sec in cell.all:
            self.value_dict[sec] = dict()
            for seg in sec:
                self.value_dict[sec][seg] = func(seg)
        self.colors = color_func_norm(self.value_dict, bounds=bounds, cmap=cmap)
        self.cmap=cmap
        self.norm=self.colors.norm

    def get_seg_color(self, seg):
        return self.colors.get_seg_color(seg)

    def __call__(self, sec, parent=None, *args, **kwargs):
        return ['error'], ['error']
