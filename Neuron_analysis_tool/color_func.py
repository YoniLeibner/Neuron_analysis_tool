##########################################################################################
#
# author: Yoni Leibner
# description: color func based on each segment
#                       color_func - by segment grouping
#                       color_func_norm - by segment values
#                       color_func_by_func - by segment values based on a function
# date of modification: 11.05.2023
#
##########################################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from Neuron_analysis_tool.utils import sec_name, seg_name


class color_func:
    """
    class that gives color to each segment baed on there part name
    """
    def __init__(self, parts_dict, color_dict, defult_name='else'):
        self.parts_dict = parts_dict
        self.parts_dict_str = dict()
        for key in parts_dict:
            self.parts_dict_str[key] = set()
            for seg in self.parts_dict[key]:
                self.parts_dict_str[key].add((sec_name(seg.sec), seg_name(seg)))

        self.color_dict = color_dict
        self.defult_name = defult_name

    def get_seg_color(self, seg):
        return self.get_seg_color_str(sec_name(seg.sec), seg_name(seg))
        # for part in self.parts_dict:
        #     if seg in self.parts_dict[part]:
        #         return [self.color_dict[part], part]
        # return [self.color_dict[self.defult_name], self.defult_name]

    def get_seg_color_str(self, sec_name_, seg_name_):
        for part in self.parts_dict:
            if (sec_name_, seg_name_) in self.parts_dict_str[part]:
                return [self.color_dict[part], part]
        return [self.color_dict[self.defult_name], self.defult_name]
    #
    # def __call__(self, sec, parent=None, *args, **kwargs):
    #     """
    #     deprecated
    #     :param sec:
    #     :param parent:
    #     :param args:
    #     :param kwargs:
    #     :return:
    #     """
    #     if type(sec) == str:
    #         if sec in self.color_dict:
    #             return [self.color_dict[sec], sec], [1]
    #         return [self.color_dict[self.defult_name], self.defult_name], [1]
    #
    #     colors=[]
    #     lengths = {self.defult_name:0}
    #     segs = list(sec)
    #     if (not sec.parentseg() is None) and (not parent == sec.parentseg().sec):
    #         segs = segs[::-1]
    #     seg_L = sec.L/sec.nseg
    #     for seg in segs:
    #         colors.append(self.get_seg_color(seg))
    #         part = colors[-1][1]
    #         if not part in lengths:
    #             lengths[part]=0
    #         lengths[part]+=seg_L
    #
    #     colors_only = np.array(colors)[:,0]
    #     indexes = np.unique(colors_only, return_index=True)[1]
    #     colors = [colors[index] for index in sorted(indexes)]
    #     assert sum(lengths.values()) - sec.L < 1e-3, str(sum(lengths.values()))+'!='+str(sec.L)
    #     lengths = [lengths[part]/sec.L for [color, part] in colors]
    #     lengths = np.array(lengths)
    #     lengths /= lengths.sum()
    #     return colors, lengths


class color_func_norm:
    """
    class that gives color to each segment baed on a givin value
    """
    def __init__(self, value_dict, bounds=None, cmap=plt.cm.coolwarm):
        """

        :param value_dict: dicinary of {sec name:{seg name: value}}
        :param bounds: the bounds for the normalaziation of the cmap
        :param cmap: the color map to use (matplotlib cmap)
        """
        self.value_dict=value_dict
        if bounds is None:
            bounds = []
            for key in value_dict:
                bounds+=list(value_dict[key].values())
        self.cmap=cmap
        self.norm = mpl.colors.Normalize(vmin=np.min(bounds), vmax=np.max(bounds))
        self.value_dict=value_dict
        self.color_dict = dict()
        for key in value_dict:
            self.color_dict[key] = dict()
            for key2 in value_dict[key]:
                self.color_dict[key][key2]= cmap(self.norm(value_dict[key][key2]))

    def get_seg_color(self, seg):
        """

        :param seg:
        :return: the segment color
        """
        return self.get_seg_color_str(sec_name(seg.sec), seg_name(seg))


    def get_seg_color_str(self, sec_name_, seg_name_):
        return [self.color_dict[sec_name_][seg_name_], '']


    def change_cmap(self, cmap):
        """
        change the cmap used
        :param cmap: the new cmap (matplotlib cmap)
        :return:
        """
        self.cmap=cmap
        self.color_dict = dict()
        for key in self.value_dict:
            self.color_dict[key] = dict()
            for key2 in self.value_dict[key]:
                self.color_dict[key][key2] = cmap(self.norm(self.value_dict[key][key2]))

    def change_bounds(self, low_bound, high_bound):
        """
        change the bound of the cmap
        :param low_bound:
        :param high_bound:
        :return:
        """
        self.norm = mpl.colors.Normalize(vmin=low_bound, vmax=high_bound)
        self.change_cmap(self.cmap)


class color_func_by_func:
    """
    class that gives color to each segment baed on a func to each segment
    """
    def __init__(self, cell, func, bounds=None, cmap=plt.cm.coolwarm):
        """

        :param cell: the cell model (Neuron model
        :param func: the function to run for each segment
        :param bounds: the bounds for the normalaziation of the cmap
        :param cmap: the color map to use (matplotlib cmap)
        """
        self.value_dict=dict()
        for sec in cell.all:
            self.value_dict[sec] = dict()
            for seg in sec:
                self.value_dict[sec][seg] = func(seg)
        self.colors = color_func_norm(self.value_dict, bounds=bounds, cmap=cmap)
        self.cmap=cmap
        self.norm=self.colors.norm

    def get_seg_color(self, seg):
        """

        :param seg:
        :return: the segment color
        """
        return self.colors.get_seg_color(seg)

    def change_cmap(self, cmap):
        """
        change the cmap used
        :param cmap: the new cmap (matplotlib cmap)
        :return:
        """
        self.colors.change_cmap(cmap)

    def change_bounds(self, low_bound, high_bound):
        """
        change the bound of the cmap
        :param low_bound:
        :param high_bound:
        :return:
        """
        self.colors.change_bounds(low_bound, high_bound)
