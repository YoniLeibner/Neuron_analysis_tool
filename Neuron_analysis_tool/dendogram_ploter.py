from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from Neuron_analysis_tool.more_conductances import more_conductances_fake
from Neuron_analysis_tool.utils import *
from neuron import h

class Dendogram():
    def __init__(self,
                 cell, # person '_' cell_name
                 seg_length_function,
                 color_func,
                 dots_loc=[],
                 more_conductances = more_conductances_fake,
                 diam_factor=None, s=10,
                 fix_diam=1.0):

        self.fix_diam=fix_diam
        self.s=s
        self.cell = cell
        self.color_func = color_func
        self.dots_loc = np.array(dots_loc)
        self.tree_dendogram_dist = dict()
        self.seg_length_function = seg_length_function
        self.more_conductances = more_conductances
        self.diam_factor = diam_factor

    def cumpute_distances_helper(self, sec, reverse=False):
        start = self.roots_dict[sec]['start']
        section_length = 0
        for seg in sec:
            section_length += self.seg_length_function(seg, self.more_conductances)
        sons = sec.children()
        self.roots_dict[sec]['length'] = [section_length, section_length]
        for son in sons:
            if son not in self.roots_dict:
                self.roots_dict[sec]['sons'].append(son)
                self.roots_dict[son] = dict(start=start+section_length, is_root=False, parent=sec, sons=[], is_back=False)
                self.cumpute_distances_helper(son, reverse=False)

        if reverse and have_parent(sec):
            if sec.parentseg().sec not in self.roots_dict:
                self.roots_dict[sec]['sons'].append(sec.parentseg().sec)
                self.roots_dict[sec.parentseg().sec] = dict(start=start+section_length, is_root=False, parent=sec, sons=[], is_back=True)
                self.cumpute_distances_helper(sec.parentseg().sec, reverse=True)

    def cumpute_distances(self, base_seg):
        h.distance(0, base_seg.x, sec=base_seg.sec)
        self.base_sec = base_seg.sec
        self.base_seg = base_seg
        self.roots_dict=dict()
        length_r=0
        length_l=0
        for seg in self.base_sec:
            if seg.x < base_seg.x:
                length_l+=self.seg_length_function(seg, self.more_conductances)
            else:
                length_r += self.seg_length_function(seg, self.more_conductances)
        self.roots_dict[self.base_sec] = dict(start=0, is_root=True, parent=None, sons=[], is_back=False, length=[length_l, length_r])
        # self.roots_dict[self.base_sec]['sons']=dict()
        for sec in base_seg.sec.children():
            self.roots_dict[self.base_sec]['sons'].append(sec)
            self.roots_dict[sec] = dict(start=length_r, is_root=False, parent=self.base_sec, sons=[], is_back=False)
            # self.roots_dict[sec] = dict(start=0, is_root=False, parent=self.base_sec, sons=[], is_back=False)
            self.cumpute_distances_helper(sec)
        if have_parent(self.base_sec):
            sec=self.base_sec.parentseg().sec
            self.roots_dict[self.base_sec]['sons'].append(sec)
            self.roots_dict[sec] = dict(start=length_l, is_root=False, parent=self.base_sec, sons=[], is_back=True)
            # self.roots_dict[sec] = dict(start=0, is_root=False, parent=self.base_sec, sons=[], is_back=True)
            self.cumpute_distances_helper(sec, reverse=True)


    def get_sub_tree_size(self, sec):
        go_over = set(self.roots_dict[sec]['sons'])
        size=1
        while len(go_over)>0:
            temp_sec = go_over.pop()
            if len(self.roots_dict[temp_sec]['sons'])==0:# terminal
                size+=1
            else:
                go_over.update(self.roots_dict[temp_sec]['sons'])
        return size

    def plot_synapse(self, sec_start, sec_end, pos, x_axis, ax, mul=1):
        [color, name], _ = self.color_func("synapse")
        ax.scatter(x_axis, (sec_start + abs(sec_end - sec_start) * float(pos))*mul, color=color, s=self.s, zorder=5)

    def plot_vertical(self, x_pos, start, end, sec, ax, mul=1):
        parent = self.roots_dict[sec]['parent']
        colors, lengths = self.color_func(sec, parent)

        start_point = start
        # lengths[-1] += 1-sum(lengths)
        for [color, name], length in zip(colors, lengths):
            l = abs(start - end)*length
            y = [start_point, start_point+l]

            y = mul*np.array(y)
            ax.plot([x_pos, x_pos], y,
                     color=color,  # plot vertical
                     linewidth=self.fix_diam if self.diam_factor is None else sec.diam * self.diam_factor)

            start_point+=l

    def get_terminals(self):
        terminals = []
        for sec in self.roots_dict:
            if self.is_terminal(sec):
                terminals.append(sec)
        return terminals

    def is_terminal(self, sec):
        return len(self.roots_dict[sec]['sons'])==0

    def plot_section(self, sec, x_pos, ax, mul=1, ignore_sections = [], num=0):
        if sec in ignore_sections:
            return x_pos, x_pos
        end = self.roots_dict[sec]['start']+self.roots_dict[sec]['length'][num]
        parent = self.roots_dict[sec]['parent']
        start = self.roots_dict[sec]['start']
        sec_name = sec.name()
        sec_name = sec_name[sec_name.find(".") + 1:]

        if self.is_terminal(sec):
            self.plot_vertical(x_pos, start, end, sec, ax=ax, mul=mul)
            for sec_n, loc in self.dots_loc:
                if sec_name == sec_n:
                    self.plot_synapse(start, end, loc, x_pos, ax=ax, mul=mul)

            return x_pos+1, x_pos

        if len(self.roots_dict[sec]['sons']) == 1:
            x_pos, mid_pos = self.plot_section(self.roots_dict[sec]['sons'][0], x_pos, ax=ax, mul=mul, ignore_sections=ignore_sections)
            self.plot_vertical(mid_pos, start, end, sec, ax=ax, mul=mul)
            for sec_n, loc in self.dots_loc:
                if sec_name == sec_n:
                    self.plot_synapse(start, end, loc, mid_pos, ax=ax, mul=mul)

            return x_pos, mid_pos

        # we have several sons
        mid_points = []
        for son_sec in self.roots_dict[sec]['sons']:
            x_pos, mid_point = self.plot_section(son_sec, x_pos, ax=ax, mul=mul, ignore_sections=ignore_sections)
            mid_points.append(mid_point)
        mid_point_out = (mid_points[-1]-mid_points[0])/2+mid_points[0]
        self.plot_vertical(mid_point_out, start, end, sec, ax=ax, mul=mul)
        colors, lengths = self.color_func(sec, parent)
        y = [end] * 2
        y=mul*np.array(y)
        ax.plot([mid_points[0], mid_points[-1]], y, color=colors[-1][0], linewidth=self.fix_diam if self.diam_factor is None else sec.diam * self.diam_factor)  # plot horizontal
        for sec_n, loc in self.dots_loc:
            if sec_name == sec_n:
                self.plot_synapse(start, end, loc, mid_point_out, ax=ax, mul=mul)
        return x_pos, mid_point_out


    def plot(self, max_y=None, ax = None, plot_legend=True, ignore_sections = []):
        if ax is None:
            plt.figure(figsize=(10, 10))
            ax = plt.gca()
        x_pos = 0.0
        self.done_section = set()
        mid_points = []
        mul=1
        colors, lengths = self.color_func(self.base_sec)
        for sec in self.roots_dict[self.base_sec]['sons']:
            if sec in ignore_sections:
                print(' skiping on ', sec)
                continue
            if self.base_sec == self.cell.soma[0]:
                mul = (sec in self.cell.apic)*2-1
            else:
                mul = -1 if self.roots_dict[sec]['is_back'] else 1

            x_pos, mid_point = self.plot_section(sec, x_pos, ax=ax, mul=mul, ignore_sections=ignore_sections, num=0 if mul ==-1 else 1)
            ax.plot([mid_point, mid_point ], [0, self.roots_dict[sec]['start']*mul], color=colors[0][0], linewidth=self.fix_diam if self.diam_factor is None else self.base_sec.diam *self.diam_factor)
            mul*=-1
            mid_points.append(mid_point)

        y = [0.0, 0.0]
        ax.plot([mid_points[0], mid_points[-1]], y, color=colors[0][0], linewidth=self.fix_diam if self.diam_factor is None else self.base_sec.diam *self.diam_factor)
        plt.scatter(np.mean(mid_points), 0, color='k', s=10)
        ax.set_xticks([])

        if plot_legend:
            legend_elements = [Line2D([0], [0], color=self.color_func.color_dict[label], lw=2, label=label) for label in self.color_func.color_dict]
            ax.legend(handles=legend_elements, loc="best")
        if max_y is None:
            max_y = ax.get_ylim()[1]
        self.done_section = set()
        return max_y, x_pos
