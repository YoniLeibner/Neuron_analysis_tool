#########################################################
#
# author: Yoni Leibner
# description: Analyzer class allows to easy visualization
#              and exploration of neurons from multiple
#              point of view, and also to create neuron
#              movies.
# date of modification: 16.11.2022
#
#########################################################

from neuron import h, gui
from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake
from Neuron_analysis_tool.color_func import color_func, color_func_norm, color_func_by_func
from Neuron_analysis_tool.morph_ploter import plot_morph
from Neuron_analysis_tool.morph_ploter import get_norm
from Neuron_analysis_tool.dendogram import plot_dendogram
from Neuron_analysis_tool.cable import get_cable
from Neuron_analysis_tool.attenuation import plot_attenuation, record_to_value
from Neuron_analysis_tool.record import record, record_all
from Neuron_analysis_tool.distance import Distance
from Neuron_analysis_tool.utils import seg_Rin_func, get_segment_length_lamda, get_segment_length_um, LAMDA, MICRO
from Neuron_analysis_tool.protocols import *
from Neuron_analysis_tool.loaders import open_morph, open_swc, open_L5PC, open_ASC, open_rall_tree, get_parts_and_colors
import matplotlib.pyplot as plt
# from moviepy.editor import VideoClip
from Neuron_analysis_tool.Video import *
from moviepy.video.io.bindings import mplfig_to_npimage
import os
import matplotlib as mpl
from tqdm import tqdm


class Analyzer():
    def __init__(self, cell=None, parts_dict=None, colors_dict=None, type='input_cell',
                 morph_path = None, Rm=10000.0, Ra=100, Cm=1, e_pas=-70,
                 more_conductances_protocol =resting_protocol):
        if cell is None:
            if type == 'Rall_tree':
                cell, parts_dict, colors_dict = open_rall_tree()
            elif type == 'ASC' and morph_path:
                cell, parts_dict, colors_dict = open_ASC(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas)
            elif type == 'swc' and morph_path:
                cell, parts_dict, colors_dict = open_swc(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas)
            elif type == 'L5PC':
                cell, parts_dict, colors_dict = open_L5PC()

        if parts_dict is None:
            cell, parts_dict, colors_dict = get_parts_and_colors(cell)
        self.type=type
        self.cell=cell
        self.parts_dict = parts_dict
        self.more_conductances = more_conductances(cell, is_resting=True,
                                                   protocol=more_conductances_protocol)
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    def get_mechanism_names(self):
        mechanisms_names = set()
        for sec in self.cell.all:
            for seg in sec:
                for mechanisms in seg:
                    mechanisms_names.add(str(mechanisms))
        return list(mechanisms_names)

    def change_color_dict(self, colors_dict):
        for part in self.parts_dict:
            assert part in colors_dict

        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=self.parts_dict, color_dict=colors_dict)


    def change_parts_dict(self, parts_dict, colors_dict):
        for part in parts_dict:
            assert part in colors_dict
        self.parts_dict = parts_dict
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=self.parts_dict, color_dict=colors_dict)


    def plot_morph(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None, sec_to_change=None,
                   ignore_sections=[],
                   theta=0, scale=0, ignore_soma=True, distance=None, electrical=False,
                   time=None, dt=1, more_conductances_=None, colors=None):
        if colors is None:
            colors = self.colors
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        fig = plt.gca().figure
        ax, lines, segs = plot_morph(self.cell, color_func=colors.get_seg_color, scatter=False, add_nums=False,
                                                     seg_to_indicate=seg_to_indicate_dict,
                                                     fig=fig, ax=ax, diam_factor=diam_factor,
                                                     sec_to_change=sec_to_change,
                                                     plot_color_bar=False,
                                                     theta=theta, ignore_sections=ignore_sections,
                                                     ignore_soma=ignore_soma,
                                                     distance=distance,
                                                     more_conductances=self.more_conductances if more_conductances_ is None else more_conductances_,
                                                     electrical=electrical,
                                                     time=time, dt=dt
                                                     )
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')
            fontsize = (y_lim[1]-y_lim[0])*0.01
            ax.text(x_lim[0]-fontsize*6, y_lim[0], str(scale)+' ('+MICRO+'m)', rotation=90, fontsize=fontsize)
            ax.set_xlim([x_lim[0]-fontsize*6, x_lim[1]])
        return ax, lines, segs

    def plot_morph_with_values(self, seg_val_dict, ax=None, seg_to_indicate_dict = {}, diam_factor=None,
                               sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None, ignore_soma=True, color_bar_idx = [0.9, 0.2, 0.02, 0.6],
                               colors=None, distance=None, electrical=False, time=None, dt=1, more_conductances_=None):
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        ax.set_aspect('equal', adjustable='box')
        if colors is None:
            colors = color_func_norm(value_dict=seg_val_dict, bounds=bounds, cmap=cmap)
        ax, lines, segs =  self.plot_morph(ax=ax, seg_to_indicate_dict = seg_to_indicate_dict, diam_factor=diam_factor, sec_to_change=sec_to_change,
                   ignore_sections=ignore_sections,
                   theta=theta, scale=scale, ignore_soma=ignore_soma, distance=distance, electrical=electrical,
                   time=time, dt=dt, more_conductances_=more_conductances_, colors=colors)

        if plot_color_bar:
            cax = fig.add_axes(color_bar_idx)
            color_bar = mpl.colorbar.ColorbarBase(cax, cmap=colors.cmap, norm=colors.norm, spacing='uniform')
        else:
            color_bar = None
            cax = None
        return ax, cax, colors, lines, segs

    def plot_morph_with_value_func(self, func=seg_Rin_func, run_time=0, theta=0.0, diam_factor=None, cmap=plt.cm.turbo,
                                   ax = None, seg_to_indicate_dict = {},
                                   sec_to_change = None, ignore_sections = [], scale = 0, plot_color_bar = True, bounds = None, colors=None):
        if colors is None:
            if run_time>0:
                h.tstop=run_time
                h.run()
            value_dict = dict()
            for part in self.parts_dict:
                for seg in tqdm(self.parts_dict[part], desc=part):
                    if seg.sec not in value_dict:
                        value_dict[seg.sec]=dict()
                    value_dict[seg.sec][seg] = func(seg)
            colors = color_func_norm(value_dict=value_dict, bounds=bounds, cmap=cmap)

        return self.plot_morph_with_values({}, theta=theta, diam_factor=diam_factor, cmap=cmap, ax=ax, seg_to_indicate_dict=seg_to_indicate_dict,
                                            sec_to_change=sec_to_change, ignore_sections=ignore_sections, scale=scale, plot_color_bar=plot_color_bar, bounds=bounds, colors=colors)

    def plot_dendogram(self, start_seg = None ,ax=None, segs_to_indecate = dict(), plot_legend=True, ignore_sections=[], electrical=False, diam_factor=None, distance=None, colors=None):
        if colors is None:
            colors = self.colors
        if start_seg is None:
            start_sec = list(self.cell.soma[0])
            start_seg = start_sec[len(start_sec)//2]
        max_y, x_pos, lines, segs = plot_dendogram(self.cell, start_seg, self.more_conductances,
                       colors, ax=ax, plot_legend=plot_legend, ignore_sections=ignore_sections,
                       segs_to_indecate=segs_to_indecate, electrical=electrical, diam_factor=diam_factor, distance=distance)
        return ax, x_pos, lines, segs

    def plot_dendogram_with_values(self, seg_val_dict, start_seg = None ,ax=None, segs_to_indecate = dict(), plot_legend=True,
                                   ignore_sections=[], electrical=False, diam_factor=None, distance=None, bounds=None, cmap=plt.cm.turbo,
                                   plot_color_bar=True, color_bar_idx = [0.9, 0.2, 0.02, 0.6], colors=None):
        if start_seg is None:
            start_sec = list(self.cell.soma[0])
            start_seg = start_sec[len(start_sec)//2]
        if ax is None:
            ax = plt.gca()
        if colors is None:
            colors = color_func_norm(value_dict=seg_val_dict, bounds=bounds, cmap=cmap)

        max_y, x_pos, lines, segs = self.plot_dendogram(start_seg = start_seg ,ax=ax, segs_to_indecate = segs_to_indecate, plot_legend=plot_legend,
                            ignore_sections=ignore_sections, electrical=electrical, diam_factor=diam_factor, distance=distance, colors=colors)
        if plot_color_bar:
            cax = ax.get_figure().add_axes(color_bar_idx)
            color_bar = mpl.colorbar.ColorbarBase(cax, cmap=colors.cmap, norm=colors.norm, spacing='uniform')
        else:
            color_bar = None
            cax = None
        return ax, x_pos, cax, colors, lines, segs


    def plot_dendogram_with_value_func(self, func, start_seg = None ,ax=None, segs_to_indecate = dict(),
                                       ignore_sections=[], electrical=False, diam_factor=None, distance=None, bounds=None,
                                       cmap=plt.cm.turbo, plot_color_bar=True, color_bar_idx = [0.9, 0.2, 0.02, 0.6], colors=None):
        if colors is None:
            value_dict = dict()
            for part in self.parts_dict:
                for seg in tqdm(self.parts_dict[part], desc=part):
                    if seg.sec not in value_dict:
                        value_dict[seg.sec]=dict()
                    value_dict[seg.sec][seg] = func(seg)
            colors = color_func_norm(value_dict=value_dict, bounds=bounds, cmap=cmap)
        return self.plot_dendogram_with_values({}, start_seg=start_seg, ax=ax, segs_to_indecate=segs_to_indecate,
                                               plot_legend=False, ignore_sections=ignore_sections, electrical=electrical,
                                               diam_factor=diam_factor, distance=distance, bounds=bounds, cmap=cmap,
                                               plot_color_bar=plot_color_bar, color_bar_idx=color_bar_idx, colors=colors)

    def plot_cable(self, start_seg = None ,ax=None,
                   factor_e_space=25, factor_m_space=10 ,
                   dots_loc_seg = [], ignore_sections=[],
                   cable_type='electric', start_loc=0, x_axis=True,
                    factor=1, dots_size=10, start_color='k', plot_legend=True, distance=None, extra=5): #'d3_2', 'dist'
        if ax is None:
            ax = plt.gca()
        if start_seg is None:
            start_seg = self.cell.soma[0](0.5)

        seg_dist_dict = dict()
        for part in ['sons', 'parent']:
            seg_dist_dict[part] = dict()
            seg_dist_dict[part][self.cell.soma[0](0.5)] = []
            for hh in dots_loc_seg:
                seg_dist_dict[part][hh] = []
        results, seg_dist, cross_dist_dict = get_cable(self.cell,
                                                        start_seg=start_seg,
                                                        factor_e_space=factor_e_space,
                                                        factor_m_space=factor_m_space,
                                                        more_conductances=self.more_conductances,
                                                        seg_dist_dict=seg_dist_dict,
                                                        part_dict=self.parts_dict,
                                                        ignore_sections=ignore_sections,
                                                        distance=distance)
        max_cable = 0
        shift = 0
        for part, direction in zip(results.keys(), [1, -1]):
            cable = results[part]['all'][cable_type].flatten()
            if max_cable < cable.max() / 2.0:
                max_cable = cable.max() / 2.0
                shift = start_loc + cable.max() / factor / 2.0
        for part, direction in zip(results, [1, -1]):
            cable = results[part]['all'][cable_type].flatten()
            if cable_type == 'd3_2':
                befor_d_3_2 = np.power(cable, 3.0 / 2.0)
            cable /= factor
            y = np.arange(0, len(cable), 1) / factor_e_space

            start_pos = -cable / 2.0 + shift
            for morph_part in self.parts_dict.keys():
                remove_start_diam = False
                part_cable = results[part][morph_part][cable_type].flatten() / factor
                if cable_type == 'd3_2':
                    part_cable = results[part][morph_part][cable_type].flatten()
                    part_cable_befor_d_3_2 = np.power(part_cable, 3.0 / 2.0)
                    part_cable = cable * (part_cable_befor_d_3_2 / befor_d_3_2)
                if part_cable[1]>0 and part_cable[0] == 0:
                    part_cable[0] = (results[part]['all'][cable_type].flatten() / factor)[0]
                    remove_start_diam=True
                    temp_=start_pos[0]
                    start_pos[0] = -cable[0] / 2.0+shift
                plot_cable = part_cable[part_cable > 0]
                start_pos_temp = start_pos[part_cable > 0]
                if x_axis:
                    ax.fill_betweenx(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable, label=morph_part, color=self.colors_dict[morph_part])
                else:
                    ax.fill_between(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable, label=morph_part, color=self.colors_dict[morph_part])
                if remove_start_diam:
                    part_cable[0]=0
                    start_pos[0] = temp_
                start_pos += part_cable
        if x_axis:
            ax.scatter(shift, 0, s=dots_size, color=start_color, label='start')
        else:
            ax.scatter(0, shift, s=dots_size, color=start_color, label='start')
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        if plot_legend:
            ax.legend(by_label.values(), by_label.keys(), loc='best')
        else:
            ax.legend([], [], loc='best')
        return ax

    def plot_attenuation(self, protocol=long_pulse_protocol, ax=None, seg_to_indicate_dict=dict(), start_seg =None,
                         record_to_value_func=record_to_value, norm=True, record_name='v', norm_by=None, electrical=True, distance=None, **kwargs):

        ax, norm_by, lines, segs = plot_attenuation(cell=self.cell, start_seg=start_seg, protocol=protocol, more_conductances=self.more_conductances, color_func=self.colors, ax=ax, record_name=record_name,
                         cut_start_ms=None, record_to_value=record_to_value_func, norm_by=norm_by, norm=norm, electrical=electrical,
                         seg_to_indicate=seg_to_indicate_dict, distance=distance, **kwargs)

        ax.set_yscale('log')
        ax.set_xlabel('distance from origin (x / '+LAMDA+')')
        if norm or (not norm_by==1.0):
            ax.set_ylabel('V(x)/V(0)')
        else:
            ax.set_ylabel('attanuation (mV)')

        return ax, norm_by, lines, segs

    def create_card(self, start_seg=None, theta=0, scale=500, factor_e_space=50, cable_type='d3_2', diam_factor=None, plot_legend=False, start_color='green', start_dots_size=50, **kwargs):
        fig, ax = plt.subplots(1, 4, figsize=(12, 3), gridspec_kw={'width_ratios': [0.5, 1.5, 1, 1]})
        plt.subplots_adjust(wspace=0.5)
        if start_seg is None:
            start_seg = self.cell.soma[0]
            start_seg = list(start_seg)
            start_seg = start_seg[len(start_seg)//2]
        seg_to_indicate_dict = {start_seg: dict(size=start_dots_size, color=start_color, alpha=1)}
        distance = Distance(self.cell, self.more_conductances)
        distance.compute(start_seg=start_seg)

        _,_,_ = self.plot_morph(ax=ax[0], theta=theta, seg_to_indicate_dict=seg_to_indicate_dict,
                                scale=scale, diam_factor=diam_factor,
                                ignore_soma=not self.type.startswith('Rall_tree'), distance=distance)
        _, x_pos, _, _ = self.plot_dendogram(start_seg=start_seg, ax=ax[1], electrical=True, plot_legend=False, segs_to_indecate=seg_to_indicate_dict, distance=distance)
        self.plot_cable(start_seg=start_seg, ax=ax[1], factor_e_space=factor_e_space, cable_type=cable_type, plot_legend=plot_legend, start_loc=x_pos+15, start_color=start_color, dots_size=start_dots_size, distance=distance)
        for a in ax[1:]:
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
        ax[1].set_axis_off()
        x_lim = ax[1].get_xlim()
        y_lim = ax[1].get_ylim()

        ax[1].plot([x_lim[0], x_lim[0]+10], [y_lim[0], y_lim[0]], color='k')
        if cable_type=='d3_2':
            ax[1].text(x_lim[0], y_lim[0]-0.25, '10 ('+MICRO+'m)')
        else:

            ax[1].text(x_lim[0], y_lim[0]-0.25, '10 ('+MICRO+'$m^{2}$)')

        ax[1].plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + 0.5], color='k')
        ax[1].text(x_lim[0]-15, y_lim[0], '0.5 ('+LAMDA+')', rotation=90)
        ax[1].set_xlim([x_lim[0]-15, x_lim[1]])
        ax[1].set_ylim([y_lim[0]-0.25, y_lim[1]])
        ax[1].set_ylabel('')

        self.plot_attenuation(start_seg=start_seg, ax=ax[2], protocol=long_pulse_protocol, seg_to_indicate_dict=seg_to_indicate_dict, distance=distance, ** kwargs)
        self.plot_attenuation(start_seg=start_seg, ax=ax[3], protocol=short_pulse_protocol, seg_to_indicate_dict=seg_to_indicate_dict, distance=distance, ** kwargs)
        ylim = ax[2].get_ylim()
        ax[3].set_xlim(ax[2].get_xlim())
        ax[3].set_ylim(ylim)
        y_ticks = np.array([10**-i for i in range(5)])
        y_ticks = y_ticks[np.logical_and(y_ticks > ylim[0], y_ticks < ylim[-1])]
        if len(y_ticks)==1:
            y_ticks = np.array([10**0, 10**-1])
        ax[2].set_yticks(y_ticks)
        # ax[2].set_yticklabels(y_ticks)
        ax[3].set_ylabel('')
        ax[3].set_yticks(y_ticks)
        ax[3].set_yticklabels(['' for i in range(len(y_ticks))])
        ax[0].set_title('morphology')
        ax[1].set_title(cable_type +' dendogram')
        ax[2].set_title('long pulse attanuation')
        ax[3].set_title('short pulse attanuation')
        return fig, ax


    def create_small_card(self, start_seg=None, theta=0, scale=500, factor_e_space=50, cable_type='d3_2', diam_factor=None, plot_legend=False, start_color='green', start_dots_size=50, **kwargs):
        fig, ax = plt.subplots(1, 3, figsize=(12, 3), gridspec_kw={'width_ratios': [0.5, 1.5, 1]})
        plt.subplots_adjust(wspace=0.5)
        if start_seg is None:
            start_seg = self.cell.soma[0]
            start_seg = list(start_seg)
            start_seg = start_seg[len(start_seg)//2]
        seg_to_indicate_dict = {start_seg: dict(size=start_dots_size, color=start_color, alpha=1)}
        distance = Distance(self.cell, self.more_conductances)
        distance.compute(start_seg=start_seg)

        self.plot_morph(ax=ax[0], theta=theta, seg_to_indicate_dict=seg_to_indicate_dict,
                        scale=scale, diam_factor=diam_factor, ignore_soma=not self.type.startswith('Rall_tree'), distance=distance)
        _, x_pos, _, _ = self.plot_dendogram(start_seg=start_seg, ax=ax[1], electrical=True, plot_legend=False, segs_to_indecate=seg_to_indicate_dict, distance=distance)
        self.plot_cable(start_seg=start_seg, ax=ax[2], factor_e_space=factor_e_space, cable_type=cable_type, plot_legend=plot_legend, start_loc=0, start_color=start_color, dots_size=start_dots_size, distance=distance)
        for a in ax[1:]:
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
        ax[1].set_axis_off()
        x_lim = ax[1].get_xlim()
        y_lim = ax[1].get_ylim()

        ax[1].plot([x_lim[0], x_lim[0]+10], [y_lim[0], y_lim[0]], color='k')
        if cable_type=='d3_2':
            ax[1].text(x_lim[0], y_lim[0]-0.25, '10 ('+MICRO+'m)')
            ax[2].set_xlabel('diameter (' + MICRO + 'm)')
        else:

            ax[1].text(x_lim[0], y_lim[0]-0.25, '10 ('+MICRO+'$m^{2}$)')
            ax[2].set_xlabel('diameter (' + MICRO + '$m^{2}$)')

        ax[1].plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + 0.5], color='k')
        ax[1].text(x_lim[0]-15, y_lim[0], '0.5 ('+LAMDA+')', rotation=90)
        ax[2].set_ylabel('distance ('+LAMDA+')')

        ax[1].set_xlim([x_lim[0]-15, x_lim[1]])
        ax[1].set_ylim([y_lim[0]-0.25, y_lim[1]])
        ax[1].set_ylabel('')

        ax[0].set_title('morphology')
        ax[1].set_title('dendogram')
        ax[2].set_title('$d^{3/2}$ equivalent cable' if cable_type=='d3_2' else cable_type)
        return fig, ax


    def record_protocol(self, protocol=spike_protocol, cut_start_ms=None, record_name='v', start_seg=None, compute_more_condunctances = False):
        if start_seg is None:
            start_seg = self.cell.soma[0](0.5)
        records = record_all(self.cell, record_name=record_name)
        if compute_more_condunctances:
            more_conductances_ = more_conductances(self.cell, is_resting=False, protocol=None)
        delay, extra = protocol(self.cell, start_seg)
        if cut_start_ms is None:
            cut_start_ms = max(delay - 50, 0)
        extraction_func = lambda x: np.array(x)[int(cut_start_ms / h.dt):]
        records.extract(extraction_func )
        if compute_more_condunctances:
            more_conductances_.extract(extraction_func)
            extra['more_conductances']=more_conductances_
        if 'draw_funcs' not in extra:
            extra['draw_funcs'] = []
        return records, extra#['draw_funcs'] if 'draw_funcs' in extra else []

    def save_movie_from_rec(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                              sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                              plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4,
                              preset='ultrafast', slow_down_factor=1, func_for_missing_frames=np.max):
        save_movie_from_rec(self, record_dict, time, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                        sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta, scale=scale, cmap=cmap,
                        plot_color_bar=plot_color_bar, save_to=save_to, clip_name=clip_name, fps=fps, threads=threads,
                        preset=preset, slow_down_factor=slow_down_factor, func_for_missing_frames=func_for_missing_frames)

    def create_movie_from_rec(self, records, seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                            plot_color_bar=True, slow_down_factor=1, func_for_missing_frames=np.max, bounds = None,
                            show_records_from=dict(), voltage_window=50, ylabel='v (mV)', xlabel='time (ms)', margin=0, draw_funcs=[],
                            base_plot_type='morph', start_seg = None, electrical=True, figsize=(5,5),
                            dancing=False, distance=None, dt=1, more_conductances_=None):

        return create_movie_from_rec(self, records, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                          sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta, scale=scale, cmap=cmap,
                          plot_color_bar=plot_color_bar, slow_down_factor=slow_down_factor, func_for_missing_frames=func_for_missing_frames, bounds=bounds,
                          show_records_from=show_records_from, voltage_window=voltage_window, ylabel=ylabel, xlabel=xlabel, margin=margin,
                          draw_funcs=draw_funcs,
                          base_plot_type=base_plot_type, start_seg=start_seg, electrical=electrical, figsize=figsize,
                          dancing=dancing, distance=distance, dt=dt, more_conductances_=more_conductances_)


    def dancing_morph(self, protocol, seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=500,
                            slow_down_factor=1,
                             figsize=(5,5)):

        return dancing_morph(self, protocol, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                  sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta, scale=scale,
                  slow_down_factor=slow_down_factor,
                  figsize=figsize)
        # import matplotlib.style as mplstyle
        # more_conductances_ = more_conductances(self.cell, is_resting=False, protocol=protocol)
        # mplstyle.use('fast')
        # time = more_conductances_.time.copy()
        # time /= 1000.0
        # time *= slow_down_factor
        #
        # fig = plt.figure(figsize=figsize)
        # ax = plt.gca()
        # ax, lines, segs = self.plot_morph(ax=ax, seg_to_indicate_dict = seg_to_indicate_dict, diam_factor=diam_factor,
        #                           sec_to_change=sec_to_change,  ignore_sections=ignore_sections,
        #                           theta=theta, scale=scale, ignore_soma=True, distance=None,
        #                           electrical=True, time=0.5, dt=1, more_conductances_=more_conductances_)
        #
        # segs = np.array(segs)
        # lines = np.array(lines)
        # lim_x = ax.get_xlim()
        # lim_y = ax.get_ylim()
        # time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
        # self.last_t = 0
        # self.to_remove = []
        #
        # def make_frame_r(distance, sec, start_point):
        #     for seg in sec:
        #         for idx, seg2 in enumerate(segs):
        #             if seg2 == seg:
        #                 x, y = lines[idx].get_data()
        #                 dx = x[:-1] - x[1:]
        #                 dy = y[:-1] - y[1:]
        #                 current_lens = (dx ** 2 + dy ** 2) ** 0.5
        #                 wanted_len = distance.get_length(seg, electrical=True)
        #                 f = wanted_len / sum(current_lens)
        #                 # fixing the shift
        #                 x += start_point['x'] - x[0]
        #                 y += start_point['y'] - y[0]
        #
        #                 # fixing the length
        #                 for point_idx in range(1, len(x), 1):
        #                     new_x = x[point_idx - 1] + (x[point_idx] - x[point_idx - 1]) * f
        #                     new_y = y[point_idx - 1] + (y[point_idx] - y[point_idx - 1]) * f
        #                     x[point_idx:] += (new_x - x[point_idx])
        #                     y[point_idx:] += (new_y - y[point_idx])
        #                 lines[idx].set_data(x,y)
        #                 start_point = dict(x=x[-1], y=y[-1])
        #     for son in sec.children():
        #         make_frame_r(distance, son, start_point)
        #
        # def make_frame(t):
        #     time_in_ms = t / slow_down_factor * 1000.0
        #     distance = Distance(self.cell, more_conductances_)
        #     distance.compute(time=time_in_ms, dt=1)
        #     for son in self.cell.soma[0].children():
        #         make_frame_r(distance, son, start_point=dict(x=0,y=0))
        #     time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
        #     return mplfig_to_npimage(fig)
        # animation = VideoClip(make_frame, duration=time[-1])
        # # animation.write_gif(os.path.join(save_to, clip_name + '.gif'), fps=int(1.0 / h.dt) if fps is None else fps/slow_down_factor)
        # return animation

