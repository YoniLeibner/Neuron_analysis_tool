from neuron import h, gui
from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake
from Neuron_analysis_tool.color_func import color_func
from Neuron_analysis_tool.morph_ploter import plot_morph
from Neuron_analysis_tool.morph_ploter import plot as plot_morph2
from Neuron_analysis_tool.morph_ploter import get_norm
from Neuron_analysis_tool.dendogram_ploter import Dendogram
from Neuron_analysis_tool.cables_bi_direction import get_cable
from Neuron_analysis_tool.attenuation_plotter import run_attenuation_ploter, attenuation
from Neuron_analysis_tool.utils import *
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import os
import matplotlib as mpl
from matplotlib import animation


def seg_Rin_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

def long_pulse_protocol(cell, initial_seg):
    delay=2000.0
    dur=1000.0
    amp=.1
    h.tstop = delay+dur+500.0
    clamp = h.IClamp(initial_seg.x, sec = initial_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return h.tstop, delay, dur, amp

def short_pulse_protocol(cell, initial_seg):
    delay=2000.0
    dur=2.0
    amp=4
    h.tstop = delay+dur+20.0
    clamp = h.IClamp(initial_seg.x, sec = initial_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return h.tstop, delay, dur, amp


class Analyzer():

    def __init__(self, cell, parts_dict, colors_dict):
        self.cell=cell
        self.parts_dict = parts_dict
        self.more_conductances = more_conductances_fake(cell)
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    def plot_morph(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None, sec_to_change=None, ignore_sections=[], theta=0, scale=0):
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        fig, ax, color_bar, points_dict, lines, segs = plot_morph(self.cell, color_func=self.colors.get_seg_color, scatter=False, add_nums=False,
                                                     seg_to_indicate=seg_to_indicate_dict,
                                                     norm_colors=False, fig=fig, ax=ax, diam_factor=diam_factor,
                                                     sec_to_change=sec_to_change, bounds=None, cmap=plt.cm.turbo,
                                                     plot_color_bar=False,
                                                     theta=theta, ignore_sections=ignore_sections)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')
        return ax

    def plot_morph_with_values(self, seg_val_dict, ax=None, seg_to_indicate_dict = {}, diam_factor=None,
                               sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None):
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure

        parts_dict = dict()
        value_dict = dict()
        for part in self.parts_dict:
            for seg in self.parts_dict[part]:
                assert seg in seg_val_dict.keys()
                key = seg.sec.name() + '_' + str(seg.x)
                parts_dict[key]=[seg]
                value_dict[key] = seg_val_dict[seg]
        colors = color_func(parts_dict=parts_dict, color_dict=value_dict)
        fig, ax, color_bar, points_dict, lines, segs =  plot_morph(self.cell, color_func=colors.get_seg_color, scatter=False, add_nums=False,
                                                     seg_to_indicate=seg_to_indicate_dict,
                                                     norm_colors=True, fig=fig, ax=ax, diam_factor=diam_factor,
                                                     sec_to_change=sec_to_change, bounds=bounds, cmap=cmap,
                                                     plot_color_bar=plot_color_bar,
                                                     theta=theta, ignore_sections=ignore_sections)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')
        return ax, color_bar, points_dict, lines, segs


    def plot_morph_with_values2(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None,
                               scale=0, cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None, points_dict=None):
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        norm=get_norm(bounds)
        ax = plot_morph2(ax, points_dict, norm=norm, cmap=cmap, add_nums=False,
                         seg_to_indicate=seg_to_indicate_dict, counter=None, diam_factor=diam_factor)
        if plot_color_bar:
            cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
            color_bar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='uniform')
        else:
            color_bar=None
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')
        return ax, color_bar

    def plot_morph_with_value_func(self, func=seg_Rin_func, run_time=0):
        if run_time>0:
            h.tstop=run_time
            h.run()
        value_dict = dict()
        for part in self.parts_dict:
            for seg in self.parts_dict[part]:
                value_dict[seg] = func(seg)
        return self.plot_morph_with_values(value_dict)[:2]

    def create_morph_movie(self, protocol=short_pulse_protocol, cut_start_ms=0, record_name='v',
                            seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                            plot_color_bar=True, save_to='', clip_name='clip', fps=None):
        mpl.use('TkAgg')
        import matplotlib.style as mplstyle
        mplstyle.use('fast')
        record_dict = dict()
        for sec in self.cell.all:
            record_dict[sec] = dict()
            for seg in sec:
                try:
                    record_dict[sec][seg.x] = h.Vector()
                    record_dict[sec][seg.x].record(getattr(sec(seg.x), '_ref_' + record_name))
                except:
                    record_dict[sec][seg.x] = None

        tstop, delay, dur, amp = protocol(self.cell, self.cell.soma[0](0.5))
        total_len_idx = int((tstop - cut_start_ms) / h.dt)
        for sec in self.cell.all:
            for seg in sec:
                if record_dict[sec][seg.x] is None:
                    record_dict[sec][seg.x] = np.zeros(total_len_idx)
                else:
                    record_dict[sec][seg.x] = np.array(record_dict[sec][seg.x])[int(cut_start_ms / h.dt):]

        min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        fig, ax = plt.subplots()
        value_dict = dict()
        for sec in self.cell.all:
            for seg in sec:
                value_dict[seg] = record_dict[sec][seg.x][0]

        ax, color_bar, points_dict, lines, segs = self.plot_morph_with_values(value_dict, ax=ax,
                                                                              seg_to_indicate_dict=seg_to_indicate_dict,
                                                                              diam_factor=diam_factor,
                                                                              sec_to_change=sec_to_change,
                                                                              ignore_sections=ignore_sections,
                                                                              theta=theta, scale=scale, cmap=cmap,
                                                                              plot_color_bar=plot_color_bar,
                                                                              bounds=[min_value, max_value])

        def make_frame(t):
            time_loc = int(t / h.dt)
            norm = get_norm([min_value, max_value])
            for line, seg in zip(lines, segs):
                line.set_color(cmap(norm(record_dict[seg.sec][seg.x][time_loc])))

            return mplfig_to_npimage(fig)

        animation = VideoClip(make_frame, duration=tstop - cut_start_ms)
        # animation.write_videofile(os.path.join(save_to, clip_name+'.mp4'), fps=int(1.0/h.dt) if fps is None else fps)
        animation.write_gif(os.path.join(save_to, clip_name + '.gif'), fps=int(1.0 / h.dt) if fps is None else fps)
        mpl.use('Qt5Agg')

    def create_morph_movie2(self, protocol=short_pulse_protocol, cut_start_ms=0, record_name='v',
                                   seg_to_indicate_dict=dict(), diam_factor=None,
                                   sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                                   plot_color_bar=True, save_to='', clip_name='clip', fps=None):
        mpl.use('TkAgg')
        import matplotlib.style as mplstyle
        mplstyle.use('fast')
        record_dict = dict()
        for sec in self.cell.all:
            record_dict[sec] = dict()
            for seg in sec:
                try:
                    record_dict[sec][seg.x] = h.Vector()
                    record_dict[sec][seg.x].record(getattr(sec(seg.x), '_ref_'+record_name))
                except:
                    record_dict[sec][seg.x]=None

        tstop, delay, dur, amp = protocol(self.cell, self.cell.soma[0](0.5))
        total_len_idx = int((tstop-cut_start_ms)/h.dt)
        for sec in self.cell.all:
            for seg in sec:
                if record_dict[sec][seg.x] is None:
                    record_dict[sec][seg.x] = np.zeros(total_len_idx)
                else:
                    record_dict[sec][seg.x] = np.array(record_dict[sec][seg.x])[int(cut_start_ms/h.dt):]

        min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        fig, ax = plt.subplots()
        value_dict = dict()
        for sec in self.cell.all:
            for seg in sec:
                value_dict[seg] = record_dict[sec][seg.x][0]

        ax, color_bar, points_dict, lines, segs = self.plot_morph_with_values(value_dict, ax=ax, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                   sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta, scale=scale, cmap=cmap,
                   plot_color_bar=plot_color_bar, bounds=[min_value, max_value])

        def make_frame(t):
            time_loc = int(t/h.dt)
            norm=get_norm([min_value, max_value])
            for line, seg in zip(lines, segs):
                line.set_color(cmap(norm(record_dict[seg.sec][seg.x][time_loc])))
                # ax.draw_artist(line)
                # fig.canvas.blit(fig.bbox)
                # fig.canvas.flush_events()
            return lines
        print(int(total_len_idx*h.dt*fps), int(1.0/fps))
        anim = animation.FuncAnimation(fig, make_frame,
                                       frames=int(total_len_idx*h.dt*fps),
                                       interval=int(1.0/fps),
                                       blit=True)
        # saving to m4 using ffmpeg writer
        anim.save(os.path.join(save_to, clip_name+'.gif'), writer='PillowWriter', fps=fps)

        # writervideo = animation.FFMpegWriter(fps=fps)
        # ani.save(os.path.join(save_to, clip_name+'.mp4'), writer=writervideo)
        plt.close()
        mpl.use('Qt5Agg')

    def plot_morph_with_value_func_after_protocol(self, func=seg_Rin_func, run_time=0):
        pass

    def plot_dendogram(self, initial_seg = None ,ax=None, dots_loc_seg = [], plot_legend=True, ignore_sections=[], electrical=False):
        if ax is None:
            ax = plt.gca()
        if initial_seg is None:
            initial_seg = self.cell.soma[0](0.5)
        dendogram = Dendogram(self.cell,
                               seg_length_function=get_segment_length_lamda if electrical else get_segment_length_um,
                               color_func=self.colors,
                               dots_loc=[[seg.sec.name().split('.')[-1], seg.x] for seg in dots_loc_seg],
                               more_conductances=self.more_conductances,
                               diam_factor=None, s=10, fix_diam=1.)
        dendogram.cumpute_distances(initial_seg)
        max_y, x_pos = dendogram.plot(ax=ax, plot_legend=plot_legend, ignore_sections=ignore_sections)
        return ax

    def plot_cable(self, initial_seg = None ,ax=None,
                   factor_e_space=25, factor_m_space=10 ,
                   dots_loc_seg = [], ignore_sections=[],
                   cable_type='electric', start_loc=0, x_axis=True,
                    factor=1, dots_size=10, start_color='k', plot_legend=True): #'d3_2', 'dist'
        if ax is None:
            ax = plt.gca()
        if initial_seg is None:
            initial_seg = self.cell.soma[0](0.5)

        seg_dist_dict = dict()
        for part in ['sons', 'parent']:
            seg_dist_dict[part] = dict()
            seg_dist_dict[part][self.cell.soma[0](0.5)] = []
            for hh in dots_loc_seg:
                seg_dist_dict[part][hh] = []
        results, seg_dist, cross_dist_dict = get_cable( self.cell,
                                                        factor_e_space=factor_e_space,
                                                        factor_m_space=factor_m_space,
                                                        start_section=initial_seg.sec,
                                                        x_start=initial_seg.x,
                                                        more_conductions=self.more_conductances,
                                                        seg_dist_dict=seg_dist_dict,
                                                        cross_dist_dict=[], # not implemented yet
                                                        part_dict=self.parts_dict,
                                                        ignore_sections=ignore_sections)
        max_cable = 0
        shift = 0
        for part, direction in zip(results.keys(), [1, -1]):
            cable = results[part]['all'][cable_type].flatten()
            if max_cable < cable.max() / 2:
                max_cable = cable.max() / 2
                shift = start_loc + cable.max() / factor / 2 + 5
        for part, direction in zip(results, [1, -1]):
            cable = results[part]['all'][cable_type].flatten()
            if cable_type == 'd3_2':
                befor_d_3_2 = np.power(cable, 3.0 / 2.0)
            cable /= factor
            y = np.arange(0, len(cable), 1) / factor_e_space

            start_pos = -cable / 2.0 + shift

            for morph_part in self.parts_dict.keys():#'basal', 'tuft', 'trunk', 'oblique']:
                part_cable = results[part][morph_part][cable_type].flatten() / factor
                if cable_type == 'd3_2':
                    part_cable = results[part][morph_part][cable_type].flatten()
                    part_cable_befor_d_3_2 = np.power(part_cable, 3.0 / 2.0)
                    part_cable = cable * (part_cable_befor_d_3_2 / befor_d_3_2)
                plot_cable = part_cable[part_cable > 0]
                start_pos_temp = start_pos[part_cable > 0]
                if x_axis:
                    ax.fill_betweenx(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable, label=morph_part, color=self.colors_dict[morph_part])
                else:
                    ax.fill_between(direction * y[part_cable > 0], start_pos_temp, start_pos_temp + plot_cable, label=morph_part, color=self.colors_dict[morph_part])
                start_pos += part_cable
        if x_axis:
            ax.scatter(shift, 0, s=dots_size, color=start_color, label='start')
        else:
            ax.scatter(0, shift, s=dots_size, color=start_color, label='start')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        if plot_legend:
            ax.legend(by_label.values(), by_label.keys(), loc='best')
        else:
            ax.legend([], [], loc='best')
        return ax

    def plot_attanuation(self, protocol=long_pulse_protocol, ax=None, seg_to_indicate=[], indication_color='orange', initial_seg =None, record_to_value_func=None, norm=True):
        if ax is None:
            ax = plt.gca()
        if initial_seg is None:
            initial_seg = self.cell.soma[0](0.5)

        att = attenuation(self.cell, color_func=self.colors, seg_length_function=get_segment_length_lamda,
                          more_conductances=self.more_conductances, param_to_record='v',
                          record_to_value_func=record_to_value_func)
        tstop, delay, dur, amp = protocol(self.cell, initial_seg)
        ax = att.plot(start_seg=initial_seg, norm=norm, cut_start=int((delay - 1) / h.dt),
                      seg_to_indicate={seg:dict(size=30, color=indication_color, alpha=1) for seg in seg_to_indicate}, ax=ax)
        ax.set_yscale('log')
        return ax

    def create_card(self):
        pass
