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
from tqdm import tqdm
from io import BytesIO
from matplotlib import animation
from matplotlib.animation import FuncAnimation


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

def spike_protocol(cell, initial_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt=spike_data.T[0][1]-spike_data.T[0][0]
    V = np.concatenate([np.zeros(int(1000.0/dt))]+[spike_data.T[1]]*10)
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)

    clamp = h.SEClamp(initial_seg.x, sec=initial_seg.sec)
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    spike_vec.play(clamp._ref_amp1, spike_data.T[0][1]-spike_data.T[0][0])

    # delay=2000.0
    # dur=2.0
    # amp=4
    h.tstop = T[-1]+20
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return h.tstop, 0, T[-1], 0


class Analyzer():

    def __init__(self, cell=None, parts_dict=None, colors_dict=None, type='input_cell'):
        if cell is None:
            if type == 'Rall_tree':
                cell, parts_dict, colors_dict = self.open_rall_tree()
        assert parts_dict is not None
        assert colors_dict is not None
        self.cell=cell
        self.parts_dict = parts_dict
        self.more_conductances = more_conductances_fake(cell)
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    def open_morph(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, nl=None):
        hoc_file_name = 'allen_model.hoc'
        h.celsius = 37
        # Create the model
        # h.load_file(model_path + hoc_file_name)
        # h.load_file("import3d.hoc")
        h.load_file("nrngui.hoc")
        h("objref cell, tobj")  # neuron object
        h.load_file('allen_model.hoc')

        h.execute("cell = new " + hoc_file_name[:-4] + "()")  # replace?
        # nl = h.Import3d_SWC_read()
        nl.input(morph_path)
        i3d = h.Import3d_GUI(nl, 0)
        i3d.instantiate(h.cell)
        cell = h.cell
        parts_dict = dict(all=list())
        colors_dict = {'all': 'k'}
        for sec in cell.all:
            sec.insert('pas')
            sec.nseg = int(sec.L / 10) + 1
            sec.e_pas = e_pas
            sec.cm = Cm
            sec.Ra = Ra
            sec.g_pas = 1.0 / Rm
            for seg in sec:
                parts_dict['all'].append(seg)
        return cell, parts_dict, colors_dict

    def open_ASC(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas = -70):
        h.load_file("import3d.hoc")
        return self.open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas = e_pas, nl=h.Import3d_Neurolucida3())

    def open_swc(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70):
        h.load_file("import3d.hoc")
        return self.open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas = e_pas, nl=h.Import3d_SWC_read())

    def open_rall_tree(self):
        morph_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/Rall_tree5.swc')
        cell, parts_dict, colors_dict = self.open_swc(morph_path)
        return cell, dict(Rall_tree = parts_dict['all']), dict(Rall_tree =colors_dict['all'])

    def change_color_dict(self, colors_dict):
        for part in self.parts_dict:
            assert part in colors_dict

        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=self.parts_dict, color_dict=colors_dict)

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
        fig, ax, color_bar, points_dict, lines, segs = plot_morph(self.cell, color_func=colors.get_seg_color, scatter=False, add_nums=False,
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

    def plot_morph_with_value_func(self, func=seg_Rin_func, run_time=0, theta=0.0):
        if run_time>0:
            h.tstop=run_time
            h.run()
        value_dict = dict()
        for part in self.parts_dict:
            for seg in self.parts_dict[part]:
                value_dict[seg] = func(seg)
        return self.plot_morph_with_values(value_dict, theta=theta)[:2]

    def record_protocol(self, protocol=spike_protocol, cut_start_ms=0, record_name='v'):
        record_dict = dict()
        for sec in self.cell.all:
            record_dict[sec] = dict()
            for seg in sec:
                try:
                    record_dict[sec][seg.x] = h.Vector()
                    record_dict[sec][seg.x].record(getattr(sec(seg.x), '_ref_' + record_name))
                except:
                    record_dict[sec][seg.x] = None
        time = h.Vector()
        time.record(h._ref_t)
        tstop, delay, dur, amp = protocol(self.cell, self.cell.soma[0](0.5))
        total_len_idx = int((tstop - cut_start_ms) / h.dt)
        for sec in self.cell.all:
            for seg in sec:
                if record_dict[sec][seg.x] is None:
                    record_dict[sec][seg.x] = np.zeros(total_len_idx)
                else:
                    record_dict[sec][seg.x] = np.array(record_dict[sec][seg.x])[int(cut_start_ms / h.dt):]
        time = np.array(time)[int(cut_start_ms / h.dt):]
        time -= time[0]
        return record_dict, time


    def create_movie_from_rec2(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                            plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset = 'ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean): # this is a good option of you have small number of segments

        mpl.use('TkAgg')
        import matplotlib.style as mplstyle
        mplstyle.use('fast')
        time /= 1000.0
        time *= slow_down_factor
        min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = min_value+10#################remove this part#################
        ax = plt.gca()
        fig = ax.get_figure()
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
        segs = np.array(segs)
        lines = np.array(lines)
        lim_x = ax.get_xlim()
        lim_y = ax.get_ylim()
        time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
        self.last_t = 0
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        norm = get_norm([min_value, max_value])
        fig1, ax1 = plt.subplots()
        cax = fig1.add_axes([0.90, 0.2, 0.02, 0.6])
        cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='uniform')

        if not scale == 0:
            ax1.plot([xlim[0], xlim[0]], [ylim[0], ylim[0] + scale], color='k')

        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.grid(False)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_axis_off()
        cb_fig = mplfig_to_npimage(fig1)
        cb_mask = ~(cb_fig.sum(axis=2)==765)

        lines_data = []
        seg_list = np.unique(segs)
        print(len(seg_list ), len(lines))
        for seg in tqdm(seg_list, desc='optimizing lines'):
            fig1, ax1 = plt.subplots()
            cax = fig1.add_axes([0.90, 0.2, 0.02, 0.6])
            for line in lines[segs==seg]:
                ax1.plot(line.get_xdata(), line.get_ydata(), lw=line.get_linewidth())
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            for a in [ax1, cax]:
                a.grid(False)
                a.set_xticks([])
                a.set_yticks([])
                a.set_axis_off()
            np_fig = mplfig_to_npimage(fig1)
            mask = ~(np_fig.sum(axis=2)==765)
            lines_data.append(dict(mask=mask, seg=seg))
        np_shape = np_fig.shape
        plt.close()

        def make_frame(t):  # get the correct index
            time_loc = np.where(time >= t)[0][0]
            prev_time_loc = np.where(time >= self.last_t)[0][0]
            self.last_t = t

            base = np.zeros(np_shape)+255
            for line_data in lines_data:
                seg = line_data['seg']
                mask = line_data['mask']
                if prev_time_loc == time_loc:
                    color = cmap(norm(record_dict[seg.sec][seg.x][time_loc]))[:3]
                else:
                    color = cmap(norm(func_for_missing_frames(record_dict[seg.sec][seg.x][prev_time_loc:time_loc])))[:3]
                base[mask] = (np.array(color) * 255).astype(int)

            base[cb_mask] = cb_fig[cb_mask]
            return base

        animation = VideoClip(make_frame, duration=time[-1])
        animation.write_videofile(os.path.join(save_to, clip_name + '.mp4'),
                                  fps=int(1000.0 / h.dt) if fps is None else fps / slow_down_factor, threads=threads,
                                  audio=False, preset=preset)
        mpl.use('Qt5Agg')
        # self.last_t = 0
        # animation.write_gif(os.path.join(save_to, clip_name + '.gif'), fps=int(1.0 / h.dt) if fps is None else fps/slow_down_factor)



    def create_movie_from_rec(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                            plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset = 'ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean):

        mpl.use('TkAgg')
        import matplotlib.style as mplstyle
        mplstyle.use('fast')
        time /= 1000.0
        time *= slow_down_factor
        min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = min_value+10#################remove this part#################
        ax = plt.gca()
        fig = ax.get_figure()
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
        segs = np.array(segs)
        lines = np.array(lines)
        lim_x = ax.get_xlim()
        lim_y = ax.get_ylim()
        time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
        self.last_t = 0

        def make_frame(t):  # get the correct index
            time_loc = np.where(time >= t)[0][0]
            prev_time_loc = np.where(time >= self.last_t)[0][0]
            self.last_t = t
            norm = get_norm([min_value, max_value])
            for line, seg in zip(lines, segs):
                if prev_time_loc == time_loc:
                    line.set_color(cmap(norm(record_dict[seg.sec][seg.x][time_loc])))
                else:
                    line.set_color(
                        cmap(norm(func_for_missing_frames(record_dict[seg.sec][seg.x][prev_time_loc:time_loc]))))
            time_text.set_text('time: ' + str(round(t / slow_down_factor * 1000, 1)) + ' (ms)')
            a=mplfig_to_npimage(fig)
            return mplfig_to_npimage(fig)

        animation = VideoClip(make_frame, duration=time[-1])
        animation.write_videofile(os.path.join(save_to, clip_name + '.mp4'),
                                  fps=int(1000.0 / h.dt) if fps is None else fps / slow_down_factor, threads=threads,
                                  audio=False, preset=preset)
        mpl.use('Qt5Agg')

    def create_morph_movie(self, protocol=spike_protocol, cut_start_ms=0, record_name='v',
                            seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                            plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset = 'ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean):

        record_dict, time = self.record_protocol(protocol=protocol, cut_start_ms=cut_start_ms, record_name=record_name)

        self.create_movie_from_rec(record_dict=record_dict, time=time, fps=fps, clip_name=clip_name,
                                       threads=threads, slow_down_factor=slow_down_factor, func_for_missing_frames=func_for_missing_frames, theta=theta)



    # def plot_morph_with_value_func_after_protocol(self, func=seg_Rin_func, run_time=0):
    #     pass

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

    def plot_attanuation(self, protocol=long_pulse_protocol, ax=None, seg_to_indicate=[], indication_color='orange', initial_seg =None, record_to_value_func=None, norm=True, norm_by=1.0, **kwargs):
        if ax is None:
            ax = plt.gca()
        if initial_seg is None:
            initial_seg = self.cell.soma[0](0.5)

        att = attenuation(self.cell, color_func=self.colors, seg_length_function=get_segment_length_lamda,
                          more_conductances=self.more_conductances, param_to_record='v',
                          record_to_value_func=record_to_value_func)
        tstop, delay, dur, amp = protocol(self.cell, initial_seg)
        ax, norm_by = att.plot(start_seg=initial_seg, norm=norm, norm_by=norm_by, cut_start=int((delay - 1) / h.dt),
                      seg_to_indicate={seg:dict(size=30, color=indication_color, alpha=1) for seg in seg_to_indicate}, ax=ax, **kwargs)
        ax.set_yscale('log')
        return ax, norm_by

    def create_card(self):
        pass
