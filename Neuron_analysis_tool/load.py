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




def seg_Rin_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

def long_pulse_protocol(cell, start_seg):
    delay=2000.0
    dur=1000.0
    amp=.1
    h.tstop = delay+dur+500.0
    clamp = h.IClamp(start_seg.x, sec = start_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return h.tstop, delay, dur, amp, {}

def short_pulse_protocol(cell, start_seg):
    delay=2000.0
    dur=2.0
    amp=0.1
    h.tstop = delay+dur+20.0
    clamp = h.IClamp(start_seg.x, sec = start_seg.sec)
    clamp.delay = delay
    clamp.dur = dur
    clamp.amp = amp
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return h.tstop, delay, dur, amp, {}

def spike_protocol(cell, start_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt=spike_data.T[0][1]-spike_data.T[0][0]
    V = np.concatenate([np.zeros(int(1000.0/dt))+spike_data.T[1][0]]+[spike_data.T[1]]*10)
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)
    h.dt=dt
    h.steps_per_ms = 1.0 / h.dt
    clamp = h.SEClamp(start_seg.x, sec=start_seg.sec)
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
    return h.tstop, 0, T[-1], 0, {}


class Analyzer():
    def __init__(self, cell=None, parts_dict=None, colors_dict=None, type='input_cell', morph_path = None, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, conductances_list=[], run_time=500):
        if cell is None:
            if type == 'Rall_tree':
                cell, parts_dict, colors_dict = self.open_rall_tree()
            elif type == 'ASC' and morph_path:
                cell, parts_dict, colors_dict = self.open_ASC(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas)
            elif type == 'swc' and morph_path:
                cell, parts_dict, colors_dict = self.open_swc(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas)
            elif type == 'L5PC':
                cell, parts_dict, colors_dict = self.open_L5PC()

        assert parts_dict is not None
        assert colors_dict is not None
        self.type=type
        self.cell=cell
        self.parts_dict = parts_dict
        self.more_conductances = more_conductances(cell, run_time=run_time, record_names=conductances_list, is_resting=len(conductances_list)==0)
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)

    def open_morph(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, nl=None):
        print('open_morph: ', morph_path)
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
        # parts_dict = dict(all=list())

        parts_dict = {'soma': [], 'basal': [], 'apical': [], 'axon': [], 'else': []}
        colors_dict = {'soma': 'k', 'basal': 'r', 'apical': 'b', 'axon': 'green', 'else': 'cyan'}
        for sec in cell.all:
            sec.insert('pas')
            sec.nseg = int(sec.L / 20) + 1
            sec.e_pas = e_pas
            sec.cm = Cm
            sec.Ra = Ra
            sec.g_pas = 1.0 / Rm
            for seg in sec:
                if sec in cell.soma:
                    parts_dict['soma'].append(seg)
                elif sec in list(cell.basal):
                    parts_dict['basal'].append(seg)
                elif sec in list(cell.apical):
                    parts_dict['apical'].append(seg)
                elif sec in list(cell.axonal):
                    parts_dict['axon'].append(seg)
                else:
                    parts_dict['else'].append(seg)
        return cell, parts_dict, colors_dict

    def open_ASC(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas = -70):
        h.load_file("import3d.hoc")
        return self.open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas = e_pas, nl=h.Import3d_Neurolucida3())

    def open_swc(self, morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70):
        h.load_file("import3d.hoc")
        return self.open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas = e_pas, nl=h.Import3d_SWC_read())

    def open_rall_tree(self):
        morph_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/Rall_tree5.swc')
        # cell, parts_dict, colors_dict = self.open_swc(morph_path)
        return self.open_swc(morph_path) #cell, dict(Rall_tree = parts_dict['all']), dict(Rall_tree =colors_dict['all'])

    def open_L5PC(self):
        h.load_file("import3d.hoc")
        morphology_file = "data/L5PC/cell1.asc"
        h.load_file("data/L5PC/L5PCbiophys3.hoc")
        h.load_file("data/L5PC/L5PCtemplate.hoc")
        cell = h.L5PCtemplate(morphology_file)

        parts_dict = {'soma': [], 'basal': [], 'apical': [], 'axon': [], 'else': []}
        colors_dict = {'soma': 'k', 'basal': 'r', 'apical': 'b', 'axon': 'green', 'else': 'cyan'}
        for sec in cell.all:
            sec.nseg = int(sec.L / 20) + 1
            for seg in sec:
                if sec in cell.soma:
                    parts_dict['soma'].append(seg)
                elif sec in cell.dend:
                    parts_dict['basal'].append(seg)
                elif sec in cell.apic:
                    parts_dict['apical'].append(seg)
                elif sec in cell.axon:
                    parts_dict['axon'].append(seg)
                else:
                    parts_dict['else'].append(seg)
        return cell, parts_dict, colors_dict

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


    def plot_morph(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None, sec_to_change=None, ignore_sections=[], theta=0, scale=0, ignore_soma=True):
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        fig = plt.gca().figure
        fig, ax, color_bar, points_dict, lines, segs = plot_morph(self.cell, color_func=self.colors.get_seg_color, scatter=False, add_nums=False,
                                                     seg_to_indicate=seg_to_indicate_dict,
                                                     norm_colors=False, fig=fig, ax=ax, diam_factor=diam_factor,
                                                     sec_to_change=sec_to_change, bounds=None, cmap=plt.cm.turbo,
                                                     plot_color_bar=False,
                                                     theta=theta, ignore_sections=ignore_sections, ignore_soma=ignore_soma)
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
        return ax

    def plot_morph_with_values(self, seg_val_dict, ax=None, seg_to_indicate_dict = {}, diam_factor=None,
                               sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None, ignore_soma=True, color_bar_idx = [0.9, 0.2, 0.02, 0.6]):
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        ax.set_aspect('equal', adjustable='box')
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
                                                     theta=theta, ignore_sections=ignore_sections, ignore_soma=ignore_soma, color_bar_idx=color_bar_idx)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')

            ax.text(x_lim[0]-20, y_lim[0], str(scale)+' (um)', rotation=90, fontsize=(y_lim[1]-y_lim[0])*0.01)
            ax.set_xlim([x_lim[0] - 20, x_lim[1]])
        return ax, color_bar, points_dict, lines, segs


    def plot_morph_with_values2(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None,
                               scale=0, cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None, points_dict=None):
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        norm=get_norm(bounds)
        ax.set_aspect('equal', adjustable='box')
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
            ax.text(x_lim[0]-20, y_lim[0], str(scale)+' (um)', rotation=90)
            ax.set_xlim([x_lim[0] - 20, x_lim[1]])
        return ax, color_bar

    def plot_morph_with_value_func(self, func=seg_Rin_func, run_time=0, theta=0.0, diam_factor=None, cmap=plt.cm.turbo,
                                   ax = None, seg_to_indicate_dict = {},
                                   sec_to_change = None, ignore_sections = [], scale = 0, plot_color_bar = True, bounds = None):
        if run_time>0:
            h.tstop=run_time
            h.run()
        value_dict = dict()
        for part in self.parts_dict:
            for seg in self.parts_dict[part]:
                value_dict[seg] = func(seg)
        return self.plot_morph_with_values(value_dict, theta=theta, diam_factor=diam_factor, cmap=cmap, ax=ax, seg_to_indicate_dict=seg_to_indicate_dict,
                                            sec_to_change=sec_to_change, ignore_sections=ignore_sections, scale=scale, plot_color_bar=plot_color_bar, bounds=bounds)[:2]

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
        tstop, delay, dur, amp, extra = protocol(self.cell, self.cell.soma[0](0.5))

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
    #
    # def create_movie_from_rec2(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
    #                         sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
    #                         plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset = 'ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean): # this is a good option of you have small number of segments
    #
    #     mpl.use('TkAgg')
    #     import matplotlib.style as mplstyle
    #     mplstyle.use('fast')
    #     time /= 1000.0
    #     time *= slow_down_factor
    #     min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
    #     max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
    #     # max_value = min_value+10#################remove this part#################
    #     ax = plt.gca()
    #     fig = ax.get_figure()
    #     value_dict = dict()
    #     for sec in self.cell.all:
    #         for seg in sec:
    #             value_dict[seg] = record_dict[sec][seg.x][0]
    #
    #     ax, color_bar, points_dict, lines, segs = self.plot_morph_with_values(value_dict, ax=ax,
    #                                                                           seg_to_indicate_dict=seg_to_indicate_dict,
    #                                                                           diam_factor=diam_factor,
    #                                                                           sec_to_change=sec_to_change,
    #                                                                           ignore_sections=ignore_sections,
    #                                                                           theta=theta, scale=scale, cmap=cmap,
    #                                                                           plot_color_bar=plot_color_bar,
    #                                                                           bounds=[min_value, max_value])
    #     segs = np.array(segs)
    #     lines = np.array(lines)
    #     lim_x = ax.get_xlim()
    #     lim_y = ax.get_ylim()
    #     time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
    #     self.last_t = 0
    #     xlim = ax.get_xlim()
    #     ylim = ax.get_ylim()
    #
    #     norm = get_norm([min_value, max_value])
    #     fig1, ax1 = plt.subplots()
    #     cax = fig1.add_axes([0.90, 0.2, 0.02, 0.6])
    #     cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='uniform')
    #
    #     if not scale == 0:
    #         ax1.plot([xlim[0], xlim[0]], [ylim[0], ylim[0] + scale], color='k')
    #
    #     ax1.set_xlim(xlim)
    #     ax1.set_ylim(ylim)
    #     ax1.grid(False)
    #     ax1.set_xticks([])
    #     ax1.set_yticks([])
    #     ax1.set_axis_off()
    #     cb_fig = mplfig_to_npimage(fig1)
    #     cb_mask = ~(cb_fig.sum(axis=2)==765)
    #
    #     lines_data = []
    #     seg_list = np.unique(segs)
    #     print(len(seg_list ), len(lines))
    #     for seg in tqdm(seg_list, desc='optimizing lines'):
    #         fig1, ax1 = plt.subplots()
    #         cax = fig1.add_axes([0.90, 0.2, 0.02, 0.6])
    #         for line in lines[segs==seg]:
    #             ax1.plot(line.get_xdata(), line.get_ydata(), lw=line.get_linewidth())
    #         ax1.set_xlim(xlim)
    #         ax1.set_ylim(ylim)
    #         for a in [ax1, cax]:
    #             a.grid(False)
    #             a.set_xticks([])
    #             a.set_yticks([])
    #             a.set_axis_off()
    #         np_fig = mplfig_to_npimage(fig1)
    #         mask = ~(np_fig.sum(axis=2)==765)
    #         lines_data.append(dict(mask=mask, seg=seg))
    #     np_shape = np_fig.shape
    #     plt.close()
    #
    #     def make_frame(t):  # get the correct index
    #         time_loc = np.where(time >= t)[0][0]
    #         prev_time_loc = np.where(time >= self.last_t)[0][0]
    #         self.last_t = t
    #
    #         base = np.zeros(np_shape)+255
    #         for line_data in lines_data:
    #             seg = line_data['seg']
    #             mask = line_data['mask']
    #             if prev_time_loc == time_loc:
    #                 color = cmap(norm(record_dict[seg.sec][seg.x][time_loc]))[:3]
    #             else:
    #                 color = cmap(norm(func_for_missing_frames(record_dict[seg.sec][seg.x][prev_time_loc:time_loc])))[:3]
    #             base[mask] = (np.array(color) * 255).astype(int)
    #
    #         base[cb_mask] = cb_fig[cb_mask]
    #         return base
    #
    #     animation = VideoClip(make_frame, duration=time[-1])
    #     animation.write_videofile(os.path.join(save_to, clip_name + '.mp4'),
    #                               fps=int(1000.0 / h.dt) if fps is None else fps / slow_down_factor, threads=threads,
    #                               audio=False, preset=preset)
    #     mpl.use('Qt5Agg')
    #     # self.last_t = 0
    #     # animation.write_gif(os.path.join(save_to, clip_name + '.gif'), fps=int(1.0 / h.dt) if fps is None else fps/slow_down_factor)

    def save_movie_from_rec(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                              sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                              plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4,
                              preset='ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean):

        animation = self.create_movie_from_rec(record_dict, time, seg_to_indicate_dict, diam_factor,
                              sec_to_change, ignore_sections, theta, scale, cmap,
                              plot_color_bar, clip_name, fps, threads,
                              preset, slow_down_factor, func_for_missing_frames)

        animation.write_videofile(os.path.join(save_to, clip_name + '.mp4'),
                              fps=int(1000.0 / h.dt) if fps is None else fps / slow_down_factor, threads=threads,
                              audio=False, preset=preset)


    def create_movie_from_rec(self, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                            plot_color_bar=True, slow_down_factor=1, func_for_missing_frames=np.mean, bounds = None,
                            show_records_from=dict(), voltage_window=50, ylabel='v (mV)', xlabel='time (ms)'):

        import matplotlib.style as mplstyle
        mplstyle.use('fast')
        time /= 1000.0
        time *= slow_down_factor
        min_value = np.min([np.min(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        max_value = np.max([np.max(record_dict[sec][pos]) for sec in record_dict for pos in record_dict[sec]])
        voltage_segs = list(show_records_from.keys())
        if len(show_records_from)==0:
            ax = plt.gca()
            fig = ax.get_figure()
            color_bar_idx = [0.8, 0.2, 0.02, 0.6]
        else:
            from matplotlib.gridspec import GridSpec
            fig = plt.figure(constrained_layout=True)
            plt.subplots_adjust(wspace=.75, hspace=0.35)
            gs = GridSpec(2, 2, figure=fig)
            ax = fig.add_subplot(gs[:, 0])
            # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
            ax2 = fig.add_subplot(gs[1, 1])
            ax3 = fig.add_subplot(gs[0, 1])
            for a in [ax2, ax3]:
                a.spines['top'].set_visible(False)
                a.spines['right'].set_visible(False)
            if len(show_records_from)==1:
                ax3.set_axis_off()


            color_bar_idx=[0.4, 0.2, 0.02, 0.6]
        value_dict = dict()
        for sec in self.cell.all:
            for seg in sec:
                value_dict[seg] = record_dict[sec][seg.x][0]
        seg_to_indicate_dict.update(show_records_from)
        ax, color_bar, points_dict, lines, segs = self.plot_morph_with_values(value_dict, ax=ax,
                                                                              seg_to_indicate_dict = seg_to_indicate_dict,
                                                                              diam_factor=diam_factor,
                                                                              sec_to_change=sec_to_change,
                                                                              ignore_sections=ignore_sections,
                                                                              theta=theta, scale=scale, cmap=cmap,
                                                                              plot_color_bar=plot_color_bar,
                                                                              bounds=[min_value, max_value], color_bar_idx = color_bar_idx)
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
            if bounds:
                norm = get_norm(bounds)
            else:
                norm = get_norm([min_value, max_value])
            for line, seg in zip(lines, segs):
                if prev_time_loc == time_loc:
                    line.set_color(cmap(norm(record_dict[seg.sec][seg.x][time_loc])))
                else:
                    line.set_color(
                        cmap(norm(func_for_missing_frames(record_dict[seg.sec][seg.x][prev_time_loc:time_loc]))))
            time_in_ms = t / slow_down_factor * 1000
            time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
            if len(show_records_from) > 0:
                try: start_idx = np.where(time >= t-(voltage_window)/1000*slow_down_factor/2)[0][0]
                except: start_idx =0
                try:end_idx = np.where(time >= t+voltage_window/1000*slow_down_factor/2)[0][0]
                except: end_idx = -1
                v1 = record_dict[voltage_segs[0].sec][voltage_segs[0].x][start_idx:end_idx]
                t1 = time[start_idx:end_idx]*1000.0/slow_down_factor
                ax2.clear()
                ax2.plot(t1,v1, color=show_records_from[voltage_segs[0]]['color'])
                ax2.set_ylim(color_bar.vmin, color_bar.vmax)
                ax2.set_ylabel(ylabel)
                ax2.set_xlabel(xlabel)
                ax2.set_title(show_records_from[voltage_segs[0]]['label'])
                ax2.axvline(t*1000.0/slow_down_factor, color='r', ls='--')
                ax2.set_xlim(xmin=time_in_ms-voltage_window/2, xmax=time_in_ms+voltage_window/2)
                if len(show_records_from) > 1:
                    ax3.clear()
                    v2 = record_dict[voltage_segs[1].sec][voltage_segs[1].x][start_idx:end_idx]
                    ax3.plot(t1, v2, color=show_records_from[voltage_segs[1]]['color'])
                    ax3.set_ylim(color_bar.vmin, color_bar.vmax)
                    ax3.axvline(t*1000.0/slow_down_factor, color='r', ls='--')
                    ax3.set_title(show_records_from[voltage_segs[1]]['label'])
                    ax3.set_ylabel(ylabel)
                    ax3.set_xticks([])
                    ax3.set_xlim(xmin=time_in_ms-voltage_window/2, xmax=time_in_ms+voltage_window/2)

            return mplfig_to_npimage(fig)
        animation = VideoClip(make_frame, duration=time[-1])
        # for a in [ax, ax2, ax3]:
        #     a.clear()
        return animation

    def create_morph_movie(self, protocol=spike_protocol, cut_start_ms=0, record_name='v',
                            seg_to_indicate_dict=dict(), diam_factor=None,
                            sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                            plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset = 'ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean):

        record_dict, time = self.record_protocol(protocol=protocol, cut_start_ms=cut_start_ms, record_name=record_name)

        self.create_movie_from_rec(record_dict=record_dict, time=time, slow_down_factor=slow_down_factor, func_for_missing_frames=func_for_missing_frames, theta=theta,
                                   scale=scale, cmap=cmap, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,sec_to_change=sec_to_change, ignore_sections=ignore_sections,
                                   plot_color_bar=plot_color_bar)



    # def plot_morph_with_value_func_after_protocol(self, func=seg_Rin_func, run_time=0):
    #     pass

    def plot_dendogram(self, start_seg = None ,ax=None, dots_loc_seg = [], plot_legend=True, ignore_sections=[], electrical=False):
        if ax is None:
            ax = plt.gca()
        if start_seg is None:
            start_seg = self.cell.soma[0](0.5)
        dendogram = Dendogram(self.cell,
                               seg_length_function=get_segment_length_lamda if electrical else get_segment_length_um,
                               color_func=self.colors,
                               dots_loc=[[seg.sec.name().split('.')[-1], seg.x] for seg in dots_loc_seg],
                               more_conductances=self.more_conductances,
                               diam_factor=None, s=10, fix_diam=1.)
        dendogram.cumpute_distances(start_seg)
        max_y, x_pos = dendogram.plot(ax=ax, plot_legend=plot_legend, ignore_sections=ignore_sections)
        return ax, x_pos

    def plot_cable(self, start_seg = None ,ax=None,
                   factor_e_space=25, factor_m_space=10 ,
                   dots_loc_seg = [], ignore_sections=[],
                   cable_type='electric', start_loc=0, x_axis=True,
                    factor=1, dots_size=10, start_color='k', plot_legend=True): #'d3_2', 'dist'
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
        results, seg_dist, cross_dist_dict = get_cable( self.cell,
                                                        factor_e_space=factor_e_space,
                                                        factor_m_space=factor_m_space,
                                                        start_section=start_seg.sec,
                                                        x_start=start_seg.x,
                                                        more_conductions=self.more_conductances,
                                                        seg_dist_dict=seg_dist_dict,
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

    def plot_attanuation(self, protocol=long_pulse_protocol, ax=None, seg_to_indicate_dict=dict(), start_seg =None, record_to_value_func=None, norm=True, norm_by=1.0, **kwargs):
        if ax is None:
            ax = plt.gca()
        if start_seg is None:
            start_seg = self.cell.soma[0](0.5)

        att = attenuation(self.cell, color_func=self.colors, seg_length_function=get_segment_length_lamda,
                          more_conductances=self.more_conductances, param_to_record='v',
                          record_to_value_func=record_to_value_func)
        tstop, delay, dur, amp, extra = protocol(self.cell, start_seg)
        ax, norm_by = att.plot(start_seg=start_seg, norm=norm, norm_by=norm_by, cut_start=int((delay - 1) / h.dt), seg_to_indicate=seg_to_indicate_dict, ax=ax, **kwargs)
        ax.set_yscale('log')
        # ax.get_yaxis().get_major_formatter().labelOnlyBase = True
        # ax.yaxis.set_major_locator(plt.MaxNLocator(number_of_lacations))
        ax.set_xlabel('distance from origin (x / '+LAMDA+')')
        if norm or (not norm_by==1.0):
            ax.set_ylabel('V(x)/V(0)')
        else:
            ax.set_ylabel('attanuation (mV)')

        return ax, norm_by

    def create_card(self, start_seg=None, theta=0, scale=500, factor_e_space=25, cable_type='d3_2', diam_factor=None, plot_legend=False, start_color='green', start_dots_size=50, **kwargs):
        fig, ax = plt.subplots(1, 4, figsize=(12, 3), gridspec_kw={'width_ratios': [0.5, 1.5, 1, 1]})
        plt.subplots_adjust(wspace=0.5)
        if start_seg is None:
            start_seg = self.cell.soma[0]
            start_seg = list(start_seg)
            start_seg = start_seg[len(start_seg)//2]
        seg_to_indicate_dict = {start_seg: dict(size=start_dots_size, color=start_color, alpha=1)}
        self.plot_morph(ax=ax[0], theta=theta, seg_to_indicate_dict=seg_to_indicate_dict, scale=scale, diam_factor=diam_factor, ignore_soma=not self.type.startswith('Rall_tree'))
        _, x_pos = self.plot_dendogram(start_seg=start_seg, ax=ax[1], electrical=True, plot_legend=False)
        self.plot_cable(start_seg=start_seg, ax=ax[1], factor_e_space=factor_e_space, cable_type=cable_type, plot_legend=plot_legend, start_loc=x_pos+10, start_color=start_color, dots_size=start_dots_size)
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

        self.plot_attanuation(start_seg=start_seg, ax=ax[2], protocol=long_pulse_protocol, seg_to_indicate_dict=seg_to_indicate_dict, ** kwargs)
        self.plot_attanuation(start_seg=start_seg, ax=ax[3], protocol=short_pulse_protocol, seg_to_indicate_dict=seg_to_indicate_dict, ** kwargs)
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
