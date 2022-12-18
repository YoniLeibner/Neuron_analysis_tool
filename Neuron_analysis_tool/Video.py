from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake
from Neuron_analysis_tool.morph_ploter import get_norm
from Neuron_analysis_tool.record import record, record_all
from Neuron_analysis_tool.distance import Distance
from Neuron_analysis_tool.protocols import *
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import os


def save_movie_from_rec(analyzer, record_dict, time, seg_to_indicate_dict=dict(), diam_factor=None,
                        sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                        plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4,
                        preset='ultrafast', slow_down_factor=1, func_for_missing_frames=np.mean):
    animation = create_movie_from_rec(analyzer, record_dict, time, seg_to_indicate_dict, diam_factor,
                                           sec_to_change, ignore_sections, theta, scale, cmap,
                                           plot_color_bar, clip_name, fps, threads,
                                           preset, slow_down_factor, func_for_missing_frames)

    animation.write_videofile(os.path.join(save_to, clip_name + '.mp4'),
                              fps=int(1000.0 / h.dt) if fps is None else fps / slow_down_factor, threads=threads,
                              audio=False, preset=preset)


def create_movie_from_rec(analyzer, records, seg_to_indicate_dict=dict(), diam_factor=None,
                          sec_to_change=None, ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo,
                          plot_color_bar=True, slow_down_factor=1, func_for_missing_frames=np.max, bounds=None,
                          show_records_from=dict(), voltage_window=50, ylabel='v (mV)', xlabel='time (ms)', margin=0,
                          draw_funcs=[],
                          base_plot_type='morph', start_seg=None, electrical=True, figsize=(5, 5),
                          dancing=False,distance=None, dt=1, more_conductances_=None):
    assert dancing == False or more_conductances_ is not None
    if base_plot_type='morph':
        assert dancing == electrical
    import matplotlib.style as mplstyle
    if start_seg is None:
        start_seg = list(analyzer.cell.soma[0])
        start_seg = start_seg[len(start_seg) // 2]
    mplstyle.use('fast')
    time = records.time.copy()
    time /= 1000.0
    time *= slow_down_factor
    min_value = records.get_min() - margin
    max_value = records.get_max() + margin
    if bounds is not None:
        min_value = bounds[0]
        max_value = bounds[1]

    voltage_segs = list(show_records_from.keys())
    if len(show_records_from) == 0:
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()
        color_bar_idx = [0.8, 0.2, 0.02, 0.6]
    else:
        from matplotlib.gridspec import GridSpec
        fig = plt.figure(constrained_layout=True, figsize=figsize)
        plt.subplots_adjust(wspace=.75, hspace=0.35)
        gs = GridSpec(len(show_records_from), 2, figure=fig)
        ax = fig.add_subplot(gs[:, 0])
        ax_v = []
        for i in range(len(show_records_from)):
            ax_v.append(fig.add_subplot(gs[i, 1]))
        # ax3 = fig.add_subplot(gs[0, 1])
        # for a in [ax2, ax3]:
            ax_v[-1].spines['top'].set_visible(False)
            ax_v[-1].spines['right'].set_visible(False)
        # if len(show_records_from) == 1:
        #     ax3.set_axis_off()

        color_bar_idx = [0.4, 0.2, 0.02, 0.6]
    value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
    seg_to_indicate_dict.update(show_records_from)
    if base_plot_type == 'morph':
        ax, cax, colors, lines, segs = analyzer.plot_morph_with_values(value_dict_by_sec, ax=ax,
                                                                   seg_to_indicate_dict=seg_to_indicate_dict,
                                                                   diam_factor=diam_factor,
                                                                   sec_to_change=sec_to_change,
                                                                   ignore_sections=ignore_sections,
                                                                   theta=theta, scale=scale,
                                                                   cmap=cmap, bounds=[min_value, max_value],
                                                                   plot_color_bar=plot_color_bar,
                                                                   color_bar_idx=color_bar_idx,
                                                                   distance=distance,
                                                                   electrical=electrical,
                                                                   time=0.5,
                                                                   dt=dt,
                                                                   more_conductances_= more_conductances_
                                                                   )
    elif base_plot_type == 'dendogram':
        ax, x_pos, cax, colors, lines, segs = analyzer.plot_dendogram_with_values(value_dict_by_sec, start_seg=start_seg,
                                                                              ax=ax,
                                                                              segs_to_indecate=seg_to_indicate_dict,
                                                                              plot_legend=False,
                                                                              ignore_sections=ignore_sections,
                                                                              electrical=electrical,
                                                                              diam_factor=diam_factor, distance=None,
                                                                              bounds=[min_value, max_value], cmap=cmap,
                                                                              plot_color_bar=plot_color_bar,
                                                                              color_bar_idx=color_bar_idx, colors=None)
    # elif base_plot_type == 'attenuation':
    #     pass

    else:
        raise Exception('base plot_type:', +str(base_plot_type) + ' not implementes')
    segs = np.array(segs)
    lines = np.array(lines)
    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    if dancing:
        x_range = abs(lim_x[1]-lim_x[0]) * 0.2
        y_range = abs(lim_y[1]-lim_y[0]) * 0.2
        lim_x = [lim_x[0]-x_range, lim_x[1]+x_range]
        lim_y = [lim_y[0]-y_range, lim_y[1]+y_range]
        ax.set_xlim(lim_x)
        ax.set_ylim(lim_y)
    time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
    analyzer.last_t = 0
    analyzer.to_remove = []

    def make_frame_dendorgam_r(distance, sec, start_point):
        sec_segs = list(sec)
        flag = False
        for seg in sec_segs:
            if flag: break
            for line in lines[segs==seg]:
                y = line.get_ydata()
                if not flag and y[1] < y[0]:
                    sec_segs=sec_segs[::-1]
                    flag = True
        for seg in sec_segs:
            for idx, seg2 in enumerate(segs):
                if seg2 == seg:
                    y = lines[idx].get_ydata()
                    if not np.all(y==y[0]):
                        y_prime = y.copy()
                        # fixing the shift
                        y += start_point['y'] - y[0]
                        dy = y[1] - y[0]
                        current_len = abs(dy) # the x is the same in dendogram
                        wanted_len = distance.get_length(seg, electrical=True)
                        f = wanted_len / current_len

                        # fixing the length
                        y[1] *= f
                        # print(y, y_prime, start_point['y'], f, sec, seg)

                        lines[idx].set_ydata(y)
                        start_point = dict(y=y[-1])
                        all_lines = lines[segs==seg]
                        if len(all_lines)>1:
                            all_lines[np.logical_not(all_lines==lines[idx])][0].set_ydata([y[-1]]*2)
        for son in sec.children():
            make_frame_dendorgam_r(distance, son, start_point)

    def make_frame_r(distance, sec, start_point):
        for seg in sec:
            for idx, seg2 in enumerate(segs):
                if seg2 == seg:
                    x, y = lines[idx].get_data()
                    # check if line is horizental
                    dx = x[:-1] - x[1:]
                    dy = y[:-1] - y[1:]
                    current_lens = (dx ** 2 + dy ** 2) ** 0.5
                    wanted_len = distance.get_length(seg, electrical=True)
                    f = wanted_len / sum(current_lens)
                    # fixing the shift
                    x += start_point['x'] - x[0]
                    y += start_point['y'] - y[0]

                    # fixing the length
                    for point_idx in range(1, len(x), 1):
                        new_x = x[point_idx - 1] + (x[point_idx] - x[point_idx - 1]) * f
                        new_y = y[point_idx - 1] + (y[point_idx] - y[point_idx - 1]) * f
                        x[point_idx:] += (new_x - x[point_idx])
                        y[point_idx:] += (new_y - y[point_idx])
                    lines[idx].set_data(x, y)
                    start_point = dict(x=x[-1], y=y[-1])
        for son in sec.children():
            make_frame_r(distance, son, start_point)


    def make_frame(t):
        for elements_to_remove in analyzer.to_remove:
            for element in elements_to_remove:
                element.remove()
        analyzer.to_remove = []
        time_in_ms = t / slow_down_factor * 1000.0
        if analyzer.last_t >= time_in_ms:
            value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
        else:
            value_dict_by_sec = records.get_vals_at_dt(t1=analyzer.last_t, t2=time_in_ms, default_res=0,
                                                       dt_func=func_for_missing_frames)
        for draw_func in draw_funcs:
            analyzer.to_remove.append(draw_func(analyzer.last_t, time_in_ms, segs, lines, ax, records))
        analyzer.last_t = time_in_ms
        if bounds:
            norm = get_norm(bounds)
        else:
            norm = get_norm([min_value, max_value])
        for line, seg in zip(lines, segs):
            line.set_color(cmap(norm(value_dict_by_sec[seg.sec][seg])))

        if dancing:
            distance = Distance(analyzer.cell, more_conductances_)
            distance.compute(time=time_in_ms, dt=1)
            for son in analyzer.cell.soma[0].children():
                if base_plot_type == 'dendogram':
                    make_frame_dendorgam_r(distance, son, start_point=dict(x=0, y=0))
                else:
                    make_frame_r(distance, son, start_point=dict(x=0, y=0))

        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
        if len(show_records_from) > 0:
            try:
                start_idx = np.where(time >= t - (voltage_window) / 1000 * slow_down_factor / 2)[0][0]
            except:
                start_idx = 0
            try:
                end_idx = np.where(time >= t + voltage_window / 1000 * slow_down_factor / 2)[0][0]
            except:
                end_idx = -1

            for i, seg_v in enumerate(voltage_segs):
                v1 = records.get_record(seg_v)[start_idx:end_idx]
                t1 = time[start_idx:end_idx] * 1000.0 / slow_down_factor
                ax_v[i].clear()
                ax_v[i].plot(t1, v1, color=show_records_from[seg_v]['color'])
                ax_v[i].set_ylim(min_value, max_value)
                ax_v[i].set_ylabel(ylabel)
                ax_v[i].set_xlabel(xlabel)
                ax_v[i].set_title(show_records_from[seg_v]['label'])
                ax_v[i].axvline(t * 1000.0 / slow_down_factor, color='r', ls='--')
                ax_v[i].set_xlim(xmin=time_in_ms - voltage_window / 2, xmax=time_in_ms + voltage_window / 2)

        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=time[-1])
    return animation


def create_morph_movie(analyzer, protocol=spike_protocol, cut_start_ms=0, record_name='v',
                       seg_to_indicate_dict=dict(), diam_factor=None,
                       sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
                       plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset='ultrafast',
                       slow_down_factor=1, func_for_missing_frames=np.mean):
    records, draw_funcs = analyzer.record_protocol(protocol=protocol, cut_start_ms=cut_start_ms, record_name=record_name)

    analyzer.create_movie_from_rec(records=records, slow_down_factor=slow_down_factor,
                               func_for_missing_frames=func_for_missing_frames, theta=theta,
                               scale=scale, cmap=cmap, seg_to_indicate_dict=seg_to_indicate_dict,
                               diam_factor=diam_factor, sec_to_change=sec_to_change, ignore_sections=ignore_sections,
                               plot_color_bar=plot_color_bar, draw_funcs=[])


def dancing_morph(analyzer, protocol, seg_to_indicate_dict=dict(), diam_factor=None,
                  sec_to_change=None, ignore_sections=[], theta=0, scale=500,
                  slow_down_factor=1,
                  figsize=(5, 5)):
    import matplotlib.style as mplstyle
    more_conductances_ = more_conductances(analyzer.cell, is_resting=False, protocol=protocol)
    mplstyle.use('fast')
    time = more_conductances_.time.copy()
    time /= 1000.0
    time *= slow_down_factor

    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    ax, lines, segs = analyzer.plot_morph(ax=ax, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                                      sec_to_change=sec_to_change, ignore_sections=ignore_sections,
                                      theta=theta, scale=scale, ignore_soma=True, distance=None,
                                      electrical=True, time=0.5, dt=1, more_conductances_=more_conductances_)

    segs = np.array(segs)
    lines = np.array(lines)
    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
    analyzer.last_t = 0
    analyzer.to_remove = []

    def make_frame_r(distance, sec, start_point):
        for seg in sec:
            for idx, seg2 in enumerate(segs):
                if seg2 == seg:
                    x, y = lines[idx].get_data()
                    dx = x[:-1] - x[1:]
                    dy = y[:-1] - y[1:]
                    current_lens = (dx ** 2 + dy ** 2) ** 0.5
                    wanted_len = distance.get_length(seg, electrical=True)
                    f = wanted_len / sum(current_lens)
                    # fixing the shift
                    x += start_point['x'] - x[0]
                    y += start_point['y'] - y[0]

                    # fixing the length
                    for point_idx in range(1, len(x), 1):
                        new_x = x[point_idx - 1] + (x[point_idx] - x[point_idx - 1]) * f
                        new_y = y[point_idx - 1] + (y[point_idx] - y[point_idx - 1]) * f
                        x[point_idx:] += (new_x - x[point_idx])
                        y[point_idx:] += (new_y - y[point_idx])
                    lines[idx].set_data(x, y)
                    start_point = dict(x=x[-1], y=y[-1])
        for son in sec.children():
            make_frame_r(distance, son, start_point)

    def make_frame(t):
        time_in_ms = t / slow_down_factor * 1000.0
        distance = Distance(analyzer.cell, more_conductances_)
        distance.compute(time=time_in_ms, dt=1)
        for son in analyzer.cell.soma[0].children():
            make_frame_r(distance, son, start_point=dict(x=0, y=0))
        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=time[-1])
    # animation.write_gif(os.path.join(save_to, clip_name + '.gif'), fps=int(1.0 / h.dt) if fps is None else fps/slow_down_factor)
    return animation

