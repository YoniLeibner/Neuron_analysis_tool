from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake
from Neuron_analysis_tool.morph_ploter import get_norm
from Neuron_analysis_tool.record import record, record_all
from Neuron_analysis_tool.distance import Distance
from Neuron_analysis_tool.protocols import *
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import os
from Neuron_analysis_tool.record import sec_name, seg_name




def plot_all_records_func(records, distance, ax, analyzer, distance_factor=1, slow_down_factor=1, plot_every=0.25):
    t1 = records.time.copy()
    for sec in analyzer.cell.all:
        for seg in sec:
            try:
                v1 = records.get_record(seg)
                start_end = distance.get_start_end(seg)
                if start_end['start']//plot_every == start_end['end']//plot_every: continue
                d = distance.get_mid_point(seg) * distance_factor
                if distance.get_part(seg) == 'parent':
                    d=-d
                color, _ = analyzer.colors.get_seg_color(seg)
                ax.plot(t1, v1+d, color=color)
            except:
                pass

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
                          base_plot_type='morph', start_seg=None, electrical=True, figsize=(8, 5),
                          dancing=False,distance=None, dt=1, more_conductances_=None,
                          plot_all_records=False, distance_factor=1, plot_every = 0.25,
                          color_bar_idx=None):

    assert dancing == False or more_conductances_ is not None
    if dancing:
        assert dancing == electrical, 'only electrical plot can dance!!!, got electrical='+str(electrical)+', and dancing='+str(dancing)
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
    time_lines = []
    voltage_segs = list(show_records_from.keys())
    if plot_all_records:
        number_of_record_plots = len(show_records_from)+1
    else:
        number_of_record_plots = len(show_records_from)
    if number_of_record_plots == 0:
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()
        if color_bar_idx  is None:
            color_bar_idx = [0.8, 0.2, 0.02, 0.6]
    else:
        from matplotlib.gridspec import GridSpec
        fig = plt.figure(constrained_layout=True, figsize=figsize)
        plt.subplots_adjust(wspace=.75, hspace=0.35)
        gs = GridSpec(number_of_record_plots, 2, figure=fig, width_ratios=[0.7, 0.3])
        ax = fig.add_subplot(gs[:, 0])
        ax_v = []
        for i in range(len(show_records_from)):
            ax_v.append(fig.add_subplot(gs[i, 1]))
            ax_v[-1].spines['top'].set_visible(False)
            ax_v[-1].spines['right'].set_visible(False)
            seg_v = voltage_segs[i]
            v1 = records.get_record(seg_v)
            t1 = time * 1000.0 / slow_down_factor
            ax_v[-1].plot(t1, v1, color=show_records_from[seg_v]['color'])
            # ax_v[-1].set_ylim(min_value, max_value)
            ax_v[-1].set_ylabel(ylabel)
            ax_v[-1].set_xlabel(xlabel)
            time_lines.append(ax_v[-1].axvline(0, color='r', ls='--'))
            if 'label' in show_records_from[seg_v]['label']:
                ax_v[-1].set_title(show_records_from[seg_v]['label'])
            else:
                ax_v[-1].set_title(str(seg_v.sec))
            ax_v[-1].set_xlim(xmin= -voltage_window / 2, xmax= voltage_window/2)

        if plot_all_records:
            if distance is None:
                distance = Distance(analyzer.cell, analyzer.more_conductances)
                distance.compute(start_seg=start_seg)
            ax_v.append(fig.add_subplot(gs[-1, 1]))
            plot_all_records_func(records, distance, ax_v[-1], analyzer, distance_factor=distance_factor, slow_down_factor=slow_down_factor, plot_every=plot_every)
            ax_v[-1].set_ylabel(ylabel)
            ax_v[-1].set_xlabel(xlabel)
            # max_distances = distance.get_max(electrical=True)
            # ax_v[-1].set_ylim(min_value-max_distances['parent']*distance_factor, max_value+max_distances['sons']*distance_factor)
            ax_v[-1].set_title('all records')
            ax_v[-1].set_xlim(xmin=-voltage_window / 2, xmax=voltage_window / 2)
            time_lines.append(ax_v[-1].axvline(0, color='r', ls='--'))
        if color_bar_idx is None:
            color_bar_idx = [0.5, 0.2, 0.02, 0.6]
    value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
    seg_to_indicate_dict.update(show_records_from)
    if base_plot_type == 'morph':
        ax, cax, colors, lines, segs = analyzer.plot_morph_with_values(value_dict_by_sec, ax=ax,
                                                                   seg_to_indicate_dict=dict() if dancing else seg_to_indicate_dict,
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
                                                                              segs_to_indecate=dict() if dancing else seg_to_indicate_dict,
                                                                              plot_legend=False,
                                                                              ignore_sections=ignore_sections,
                                                                              electrical=electrical,
                                                                              diam_factor=diam_factor, distance=None,
                                                                              bounds=[min_value, max_value], cmap=cmap,
                                                                              plot_color_bar=plot_color_bar,
                                                                              color_bar_idx=color_bar_idx, colors=None)

    elif base_plot_type == 'attenuation':
        def record_to_value(rec):
            return rec.max()
        ax, norm_by, lines, segs, records = analyzer.plot_attenuation(protocol=None,
                                                                      ax=ax,
                                                                      seg_to_indicate_dict=seg_to_indicate_dict,
                                                                      start_seg =start_seg,
                                                                      record_to_value_func=record_to_value,
                                                                      norm=False,
                                                                      record_name=records.record_name,
                                                                      norm_by=None,
                                                                      electrical=electrical,
                                                                      distance=None,
                                                                      records=records,
                                                                      ignore_sections=ignore_sections)
        ax.set_yscale('linear')
        ax.set_ylim([records.get_min()-0.5, records.get_max()+0.5])
            # .plot_dendogram_with_values(value_dict_by_sec, start_seg=start_seg,
            #                                                                   ax=ax,
            #                                                                   segs_to_indecate=dict() if dancing else seg_to_indicate_dict,
            #                                                                   plot_legend=False,
            #                                                                   ignore_sections=ignore_sections,
            #                                                                   electrical=electrical,
            #                                                                   diam_factor=diam_factor, distance=None,
            #                                                                   bounds=[min_value, max_value], cmap=cmap,
            #                                                                   plot_color_bar=plot_color_bar,
            #                                                                   color_bar_idx=color_bar_idx, colors=None)
    # elif base_plot_type == 'attenuation':
    #     pass

    else:
        raise Exception('base plot_type:', +str(base_plot_type) + ' not implementes')
    segs = np.array(segs)
    lines = np.array(lines)
    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    if dancing:
        x_range = abs(lim_x[1]-lim_x[0]) * 0.4
        y_range = abs(lim_y[1]-lim_y[0]) * 0.4
        lim_y = [lim_y[0] - y_range, lim_y[1] + y_range]
        ax.set_ylim(lim_y)

        if base_plot_type == 'morph':
            lim_x = [lim_x[0] - x_range, lim_x[1] + x_range]
            ax.set_xlim(lim_x)

    time_text = ax.text(lim_x[0], lim_y[0], 'time: 0.0 (ms)')
    analyzer.last_t = 0
    analyzer.to_remove = []

    def make_frame_dendorgam_r(distance, seg_to_indicate_dict=dict()):
        for seg, line in zip(segs, lines):
            start_end = distance.get_start_end(seg, electrical=True)
            y = line.get_ydata()
            if y[0] == y[1] == 0:
                continue
            else:
                if y[0]<0 or y[1]<0:
                    mul=-1
                else: mul=1
                if y[0] == y[1]: #horizontal_line
                    line.set_ydata([start_end['end']*mul]*2)
                else:
                    line.set_ydata([start_end['start']*mul, start_end['end']*mul])
                    if seg in seg_to_indicate_dict.keys():
                        analyzer.to_remove.append([ax.scatter(np.mean(line.get_xdata()), np.mean(line.get_ydata()),
                               color=seg_to_indicate_dict[seg]['color'], s=seg_to_indicate_dict[seg]['size'],
                               alpha=seg_to_indicate_dict[seg]['alpha'], zorder=3)])

    def make_frame_attenuation_r(distance, start_time=0, end_time=1000):
        if start_time == end_time: return
        # norm_by = records.get_record_at_dt(start_seg, start_time, end_time, dt_func=record_to_value)
        for seg, line in zip(segs, lines):
            # start_end = distance.get_start_end(seg, electrical=electrical)
            y_ = line.get_ydata()
            if y_[0] == y_[1] == 0:
                continue
            else:
                parent_seg = distance.get_seg_parent(seg)
                y = [records.get_record_at_dt(parent_seg, start_time, end_time, dt_func=record_to_value) / norm_by,
                     records.get_record_at_dt(seg, start_time, end_time, dt_func=record_to_value) / norm_by]

                line.set_ydata(y)


    def make_frame_r(distance, sec, start_point, seg_to_indicate_dict=dict()):
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
                    if seg in seg_to_indicate_dict.keys():
                        analyzer.to_remove.append([ax.scatter(x[len(x)//2], y[len(y)//2],
                               color=seg_to_indicate_dict[seg]['color'], s=seg_to_indicate_dict[seg]['size'],
                               alpha=seg_to_indicate_dict[seg]['alpha'], zorder=3)])
                    start_point = dict(x=x[-1], y=y[-1])
        for son in sec.children():
            make_frame_r(distance, son, start_point, seg_to_indicate_dict=seg_to_indicate_dict)


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

        if bounds:
            norm = get_norm(bounds)
        else:
            norm = get_norm([min_value, max_value])
        for line, seg in zip(lines, segs):
            line.set_color(cmap(norm(value_dict_by_sec[sec_name(seg.sec)][seg_name(seg)])))
        if base_plot_type == 'attenuation':
            distance = Distance(analyzer.cell, analyzer.more_conductances)
            distance.compute(start_seg=start_seg)

            make_frame_attenuation_r(distance, start_time=analyzer.last_t, end_time=time_in_ms)
        elif dancing:
            distance = Distance(analyzer.cell, more_conductances_)
            distance.compute(start_seg=start_seg, time=time_in_ms, dt=1)
            if base_plot_type == 'dendogram':
                make_frame_dendorgam_r(distance, seg_to_indicate_dict=seg_to_indicate_dict)

            else:
                for son in analyzer.cell.soma[0].children():
                    make_frame_r(distance, son, start_point=dict(x=0, y=0), seg_to_indicate_dict=seg_to_indicate_dict)
        analyzer.last_t = time_in_ms
        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
        if number_of_record_plots > 0:
            for i, t_line in enumerate(time_lines):
                ax_v[i].set_xlim(xmin=time_in_ms - voltage_window / 2, xmax=time_in_ms + voltage_window / 2)
                t_line.set_xdata([t * 1000.0 / slow_down_factor]*2)

        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=time[-1])
    return animation


# def create_morph_movie(analyzer, protocol=spike_protocol, cut_start_ms=0, record_name='v',
#                        seg_to_indicate_dict=dict(), diam_factor=None,
#                        sec_to_change=None, ignore_sections=[], theta=0, scale=0, cmap=plt.cm.turbo,
#                        plot_color_bar=True, save_to='', clip_name='clip', fps=None, threads=4, preset='ultrafast',
#                        slow_down_factor=1, func_for_missing_frames=np.mean):
#     records, draw_funcs = analyzer.record_protocol(protocol=protocol, cut_start_ms=cut_start_ms, record_name=record_name)
#
#     analyzer.create_movie_from_rec(records=records, slow_down_factor=slow_down_factor,
#                                func_for_missing_frames=func_for_missing_frames, theta=theta,
#                                scale=scale, cmap=cmap, seg_to_indicate_dict=seg_to_indicate_dict,
#                                diam_factor=diam_factor, sec_to_change=sec_to_change, ignore_sections=ignore_sections,
#                                plot_color_bar=plot_color_bar, draw_funcs=[])


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

