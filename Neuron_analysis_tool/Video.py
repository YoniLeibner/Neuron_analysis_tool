from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake
from Neuron_analysis_tool.morph_ploter import get_norm
from Neuron_analysis_tool.record import record, record_all
from Neuron_analysis_tool.distance import Distance
from Neuron_analysis_tool.protocols import *
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import os
from Neuron_analysis_tool.utils import sec_name, seg_name, LAMDA, MICRO




def plot_all_records_func(records, distance, ax, analyzer, distance_factor=1, plot_every=0.25):
    t1 = records.time.copy()
    soma_sec = list(analyzer.cell.soma[0])
    soma_seg = soma_sec[len(soma_sec)//2]
    for sec in analyzer.cell.all:
        for seg in sec:
            try:
                v1 = records.get_record(seg)
                start_end = distance.get_start_end(seg)
                if (start_end['start']//plot_every == start_end['end']//plot_every) and (not seg == soma_seg): # show the soma in all cases!!
                    continue
                d = distance.get_mid_point(seg) * distance_factor
                if distance.get_part(seg) == 'parent':
                    d=-d
                color, _ = analyzer.colors.get_seg_color(seg)
                ax.plot(t1, v1+d, color=color)
            except:
                pass


def initiate_ax_single_record(ax, analyzer, seg, records, slow_down_factor, color, voltage_window=50, xlabel='time (ms)', ylabel=None, title=None):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylabel is None:
        ylabel = records.record_name
    time = records.time.copy()
    time /= 1000.0
    time *= slow_down_factor
    v1 = records.get_record(seg)
    t1 = time * 1000.0 / slow_down_factor
    ax.plot(t1, v1, color=color)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    time_lines = ax.axvline(0, color='r', ls='--')

    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(str(seg.sec.name()))
    ax.set_xlim(xmin=-voltage_window / 2, xmax=voltage_window / 2)

    def update(time_in_ms):

        ax.set_xlim(xmin=time_in_ms - voltage_window / 2, xmax=time_in_ms + voltage_window / 2)
        time_lines.set_xdata([time_in_ms] * 2)

    return update

def initiate_ax_all_records(ax, analyzer, records, start_seg, distance=None, distance_factor=0, plot_every = 0.25,
                            voltage_window=50, xlabel='time (ms)', ylabel=None, dt_func= lambda x: np.mean(x)):

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylabel is None:
        ylabel = records.record_name
    if distance is None:
        distance = Distance(analyzer.cell, analyzer.more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)
    plot_all_records_func(records, distance, ax, analyzer, distance_factor=distance_factor, plot_every=plot_every)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title('all records')
    ax.set_xlim(xmin=-voltage_window / 2, xmax=voltage_window / 2)
    time_lines=ax.axvline(0, color='r', ls='--')

    def update(time_in_ms):
        ax.set_xlim(xmin=time_in_ms - voltage_window / 2, xmax=time_in_ms + voltage_window / 2)
        time_lines.set_xdata([time_in_ms] * 2)

    return update

def initiate_ax_morph(ax, analyzer, records, bounds=None, margin=0, dancing=False, seg_to_indicate_dict=dict(), diam_factor=None, sec_to_change=None,
                      ignore_sections=[], theta=0, scale=500, cmap=plt.cm.turbo, plot_color_bar=True, color_bar_kwarts = dict(shrink=0.6), distance=None,
                      electrical=False, more_conductances_=None, draw_funcs=[], func_for_missing_frames=np.max, dt_func= lambda x: np.mean(x)):
    assert dancing == False or more_conductances_ is not None
    value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
    min_value = records.get_min() - margin
    max_value = records.get_max() + margin
    if bounds is not None:
        min_value = bounds[0]
        max_value = bounds[1]
    ax, cax, colors, lines, segs = analyzer.plot_morph_with_values(value_dict_by_sec, ax=ax,
                                                                   seg_to_indicate_dict=dict() if dancing else seg_to_indicate_dict,
                                                                   diam_factor=diam_factor,
                                                                   sec_to_change=sec_to_change,
                                                                   ignore_sections=ignore_sections,
                                                                   theta=theta, scale=scale,
                                                                   cmap=cmap, bounds=[min_value, max_value],
                                                                   plot_color_bar=plot_color_bar,
                                                                   color_bar_kwarts=color_bar_kwarts,
                                                                   distance=distance,
                                                                   electrical=electrical,
                                                                   time=0.5,
                                                                   dt=1,
                                                                   more_conductances_=more_conductances_,
                                                                   dt_func=dt_func)
    segs=np.array(segs)
    lines=np.array(lines)
    if bounds:
        norm = get_norm(bounds)
    else:
        norm = get_norm([min_value, max_value])

    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    if dancing:

        x_range = abs(lim_x[1] - lim_x[0]) * 0.4
        lim_x = [lim_x[0] - x_range, lim_x[1] + x_range]

        y_range = abs(lim_y[1] - lim_y[0]) * 0.4
        lim_y = [lim_y[0] - y_range, lim_y[1] + y_range]

        ax.set_xlim(lim_x)
        ax.set_ylim(lim_y)


    time_text = ax.annotate('time: 0.0 (ms)'.format(lim_y[1]), ha='center',
            size=10, xy=(np.mean(lim_x), 1), xytext=(0, 5),
            textcoords='offset points',
            xycoords=('data', 'axes fraction'))

    # time_text = ax.text(lim_x[1], lim_y[1], 'time: 0.0 (ms)')


    def dancing_morph_r(distance, sec, start_point, seg_to_indicate_dict=dict()):
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
            dancing_morph_r(distance, son, start_point, seg_to_indicate_dict=seg_to_indicate_dict)

    def update(time_in_ms):

        if analyzer.last_t >= time_in_ms:
            value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
        else:
            value_dict_by_sec = records.get_vals_at_dt(t1=analyzer.last_t, t2=time_in_ms, default_res=0,
                                                       dt_func=func_for_missing_frames)
        for draw_func in draw_funcs:
            analyzer.to_remove.append(draw_func(analyzer.last_t, time_in_ms, segs, lines, ax, records))

        for line, seg in zip(lines, segs):
            line.set_color(cmap(norm(value_dict_by_sec[sec_name(seg.sec)][seg_name(seg)])))

        if dancing:
            distance = Distance(analyzer.cell, more_conductances_, dt_func=dt_func)
            distance.compute(start_seg=None, time=time_in_ms, dt=1)
            for son in analyzer.cell.soma[0].children():
                dancing_morph_r(distance, son, start_point=dict(x=0, y=0), seg_to_indicate_dict=seg_to_indicate_dict)

        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')

    return update



def initiate_ax_dendogram(ax, analyzer, records, start_seg=None, seg_to_indicate_dict=dict(), ignore_sections=[], electrical=True,
                          diam_factor=None, margin=0, bounds=None, cmap=plt.cm.turbo, plot_color_bar=True, color_bar_kwarts = dict(shrink=0.6),
                          dancing=False, more_conductances_=None, draw_funcs=[], func_for_missing_frames=np.max, dt_func= lambda x: np.mean(x)):
    assert dancing == False or more_conductances_ is not None
    if start_seg is None:
        start_sec = list(analyzer.cell.soma[0])
        start_seg = start_sec[len(start_sec) // 2]
    value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
    min_value = records.get_min() - margin
    max_value = records.get_max() + margin
    if bounds is not None:
        min_value = bounds[0]
        max_value = bounds[1]
    ax, x_pos, cax, colors, lines, segs = analyzer.plot_dendogram_with_values(value_dict_by_sec, start_seg=start_seg,
                                                                              ax=ax,
                                                                              segs_to_indecate=dict() if dancing else seg_to_indicate_dict,
                                                                              plot_legend=False,
                                                                              ignore_sections=ignore_sections,
                                                                              electrical=electrical,
                                                                              diam_factor=diam_factor, distance=None,
                                                                              bounds=[min_value, max_value], cmap=cmap,
                                                                              plot_color_bar=plot_color_bar,
                                                                              color_bar_kwarts=color_bar_kwarts, colors=None, dt_func=dt_func)
    # cax.set_aspect(5, anchor='C')
    segs=np.array(segs)
    lines=np.array(lines)
    if bounds:
        norm = get_norm(bounds)
    else:
        norm = get_norm([min_value, max_value])

    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    if dancing:
        x_range = abs(lim_x[1] - lim_x[0]) * 0.4
        lim_x = [lim_x[0] - x_range, lim_x[1] + x_range]
        y_range = abs(lim_y[1] - lim_y[0]) * 0.4
        lim_y = [lim_y[0] - y_range, lim_y[1] + y_range]

        ax.set_xlim(lim_x)
        ax.set_ylim(lim_y)
    time_text = ax.annotate('time: 0.0 (ms)'.format(lim_y[1]), ha='center',
            size=10, xy=(np.mean(lim_x), 1), xytext=(0, 5),
            textcoords='offset points',
            xycoords=('data', 'axes fraction'))
    # time_text = ax.text(lim_x[1], lim_y[1], 'time: 0.0 (ms)')


    def dancing_dendogram_r(distance):
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

    def update(time_in_ms):

        if analyzer.last_t >= time_in_ms:
            value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
        else:
            value_dict_by_sec = records.get_vals_at_dt(t1=analyzer.last_t, t2=time_in_ms, default_res=0,
                                                       dt_func=func_for_missing_frames)
        for draw_func in draw_funcs:
            analyzer.to_remove.append(draw_func(analyzer.last_t, time_in_ms, segs, lines, ax, records))

        for line, seg in zip(lines, segs):
            line.set_color(cmap(norm(value_dict_by_sec[sec_name(seg.sec)][seg_name(seg)])))

        if dancing:
            distance = Distance(analyzer.cell, more_conductances_, dt_func=dt_func)
            distance.compute(start_seg=start_seg, time=time_in_ms, dt=1)
            dancing_dendogram_r(distance)

        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
    return update

def initiate_ax_attenuation(ax, analyzer, records, electrical=True, start_seg=None, seg_to_indicate_dict=dict(), ignore_sections=[],
                            draw_funcs=[], cmap=plt.cm.turbo, func_for_missing_frames=np.max, bounds=None, margin=0.5, dt_func= lambda x: np.mean(x)):
    def record_to_value(rec):
        return rec.max()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    min_value = records.get_min()
    max_value = records.get_max()
    if bounds is not None:
        min_value = bounds[0]
        max_value = bounds[1]
    ax, norm_by, lines, segs, records = analyzer.plot_attenuation(protocol=None,
                                                                  ax=ax,
                                                                  seg_to_indicate_dict=seg_to_indicate_dict,
                                                                  start_seg=start_seg,
                                                                  record_to_value_func=record_to_value,
                                                                  norm=False,
                                                                  record_name=records.record_name,
                                                                  norm_by=None,
                                                                  electrical=electrical,
                                                                  distance=None,
                                                                  records=records,
                                                                  ignore_sections=ignore_sections,
                                                                  dt_func=dt_func)
    segs=np.array(segs)
    lines=np.array(lines)
    ax.set_yscale('linear')
    ax.set_ylabel(records.record_name)
    ax.set_ylim([records.get_min() - margin, records.get_max() + margin])
    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    time_text = ax.annotate('time: 0.0 (ms)'.format(lim_y[1]), ha='center',
            size=10, xy=(np.mean(lim_x), 1), xytext=(0, 5),
            textcoords='offset points',
            xycoords=('data', 'axes fraction'))
    # time_text = ax.text(lim_x[1], lim_y[1], 'time: 0.0 (ms)')
    if bounds:
        norm = get_norm(bounds)
    else:
        norm = get_norm([min_value, max_value])


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

    def update(time_in_ms):

        if analyzer.last_t >= time_in_ms:
            value_dict_by_sec = records.get_vals_at_t(t=0, default_res=0)
        else:
            value_dict_by_sec = records.get_vals_at_dt(t1=analyzer.last_t, t2=time_in_ms, default_res=0,
                                                       dt_func=func_for_missing_frames)
        for draw_func in draw_funcs:
            analyzer.to_remove.append(draw_func(analyzer.last_t, time_in_ms, segs, lines, ax, records))

        for line, seg in zip(lines, segs):
            line.set_color(cmap(norm(value_dict_by_sec[sec_name(seg.sec)][seg_name(seg)])))

        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')
        distance = Distance(analyzer.cell, analyzer.more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)
        make_frame_attenuation_r(distance, start_time=analyzer.last_t, end_time=time_in_ms)

    return update

def initiate_ax_cable(ax, analyzer, more_conductances, cable_type='d3_2', start_seg=None, seg_to_indicate_dict=dict(),
                      factor_e_space=25, factor_m_space=10, ignore_sections=[], start_loc=0, shift=0, x_axis = True, factor=1,
                      dots_size=10, start_color='k', plot_legend = True,  extra = 5, cable_factor = 1, labal_start = None,
                      scales=dict(x=10, y=1), dt_func= lambda x: np.mean(x), title=None):

    distance = Distance(analyzer.cell, more_conductances, dt_func=dt_func)
    distance.compute(start_seg=start_seg, time=0, dt=1)
    analyzer.plot_cable(start_seg = start_seg, ax = ax,
                        factor_e_space = factor_e_space, factor_m_space = factor_m_space, segs_to_indecate = seg_to_indicate_dict,
                        ignore_sections = ignore_sections, cable_type = cable_type, start_loc = start_loc, shift=shift, x_axis = x_axis,
                        factor = factor, dots_size = dots_size, start_color = start_color,
                        plot_legend = plot_legend, distance = distance, extra = extra, cable_factor = cable_factor, labal_start = labal_start,
                        return_shift = False, dt_func=dt_func)

    lim_x = ax.get_xlim()
    lim_y = ax.get_ylim()
    x_range = abs(lim_x[1] - lim_x[0]) * 0.4
    lim_x = [lim_x[0] - x_range, lim_x[1] + x_range]
    ax.set_xlim(lim_x)
    y_range = abs(lim_y[1] - lim_y[0]) * 0.4
    lim_y = [lim_y[0] - y_range, lim_y[1] + y_range]
    ax.set_ylim(lim_y)
    if (scales is not None) and ('x' in scales) and ('y' in scales):
        ax.set_axis_off()
        ax.plot([lim_x[0], lim_x[0]+scales['x']], [lim_y[0]]*2, color='k')
        ax.plot([lim_x[0]]*2, [lim_y[0], lim_y[0]+scales['y']], color='k')
        y_scale_text = str(scales['y']) + ' ' + LAMDA
        x_scale_text = str(scales['x']) + ' ' + MICRO + 'm'
        y_pos = lim_y[0] + scales['y'] / 2
        x_pos = lim_x[0] + scales['x'] / 2
        if 'scale_text_size' in scales:
            s = scales['scale_text_size']
        else:
            s = 10
        ax.annotate(y_scale_text,
                    xy=(lim_x[0], y_pos), xycoords='data', size=s,
                    xytext=(-s - 2, -s / 2), textcoords='offset points', rotation=90)

        ax.annotate(x_scale_text,
                    xy=(x_pos, lim_y[0]), xycoords='data', size=s,
                    xytext=(-s, -s - 2), textcoords='offset points')

    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    time_text = ax.annotate('time: 0.0 (ms)'.format(lim_y[1]), ha='center',
                            size=10, xy=(np.mean(lim_x), 1), xytext=(0, 5),
                            textcoords='offset points',
                            xycoords=('data', 'axes fraction'))
    if title is not None:
        ax.annotate(title.format(ax.get_ylim()[1]), ha='center',
                                   size=10, xy=(np.mean(ax.get_xlim()), 1), xytext=(0, 17),
                                   textcoords='offset points',
                                   xycoords=('data', 'axes fraction'))
    def update(time_in_ms):
        distance = Distance(analyzer.cell, more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg, time=time_in_ms, dt=1)
        ax.clear()
        analyzer.plot_cable(start_seg=start_seg, ax=ax,
                            factor_e_space=factor_e_space, factor_m_space=factor_m_space,
                            segs_to_indecate=seg_to_indicate_dict,
                            ignore_sections=ignore_sections, cable_type=cable_type, start_loc=start_loc, shift=shift, x_axis=x_axis,
                            factor=factor, dots_size=dots_size, start_color=start_color,
                            plot_legend=plot_legend, distance=distance, extra=extra, cable_factor=cable_factor,
                            labal_start=labal_start, return_shift=False)
        ax.set_xlim(lim_x)
        ax.set_ylim(lim_y)
        if (scales is not None) and ('x' in scales) and ('y' in scales):
            ax.set_axis_off()
            ax.plot([lim_x[0], lim_x[0] + scales['x']], [lim_y[0]] * 2, color='k')
            ax.plot([lim_x[0]] * 2, [lim_y[0], lim_y[0] + scales['y']], color='k')
            y_scale_text = str(scales['y']) + ' ' + LAMDA
            x_scale_text = str(scales['x']) + ' ' + MICRO + 'm'
            y_pos = lim_y[0] + scales['y'] / 2
            x_pos = lim_x[0] + scales['x'] / 2
            if 'scale_text_size' in scales:
                s = scales['scale_text_size']
            else:
                s = 10
            ax.annotate(y_scale_text,
                        xy=(lim_x[0], y_pos), xycoords='data', size=s,
                        xytext=(-s - 2, -s / 2), textcoords='offset points', rotation=90)

            ax.annotate(x_scale_text,
                        xy=(x_pos, lim_y[0]), xycoords='data', size=s,
                        xytext=(-s, -s - 2), textcoords='offset points')

        else:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        time_text = ax.annotate('time: 0.0 (ms)'.format(lim_y[1]), ha='center',
                                size=10, xy=(np.mean(lim_x), 1), xytext=(0, 5),
                                textcoords='offset points',
                                xycoords=('data', 'axes fraction'))
        time_text.set_text('time: ' + str(round(time_in_ms, 1)) + ' (ms)')

        if title is not None:
            ax.annotate(title.format(ax.get_ylim()[1]), ha='center',
                        size=10, xy=(np.mean(ax.get_xlim()), 1), xytext=(0, 17),
                        textcoords='offset points',
                        xycoords=('data', 'axes fraction'))

    return update

def initiate_ax(ax, analyzer, plot_type, records=None, slow_down_factor=1,  electrical=False, seg=None, dancing=False, more_conductances_=None,
                color='k', voltage_window=50, xlabel='time (ms)', ylabel=None, title=None, # for voltage
                distance_factor=1, plot_every=0.25, # for voltage_all
                seg_to_indicate_dict=dict(), ignore_sections=[], draw_funcs=[], cmap=plt.cm.turbo, func_for_missing_frames=np.max, bounds=None, # for attenuation
                margin=0, diam_factor=None, sec_to_change=None, theta=0, scale=0, plot_color_bar=True, color_bar_kwarts = dict(shrink=0.6),# for morph
                cable_type='d3_2', factor_e_space=25, factor_m_space=10, start_loc=0, shift=None, x_axis=True, factor=1,
                dots_size=10, start_color='k', plot_legend=True, extra=5, cable_factor=1, labal_start=None, scales=dict(x=10, y=2),
                dt_func=lambda x: np.mean(x)):

    if seg is None:
        seg = list(analyzer.cell.soma[0])
        seg = seg[len(seg) // 2]

    if plot_type == 'single_record':
        return initiate_ax_single_record(ax, analyzer, seg, records, slow_down_factor, color, voltage_window=voltage_window, xlabel=xlabel, ylabel=ylabel, title=title)
    elif plot_type == 'all_records':
        return initiate_ax_all_records(ax, analyzer, records, seg, distance=None, distance_factor=distance_factor, plot_every = plot_every,
                            voltage_window=voltage_window, xlabel=xlabel, ylabel=ylabel, dt_func=dt_func)
    elif plot_type == 'attenuation':
        return initiate_ax_attenuation(ax, analyzer, records, electrical=electrical, start_seg=seg, seg_to_indicate_dict=seg_to_indicate_dict, ignore_sections=ignore_sections,
                            draw_funcs=draw_funcs, cmap=cmap, func_for_missing_frames=func_for_missing_frames, bounds=bounds, margin=margin, dt_func=dt_func)

    if dancing:
        assert dancing == electrical, 'only electrical plot can dance!!!, got electrical=' + str(
            electrical) + ', and dancing=' + str(dancing)
        assert more_conductances_ is not None
        dancing_text = ax.annotate('dancing'.format(ax.get_ylim()[1]), ha='center',
                                size=10, xy=(np.mean(ax.get_xlim()), 1), xytext=(0, 17),
                                textcoords='offset points',
                                xycoords=('data', 'axes fraction'))

    if plot_type == 'morph':
        return initiate_ax_morph(ax, analyzer, records, bounds=bounds, margin=margin, dancing=dancing, seg_to_indicate_dict=seg_to_indicate_dict, diam_factor=diam_factor,
                                 sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta, scale=scale, cmap=cmap, plot_color_bar=plot_color_bar, color_bar_kwarts = color_bar_kwarts ,
                                 distance=None, electrical=electrical, more_conductances_=more_conductances_, draw_funcs=draw_funcs, func_for_missing_frames=func_for_missing_frames, dt_func=dt_func)
    elif plot_type == 'dendogram':
        return initiate_ax_dendogram(ax, analyzer, records, start_seg=seg, seg_to_indicate_dict=seg_to_indicate_dict, ignore_sections=ignore_sections, electrical=electrical,
                          diam_factor=diam_factor, margin=margin, bounds=bounds, cmap=cmap, plot_color_bar=plot_color_bar, color_bar_kwarts = color_bar_kwarts ,
                          dancing=dancing, more_conductances_=more_conductances_, draw_funcs=draw_funcs, func_for_missing_frames=func_for_missing_frames, dt_func=dt_func)
    elif plot_type == 'cable':
        assert more_conductances_ is not None
        return initiate_ax_cable(ax, analyzer, more_conductances =more_conductances_, cable_type=cable_type, start_seg=seg,
                                 seg_to_indicate_dict=seg_to_indicate_dict, factor_e_space=factor_e_space, factor_m_space=factor_m_space,
                                 ignore_sections=ignore_sections, start_loc=start_loc, shift=shift, x_axis=x_axis, factor=factor, dots_size=dots_size,
                                 start_color=start_color, plot_legend=plot_legend, extra=extra, cable_factor=cable_factor, labal_start=labal_start, scales=scales, title=title, dt_func=dt_func)

    else:
        raise Exception('not a valid plot_type: '+plot_type +', not implemented')


def create_movie_from_rec(analyzer, fig, slow_down_factor=1, plot_kwargs=[], func_before_run=[], func_during_run=[]):
    analyzer.last_t = 0
    analyzer.to_remove = []
    update_funcs = []
    for kwargs in plot_kwargs:
        update_funcs.append(initiate_ax(analyzer=analyzer, slow_down_factor=slow_down_factor, **kwargs))
    for func in func_before_run:
        func()

    for func in func_during_run:
        func(plot_kwargs)
    def make_frame(t):
        for elements_to_remove in analyzer.to_remove:
            for element in elements_to_remove:
                element.remove()
        analyzer.to_remove = []
        time_in_ms = t / slow_down_factor * 1000.0
        if time_in_ms  == 0.0:
            analyzer.last_t = 0
        for update_fun in update_funcs:
            update_fun(time_in_ms)
        analyzer.last_t = time_in_ms
        for func in func_during_run:
            func(plot_kwargs)
        return mplfig_to_npimage(fig)
    print('duration=', plot_kwargs[0]['records'].time[-1]*slow_down_factor/1000.0)
    animation = VideoClip(make_frame, duration=plot_kwargs[0]['records'].time[-1]*slow_down_factor/1000.0)
    return animation

def save_movie_from_rec(analyzer, fig, slow_down_factor=1, plot_kwargs=[], func_before_run=[], f2=[], save_to='', clip_name='clip', fps=None, threads=4, preset='medium'): #ultrafast
    
    animation = create_movie_from_rec(analyzer=analyzer, fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs, func_before_run=func_before_run, f2=f2)
    if not clip_name.endswith('.mp4'):
        clip_name+='.mp4'
    animation.write_videofile(os.path.join(save_to, clip_name),
                              fps=int(1000.0 / h.dt) if fps is None else fps, threads=threads,
                              audio=False, preset=preset)
