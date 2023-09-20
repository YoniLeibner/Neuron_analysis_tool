#########################################################
#
# author: Yoni Leibner
# description: Analyzer utils
# date of modification: 16.11.2022
#
#########################################################

import numpy as np
from neuron import h
import matplotlib as mpl

# chars of special characters lamda and micro
LAMDA = '\u03BB'
MICRO = '\u03BC'

def sec_name(sec):
    """
    change section into string for the section name
    :param sec:
    :return:
    """
    sec_full_name = sec.name()
    return sec_full_name.split('.')[-1]

def seg_name(seg):
    """
    change segment into string for the segment name
    :param seg:
    :return:
    """
    return str(seg.x)

def get_segment_length_lamda(seg, more_conductances, time=None, dt=1, dt_func= lambda x: np.mean(x)):
    """
    return the segment  e_length
    :param seg_len:
    :param RM:
    :param RA:
    :return:
    """
    sec = seg.sec
    seg_len = sec.L/sec.nseg

    # # old calc
    # R_total2 = (1.0/seg.g_pas)/(np.pi*seg.diam)
    # R_total2 = (1.0/(seg.g_pas+seg.Ih_human_linear.gIh))
    # lamda2 = np.sqrt(R_total2/sec.Ra*seg.diam*10**4/4)
    # return seg_len/lamda2

    R_total = more_conductances.cumpute(seg, time=time, dt=dt, dt_func=dt_func)*10**-4 # this in u ohm*um change to ohm*cm
    ri = sec.Ra * 4 / (np.pi * (seg.diam*10**-4) ** 2) # Ra is in ohm*cm we change ohm/cm
    lamda = (R_total / ri) ** 0.5 *10**2 #result i non units, this is in cm (10^2)
    return (float(seg_len)*10**-5) / lamda # change the length to cm to match with lamda

def get_segment_length_um(seg):
    """
    return the segment  e_length
    :param seg_len:
    :param RM:
    :param RA:
    :return:
    """
    sec = seg.sec
    return sec.L/sec.nseg

def have_parent(sec):
    """
    check if a section have parent seg or its a source
    :param sec:
    :return:
    """
    return not sec.parentseg() is None

def seg_Rin_func(seg):
    """
    function that get the input resistance of a segment
    :param seg:
    :return:
    """
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)

def video_player(local_path, video, mtype="video/mp4"):
    """ Displays mp4 video in Jupyter cell. Jupyter requires video
    in the same directory as the calling notebook. An assertion error
    will be thrown if this is not true.

    Parameters
    ----------
    video (str): the filename of the video. Example: "generating_bootstrap_replicates.mp4"
    mtype (str): the mime type of the video. Example: "video/mp4"

    """

    from IPython.display import HTML, display

    cwd = local_path

    # assert video in [file.name for file in list(cwd.glob('*'))], \
    #     f'{video} must be in local directory: {cwd}'

    call = """
    <video width=100% controls autoplay loop>
        <source src="{}" type="{}">
    </video>""".format(video+'?t='+str(np.random.rand(1)[0]), mtype)

    display(HTML(call))

def get_norm(all_vals):
    """
    getting the norm for ploting (allow working with colormaps
    :param all_vals:
    :return:
    """
    norm = mpl.colors.Normalize(vmin=np.min(all_vals), vmax=np.max(all_vals))
    return norm

def plot_shring_axes(plot_kwargs, share_idxs, scales=dict(x=10, y=2, x_text=MICRO + 'm', y_text=LAMDA)):
    def cable_limets_func(plot_kwargs):
        x_lims = [None, None]
        y_lims = [None, None]
        for i, kwargs in enumerate(plot_kwargs):
            if i not in share_idxs: continue
            xlim = kwargs['ax'].get_xlim()
            ylim = kwargs['ax'].get_ylim()
            if x_lims[0] is None:
                x_lims[0] = xlim[0]
                x_lims[1] = xlim[1]
                y_lims[0] = ylim[0]
                y_lims[1] = ylim[1]
            else:
                if x_lims[0] > xlim[0]: x_lims[0] = xlim[0]
                if x_lims[1] < xlim[1]: x_lims[1] = xlim[1]
                if y_lims[0] > ylim[0]: y_lims[0] = ylim[0]
                if y_lims[1] < ylim[1]: y_lims[1] = ylim[1]

        for i, kwargs in enumerate(plot_kwargs):
            if i not in share_idxs: continue
            ax=kwargs['ax']
            ax.set_xlim(x_lims)
            ax.set_ylim(y_lims)

            if scales is not None and 'x' in scales and 'y'in scales:
                ax.set_axis_off()
                ax.plot([x_lims[0], x_lims[0] + scales['x']], [y_lims[0]] * 2, color='k')
                ax.plot([x_lims[0]] * 2, [y_lims[0], y_lims[0] + scales['y']], color='k')
                y_scale_text = str(scales['y']) + ' ' + scales['y_text']
                x_scale_text = str(scales['x']) + ' ' + scales['x_text']
                y_pos = y_lims[0] + scales['y'] / 2
                x_pos = x_lims[0] + scales['x'] / 2

                if 'scale_text_size' in scales:
                    s = scales['scale_text_size']
                else:
                    s = 10
                ax.annotate(y_scale_text,
                            xy=(x_lims[0], y_pos), xycoords='data', size=s,
                            xytext=(-s - 2, -s / 2), textcoords='offset points', rotation=90)

                ax.annotate(x_scale_text,
                            xy=(x_pos, y_lims[0]), xycoords='data', size=s,
                            xytext=(-s, -s - 2), textcoords='offset points')
    return cable_limets_func
