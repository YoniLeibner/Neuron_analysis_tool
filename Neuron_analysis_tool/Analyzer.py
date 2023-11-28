#########################################################
#
# author: Yoni Leibner
# description: Analyzer class allows to easy visualization
#              and exploration of neurons from multiple
#              point of view, and also to create neuron
#              movies.
# date of modification:  11.05.2023
#
#########################################################

# from neuron import h, gui
# from Neuron_analysis_tool.distance import Distance
# from Neuron_analysis_tool.morph_ploter import get_norm
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import os
# import matplotlib as mpl
# import numbers

from Neuron_analysis_tool.protocols import *
from Neuron_analysis_tool.color_func import color_func, color_func_norm #, color_func_by_func
from Neuron_analysis_tool.morph_ploter import plot_morph
from Neuron_analysis_tool.dendogram import plot_dendogram
from Neuron_analysis_tool.cable import get_cable, plot_cable
from Neuron_analysis_tool.attenuation import plot_attenuation, record_to_value
from Neuron_analysis_tool.record import multi_record_all #record, record_all,
from Neuron_analysis_tool.utils import seg_Rin_func #, get_segment_length_lamda, get_segment_length_um, LAMDA, MICRO, sec_name, seg_name
from Neuron_analysis_tool.loaders import open_swc, open_L5PC, open_ASC, open_rall_tree, get_parts_and_colors, insert_g_total #, open_morph
from Neuron_analysis_tool.Video import *
from Neuron_analysis_tool.more_conductances import more_conductances, more_conductances_fake, callback
from tqdm import tqdm
from neuron import h
import neuron
# import os


import sys
if hasattr(sys, 'getwindowsversion'):
    print('error in loading mechanism')
    # h.load_dll(os.path.dirname(os.path.realpath(__file__))+'/')
else:
    try:
        neuron.load_mechanisms(os.path.dirname(os.path.realpath(__file__)))
    except:
        neuron.load_mechanisms(os.path.dirname(os.path.realpath(__file__))+'/x86_64/special')

class Analyzer:
    """
    class that help to analyze neuron model from every location on the morphology and a givin protocol
    main components:
    distance - compute both the phisical distance in (um) and in electrical units (lamda(0)) from a givin segment
    records - record the activity of the neuron (record a givin values (such as voltage/conductances .etc)

    you can generate the folowing plots:
        morphology - in 2d
        dendogram from a givin starting point
        cable from a givin starting point
        attenuation of a recorded value along distance
        plot a record
        plot all the records

    you can also create a vidio from this protocol and see the resolting record as a color map on the plots
    (morphology, dendorgam, attenuation)
    note that the morphology, dendorgam and cable plots are dancing!!
    """
    def __init__(self, cell=None, parts_dict=None, colors_dict=None, type='input_cell',
                 morph_path = None, Rm=10000.0, Ra=100, Cm=1, e_pas=-70,
                 more_conductances_protocol =resting_protocol, seg_every=20):
        """
        you can initiate this tool with
            type='input_cell' -> with a givin cell model
            type == 'Rall_tree' -> will give you a Rall tree model
            type == 'ASC'/ type == 'swc' -> will load a morphology file from morph_path
            type == 'L5PC' -> will load Itay Hay L5PC model

        the init function can take a short time since its runing a resting membrane protocol in order to get the
        resting conductances of the model (you can change it by passing a diffrent more_conductances_protocol

        :param cell: a neuron cell model.
        :param parts_dict: split os the model into parts, this is a dictionary of {part_name: [seg list]}. default is soma, basal, apical, axon.
        :param colors_dict: a dictinary of colors to use for each part {part_name: color}. default is soma, basal, apical, axon.
        :param type: the type of open from: input_cell, Rall_tree, ASC, swc, L5PC
        :param morph_path: the path to the ASC/swc file used only if type == 'ASC'/ type == 'swc'
        :param Rm: the Rm value type == 'ASC'/ type == 'swc'
        :param Ra: the Ra value type == 'ASC'/ type == 'swc'
        :param Cm: the Cm value type == 'ASC'/ type == 'swc'
        :param e_pas: the e_pas value type == 'ASC'/ type == 'swc'
        :param more_conductances_protocol: protocol to run to get the conductances for the distance. default is 500 ms of nothing
        :param seg_every: how many segments there are in every um
        """
        if cell is None:
            if type == 'Rall_tree':
                cell, parts_dict, colors_dict = open_rall_tree(seg_every=seg_every)
            elif type == 'ASC' and morph_path:
                cell, parts_dict, colors_dict = open_ASC(morph_path, Rm=Rm, Ra=Ra,
                                                         Cm=Cm, e_pas=e_pas, seg_every=seg_every)
            elif type == 'swc' and morph_path:
                cell, parts_dict, colors_dict = open_swc(morph_path, Rm=Rm, Ra=Ra,
                                                         Cm=Cm, e_pas=e_pas, seg_every=seg_every)
            elif type == 'L5PC':
                cell, parts_dict, colors_dict = open_L5PC(seg_every=seg_every)

        if parts_dict is None:
            cell, parts_dict, colors_dict = get_parts_and_colors(cell)
        insert_g_total(cell)
        self.type=type
        self.cell=cell
        self.parts_dict = parts_dict
        #
        soma_segs = list(cell.soma[0])
        self.bscallback = h.beforestep_callback(soma_segs[len(soma_segs)//2])  # starting the ref of g_tatals (this asume you have less then 19 mechanisms and les then 39 synapses per segment!!!!!!!
        self.bscallback.set_callback((callback, self))
        self.more_conductances = more_conductances(cell, is_resting=True,
                                                   protocol=more_conductances_protocol)
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=parts_dict, color_dict=colors_dict)


    def get_mechanism_names(self):
        """

        :return: the names of the mechanism in the model
        """
        mechanisms_names = set()
        for sec in self.cell.all:
            for seg in sec:
                for mechanisms in seg:
                    if not mechanisms.is_ion():
                        mechanisms_names.add(str(mechanisms))
        return list(mechanisms_names)

    def get_synapses(self):
        """

        :return: the names of the mechanism in the model
        """
        synapses = dict()
        for sec in self.cell.all:
            for seg in sec:
                for syn in seg.point_processes():
                    if sec_name(sec) not in synapses:
                        synapses[sec_name(sec)] = []
                    syn_type = syn.hname().split('[')[0]
                    if hasattr(syn, 'g'):
                            synapses[sec_name(sec)].append(dict(syn=syn, g_name='g'))
                    else:
                        synapses[sec_name(sec)].append(dict(syn=syn, g_name='g_'+syn_type))
        return synapses

    def change_color_dict(self, colors_dict):
        """
        change the colors in the plots
        :param colors_dict: a dictinary of colors to use for each part {part_name: color}.
        :return:
        """
        for part in self.parts_dict:
            assert part in colors_dict

        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=self.parts_dict, color_dict=colors_dict)

    def change_parts_dict(self, parts_dict, colors_dict):
        """
        change the naming parts and color in the model
        :param parts_dict: split os the model into parts, this is a dictionary of {part_name: [seg list]}.
        :param colors_dict:  a dictinary of colors to use for each part {part_name: color}.
        :return:
        """
        for part in parts_dict:
            assert part in colors_dict
        self.parts_dict = parts_dict
        self.colors_dict = colors_dict
        self.colors = color_func(parts_dict=self.parts_dict, color_dict=colors_dict)

    def save_morph_to_swc_helper(self, f, all_points_3d, sec, counter=-1, parent_counter=-1):
        if sec not in all_points_3d: return counter
        if sec in self.cell.soma:
            type_num = '1'
            sec_data = all_points_3d[sec]
            total_x = [sec_data[0]['x'].min(), sec_data[0]['x'].max()]
            total_y = [sec_data[0]['y'].min(), sec_data[0]['y'].max()]
            total_z = [sec_data[0]['z'].min(), sec_data[0]['z'].max()]
            for seg_data in sec_data:
                total_x = [min(seg_data['x'].min(), total_x[0]),
                           max(seg_data['x'].max(), total_x[1])]
                total_y = [min(seg_data['y'].min(), total_y[0]),
                           max(seg_data['y'].max(), total_y[1])]
                total_z = [min(seg_data['z'].min(), total_z[0]),
                           max(seg_data['z'].max(), total_z[1])]
            f.write(' '.join([str(counter),
                              type_num,
                              str(total_x[0]),
                              str(total_y[0]),
                              str(total_z[0]),
                              str(sec.diam / 2.0),
                              str( )
                              ]) + '\n')
            parent_counter=counter
            counter+=1
            f.write(' '.join([str(counter),
                              type_num,
                              str(total_x[1]),
                              str(total_y[1]),
                              str(total_z[1]),
                              str(sec.diam / 2.0),
                              str(parent_counter)
                              ]) + '\n')
            parent_counter=counter
            counter+=1
        else:
            sec_data = all_points_3d[sec]
            for seg in sec:
                for seg_data in sec_data:
                    if seg_data['seg']==seg:
                        for i in range(seg_data['x'].shape[0]):
                            type_num = '0'
                            if len(self.cell.axon) > 1 and sec in self.cell.axon:
                                type_num = '2'
                            elif len(self.cell.dend) > 1 and sec in self.cell.dend:
                                type_num = '3'
                            elif len(self.cell.apic) > 1 and sec in self.cell.apic:
                                type_num = '4'
                            f.write(' '.join([str(counter),
                                              type_num,
                                              str(seg_data['x'][i]),
                                              str(seg_data['y'][i]),
                                              str(seg_data['z'][i]),
                                              str(seg.diam/2.0),
                                              str(parent_counter)
                                              ])+'\n')
                            parent_counter=counter
                            counter+=1

        for son in sec.children():
            counter = self.save_morph_to_swc_helper(f, all_points_3d, son, counter=counter, parent_counter=parent_counter)
        return counter

    def save_morph_to_swc(self, save_name, electrical=False, distance=None, more_conductances=None, time=None, dt=1,
                          dt_func= lambda x: np.mean(x)):
        if more_conductances is None:
            more_conductances=self.more_conductances
        from .morph_ploter import get_dots
        all_points_3d = get_dots(self.cell, electrical=electrical, distance=distance, more_conductances=more_conductances,
                               time=time, dt=dt, dt_func=dt_func)
        start_lines = ['# generated by Neuron_analysis_tool and Yoni Leibner',
                       '# please note this tool is prilimanary and needed to be checked',
                       '# values in each line are:',
                       '# id,type,x,y,z,r,pid']
        swc_out = open(save_name, 'w')
        for line in start_lines:
            swc_out.write(line + '\n')

        self.save_morph_to_swc_helper(f=swc_out, all_points_3d=all_points_3d, sec=self.cell.soma[0], counter=0, parent_counter=-1)

    def plot_morph(self, ax=None, seg_to_indicate_dict = {}, diam_factor=None, diam_const=1, diam_min=0, sec_to_change=None,
                   ignore_sections=[], theta=0, scale=0, scale_text=True, ignore_soma=True, distance=None,
                   electrical=False, time=None, dt=1, more_conductances_=None, colors=None,
                   dt_func= lambda x: np.mean(x)):
        """
        plot the morphology
        :param ax: the ax to plot on
        :param seg_to_indicate_dict: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param sec_to_change:sec to change there location. the cange is only to the plot and not the model
        :param ignore_sections: section to not plot
        :param theta:the rotation to the plot in degrees [0-360]
        :param scale:if this is>0 it will add a scale bar with length scale
        :param scale_text: the text to put next to the scalebar
        :param ignore_soma: True in order to plot the soma as a circle and not a line
        :param distance: a distance for the segments
                        (if not givin its uses the default distance from the init more_conductances function)
        :param electrical: if True it will change the distance units from um to lamda
        :param time: time to plot from the more_conductances
        :param dt: the dt to use in the more_conductances
        :param more_conductances_: recorded more_conductances that are diffrent from the init
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        if colors is None:
            colors = self.colors
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        fig = plt.gca().figure
        ax, lines, segs = plot_morph(self.cell, color_func=colors.get_seg_color, add_nums=False,
                                                     seg_to_indicate=seg_to_indicate_dict,
                                                     ax=ax, diam_factor=diam_factor, diam_const=diam_const,
                                                     diam_min=diam_min, sec_to_change=sec_to_change,
                                                     theta=theta, ignore_sections=ignore_sections,
                                                     ignore_soma=ignore_soma,
                                                     distance=distance,
                                                     more_conductances=self.more_conductances if more_conductances_ is None else more_conductances_,
                                                     electrical=electrical,
                                                     time=time, dt=dt, dt_func=dt_func
                                                     )
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        if not scale == 0:
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + scale], color='k')
            if scale_text:
                text_scale = str(scale)+' '+MICRO+'m'
                y_pos = y_lim[0]+scale/3.0
                text_scale_size=10
                ax.annotate(text_scale,
                            xy=(x_lim[0], y_pos), xycoords='data', size=text_scale_size,
                            xytext=(-text_scale_size - 2, -text_scale_size / 2), textcoords='offset points', rotation=90)


        return ax, lines, segs

    def plot_morph_with_values(self, seg_val_dict, ax=None, seg_to_indicate_dict={}, diam_factor=None,
                               sec_to_change=None, ignore_sections=[], theta=0, scale=500, scale_text=True,
                               cmap=plt.cm.turbo,
                               plot_color_bar=True, bounds=None, ignore_soma=True,
                               color_bar_kwarts=dict(shrink=0.6),
                               colors=None, distance=None, electrical=False, time=None, dt=1,
                               more_conductances_=None, dt_func=lambda x: np.mean(x), bar_name=''):
        """
        plot the morph with colorcode
        :param seg_val_dict: dictionary of {section name: {segment name: numarical value}}
        :param ax: the ax to plot on
        :param seg_to_indicate_dict: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param sec_to_change:sec to change there location. the cange is only to the plot and not the model
        :param ignore_sections: section to not plot
        :param theta:the rotation to the plot in degrees [0-360]
        :param scale:if this is>0 it will add a scale bar with length scale
        :param scale_text: the text to put next to the scalebar
        :param ignore_soma: True in order to plot the soma as a circle and not a line
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param electrical: if True it will change the distance units from um to lamda
        :param time: time to plot from the more_conductances
        :param dt: the dt to use in the more_conductances
        :param more_conductances_: recorded more_conductances that are diffrent from the init
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param dt_func: function for dt in the more_conductances
        :param plot_color_bar: if to plot a color bar
        :param bounds: the bounds for the color map. default is the max and min value givin in seg_val_dict.
        :param cmap: cmap to use, a matplotlib cmap.
        :param color_bar_kwarts: kwarts for the color bar
        :param bar_name: text for the color-bar
        :return:
        """
        if self.type.startswith('Rall_tree'):
            ignore_soma=False
        if ax is None:
            ax = plt.gca()
        fig = plt.gca().figure
        ax.set_aspect('equal', adjustable='box')
        if colors is None:
            colors = color_func_norm(value_dict=seg_val_dict, bounds=bounds, cmap=cmap)

        ax, lines, segs =  self.plot_morph(ax=ax, seg_to_indicate_dict = seg_to_indicate_dict, diam_factor=diam_factor,
                                           sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta,
                                           scale=scale, scale_text=scale_text, ignore_soma=ignore_soma, distance=distance,
                                           electrical=electrical, time=time, dt=dt,
                                           more_conductances_=more_conductances_, colors=colors, dt_func=dt_func)

        if plot_color_bar:
            im = plt.cm.ScalarMappable(norm=colors.norm, cmap=colors.cmap)
            color_bar = plt.colorbar(im, ax=ax, **color_bar_kwarts )
            if bounds is not None:
                ticks = color_bar.get_ticks()
                ticks[0]=bounds[0]
                ticks[-1]=bounds[-1]
                color_bar.set_ticks(ticks)
            cax = color_bar.ax
            fig.add_axes(cax)
            cax.set_title(bar_name)
        else:
            color_bar = None
            cax = None
        return ax, cax, colors, lines, segs

    def plot_morph_with_value_func(self, func=seg_Rin_func, run_time=0, ax=None, seg_to_indicate_dict={},
                                   diam_factor=None,
                                   sec_to_change=None, ignore_sections=[], theta=0, scale=0, scale_text=True,
                                   cmap=plt.cm.turbo,
                                   plot_color_bar=True, bounds=None, ignore_soma=True,
                                   color_bar_kwarts=dict(shrink=0.6),
                                   colors=None, distance=None, electrical=False, time=None, dt=1,
                                   more_conductances_=None, dt_func=lambda x: np.mean(x), bar_name=''):
        """
        plot the morph with color code resulting of a function for each segment
        :param func: the function to run on each segment
        :param run_time: time to run befor the function
        :param seg_val_dict: dictionary of {section name: {segment name: numarical value}}
        :param ax: the ax to plot on
        :param seg_to_indicate_dict: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param sec_to_change:sec to change there location. the cange is only to the plot and not the model
        :param ignore_sections: section to not plot
        :param theta:the rotation to the plot in degrees [0-360]
        :param scale:if this is>0 it will add a scale bar with length scale
        :param scale_text: the text to put next to the scalebar
        :param ignore_soma: True in order to plot the soma as a circle and not a line
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param electrical: if True it will change the distance units from um to lamda
        :param time: time to plot from the more_conductances
        :param dt: the dt to use in the more_conductances
        :param more_conductances_: recorded more_conductances that are diffrent from the init
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param dt_func: function for dt in the more_conductances
        :param plot_color_bar: if to plot a color bar
        :param bounds: the bounds for the color map. default is the max and min value givin in seg_val_dict.
        :param cmap: cmap to use, a matplotlib cmap.
        :param color_bar_kwarts: kwarts for the color bar
        :param bar_name: text for the color-bar
        :return:
        """
        if colors is None:
            if run_time>0:
                h.tstop=run_time
                h.run()
            value_dict = dict()
            for part in self.parts_dict:
                for seg in tqdm(self.parts_dict[part], desc=part):
                    if sec_name(seg.sec) not in value_dict:
                        value_dict[sec_name(seg.sec)]=dict()
                    value_dict[sec_name(seg.sec)][seg_name(seg)] = func(seg)
            colors = color_func_norm(value_dict=value_dict, bounds=bounds, cmap=cmap)

        return self.plot_morph_with_values({}, ax=ax, seg_to_indicate_dict=seg_to_indicate_dict,
                                           diam_factor=diam_factor,
                                           sec_to_change=sec_to_change, ignore_sections=ignore_sections, theta=theta,
                                           scale=scale, scale_text=scale_text, cmap=cmap,
                                           plot_color_bar=plot_color_bar, bounds=bounds, ignore_soma=ignore_soma,
                                           color_bar_kwarts=color_bar_kwarts,
                                           colors=colors, distance=distance, electrical=electrical, time=time, dt=dt,
                                           more_conductances_=more_conductances_, dt_func=dt_func, bar_name=bar_name)

    def plot_dendogram(self, start_seg = None ,ax=None, segs_to_indecate = dict(), plot_legend=True, ignore_sections=[],
                       electrical=False, diam_factor=None, distance=None, colors=None,
                       BRANCH_SPACE_=None, dt_func= lambda x: np.mean(x)):
        """
        plot dendorgam of the tree from a givin start point
        :param start_seg: the seg to start the dendogram from
        :param ax: the axes to plot on
        :param segs_to_indecate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param plot_legend: if True it will dd a legend with the parts_dict keys and colors
        :param ignore_sections: section not to plot
        :param electrical: if True it will plot the dendogram in electrical units
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param BRANCH_SPACE_: spacing between the terminal branches (makes the dendogram more/less sparse)
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        if colors is None:
            colors = self.colors
        if start_seg is None:
            start_sec = list(self.cell.soma[0])
            start_seg = start_sec[len(start_sec)//2]
        max_y, x_pos, lines, segs = plot_dendogram(self.cell, start_seg, self.more_conductances, colors, ax=ax,
                                                   plot_legend=plot_legend, ignore_sections=ignore_sections,
                                                   segs_to_indecate=segs_to_indecate, electrical=electrical,
                                                   diam_factor=diam_factor, distance=distance,
                                                   BRANCH_SPACE_=BRANCH_SPACE_, dt_func=dt_func)
        return ax, x_pos, lines, segs

    def plot_dendogram_with_values(self, seg_val_dict, start_seg = None ,ax=None, segs_to_indecate = dict(),
                                   ignore_sections=[], electrical=False, diam_factor=None, distance=None, bounds=None,
                                   cmap=plt.cm.turbo, plot_color_bar=True, color_bar_kwarts=dict(shrink=0.6),
                                   colors=None, BRANCH_SPACE_=None, dt_func= lambda x: np.mean(x), bar_name=''):
        """
        plot the dendorgam with colorcode
        :param seg_val_dict:
        :param start_seg: the seg to start the dendogram from
        :param ax: the axes to plot on
        :param segs_to_indecate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param ignore_sections: section not to plot
        :param electrical: if True it will plot the dendogram in electrical units
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param BRANCH_SPACE_: spacing between the terminal branches (makes the dendogram more/less sparse)
        :param dt_func: function for dt in the more_conductances
        :param bounds: the bounds for the color map. default is the max and min value givin in seg_val_dict.
        :param cmap: cmap to use, a matplotlib cmap.
        :param plot_color_bar: if to plot a color bar
        :param color_bar_kwarts: kwarts for the color bar
        :return:
        """
        if start_seg is None:
            start_sec = list(self.cell.soma[0])
            start_seg = start_sec[len(start_sec)//2]
        if ax is None:
            ax = plt.gca()
        if colors is None:
            colors = color_func_norm(value_dict=seg_val_dict, bounds=bounds, cmap=cmap)

        max_y, x_pos, lines, segs = self.plot_dendogram(start_seg = start_seg ,ax=ax, plot_legend=False,
                                                        segs_to_indecate = segs_to_indecate,
                                                        ignore_sections=ignore_sections, electrical=electrical,
                                                        diam_factor=diam_factor, distance=distance, colors=colors,
                                                        BRANCH_SPACE_=BRANCH_SPACE_, dt_func=dt_func)
        if plot_color_bar:
            fig = ax.figure
            im = plt.cm.ScalarMappable(norm=colors.norm, cmap=colors.cmap)
            color_bar = plt.colorbar(im, ax=ax, **color_bar_kwarts)
            if bounds is not None:
                ticks = color_bar.get_ticks()
                ticks[0]=bounds[0]
                ticks[-1]=bounds[-1]
                color_bar.set_ticks(ticks)
            cax = color_bar.ax
            fig.add_axes(cax)
            cax.set_title(bar_name)
        else:
            color_bar = None
            cax = None

        return ax, x_pos, cax, colors, lines, segs

    def plot_dendogram_with_value_func(self, func, start_seg=None, ax=None, segs_to_indecate=dict(),
                                       ignore_sections=[], electrical=False, diam_factor=None, distance=None,
                                       bounds=None,
                                       cmap=plt.cm.turbo, plot_color_bar=True, color_bar_kwarts=dict(shrink=0.6),
                                       colors=None,
                                       BRANCH_SPACE_=None, dt_func=lambda x: np.mean(x)):
        """
        plot the morph with color code resulting of a function for each segment
        :param func: the function to run on each segment
        :param start_seg: the seg to start the dendogram from
        :param ax: the axes to plot on
        :param segs_to_indecate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param ignore_sections: section not to plot
        :param electrical: if True it will plot the dendogram in electrical units
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param colors: a givin color class that gives color to each segment via get_seg_color func
        :param BRANCH_SPACE_: spacing between the terminal branches (makes the dendogram more/less sparse)
        :param dt_func: function for dt in the more_conductances
        :param bounds: the bounds for the color map. default is the max and min value givin in seg_val_dict.
        :param cmap: cmap to use, a matplotlib cmap.
        :param plot_color_bar: if to plot a color bar
        :param color_bar_kwarts: kwarts for the color bar
        :return:
        """
        if colors is None:
            value_dict = dict()
            for part in self.parts_dict:
                for seg in tqdm(self.parts_dict[part], desc=part):
                    if seg.sec not in value_dict:
                        value_dict[sec_name(seg.sec)]=dict()
                    value_dict[sec_name(seg.sec)][seg_name(seg)] = func(seg)
            colors = color_func_norm(value_dict=value_dict, bounds=bounds, cmap=cmap)
        return self.plot_dendogram_with_values({}, start_seg=start_seg, ax=ax, segs_to_indecate=segs_to_indecate,
                                               ignore_sections=ignore_sections, electrical=electrical,
                                               diam_factor=diam_factor, distance=distance, bounds=bounds, cmap=cmap,
                                               plot_color_bar=plot_color_bar, color_bar_kwarts=color_bar_kwarts,
                                               colors=colors, BRANCH_SPACE_=BRANCH_SPACE_, dt_func=dt_func)

    def get_cable(self, start_seg=None, type='d3_2', factor_e_space=100, factor_m_space=10,
                  ignore_sections = [], distance=None, dt_func=lambda x: np.mean(x)):
        """
        calculate the cables.
        for each direction in 'sons', 'parent' and for each part from part dict and all
        :param cell: the cell neuron model
        :param start_seg: the seg to start the cable from
        :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
        :param factor_m_space: the bin size in physical units, every factor_m_space um is a bin
        :param more_conductances: more_conductances to initiate the distance if distance is None
        :param seg_dist_dict: dict in the form {direction: {seg: []}} (segs to get there electric distane
        :param part_dict: the part dict as {part name: list of segments}
        :param ignore_sections: section to skip in the cable
        :param distance: a distance for the segments
                        (if not givin its uses the default distance from the init more_conductances function)
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        if start_seg is None:
            start_sec = list(self.cell.soma[0])
            start_seg = start_sec[len(start_sec)//2]
        cable = get_cable(self.cell,
                         start_seg=start_seg,
                         factor_e_space=factor_e_space,
                         factor_m_space=factor_m_space,
                         more_conductances=self.more_conductances,
                         seg_dist_dict={},
                         part_dict=self.parts_dict,
                         ignore_sections = ignore_sections,
                         distance=distance,
                         dt_func=dt_func)[0]
        return cable['parent']['all'][type], cable['sons']['all'][type]

    def plot_cable(self, start_seg = None ,ax=None, factor_e_space=25, factor_m_space=10, segs_to_indecate= dict(),
                   ignore_sections=[], cable_type='d3_2', start_loc=0, shift=None, vertical=True, dots_size=10,
                   start_color='k', plot_legend=True, distance=None, cable_factor=1, labal_start=None,
                   return_shift=False, dt_func= lambda x: np.mean(x)):
        """
        plot a cable from a givin start point
        :param start_seg: the seg to start the cable from
        :param ax: the axes to plot on
        :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
        :param factor_m_space: the bin size in physical units, every factor_m_space um is a bin
        :param segs_to_indecate: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param ignore_sections: section to skip in the cable
        :param cable_type: type of the cable from: d3_2/electric/dist, d3_2 is the equivelent cable, electricis dA/dX and dist is dA/dx
        :param start_loc: the start location on the axes
        :param shift: a constant shift or a chosen from the data
        :param vertical: if to plot the cable vertically or horizontally
        :param dots_size: the size of the dot at the origin
        :param start_color: the color of the dot in the origin
        :param plot_legend: if True it will dd a legend with the parts_dict keys and colors
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param cable_factor: factor to add to the cable to match the axes
        :param labal_start: the label for the legend of the start seg
        :param return_shift: if return the result shift
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        if ax is None:
            ax = plt.gca()
        if start_seg is None:
            start_sec_list = list(self.cell.soma[0])
            start_seg = start_sec_list[len(start_sec_list)//2]

        return plot_cable(self.cell, start_seg=start_seg, ax=ax, cable_type=cable_type, factor_e_space=factor_e_space,
                          factor_m_space=factor_m_space, more_conductances=self.more_conductances,
                          part_dict=self.parts_dict, colors_dict=self.colors_dict, ignore_sections=ignore_sections,
                          distance=distance, segs_to_indecate=segs_to_indecate, start_loc=start_loc, vertical=vertical,
                          dots_size=dots_size, start_color=start_color, shift=shift, plot_legend=plot_legend,
                          cable_factor=cable_factor, labal_start=labal_start, return_shift=return_shift, dt_func=dt_func)


    def plot_attenuation(self, protocol=long_pulse_protocol, ax=None, seg_to_indicate_dict=dict(), start_seg =None,
                         record_to_value_func=record_to_value, norm=True, record_name='v', norm_by=None,
                         electrical=True, distance=None, records=None, ignore_sections=[],
                         dt_func= lambda x: np.mean(x), direction_dist_factors=dict(sons=1, parent=1), **kwargs):
        """
        plot the record along the distance
        :param protocol: protocol to run (if records are not givin)
        :param ax: the axes to plot on
        :param seg_to_indicate_dict: dictinary {seg: dict(size=100, color='b', alpha=0.75)}
        :param start_seg: the seg to start the plot from
        :param record_to_value_func: function to do after protocol is done
        :param norm: if to norm the records
        :param record_name: the name of the record
        :param norm_by: value for norm, if not givin its will be the value of the record at start_seg
        :param electrical: if True it will plot the dendogram in electrical units
        :param distance: a distance for the segments (if not givin its uses the default distance from the init more_conductances function)
        :param records: the records from a pre-calc protocol
        :param ignore_sections: section not to plot
        :param dt_func: function for dt in the more_conductances
        :param kwargs: extra kwargs for the each line in matplotlib such lw, ls etc.
        :return:
        """
        ax, norm_by, lines, segs, records = plot_attenuation(cell=self.cell, start_seg=start_seg, protocol=protocol,
                                                             more_conductances=self.more_conductances,
                                                             color_func=self.colors, ax=ax, record_name=record_name,
                                                             cut_start_ms=None,  record_to_value=record_to_value_func,
                                                             norm_by=norm_by,  norm=norm, electrical=electrical,
                                                             seg_to_indicate=seg_to_indicate_dict, distance=distance,
                                                             records=records, ignore_sections=ignore_sections,
                                                             dt_func=dt_func, direction_dist_factors=direction_dist_factors,
                                                             **kwargs)

        ax.set_yscale('log')
        ax.set_xlabel('distance from origin (x / '+LAMDA+')')
        if norm or (not norm_by==1.0):
            ax.set_ylabel('V(x)/V(0)')
        else:
            ax.set_ylabel('attanuation (mV)')

        return ax, norm_by, lines, segs, records

    def create_card(self, start_seg=None, theta=0, scale=500, factor_e_space=50, cable_type='d3_2', diam_factor=None,
                    plot_legend=False, start_color='green', start_dots_size=50, cable_factor=1, dt_func= lambda x: np.mean(x),
                    protocol1=long_pulse_protocol, protocol2=short_pulse_protocol, **kwargs):
        """
        create a viewing card from a segment prospective
        the card incluses:
            morphology with the segment indication
            dendogram from the segment prespective
            attenuation in response to 2 protocols
        :param start_seg: the seg to start the plot from
        :param theta: the rotation to the plot in degrees [0-360]
        :param scale: if this is>0 it will add a scale bar with length scale
        :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
        :param cable_type: type of the cable from: d3_2/electric/dist, d3_2 is the equivelent cable, electricis dA/dX and dist is dA/dx
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param plot_legend: if True it will dd a legend with the parts_dict keys and colors
        :param start_dots_size: the size of the dot at the origin
        :param start_color: the color of the dot in the origin
        :param cable_factor: factor for the cable width
        :param dt_func: function for dt in the more_conductances
        :param protocol1: the first protocol for attenuation plot. default=long_pulse_protocol (1 sec current ingection of 0.1 pA)
        :param protocol2: the second protocol for attenuation plot. default=short_pulse_protocol (2 ms current ingection of 0.1 pA)
        :param kwargs: extra kwargs for attenuation (for the each line in matplotlib such lw, ls etc.)
        :return:
        """
        fig, ax = plt.subplots(1, 4, figsize=(12, 3), gridspec_kw={'width_ratios': [0.5, 1.5, 1, 1]})
        plt.subplots_adjust(wspace=0.5)
        if start_seg is None:
            start_seg = self.cell.soma[0]
            start_seg = list(start_seg)
            start_seg = start_seg[len(start_seg)//2]
        seg_to_indicate_dict = {start_seg: dict(size=start_dots_size, color=start_color, alpha=1)}
        distance = Distance(self.cell, self.more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)

        _,_,_ = self.plot_morph(ax=ax[0], theta=theta, seg_to_indicate_dict=seg_to_indicate_dict, scale=scale,
                                diam_factor=diam_factor, ignore_soma=not self.type.startswith('Rall_tree'),
                                distance=distance)
        _, x_pos, _, _ = self.plot_dendogram(start_seg=start_seg, ax=ax[1], electrical=True, plot_legend=False,
                                             segs_to_indecate=seg_to_indicate_dict, distance=distance)
        _, shift = self.plot_cable(start_seg=start_seg, ax=ax[1], factor_e_space=factor_e_space, cable_type=cable_type,
                        plot_legend=plot_legend, start_loc=x_pos+15, start_color=start_color, dots_size=start_dots_size,
                        distance=distance, cable_factor=cable_factor, return_shift=True)
        for a in ax[1:]:
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
        ax[1].set_axis_off()
        x_lim = ax[1].get_xlim()
        y_lim = ax[1].get_ylim()

        ax[1].plot([shift-5*cable_factor, shift+5*cable_factor], [y_lim[0], y_lim[0]], color='k')
        text_size = 10
        if cable_type=='d3_2':
            x_scale_text = '10 '+MICRO+'m'
            ax[1].annotate(x_scale_text,
                        xy=(shift, y_lim[0]), xycoords='data', size=text_size,
                        xytext=(-text_size, -text_size - 2), textcoords='offset points')
        else:
            x_scale_text = '10 '+MICRO+'$m^{2}$'
            ax[1].annotate(x_scale_text,
                        xy=(shift, y_lim[0]), xycoords='data', size=text_size,
                        xytext=(-text_size, -text_size - 2), textcoords='offset points')

        ax[1].plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + 0.5], color='k')
        y_scale_text = '0.5 '+LAMDA+''
        y_pos = y_lim[0] + 0.5 / 3.0
        ax[1].annotate(y_scale_text,
                    xy=(x_lim[0], y_pos), xycoords='data', size=text_size,
                    xytext=(-text_size - 2, -text_size / 2), textcoords='offset points', rotation=90)
        ax[1].set_ylabel('')

        self.plot_attenuation(start_seg=start_seg, ax=ax[2], protocol=protocol1,
                              seg_to_indicate_dict=seg_to_indicate_dict, distance=distance, ** kwargs)
        self.plot_attenuation(start_seg=start_seg, ax=ax[3], protocol=protocol2,
                              seg_to_indicate_dict=seg_to_indicate_dict, distance=distance, ** kwargs)
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

    def create_morpho_card(self, start_seg=None, theta=0, scale=500, factor_e_space=50, cable_type='d3_2',
                           diam_factor=None, plot_legend=False, start_color='green', start_dots_size=50,
                           dt_func= lambda x: np.mean(x)):
        """
        create a morphology viewing card from a segment prospective
        the card incluses:
            morphology with the segment indication
            dendogram from the segment prespective
        :param start_seg: the seg to start the plot from
        :param theta: the rotation to the plot in degrees [0-360]
        :param scale: if this is>0 it will add a scale bar with length scale
        :param factor_e_space: the bin size in electrical units, every 1/factor_e_space lamda is a bin
        :param cable_type: type of the cable from: d3_2/electric/dist, d3_2 is the equivelent cable, electricis dA/dX and dist is dA/dx
        :param diam_factor: factor to change the ploting diam so the line width will be seg.diam*diam_factor. default=1
        :param plot_legend: if True it will dd a legend with the parts_dict keys and colors
        :param start_dots_size: the size of the dot at the origin
        :param start_color: the color of the dot in the origin
        :param dt_func: function for dt in the more_conductances
        :return:
        """
        fig, ax = plt.subplots(1, 3, figsize=(12, 3), gridspec_kw={'width_ratios': [0.5, 1.5, 1]})
        plt.subplots_adjust(wspace=0.5)
        if start_seg is None:
            start_seg = self.cell.soma[0]
            start_seg = list(start_seg)
            start_seg = start_seg[len(start_seg)//2]
        seg_to_indicate_dict = {start_seg: dict(size=start_dots_size, color=start_color, alpha=1)}
        distance = Distance(self.cell, self.more_conductances, dt_func=dt_func)
        distance.compute(start_seg=start_seg)

        self.plot_morph(ax=ax[0], theta=theta, seg_to_indicate_dict=seg_to_indicate_dict, scale=scale,
                        diam_factor=diam_factor, ignore_soma=not self.type.startswith('Rall_tree'), distance=distance)
        _, x_pos, _, _ = self.plot_dendogram(start_seg=start_seg, ax=ax[1], electrical=True, plot_legend=False,
                                             segs_to_indecate=seg_to_indicate_dict, distance=distance)
        self.plot_cable(start_seg=start_seg, ax=ax[2], factor_e_space=factor_e_space, cable_type=cable_type,
                        plot_legend=plot_legend, start_loc=0, start_color=start_color, dots_size=start_dots_size,
                        distance=distance)
        for a in ax[1:]:
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
        ax[1].set_axis_off()
        x_lim = ax[1].get_xlim()
        y_lim = ax[1].get_ylim()

        text_size=10
        ax[1].plot([x_lim[0], x_lim[0]], [y_lim[0], y_lim[0] + 0.5], color='k')
        y_scale_text = '0.5 '+LAMDA+''
        y_pos = y_lim[0] + 0.5 / 3.0
        ax[1].annotate(y_scale_text,
                    xy=(x_lim[0], y_pos), xycoords='data', size=text_size,
                    xytext=(-text_size - 2, -text_size / 2), textcoords='offset points', rotation=90)
        ax[1].set_ylabel('')

        ax[2].set_ylabel('distance ('+LAMDA+')')
        ax[1].set_ylabel('')

        ax[0].set_title('morphology')
        ax[1].set_title('dendogram')
        ax[2].set_title('$d^{3/2}$ equivalent cable' if cable_type=='d3_2' else cable_type)
        return fig, ax

    def record_protocol(self, protocol=spike_protocol, cut_start_ms=None, record_names=['v'], start_seg=None,
                        compute_more_condunctances = False):
        """
        run and record a protocol
        :param protocol: protocol to run
        :param cut_start_ms: the start time to cut from the starting of the protocol (in ms)
        :param record_names: the names of variable to record
        :param start_seg: the start_seg to pass for the rotocol
        :param compute_more_condunctances: if to record all the conductances in the neuron model whie runing
        the protocol(needed for dancing video plots)
        :return: the multi records, and extra stuff (more_conductances/draw_funcs etc.)
        """
        if start_seg is None:
            start_seg = self.cell.soma[0](0.5)
        records = multi_record_all(self.cell, record_names=record_names)
        if compute_more_condunctances:
            more_conductances_ = more_conductances(self.cell, is_resting=False, protocol=None)
        delay, extra = protocol(self.cell, start_seg)
        if cut_start_ms is None:
            cut_start_ms = max(delay, 0)
        extraction_func = lambda x: np.array(x)[int(cut_start_ms / h.dt):]
        records.extract(extraction_func )
        if compute_more_condunctances:
            more_conductances_.extract(extraction_func)
            extra['more_conductances']=more_conductances_
        if 'draw_funcs' not in extra:
            extra['draw_funcs'] = []
        return records, extra#['draw_funcs'] if 'draw_funcs' in extra else []

    def save_movie_from_rec(self, fig,slow_down_factor=1, plot_kwargs=[], func_before_run=[lambda: plt.tight_layout()],
                            func_during_run=[], save_to='', clip_name='clip', fps=None, threads=16, preset='medium', auto_del=False): # ultrafast
        """
        create and save a video
        :param fig: the figure of the movie
        :param slow_down_factor: slowdown factor in the movie
        :param plot_kwargs: list of dictionarys each is used to create a single axes in the movie (see Video script)
        :param func_before_run: function to run after axes are initialized before the video creation start
        :param func_during_run: function to run on every dt in the video
        :param save_to: the folder to save the video to
        :param clip_name: the clip name for saving
        :param fps: the fsp of the video, note that fps*slow_down_factor is the number of frames that will be generated
        :param threads: number of threds to use in saving
        :param preset: the preset to use in frame generation
        :return:
        """
        if not clip_name.endswith('.mp4'):
            clip_name += '.mp4'
        if os.path.isfile(os.path.join(save_to, clip_name)):
            if auto_del:
                os.remove(os.path.join(save_to, clip_name))
            else:
                to_del=input('the file name is taken, delete the file? (Y, N)')
                if to_del=='Y':
                    os.remove(os.path.join(save_to, clip_name))
                else:
                    print('exit without running')
                    return
        save_movie_from_rec(analyzer=self, fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs,
                            func_before_run=func_before_run, func_during_run=func_during_run, save_to=save_to,
                            clip_name=clip_name, fps=fps, threads=threads, preset=preset)

    def create_movie_from_rec(self, fig, slow_down_factor=1, plot_kwargs=[],
                              func_before_run=[lambda: plt.tight_layout()], func_during_run=[]):
        """

        :param fig: the figure of the movie
        :param slow_down_factor: slowdown factor in the movie
        :param plot_kwargs: list of dictionarys each is used to create a single axes in the movie (see Video script)
        :param func_before_run: function to run after axes are initialized before the video creation start
        :param func_during_run: function to run on every dt in the video
        :return:
        """
        return create_movie_from_rec(analyzer=self, fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs,
                                     func_before_run=func_before_run, func_during_run=func_during_run)

