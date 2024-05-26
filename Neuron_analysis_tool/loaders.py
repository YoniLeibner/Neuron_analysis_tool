#########################################################
#
# author: Yoni Leibner
# description: cell loaders, we seport:
#               Hay. et al neuron
#               any ASC/swc file (you can give RM Ra and Cm
#               Rall model of degree of 5 as was created by Rall_tree.py
# date of modification: 11.05.2023
#
#########################################################
from neuron import h
import os

def open_morph(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, nl=None, seg_every=20):
    """
    load a model from a givin morphology
    :param morph_path: the morphology file path
    :param Rm: Rm for the model
    :param Ra: Ra for the model
    :param Cm: Cm for the model
    :param e_pas: e_pas for the model
    :param nl: loading case, diffrent for ASC and swc
    :param seg_every: control how many segment each section will have, result with sec, nseg=L/seg_every+1
    :return: Neuron model, parts dict, color dict
    """
    # print('open_morph: ', morph_path)
    hoc_file_name = 'allen_model.hoc'
    h.celsius = 37
    h.load_file("nrngui.hoc")
    h("objref cell, tobj")  # neuron object
    h.load_file(os.path.realpath(__file__)[:-10]+'allen_model.hoc')

    h.execute("cell = new " + hoc_file_name[:-4] + "()")  # replace?
    nl.input(morph_path)
    i3d = h.Import3d_GUI(nl, 0)
    i3d.instantiate(h.cell)
    cell = h.cell
    # parts_dict = dict(all=list())
    for sec in cell.all:
        sec.insert('pas')
        sec.nseg = int(sec.L / seg_every) + 1
        sec.e_pas = e_pas
        sec.cm = Cm
        sec.Ra = Ra
        sec.g_pas = 1.0 / Rm
    return get_parts_and_colors(cell)


def open_ASC(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, seg_every=20):
    """
    load a model from a givin ASC morphology
    :param morph_path: the morphology file path
    :param Rm: Rm for the model
    :param Ra: Ra for the model
    :param Cm: Cm for the model
    :param e_pas: e_pas for the model
    :param seg_every: control how many segment each section will have, result with sec, nseg=L/seg_every+1
    :return: Neuron model, parts dict, color dict
    :return:
    """
    h.load_file("import3d.hoc")
    return open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas, nl=h.Import3d_Neurolucida3(), seg_every=seg_every)


def open_swc(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, seg_every=20):
    """
    load a model from a givin swc morphology
    :param morph_path: the morphology file path
    :param Rm: Rm for the model
    :param Cm: Cm for the model
    :param e_pas: e_pas for the model
    :param seg_every: control how many segment each section will have, result with sec, nseg=L/seg_every+1
    :return: Neuron model, parts dict, color dict
    :return:
    """
    h.load_file("import3d.hoc")
    return open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas, nl=h.Import3d_SWC_read(), seg_every=seg_every)


def open_rall_tree(seg_every=20):
    """
    load a rall_tree (with 5 degree)
    :return: Neuron model, parts dict, color dict
    """
    morph_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/Rall_tree5.swc')
    cell, parts_dict, colors_dict = open_swc(morph_path)
    for sec in cell.soma[0].children():
        h.disconnect(sec=sec)
        sec.connect(cell.soma[0](1))
    for sec in cell.all:
        sec.nseg = int(sec.L / seg_every) + 1
    return cell, parts_dict, colors_dict

def open_L5PC(seg_every=20):
    """
    open Itay Hay L5PC model
    :return:  Neuron model, parts dict, color dict
    """
    h.load_file("import3d.hoc")
    f = os.path.dirname(os.path.realpath(__file__))
    morphology_file = os.path.join(f,"data/L5PC/cell1.asc")
    h.load_file(os.path.join(f,"data/L5PC/L5PCbiophys3.hoc"))
    h.load_file(os.path.join(f,"data/L5PC/L5PCtemplate.hoc"))
    cell = h.L5PCtemplate(morphology_file)
    for sec in cell.all:
        sec.nseg = int(sec.L / seg_every) + 1
    return get_parts_and_colors(cell)


def get_parts_and_colors(cell):
    """
    defult part for a neuron into:'soma','basal', 'apical', 'axon', 'else'
    :param cell: Neuron model
    :return: Neuron model, parts dict, color dict
    """
    parts_dict = {'soma': [], 'basal': [], 'apical': [], 'axon': [], 'else': []}
    colors_dict = {'soma': 'k', 'basal': 'r', 'apical': 'b', 'axon': 'green', 'else': 'cyan'}
    for sec in cell.all:
        if sec.nseg <= 0: continue
        sec.nseg = int(sec.L / 20) + 1
        for seg in sec:
            if sec in cell.soma:
                parts_dict['soma'].append(seg)
            elif len(cell.dend)>1 and sec in cell.dend:
                parts_dict['basal'].append(seg)
            elif len(cell.apic)>1 and sec in cell.apic:
                parts_dict['apical'].append(seg)
            elif len(cell.axon)>1 and sec in cell.axon:
                parts_dict['axon'].append(seg)
            else:
                parts_dict['else'].append(seg)
    return cell, parts_dict, colors_dict

def insert_g_total(cell):
    """
    defult part for a neuron into:'soma','basal', 'apical', 'axon', 'else'
    :param cell: Neuron model
    :return: Neuron model, parts dict, color dict
    """
    for sec in cell.all:
        sec.insert('g_total')
