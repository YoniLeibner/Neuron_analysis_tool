#########################################################
#
# author: Yoni Leibner
# description: cell loaders, we seport:
#               Hay. et al neuron
#               any ASC/swc file (you can give RM Ra and Cm
#               Rall model of degree of 5 as was created by Rall_tree.py
# date of modification: 16.11.2022
#
#########################################################

from neuron import h
import os

def open_morph(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70, nl=None):
    print('open_morph: ', morph_path)
    hoc_file_name = 'allen_model.hoc'
    h.celsius = 37
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
        sec.nseg = int(sec.L / 10) + 1
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


def open_ASC(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70):
    h.load_file("import3d.hoc")
    return open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas, nl=h.Import3d_Neurolucida3())


def open_swc(morph_path, Rm=10000.0, Ra=100, Cm=1, e_pas=-70):
    h.load_file("import3d.hoc")
    return open_morph(morph_path, Rm=Rm, Ra=Ra, Cm=Cm, e_pas=e_pas, nl=h.Import3d_SWC_read())


def open_rall_tree():
    morph_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/Rall_tree5.swc')
    # cell, parts_dict, colors_dict = self.open_swc(morph_path)
    return open_swc(morph_path)  # cell, dict(Rall_tree = parts_dict['all']), dict(Rall_tree =colors_dict['all'])


def open_L5PC():
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
