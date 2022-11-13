import numpy as np
from neuron import h


LAMDA = '\u03BB'
MICRO = '\u03BC'

def get_segment_length_lamda(seg, more_conductances):
    """
    return the segment  e_length
    :param seg_len:
    :param RM:
    :param RA:
    :return:
    """
    sec = seg.sec
    seg_len = sec.L/sec.nseg
    d = seg.diam
    R_total = more_conductances.cumpute(seg)
    lamda = np.sqrt((R_total / sec.Ra) * (d / 10000.0) / 4.0)
    return (float(seg_len) / 10000.0) / lamda

def get_segment_length_um(seg, more_conductances):
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
    return not sec.parentseg() is None

def seg_Rin_func(seg):
    imp = h.Impedance(seg.x, sec=seg.sec)
    imp.loc(seg.x, sec=seg.sec)
    imp.compute(0, 1)
    return imp.input(seg.x, sec=seg.sec)
