#########################################################
#
# author: Yoni Leibner
# description: preset protocols for analyzer
#                    somatic  long pulse
#                    somatic short pulse
#                    somatic simulated spikes (10 spikes)
# date of modification: 16.11.2022
#
#########################################################

from neuron import h, gui
import os
import numpy as np


def resting_protocol(cell, start_seg=None):
    h.tstop = 500.0
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 400, {}

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
    return delay, {}

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
    return delay, {},

def spike_protocol(cell, start_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt=spike_data.T[0][1]-spike_data.T[0][0]
    v = spike_data.T[1]
    V = np.concatenate([np.zeros(int(1000.0/dt))+spike_data.T[1][0]]+[v]*10)
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)
    h.dt=dt
    h.steps_per_ms = 1.0 / h.dt
    clamp = h.SEClamp(start_seg.x, sec=start_seg.sec)
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    spike_vec.play(clamp._ref_amp1, spike_data.T[0][1]-spike_data.T[0][0])
    h.tstop = T[-1]+20
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}


def spike_protocol2(cell, start_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt=spike_data.T[0][1]-spike_data.T[0][0]
    v = spike_data.T[1][750:1200]
    V = np.concatenate([np.zeros(int(1000.0/dt))+spike_data.T[1][0]]+[v]*3)
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)
    h.dt=dt
    h.steps_per_ms = 1.0 / h.dt
    clamp = h.SEClamp(start_seg.x, sec=start_seg.sec)
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    spike_vec.play(clamp._ref_amp1, spike_data.T[0][1]-spike_data.T[0][0])
    h.tstop = T[-1]+20
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}


def spike_protocol3(cell, start_seg):
    spike_data = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/spike.txt'))
    dt=spike_data.T[0][1]-spike_data.T[0][0]
    v = spike_data.T[1][800:900]
    V = np.concatenate([np.zeros(int(1000.0/dt))+spike_data.T[1][0]]+[v]*3)
    T = np.arange(0, len(V), 1) * dt
    spike_vec = h.Vector(V)
    h.dt=dt
    h.steps_per_ms = 1.0 / h.dt
    clamp = h.SEClamp(start_seg.x, sec=start_seg.sec)
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    spike_vec.play(clamp._ref_amp1, spike_data.T[0][1]-spike_data.T[0][0])
    h.tstop = T[-1]+20
    h.v_init = cell.soma[0].e_pas
    h.celsius = 37
    h.run()
    return 0, {}
