from neuron import h, gui
from Neuron_analysis_tool.Analyzer import Analyzer
import matplotlib.pyplot as plt
from Neuron_analysis_tool.record import record_all
import numpy as np


h.load_file("stdrun.hoc")

def add_syn(seg, g_AMPA=0.0004, g_NMDA=0.0004):
    netstim = h.NetStim()
    netstim.interval = 10 # fr of 1
    netstim.start = 100
    netstim.noise = 0.0
    netstim.number = 1
    # AMPA part
    AMPA = h.Exp2Syn(seg.x, sec=seg.sec)
    AMPA_con= h.NetCon(netstim, AMPA)
    AMPA.e = 0
    AMPA.tau1 = 0.3
    AMPA.tau2 = 1.5
    AMPA_con.weight[0] = g_AMPA
    AMPA_con.delay = 0

    # NMDA part
    NMDA=h.NMDA(seg.x, sec=seg.sec)
    NMDA_con = h.NetCon(netstim, NMDA)
    NMDA.e = 0
    NMDA.tau_r_NMDA = 8
    NMDA.tau_d_NMDA = 35
    NMDA.n_NMDA = 0.27
    NMDA.gama_NMDA = 0.076
    NMDA_con.weight[0] = g_NMDA
    NMDA_con.delay = 0
    return [AMPA, AMPA_con], [NMDA, NMDA_con], netstim

def short_pulse_protocol(cell, start_seg):
    delay=200.0
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
    return delay, {}

def get_condactance_name(mechanisms):
    try:
        valids = []
        for param in list(mechanisms):
            if param.name()[0] == 'g' and 'bar' not in param.name():
                valids.append(param.name())
        assert len(valids)==1, 'check the more_conductances, and change the get_condactance func to seport your case'
        return valids[0]
        # return getattr(mechanisms, 'g_'+str(mechanisms))
    except:
        return ''

def get_condactance_point_prosses(point_prosses):
    if hasattr(point_prosses, 'g'):
        return True, getattr(point_prosses, 'g')
    name = point_prosses.hname().split('[')[0]
    if hasattr(point_prosses, 'g_'+name):
        return True, getattr(point_prosses, 'g_'+name)
    return False, 0

def callback(cell):
    for sec in cell.all:
        for seg in sec:
            g_total=0
            for mechanisms in seg:
                if mechanisms.name() == 'g_total': continue
                if mechanisms.is_ion(): continue
                if mechanisms.name() in ['CaDynamics_E2']: continue
                record_name = get_condactance_name(mechanisms)
                g_total += getattr(seg, record_name)
            for point_prosses in seg.point_processes():
                check, g = get_condactance_point_prosses(point_prosses)
                if check:
                    g_total+=g
            seg.g_total.g_total = g_total

# analyser = Analyzer(type='Rall_tree')
analyser = Analyzer(type='L5PC')
# colors_dict  = analyser.colors_dict
s = list(analyser.cell.apic[29])[-1]
syn = add_syn(s, g_AMPA=0.004, g_NMDA=0.004)
# Need a home for bscallback. Location is important for threads
# and local variable step method. Not important otherwise
bscallback = h.beforestep_callback(s)



bscallback.set_callback((callback, analyser.cell))
record_dict = record_all(analyser.cell, record_name='g_total_g_total')
t=h.Vector()
v=h.Vector()
g=h.Vector()
g1=h.Vector()
g2=h.Vector()
t.record(h._ref_t)
v.record(s._ref_v)
g.record(s._ref_g_total_g_total)
g1.record(syn[0][0]._ref_g)
g2.record(syn[1][0]._ref_g_NMDA)
short_pulse_protocol(cell=analyser.cell, start_seg=list(analyser.cell.soma[0])[0])


# plt.plot(t, v, label='v')
plt.plot(t, g, label='g_total')
plt.plot(t, g1, label='g_AMPA')
plt.plot(t, g2, label='g_NMDA')
plt.plot(t, np.array(g2)+np.array(g1)+s.g_pas, label='g_check')
plt.plot(t, np.array(g1)+s.g_pas, label='g_check2')
plt.legend()
plt.show()