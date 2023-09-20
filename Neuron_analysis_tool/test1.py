from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from Neuron_analysis_tool.utils import video_player
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt
from copy import deepcopy
from pathlib import Path

def efun(z):
    res = z / (np.exp(z) - 1)
    res2 = 1 - z / 2
    res[np.abs(z) < 1e-4] = res2[np.abs(z) < 1e-4]

    return res

def a_b(v, R, qa=9, tha=25):
    return  R * qa * efun(-(v - tha) / qa)

def new_temp(celsius=34, temp=23):
    return 2.3**((celsius- temp)/10)

def tau(v):
    return 1/new_temp()/(a_b(v, 0.02)+a_b(v, 0.002))

def inf(v):
    a=a_b(v, 0.02)+a_b(v, 0.02)
    b=a_b(v, 0.02)+a_b(v, 0.002)
    return a/(a+b)

v = np.arange(-100, 100, 0.01)

fig, ax = plt.subplots(2)
ax[0].plot(v, tau(v))
ax[1].plot(v, inf(v))
plt.show()













analyser = Analyzer(type='Rall_tree')
colors_dict  = analyser.colors_dict
colors_dict['soma']='r'
colors_dict['basal']='pink'

analyser.change_color_dict(colors_dict)
records, extra = analyser.record_protocol(cut_start_ms=1000.0, record_names=['v'])

start_seg = list(analyser.cell.soma[0])[0]
fig, ax = plt.subplots(1,2)
slow_down_factor=50
plot_kwargs=[dict(ax=ax[0], seg = start_seg,
                  records=records.all_records['v'],
                  electrical=False,
                  plot_type='morph',
                  seg_to_indicate_dict=dict(),
                  plot_color_bar=True,
                  theta=0,
                  diam_factor=1,)]

plot_kwargs.append(dict(ax=ax[1], seg = start_seg, records=records.all_records['v'],
                        electrical=True, plot_type='attenuation'))
f = lambda: plt.tight_layout()

videos_folder = 'videos/rall_model/'
video_name = 'spiking_on_rall_tree_with_attanuation'
os.makedirs(videos_folder, exist_ok=True)
analyser.save_movie_from_rec(fig=fig, slow_down_factor=slow_down_factor, plot_kwargs=plot_kwargs,
                             func_before_run=[f],
                             save_to=videos_folder, clip_name=video_name, fps=10,
                             threads=16, preset='ultrafast')

# video_player(Path.cwd(), videos_folder+video_name+'.mp4')
# video_player(Path.cwd(), videos_folder + video_name)