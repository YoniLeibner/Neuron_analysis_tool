from Neuron_analysis_tool.Analyzer import Analyzer, long_pulse_protocol
from Neuron_analysis_tool.utils import video_player
from pathlib import Path
from neuron import gui, h
import numpy as np
import os
import matplotlib.pyplot as plt

analyser = Analyzer(type='swc', morph_path='/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/temp/temp.swc')

analyser.save_morph_to_swc('/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/temp/temp_e2.swc',
                           electrical=False, distance=None, more_conductances=None, time=None, dt=1,dt_func= lambda x: np.mean(x))
