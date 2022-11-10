import numpy as np
from Neuron_analysis_tool.distance import Distance


def get_part(part_dict, seg):
    for part in part_dict.keys():
        if seg in part_dict[part]:
            return part
    raise Exception('seg have no part, '+str(seg)+', '+str(seg.sec))#this shold not bereached

def get_cable(cell,
              start_seg,
              factor_e_space=100,
              factor_m_space=10,
              more_conductances={},
              seg_dist_dict={},
              part_dict=dict(),
              ignore_sections = [],
              distance=None):

    if (distance is None) or (not distance.start_seg == start_seg):
        distance = Distance(cell, more_conductances)
        distance.compute(start_seg=start_seg)

    cross_dist_dict = dict(sons=[], parent=[])
    total_res = dict()
    for direction in seg_dist_dict:
        for seg in seg_dist_dict[direction]:
            if distance.get_part(seg) == direction:
                distance.get_mid_point(seg, electrical=True)
                seg_dist_dict[direction][seg].append(dict(dist_m=distance.get_mid_point(seg, electrical=False), dist_e=distance.get_mid_point(seg, electrical=True)))

    for direction in ['sons', 'parent']:
        total_res[direction] = dict()
        total_res[direction]['all'] = dict(dist=np.zeros(1000), electric= np.zeros(1000), d3_2 = np.zeros(1000))
        for part in part_dict.keys():
            total_res[direction][part] = dict(dist=np.zeros(1000), electric= np.zeros(1000), d3_2 = np.zeros(1000))
            total_res[direction][part]['dist'][0] = total_res[direction][part]['electric'][0] = start_seg.area()
    e_threshs = np.arange(0, 1000, 1)/factor_e_space
    m_threshs = np.arange(0, 1000, 1)*factor_m_space

    for sec in distance.distance_dict:
        if sec in ignore_sections: continue
        for seg in distance.distance_dict[sec]['segs']:
            e_dist = distance.get_mid_point(seg, electrical=True)
            m_dist = distance.get_mid_point(seg, electrical=False)
            m_idx = np.where(m_threshs<=m_dist)[0][-1]
            e_idx = np.where(e_threshs<=e_dist)[0][-1]
            direction = distance.get_part(seg)
            start_end = distance.get_start_end(seg, electrical=True)

            part = get_part(part_dict, seg)
            total_res[direction]['all']['dist'][m_idx] += seg.area()
            total_res[direction][part]['dist'][m_idx] += seg.area()
            total_res[direction]['all']['electric'][e_idx]+=seg.area()
            total_res[direction][part]['electric'][e_idx]+=seg.area()
            if not (sum(e_threshs<=start_end['start']) == sum(e_threshs<start_end['end'])): #we croosed a threshold
                e_idx_start = np.where(e_threshs <= start_end['start'])[0][-1]
                total_res[direction]['all']['d3_2'][e_idx_start] += seg.diam**1.5
                total_res[direction][part]['d3_2'][e_idx_start] += seg.diam**1.5

    for direction in total_res.keys():
        for part in total_res[direction].keys():
            total_res[direction][part]['d3_2'] = np.power(total_res[direction][part]['d3_2'], 2.0/3.0)
    return total_res, seg_dist_dict, dict()
