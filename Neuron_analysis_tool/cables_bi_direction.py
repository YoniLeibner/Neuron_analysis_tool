import numpy as np
from neuron import h
import math
from copy import deepcopy

def DADX(cell, sec, cur_dist_e,
         results,
         factor_e_space=100,
         factor_m_space=10,
         done_sections=set(),
         go_to_parent=True,
         more_conductions={},
         seg_dist={},
         cross_dist_dict={},
         part_dict=dict()):

    if sec in done_sections:
        return seg_dist, cross_dist_dict, results
    length = sec.L / 10000.0 / sec.nseg

    for seg in sec:
        RM_total = more_conductions.cumpute(seg)
        lamda = np.sqrt((RM_total / sec.Ra) * ((seg.diam / 10000.0) / 4.0))
        for dist_dict in cross_dist_dict:
            if (dist_dict['thresh'] < cur_dist_e) and (dist_dict['thresh'] > cur_dist_e + length / lamda):
                dist_dict['segs'].append(seg)
        cur_dist_e += length / lamda
        dist_e = int(math.ceil(cur_dist_e * factor_e_space))
        dist_m = int(math.ceil(h.distance(sec(seg.x)) / factor_m_space))

        if seg in seg_dist.keys():
            seg_dist[seg].append({'dist_m': h.distance(sec(seg.x)), 'dist_e': cur_dist_e})
        results['all']['dist'][0][dist_m] += seg.area()
        results['all']['electric'][0][dist_e] += seg.area()
        results['all']['d3_2'][0][dist_e] += seg.diam**(3.0 / 2.0)  # Rall 3/2 pow roll
        for part in part_dict.keys():
            if seg in part_dict[part]:
                results[part]['dist'][0][dist_m] += seg.area()
                results[part]['electric'][0][dist_e] += seg.area()
                results[part]['d3_2'][0][dist_e] += seg.diam**(3.0 / 2.0)   # Rall 3/2 pow roll
        h.pop_section()
    done_sections.add(sec)
    if not h.SectionRef(sec=sec).nchild() == 0:
        for i in range(int(h.SectionRef(sec=sec).nchild())):
            seg_dist, cross_dist_dict, results = DADX(cell,
                                             h.SectionRef(sec=sec).child[i],
                                             cur_dist_e,
                                             results,
                                             factor_e_space=factor_e_space,
                                             factor_m_space=factor_m_space,
                                             done_sections=done_sections,
                                             go_to_parent=go_to_parent,
                                             more_conductions=more_conductions,
                                             seg_dist=seg_dist,
                                             cross_dist_dict=cross_dist_dict,
                                             part_dict=part_dict)
    if go_to_parent:
        s_ref = h.SectionRef(sec=sec)
        if s_ref.has_parent():
            parent = s_ref.parent
            seg_dist, cross_dist_dict, results = DADX(cell,
                                             parent,
                                             cur_dist_e,
                                             results,
                                             factor_e_space=factor_e_space,
                                             factor_m_space=factor_m_space,
                                             done_sections=done_sections,
                                             go_to_parent=go_to_parent,
                                             more_conductions=more_conductions,
                                             seg_dist=seg_dist,
                                             cross_dist_dict=cross_dist_dict,
                                             part_dict=part_dict)

    return seg_dist, cross_dist_dict, results


def start_DADX(cell, sec, results,
               factor_e_space=100,
               factor_m_space=10,
               x_start=0.5,
               more_conductions={},
               seg_dist={},
               cross_dist_dict=[],
               part_dict=dict(),
               do_sons = True):
    length = sec.L / 10000.0 / sec.nseg
    # cur_dist_e = [0, 0]  # for both directions
    cur_dist_e = 0
    segs = [seg for seg in sec]
    if do_sons :
        for seg in segs:
            if seg.x > x_start:
                RM_total = more_conductions.cumpute(seg)
                lamda = np.sqrt((RM_total / sec.Ra) * ((seg.diam / 10000.0) / 4.0))
                for dist_dict in cross_dist_dict:
                    if (dist_dict['thresh'] < cur_dist_e) and (dist_dict['thresh'] > cur_dist_e + length / lamda):
                        dist_dict['segs'].append(seg)
                cur_dist_e += length / lamda
                dist_e = int(math.ceil(cur_dist_e * factor_e_space))
                dist_m = int(math.ceil(h.distance(sec(seg.x)) / factor_m_space))
                if seg in seg_dist.keys():
                    seg_dist[seg].append({'dist_m': h.distance(sec(seg.x)), 'dist_e': cur_dist_e})

                results['all']['dist'][0][dist_m] += seg.area()
                results['all']['electric'][0][dist_e] += seg.area()
                results['all']['d3_2'][0][dist_e] += seg.diam**(3.0 / 2.0)   # Rall 3/2 pow roll
                for part in part_dict.keys():
                    if seg in part_dict[part]:
                        results[part]['dist'][0][dist_m] += seg.area()
                        results[part]['electric'][0][dist_e] += seg.area()
                        results[part]['d3_2'][0][dist_e] += seg.diam**(3.0 / 2.0)  # Rall 3/2 pow roll

        h.pop_section()
        sons = list(sec.children())
        h.pop_section()
    else:
        for seg in segs[::-1]:
            if seg.x < x_start:
                RM_total = more_conductions.cumpute(seg)
                lamda = np.sqrt((RM_total / sec.Ra) * ((seg.diam / 10000.0) / 4.0))
                for dist_dict in cross_dist_dict:
                    if (dist_dict['thresh'] < cur_dist_e) and (dist_dict['thresh'] > cur_dist_e + length / lamda):
                        dist_dict['segs'].append(seg)
                cur_dist_e += length / lamda
                dist_e = int(math.ceil(cur_dist_e * factor_e_space))
                dist_m = int(math.ceil(h.distance(sec(seg.x)) / factor_m_space))
                if seg in seg_dist.keys():
                    seg_dist[seg].append({'dist_m': h.distance(sec(seg.x)), 'dist_e': cur_dist_e})
                results['all']['dist'][0][dist_m] += seg.area()
                results['all']['electric'][0][dist_e] += seg.area()
                results['all']['d3_2'][0][dist_e] += pow(seg.diam, 3.0 / 2.0)  # Rall 3/2 pow roll
                for part in part_dict.keys():
                    if seg in part_dict[part]:
                        results[part]['dist'][0][dist_m] += seg.area()
                        results[part]['electric'][0][dist_e] += seg.area()
                        results[part]['d3_2'][0][dist_e] += pow(seg.diam, 3.0 / 2.0)  # Rall 3/2 pow roll
            h.pop_section()
        # s_ref = h.SectionRef(sec=sec)
        sons = []
        if not sec.parentseg() is None:
            sons.append(sec.parentseg().sec)
        h.pop_section()

    return cur_dist_e, sons, seg_dist, cross_dist_dict, results


def get_cable(cell,
              factor_e_space=100,
              factor_m_space=10,
              start_section=None,
              x_start=0.5,
              more_conductions={},
              seg_dist_dict={},
              cross_dist_dict=[{'thresh': 0.25, 'segs': []}],
              part_dict=dict(),
              ignore_sections = []):

    total_res = dict()
    # seg_dist_dict = {'sons':seg_dist, 'parent':deepcopy(seg_dist)}
    for direction in ['sons', 'parent']:
        results = dict(all=dict(dist=np.zeros([1, 1000]), electric= np.zeros([1, 1000]), d3_2=np.zeros([1, 1000])))
        for part in part_dict.keys():
            results[part] = dict(dist=np.zeros([1, 1000]), electric= np.zeros([1, 1000]), d3_2=np.zeros([1, 1000]))
        done_sections = set()
        if start_section is None:
            start_section = cell.soma[0]
            go_to_parent = False
            h.distance(0, 0.5, sec=start_section)
            cur_e_dist=0
            section_to_continue = list(h.SectionRef(sec=start_section).child)
        else:
            go_to_parent = True
            h.distance(0, x_start, sec=start_section)
            cur_e_dist, section_to_continue, seg_dist_dict[direction], cross_dist_dict, results = start_DADX(cell, start_section,
                                                                                 results,
                                                                                 factor_e_space=factor_e_space,
                                                                                 factor_m_space=factor_m_space,
                                                                                 x_start=x_start,
                                                                                 more_conductions=more_conductions,
                                                                                 seg_dist=seg_dist_dict[direction],
                                                                                 cross_dist_dict=cross_dist_dict,
                                                                                 part_dict=part_dict,
                                                                                 do_sons=direction == 'sons')
        if start_section == cell.soma[0]:
            go_to_parent=False
            section_to_continue = []
            for sec in start_section.children():
                if direction == 'sons' and sec in cell.apic:
                    section_to_continue.append(sec)
                if direction == 'parent' and sec not in cell.apic:
                    section_to_continue.append(sec)


        done_sections.add(start_section)
        for current_sec in section_to_continue:
                # for current_sec in current_sections:
                if current_sec in ignore_sections:
                    print('skiping section ', current_sec)
                    continue

                # if len(cell.axon)>0 and not current_sec in cell.axon:
                else:
                    seg_dist_dict[direction], cross_dist_dict, results = DADX(cell, current_sec, cur_e_dist,
                                                     results,
                                                     factor_e_space=factor_e_space,
                                                     factor_m_space=factor_m_space,
                                                     done_sections=done_sections,
                                                     go_to_parent=go_to_parent,
                                                     more_conductions=more_conductions,
                                                     seg_dist=seg_dist_dict[direction],
                                                     cross_dist_dict=cross_dist_dict,
                                                     part_dict=part_dict)
        for key in results:
            results[key]['d3_2'] = np.power(results[key]['d3_2'], 2.0 / 3.0)
        total_res[direction] = deepcopy(results)
    return total_res, seg_dist_dict, cross_dist_dict

