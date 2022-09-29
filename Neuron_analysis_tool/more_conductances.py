
#todo we need to generalize to all channels and allow to add protocol run and record all channels
class more_conductances():
    def __init__(self, cell, run_time=3000, record_name='gIhbar_Ih_human_linear'):
        self.name = 'Ih_check'
        self.cell=cell
        self.run_time = run_time
        self.record_name = '_ref_'+record_name
        self.run_resting()

    def record_Ih(self):
        record_dict = dict()
        for sec in self.cell.all:
            record_dict[sec] = dict()
            for i, seg in enumerate(sec):
                try:
                    record_dict[sec][seg] = h.Vector()
                    # record_dict[sec][seg].record(sec(seg.x)._ref_gIh_Ih_human_shifts_mul_add)
                    record_dict[sec][seg].record(getattr(sec(seg.x), self.record_name))
                except:
                    record_dict[sec][seg]=0 # no Ih hare
        return record_dict

    def run_resting(self):
        self.g_Ih_record_dict = self.record_Ih()
        h.tstop = self.run_time
        h.run()
        for sec in self.g_Ih_record_dict.keys():
            for seg in self.g_Ih_record_dict[sec].keys():
                if self.g_Ih_record_dict[sec][seg] == 0: continue
                self.g_Ih_record_dict[sec][seg] = np.array(self.g_Ih_record_dict[sec][seg])[-1] # stady state opening

    def cumpute(self, seg):
        sec= seg.sec
        g_total = seg.g_pas + self.g_Ih_record_dict[sec][seg]
        return 1.0/g_total


class more_conductances_fake():
    def __init__(self, cell):
        self.name='fake'

    def cumpute(self, seg):
        return 1.0/seg.g_pas
