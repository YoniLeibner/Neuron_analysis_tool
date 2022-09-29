import numpy as np

class color_func:
    def __init__(self, parts_dict, color_dict, defult_name='else'):
        self.parts_dict=parts_dict
        self.color_dict=color_dict
        self.defult_name=defult_name

    def get_seg_color(self, seg):
        for part in self.parts_dict:
            if seg in self.parts_dict[part]:
                return [self.color_dict[part], part]
        return [self.color_dict[self.defult_name], self.defult_name]

    def __call__(self, sec, parent=None, *args, **kwargs):
        if type(sec) == str:
            if sec in self.color_dict:
                return [self.color_dict[sec], sec], [1]
            return [self.color_dict[self.defult_name], self.defult_name], [1]

        colors=[]
        lengths = {self.defult_name:0}
        segs = list(sec)
        if (not sec.parentseg() is None) and (not parent == sec.parentseg().sec):
            segs = segs[::-1]
        seg_L = sec.L/sec.nseg
        for seg in segs:
            colors.append(self.get_seg_color(seg))
            part = colors[-1][1]
            if not part in lengths:
                lengths[part]=0
            lengths[part]+=seg_L

        colors_only = np.array(colors)[:,0]
        indexes = np.unique(colors_only, return_index=True)[1]
        colors = [colors[index] for index in sorted(indexes)]
        assert sum(lengths.values()) - sec.L < 1e-3, str(sum(lengths.values()))+'!='+str(sec.L)
        lengths = [lengths[part]/sec.L for [color, part] in colors]
        lengths = np.array(lengths)
        lengths /= lengths.sum()
        return colors, lengths