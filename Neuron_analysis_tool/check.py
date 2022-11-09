import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar


from matplotlib_scalebar.dimension import _Dimension, _PREFIXES_FACTORS, _LATEX_MU
LAMDA = '\u03BB'
MICRO = '\u03BC'
class LamdaSpaceDimension(_Dimension):
    def __init__(self):
        super().__init__(LAMDA)
        # for prefix, factor in _PREFIXES_FACTORS.items():
        # latexrepr = None
        #     # if prefix == "\u00b5" or prefix == "u":
        #     #     latexrepr = _LATEX_MU + "s"
        #     # self.add_units(prefix + "s", factor, latexrepr)
        # self.add_units(MICRO, 0.1, latexrepr)
        # self.add_units(LAMDA, 0.1, LAMDA)
        a=1



fig, ax = plt.subplots(1)
# ax.plot([1, 1000], [1,1000])
# ax.add_artist(
#     ScaleBar(1.0,units=LAMDA, dimension=LamdaSpaceDimension(), location="lower right", fixed_value=0.5)
# )
plt.plot([0, 0], [0, 10])
ax.annotate('Peak', xy=(2, 8), xytext=(4, 10), fontsize=12,
            arrowprops=dict(facecolor='green', shrink=0.05))
ax.add_artist(ScaleBar(10**-6, units='m', label='check'))
plt.show()