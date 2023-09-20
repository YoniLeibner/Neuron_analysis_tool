#########################################################
#
# author: Yoni Leibner
# description:not working
# date of modification:
#
#########################################################

import matplotlib.pyplot as plt

class matplotlib_ploter:

    def __init__(self):
        pass

    def plot_line_2d(self, ax, x, y, color='k', alpha=1, ls='-', lw=1, zorder=1):
        ax.plot(x, y, color=color, ls=ls, lw=lw, alpha=alpha, zorder=zorder)

    def scatter_2d(self, ax, x, y, color='k', alpha=1, s=10, marker='o', zorder=1):
        ax.scatter(x, y, color=color, s=s, alpha=alpha, zorder=zorder, marker=marker)

    def plot_line_3d(self, ax, x, y, z, color='k', alpha=1, ls='-', lw=1, zorder=1):
        ax.plot(x, y, z, color=color, ls=ls, lw=lw, alpha=alpha, zorder=zorder)

    def scatter_3d(self, ax, x, y, z, color='k', alpha=1, s=10, marker='o', zorder=1):
        ax.scatter(x, y, z, color=color, s=s, alpha=alpha, zorder=zorder, marker=marker)

    def text(self, ax, x, y, text, text_size=10, rotation=0):
        ax.annotate(text,
                    xy=(x, y), xycoords='data', size=text_size,
                    xytext=(-text_size - 2, -text_size / 2), textcoords='offset points', rotation=rotation)
