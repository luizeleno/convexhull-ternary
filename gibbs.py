'''
    Auxiliary functions for the Gibbs triangle 
    and tooltips for the plots
'''

import numpy as np
import mplcursors

sin60, cos60 = .5 * np.sqrt(3), .5

# transform points to Gibbs triangle
def transform_coords(x1, x2):
    xa, xb = np.array(x1), np.array(x2)
    xt = xa + xb * cos60
    yt = xb * sin60
    return xt, yt

# DRAWING GIBBS TRIANGLE
def triangle(ax, A='A', B='B', C='C'):

    # plt.axis('off')
    ax.set_axis_off()

    ax.plot([0, 1, cos60, 0], [0, 0, sin60, 0], color='k')

    ax.text(0, -.04, A, fontsize=20, ha='center', transform=ax.transAxes)
    ax.text(1, -.04, B, fontsize=20, ha='center', transform=ax.transAxes)
    ax.text(cos60, .12+sin60, C, fontsize=20, ha='center', transform=ax.transAxes)

    ax.text(.5, -.05, f"$x_\\mathrm{{{B}}}$", fontsize=14, ha='center', transform=ax.transAxes)
    ax.text(.18, .5, f"$x_\\mathrm{{{C}}}$", rotation=60, fontsize=14, ha='center', transform=ax.transAxes)

    for i in range(1, 10):
        ax.text(i/10, -.05, f'0.{i}', fontsize=12, ha='center')

        linev = ax.plot([i/10,i/10], [0, -.01], color='k')

        xc, yc = transform_coords([0, -.012], [i/10, i/10])
        ax.text(xc[1]-.03, yc[1], f'0.{i}', fontsize=12, ha='center', va='center')
        lineh = ax.plot(xc, yc, color='k')

def annotate(data, plot):
    
    def cursor_annotations(sel):
        text = f'ID: {int(data[sel.index, 4])}=({data[sel.index, 0]:.3f}, {data[sel.index, 1]:.3f})\n'
        text += f'in hull: {bool(data[sel.index, -1])}\nDistance to hull: {data[sel.index, 3]:.3f} meV'
        sel.annotation.set_text(text)
        sel.annotation.get_bbox_patch().set(fc="lightsalmon", alpha=0.9)

    crs = mplcursors.cursor(plot, hover=True)
    crs.connect("add", cursor_annotations)
