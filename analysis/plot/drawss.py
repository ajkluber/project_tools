""" Plot secondary structure icons along axes of contact map.

Description:

    Uses matplotlib's artist utilities, specifically 'patches' (i.e. patches of
color) to plot little icons to represent protein secondary structure on the
borders of the figure. This uses the Rectangle and Arrow patch.

See documentation for matplotlib.patches.Rectangle & matplotlib.patches.Arrow

Example:
To create Rectangle object with bottom left corner at data coordinates 
(xpos,ypos) with dimensions (xdim,ydim).
Retangle( (xpos,ypos), xdim, ydim,...)

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def add_secondary_struct_icons(ax,elements,seqbounds):
    """ Adds a rectangle (arrow) for beta sheet (alpha helix) along each axis. """
    for i in range(len(elements)):
        element = elements[i]
        bounds = seqbounds[i]

        if element.startswith("beta") or element.startswith("alpha"):
            struct = element.split("_")[0]
            number = element.split("_")[1]
            
            if struct == "beta":
                xpos = bounds[0]
                #ypos = -7
                ypos = -3
                xdim = bounds[1] - bounds[0]
                ydim = 2

                ## Add to x-axis
                patch1 = mpatches.Rectangle((xpos,ypos), xdim, ydim, facecolor="blue",alpha=0.7,ec='k')
                #ax.text(float(xpos)+0.3*float(xdim), ypos-4,"$\\beta %s$" % number, size=18)
                ax.text(float(xpos)+0.3*float(xdim), ypos-3,"$\\beta %s$" % number, size=18)

                ## Add to y-axis
                patch2 = mpatches.Rectangle((ypos,xpos), ydim, xdim, facecolor="blue",alpha=0.7,ec='k')
                #ax.text(ypos-7, float(xpos)+0.2*float(xdim),"$\\beta %s$" % number, size=18)
                ax.text(ypos-3, float(xpos)+0.2*float(xdim),"$\\beta %s$" % number, size=18)
            elif struct == "alpha":
                xpos = bounds[0]
                #ypos = -6
                ypos = -3
                xdim = bounds[1] - bounds[0]

                ## Add to x-axis
                patch1 = mpatches.Arrow(xpos, ypos, xdim,0, width=8.2,color='red',alpha=0.7,ec='k')
                #ax.text(float(xpos)+0.3*float(xdim), ypos-5,"$\\alpha %s$" % number, size=18)
                ax.text(float(xpos)+0.3*float(xdim), ypos-5,"$\\alpha %s$" % number, size=18)

                ## Add to y-axis
                patch2 = mpatches.Arrow(ypos, xpos, 0,xdim, width=8.2,color='red',alpha=0.7,ec='k')
                #ax.text(ypos-8, float(xpos)+0.2*float(xdim),"$\\alpha %s$" % number, size=18)
                ax.text(ypos-8, float(xpos)+0.2*float(xdim),"$\\alpha %s$" % number, size=18)

        ## plot patches to the axis
        ## prevent it being covered by figure background.
        temp = ax.add_patch(patch1)
        temp.set_clip_on(False)
        temp = ax.add_patch(patch2)
        temp.set_clip_on(False)


if __name__ == "__main__":
    Qref = np.loadtxt("1RIS/Qref_cryst.dat")

    elements = ["beta_1","alpha_1","beta_2","beta_3","alpha_2","beta_4"]
    seqbounds = [(2,10),(16,32),(36,52),(55,67),(69,76),(85,92)]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.pcolor(Qref.T,cmap='binary')
    ax.set_xticks(range(0,len(Qref),10))
    ax.set_yticks(range(0,len(Qref),10))
    ax.set_title("1RIS")
    ax.grid(True)

    add_secondary_struct_icons(ax,elements,seqbounds)

    plt.savefig("1RIS/map_ss.pdf")
    plt.show()

