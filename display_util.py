#!Python
# ===================================================================
#ScriptName     Image Display utilities
# Purpose:      Use matplotlib for diagonosis 
# Author:       Anchi Cheng
# ===================================================================
import numpy as np
import math
from scipy import ndimage
image_buffer = []

def addImage(arr, peaks=[]):
    """
    Add image to image_buffer for display
    # peak is correlation shift with 0,0 unshifted in numpy convention
    # when peaks are specified. arr is wrapped correlation image with unshifted
    # at the center of the correlation image
    """
    global image_buffer
    if np.iscomplexobj(arr):
        arr = np.log1p(np.abs(arr).copy())
        vmin, vmax = np.percentile(arr, [1, 99])
        arr = 255*((arr -vmin)/(vmax-vmin))
    else:
        vmin, vmax = np.percentile(arr,[1,99])
        arr = 255*((arr -vmin)/(vmax-vmin))
    image_buffer.append(arr.copy())

    for peak in peaks:
        print('peak (y,x)',peak)
        arr_min_shape = min(arr.shape)
        a = 0.01 * arr_min_shape
        b = 0.04 * arr_min_shape
        c = np.array(arr.shape)//2
        fill = (arr.max()-arr.min())*2 + arr.max()
        #show peak
        image_buffer[-1][int(-peak[0]+c[0]-b):int(-peak[0]+c[0]+b),int(-peak[1]+c[1]-b):int(-peak[1]+c[1]+b)] = fill
        #show center
        image_buffer[-1][int(c[0]-a):int(c[0]+a),int(c[1]-a):int(c[1]+a)] = fill*0.8

def showImages(ncols=8,panel_width=2,cmap=None,title='',savefig=False):
    import matplotlib.pyplot as plt
    # ncols: number of columns of subplots
    # panel_width: width of the each subplot panel in inches
    # cmap: colormap. None uses default in matplotlib
    number_of_buffer_images = len(image_buffer)
    nrows = int(math.ceil(number_of_buffer_images / ncols))
    # layout
    aspect = 1
    fig, ax = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(ncols * panel_width, nrows * panel_width * aspect), squeeze=False)
    for i in range(number_of_buffer_images):
        col = i % ncols
        row = i // ncols
        if i == 0:
            shape0 = image_buffer[i].shape
            display_arr = image_buffer[i]
        else:
            my_shape = image_buffer[i].shape
            factor = (shape0[0]/my_shape[0])
            display_arr = ndimage.zoom(image_buffer[i], factor, order=1)

        if cmap is not None:
            ax[row][col].imshow(display_arr, cmap=cmap)
        else:
            ax[row][col].imshow(display_arr)
    if title:
        plt.title(title)
    plt.tight_layout()
    if savefig and title:
        plt.savefig(title+'.png', dpi=300, bbox_inches="tight")
    plt.show()
