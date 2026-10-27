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

Fit_X = []
Fit_UV = []
Fit_T = []
Fit_SubTitle = []

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
        if vmin == vmax:
            print('vmin==vmax')
            vmin = 0.99*vmax
        arr = 255*((arr -vmin)/(vmax-vmin))
    else:
        vmin, vmax = np.percentile(arr,[1,99])
        vmean = arr.mean()
        # use mean std method if needed
        vstd_max = vmean+4*arr.std()
        if vmin==vmax or vstd_max < (vmax+vmin)/2:
            vmin = vmean-4*arr.std()
            vmax = vstd_max
        arr = 255*((arr -vmin)/(vmax-vmin))
    image_buffer.append(arr.copy())

    for peak in peaks:
        print('peak (y,x)',peak)
        arr_min_shape = min(arr.shape)
        a = 0.01 * arr_min_shape
        b = 0.035 * arr_min_shape
        b1 = 0.04 * arr_min_shape
        c = np.array(arr.shape)//2
        fill = (arr.max()-arr.min())*0.5 + arr.min()
        #show peak
        image_buffer[-1][int(-peak[0]+c[0]-b1):int(-peak[0]+c[0]+b1),int(-peak[1]+c[1]-b1):int(-peak[1]+c[1]+b1)] = arr.max()
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
        fig.suptitle(title)
    plt.tight_layout()
    if savefig and title:
        plt.savefig(title+'.png', dpi=300, bbox_inches="tight")
    plt.show()

def addTransformFitData(new_transform,old_cal_xts,xt_residuals_from_old,is_shifts, subtitle=''):
    xt_residuals_from_new = is_shifts @ np.linalg.inv(new_transform)
    
    global Fit_X,Fit_UV,Fit_T,Fit_SubTitle
    Fit_X.append(is_shifts)
    Fit_UV.append([xt_residuals_from_old, xt_residuals_from_new])
    Fit_T.append(new_transform)
    Fit_SubTitle.append(subtitle)

def resetTransformFitData():
    global Fit_X,Fit_UV,Fit_T,Fit_SubTitle
    Fit_X = []
    Fit_UV = []
    Fit_T = []
    Fit_SubTitle = []

def showTransformFit():
    """
    Create graph that shows residual vectors at each calibration applied shift and resulting
    transform matrix fit.
    For refining is_xt_matrix
        Fit_X: applied image shifts
        Fit_UV: xt_residuals_from_old matrix, xt_residuals_from new matrix
        Fit_T: 2x2 transform matrix
    """
    import matplotlib.pyplot as plt
    # ncols: number of columns of subplots
    # panel_width: width of the each subplot panel in inches
    # cmap: colormap. None uses default in matplotlib
    ncols = 3
    number_of_images = len(Fit_T)
    print(number_of_images)
    nrows = int(math.ceil(number_of_images / ncols))
    # layout
    aspect = 1
    panel_width = 3
    fig, ax = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(ncols * panel_width, nrows * panel_width * aspect), squeeze=False)
    fig.suptitle("Calibration residual normalized by image shift magnitude", fontweight="bold")

    U1 = []
    U2 = []
    V1 = []
    V2 = []
    for n in range(number_of_images):
        col = n % ncols
        row = n // ncols
        X = Fit_X[n][:,0]
        Y = Fit_X[n][:,1]
        # Define vector components at each (x, y) point
        scale = 0.5e6
        U1.append(Fit_UV[n][0][:,0] * scale)
        V1.append(Fit_UV[n][0][:,1] * scale)
        U2.append(Fit_UV[n][1][:,0] * scale)
        V2.append(Fit_UV[n][1][:,1] * scale)
        # Compute magnitude for coloring (optional)
        magnitude = np.sqrt(U1[0]**2 + V1[0]**2)

        # Plot the vectors
        q1 = ax[row][col].quiver(X, Y, U1[n], V1[n], label='norm_res', color='orange', scale=50)
        p1 = ax[row][col].scatter(X, Y, color='orange', s=8)
        #q2 = ax[row][col].quiver(X, Y, U2[n], V2[n], label='LAFIS new', color='green', scale=50)


    X_stack = np.vstack(Fit_X)
    for n in range(number_of_images):
        col = n % ncols
        row = n // ncols
        ax[row][col].set_xlim(X_stack[:,0].min() - 2, X_stack[:,1].max() + 2)
        ax[row][col].set_ylim(X_stack[:,1].min() - 2, X_stack[:,1].max() + 2)

        ax[row][col].set_xlabel('Image Shift X (um)')
        ax[row][col].set_ylabel('Image Shift Y (um)')
        ax[row][col].set_title(Fit_SubTitle[n])
        ax[row][col].set_aspect('equal')

        ax[row][col].legend(loc='lower right', fontsize=11, bbox_to_anchor=(0.995, 0.035))
    plt.tight_layout()
    plt.show()

