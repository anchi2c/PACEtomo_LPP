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

Scatter_X = []
Scatter_Y = []
Scatter_Fit_Slope = []
Scatter_Fit_Intercept = []

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

def addScatterFit(xs,ys,slope,intercept):
    global Scatter_X,Scatter_Y,Scatter_Fit_Slope,Scatter_Fit_Intercept
    Scatter_X.append(xs)
    Scatter_Y.append(ys)
    Scatter_Fit_Slope.append(slope)
    Scatter_Fit_Intercept.append(intercept)

def resetScatterFit():
    global Scatter_X,Scatter_Y,Scatter_Fit_Slope,Scatter_Fit_Intercept
    Scatter_X = []
    Scatter_Y = []
    Scatter_Fit_Slope = []
    Scatter_Fit_Intercept = []

def showScatterFit():
    import matplotlib.pyplot as plt
    number_of_scatter_fit = len(Scatter_X)

    for i in range(number_of_scatter_fit):
        x = Scatter_X[i]
        y = Scatter_Y[i]
        xf = np.array([x.min(), x.max()])
        slope = Scatter_Fit_Slope[i]
        intercept = Scatter_Fit_Intercept[i]
        plt.scatter(x, y, s=15, label=f'data{i+1}')
        plt.plot(xf, slope * xf + intercept, "r-",
             label=f"fit {i+1}: y = {slope:.2f}x + {intercept:.2f}")
        plt.legend()
        plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

Contrast_Phases = []
Contrast_Std = []
Contrast_Fit = []
Contrast_SubTitle = []
Contrast_Model = []

def addContrastGrid(phases, std_grid, fit_params, u_max, v_max, subtitle='', model='cosine'):
    """
    Add x-tilt contrast scan for display.
        phases: phase coordinates (rad) of the scan along each ronchiCorrMatrix vector
        std_grid: image standard deviation, indexed [u, v]
        fit_params: (c, a_u, f_u, p_u, a_v, f_v, p_v) of the 2D cosine fit, or
                    (a, b, c, d, e, f) of a + b*u + c*v + d*u^2 + e*u*v + f*v^2
        u_max, v_max: fitted maximum contrast phase coordinates
        model: 'cosine' or 'quadratic'
    """
    global Contrast_Phases, Contrast_Std, Contrast_Fit, Contrast_SubTitle, Contrast_Model
    Contrast_Phases.append(np.asarray(phases, dtype=float))
    Contrast_Std.append(np.asarray(std_grid, dtype=float))
    Contrast_Fit.append((tuple(fit_params), u_max, v_max))
    Contrast_SubTitle.append(subtitle)
    Contrast_Model.append(model)

def resetContrastGrid():
    global Contrast_Phases, Contrast_Std, Contrast_Fit, Contrast_SubTitle, Contrast_Model
    Contrast_Phases = []
    Contrast_Std = []
    Contrast_Fit = []
    Contrast_SubTitle = []
    Contrast_Model = []

def _cosine_2d_model(u, v, fit_params):
    c, a_u, f_u, p_u, a_v, f_v, p_v = fit_params
    return c + a_u * np.cos(f_u * u - p_u) + a_v * np.cos(f_v * v - p_v)

def _quadratic_2d_model(u, v, fit_params):
    a, b, c, d, e, f = fit_params
    return a + b * u + c * v + d * u * u + e * u * v + f * v * v

def showContrastGrid(panel_width=3.5, cmap='viridis', title='', savefig=False):
    """
    Heatmaps of x-tilt contrast scans, one row per scan:
    measured std, 2D fit (cosine or quadratic), and residual (measured - fit).
    Axes are phase offsets (rad) along the two ronchiCorrMatrix vectors
    from the scan center x-tilt (+). The fitted maximum contrast is marked (x).
    """
    import matplotlib.pyplot as plt
    number_of_scans = len(Contrast_Std)
    ncols = 3
    nrows = number_of_scans
    fig, ax = plt.subplots(nrows, ncols, figsize=(ncols * panel_width, nrows * panel_width), squeeze=False)
    fig.suptitle(title or "Image std vs x-tilt along ronchiCorrMatrix vectors", fontweight="bold")
    for n in range(number_of_scans):
        phases = Contrast_Phases[n]
        std_grid = Contrast_Std[n]
        fit_params, u_max, v_max = Contrast_Fit[n]
        step = phases[1] - phases[0] if len(phases) > 1 else 1.0
        extent = [phases[0] - step / 2, phases[-1] + step / 2, phases[0] - step / 2, phases[-1] + step / 2]
        model = Contrast_Model[n]
        model_fn = _quadratic_2d_model if model == 'quadratic' else _cosine_2d_model
        u_grid, v_grid = np.meshgrid(phases, phases, indexing='ij')
        fit_grid = model_fn(u_grid, v_grid, fit_params)
        # smooth fit on a finer grid for display
        fine = np.linspace(extent[0], extent[1], 101)
        u_fine, v_fine = np.meshgrid(fine, fine, indexing='ij')
        fit_fine = model_fn(u_fine, v_fine, fit_params)
        vmin = min(std_grid.min(), fit_grid.min())
        vmax = max(std_grid.max(), fit_grid.max())
        residual = std_grid - fit_grid
        rmax = np.abs(residual).max() or 1.0
        rms = np.sqrt(np.mean(residual ** 2))
        subtitle = Contrast_SubTitle[n]
        panels = (
            (std_grid, 'measured std', cmap, vmin, vmax),
            (fit_fine, f'{model} fit', cmap, vmin, vmax),
            (residual, f'residual (rms {rms:.3g})', 'RdBu_r', -rmax, rmax),
        )
        for col, (data, label, panel_cmap, lo, hi) in enumerate(panels):
            a = ax[n][col]
            # transpose so u is on the horizontal axis
            im = a.imshow(data.T, origin='lower', extent=extent, cmap=panel_cmap, vmin=lo, vmax=hi,
                          interpolation='nearest')
            fig.colorbar(im, ax=a, fraction=0.046, pad=0.04)
            a.plot(0, 0, '+', color='white', markersize=10, markeredgewidth=2, label='center')
            a.plot(u_max, v_max, 'x', color='red', markersize=10, markeredgewidth=2, label='max fit')
            a.set_xlabel('phase u along corr row 0 (rad)')
            a.set_ylabel('phase v along corr row 1 (rad)')
            a.set_title(f'{subtitle} {label}'.strip())
            a.set_aspect('equal')
        ax[n][0].legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    if savefig:
        plt.savefig((title or 'contrast_grid') + '.png', dpi=300, bbox_inches="tight")
    plt.show()

Quadrant_Stats = []
Quadrant_Selected = []
Quadrant_Names = []
Quadrant_SubTitle = []

def addQuadrantStats(quadrant_stats, selected, quadrant_names, subtitle=''):
    """
    Add quadrant pixel statistics of one x-tilt scan for display.
        quadrant_stats: {stat: array (n, n, 4)} indexed [u, v, quadrant], e.g. mean, min, max, std
        selected: indices of the quadrants used as contrast, or None if the whole box is used
        quadrant_names: names of the 4 quadrants
    """
    global Quadrant_Stats, Quadrant_Selected, Quadrant_Names, Quadrant_SubTitle
    Quadrant_Stats.append({stat: np.asarray(values, dtype=float) for stat, values in quadrant_stats.items()})
    Quadrant_Selected.append(None if selected is None else list(selected))
    Quadrant_Names.append(tuple(quadrant_names))
    Quadrant_SubTitle.append(subtitle)

def resetQuadrantStats():
    global Quadrant_Stats, Quadrant_Selected, Quadrant_Names, Quadrant_SubTitle
    Quadrant_Stats = []
    Quadrant_Selected = []
    Quadrant_Names = []
    Quadrant_SubTitle = []

def showQuadrantStats(panel_width=3.5, title='', savefig=False):
    """
    Line plots of each quadrant's pixel statistics against acquisition order,
    one row per scan and one column per statistic. Quadrants used as contrast
    are solid lines; excluded quadrants are dashed. Dotted vertical lines mark
    the start of each row (u value) of the scan grid.
    """
    import matplotlib.pyplot as plt
    number_of_scans = len(Quadrant_Stats)
    stats = list(Quadrant_Stats[0].keys())
    ncols = len(stats)
    nrows = number_of_scans
    fig, ax = plt.subplots(nrows, ncols, figsize=(ncols * panel_width, nrows * panel_width * 0.8), squeeze=False)
    fig.suptitle(title or "Box quadrant statistics in acquisition order", fontweight="bold")
    colors = ('tab:blue', 'tab:orange', 'tab:green', 'tab:red')
    for n in range(number_of_scans):
        selected = Quadrant_Selected[n]
        names = Quadrant_Names[n]
        for col, stat in enumerate(stats):
            values = Quadrant_Stats[n][stat]
            grid_size = values.shape[0]
            # acquisition order: v changes fastest, as in _scan_contrast_grid
            values = values.reshape(-1, values.shape[-1])
            index = np.arange(1, values.shape[0] + 1)
            a = ax[n][col]
            for q in range(values.shape[1]):
                used = selected is None or q in selected
                a.plot(index, values[:, q], color=colors[q % len(colors)], linestyle='-' if used else '--',
                       linewidth=1.8 if used else 1.0, marker='o', markersize=3,
                       label=names[q] + ('' if used else ' (excluded)'))
            for row_start in range(grid_size + 1, values.shape[0] + 1, grid_size):
                a.axvline(row_start - 0.5, color='gray', linestyle=':', linewidth=0.8)
            a.set_xlabel('image number')
            a.set_ylabel(stat)
            a.set_title(f'{Quadrant_SubTitle[n]} quadrant {stat}'.strip())
        ax[n][0].legend(fontsize=7)
    plt.tight_layout()
    if savefig:
        plt.savefig((title or 'quadrant_stats') + '.png', dpi=300, bbox_inches="tight")
    plt.show()
