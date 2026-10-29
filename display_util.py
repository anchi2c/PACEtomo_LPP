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
        std_grid: image contrast (normalized std = std / mean), indexed [u, v]
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
    measured contrast (std / mean), 2D fit (cosine or quadratic), and residual (measured - fit).
    Axes are phase offsets (rad) along the two ronchiCorrMatrix vectors
    from the scan center x-tilt (+). The fitted maximum contrast is marked (x).
    """
    import matplotlib.pyplot as plt
    number_of_scans = len(Contrast_Std)
    ncols = 3
    nrows = number_of_scans
    fig, ax = plt.subplots(nrows, ncols, figsize=(ncols * panel_width, nrows * panel_width), squeeze=False)
    fig.suptitle(title or "Image std/mean vs x-tilt along ronchiCorrMatrix vectors", fontweight="bold")
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
            (std_grid, 'measured std/mean', cmap, vmin, vmax),
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

C3Scan_Offsets = []
C3Scan_RegionStds = []
C3Scan_RegionNames = []
C3Scan_Best = []
C3Scan_Start = []

def addC3ContrastScan(offsets, region_stds, region_names, best_offset, start_offset):
    """
    Add a C3 (ImageDistanceOffset) contrast scan for display.
        offsets: scanned ImageDistanceOffset values
        region_stds: std/mean, array (number of offsets, number of regions); region 0 is the center
        region_names: names of the regions
    """
    global C3Scan_Offsets, C3Scan_RegionStds, C3Scan_RegionNames, C3Scan_Best, C3Scan_Start
    C3Scan_Offsets.append(np.asarray(offsets, dtype=float))
    C3Scan_RegionStds.append(np.asarray(region_stds, dtype=float))
    C3Scan_RegionNames.append(tuple(region_names))
    C3Scan_Best.append(best_offset)
    C3Scan_Start.append(start_offset)

def resetC3ContrastScan():
    global C3Scan_Offsets, C3Scan_RegionStds, C3Scan_RegionNames, C3Scan_Best, C3Scan_Start
    C3Scan_Offsets = []
    C3Scan_RegionStds = []
    C3Scan_RegionNames = []
    C3Scan_Best = []
    C3Scan_Start = []

def showC3ContrastScan(panel_width=4.5, title='', savefig=False):
    """
    One row per C3 scan: std/mean of each region against ImageDistanceOffset,
    and the difference center - mean of corners. The start (gray dotted) and
    picked best offset (red dashed) are marked.
    """
    import matplotlib.pyplot as plt
    number_of_scans = len(C3Scan_Offsets)
    fig, ax = plt.subplots(number_of_scans, 2, figsize=(2 * panel_width, number_of_scans * panel_width * 0.75), squeeze=False)
    fig.suptitle(title or "Image std/mean vs ImageDistanceOffset", fontweight="bold")
    for n in range(number_of_scans):
        offsets = C3Scan_Offsets[n]
        region_stds = C3Scan_RegionStds[n]
        names = C3Scan_RegionNames[n]
        diffs = region_stds[:, 0] - region_stds[:, 1:].mean(axis=1)
        a = ax[n][0]
        for r in range(region_stds.shape[1]):
            a.plot(offsets, region_stds[:, r], marker='o', markersize=3,
                   linewidth=2.0 if r == 0 else 1.0, color='black' if r == 0 else None, label=names[r])
        a.set_ylabel('std/mean')
        a.set_title('regions')
        a.legend(fontsize=7)
        b = ax[n][1]
        b.plot(offsets, diffs, marker='o', markersize=3, color='tab:purple')
        b.axhline(0, color='gray', linewidth=0.8)
        b.set_ylabel('center - corners std/mean')
        b.set_title('difference')
        for c in (a, b):
            c.axvline(C3Scan_Start[n], color='gray', linestyle=':', label='start')
            c.axvline(C3Scan_Best[n], color='red', linestyle='--', label='best')
            c.set_xlabel('ImageDistanceOffset')
        b.legend(fontsize=7)
    plt.tight_layout()
    if savefig:
        plt.savefig((title or 'c3_contrast_scan') + '.png', dpi=300, bbox_inches="tight")
    plt.show()

Cross_Images = []
Cross_Results = []
Cross_SubTitle = []
Cross_Display_Binning = 8

def addCrossFit(image, result, subtitle=''):
    """
    Add one on-cross scan for display.
        image: ronchigram (unbinned)
        result: dict from xlpp_on_cross.find_cross_center with center (row, col),
                coeffs of a + b*row + c*col + d*row^2 + e*row*col + f*col^2,
                positions, is_node, contrast and used of the lattice points
    The image is kept binned by Cross_Display_Binning to save memory.
    """
    global Cross_Images, Cross_Results, Cross_SubTitle
    b = Cross_Display_Binning
    arr = np.asarray(image, dtype=float)
    arr = arr[:arr.shape[0] // b * b, :arr.shape[1] // b * b]
    Cross_Images.append(arr.reshape(arr.shape[0] // b, b, arr.shape[1] // b, b).mean(axis=(1, 3)))
    Cross_Results.append(result)
    Cross_SubTitle.append(subtitle)

def resetCrossFit():
    global Cross_Images, Cross_Results, Cross_SubTitle
    Cross_Images = []
    Cross_Results = []
    Cross_SubTitle = []

def showCrossFit(panel_width=4.5, title='', savefig=False):
    """
    One row per scan:
    left: ronchigram with nodes (o) and antinodes (s) colored by local fringe
    contrast; points used in the fit have a black edge.
    right: contrast of the fitted points with contours of the 2D quadratic fit.
    The fitted cross center is a red +, the image center a cyan x.
    """
    import matplotlib.pyplot as plt
    number_of_scans = len(Cross_Images)
    b = Cross_Display_Binning
    fig, ax = plt.subplots(number_of_scans, 2, figsize=(2 * panel_width, number_of_scans * panel_width * 0.8), squeeze=False)
    fig.suptitle(title or "Laser fringe cross center fit", fontweight="bold")
    for n in range(number_of_scans):
        image = Cross_Images[n]
        r = Cross_Results[n]
        pos = np.asarray(r['positions'])
        contrast = np.asarray(r['contrast'])
        valid = ~np.isnan(contrast)
        used = np.asarray(r['used'])
        is_node = np.asarray(r['is_node'])
        center = r['center']
        rows, cols = r['image_shape']
        vmax = np.percentile(contrast[valid], 98) if valid.any() else 1.0
        # left: image with lattice points (display pixels = unbinned / b)
        a = ax[n][0]
        lo, hi = np.percentile(image, [1, 99])
        a.imshow(image, cmap='gray', vmin=lo, vmax=hi)
        for marker, mask in (('o', is_node), ('s', ~is_node)):
            m = mask & valid
            a.scatter(pos[m, 1] / b, pos[m, 0] / b, c=contrast[m], cmap='viridis', vmin=0, vmax=vmax,
                      marker=marker, s=14, edgecolors=np.where(used[m], 'black', 'none'), linewidths=0.6)
        a.plot(center[1] / b, center[0] / b, '+', color='red', markersize=14, markeredgewidth=2.5)
        a.plot(cols / 2 / b, rows / 2 / b, 'x', color='cyan', markersize=10, markeredgewidth=2)
        a.set_title(f'{Cross_SubTitle[n]} nodes (o) antinodes (s)'.strip(), fontsize=9)
        a.set_xlabel(f'x (px / {b})')
        a.set_ylabel(f'y (px / {b})')
        # right: fitted points and quadratic fit contours
        c = ax[n][1]
        sc = c.scatter(pos[used, 1], pos[used, 0], c=contrast[used], cmap='viridis', vmin=0, vmax=vmax, s=16)
        fig.colorbar(sc, ax=c, fraction=0.046, pad=0.04, label='contrast')
        rr, cc = np.meshgrid(np.linspace(0, rows, 80), np.linspace(0, cols, 80), indexing='ij')
        a0, a1, a2, a3, a4, a5 = r['coeffs']
        fit = a0 + a1 * rr + a2 * cc + a3 * rr * rr + a4 * rr * cc + a5 * cc * cc
        cs = c.contour(cc, rr, fit, levels=np.linspace(0, vmax, 8), colors='gray', linewidths=0.8)
        c.clabel(cs, fontsize=6, fmt='%.0f')
        c.plot(center[1], center[0], '+', color='red', markersize=14, markeredgewidth=2.5)
        c.plot(cols / 2, rows / 2, 'x', color='cyan', markersize=10, markeredgewidth=2)
        c.set_xlim(0, cols)
        c.set_ylim(rows, 0)
        c.set_aspect('equal')
        c.set_title(f'{Cross_SubTitle[n]} fit, cross (x, y) = ({center[1]:.0f}, {center[0]:.0f})'.strip(), fontsize=9)
        c.set_xlabel('x (px)')
        c.set_ylabel('y (px)')
    plt.tight_layout()
    if savefig:
        plt.savefig((title or 'cross_fit') + '.png', dpi=300, bbox_inches="tight")
    plt.show()

CrossSteps_C3 = []
CrossSteps_Centers = []
CrossSteps_Target = []
CrossSteps_XT = []
CrossSteps_Final = []
CrossSteps_Applied = []

def addCrossSteps(c3_offsets, centers_xy, target_xy, xt_values, final_center_xy, applied=None):
    """
    Add an on-cross correction run for display.
        c3_offsets: C3 offset of each step
        centers_xy: cross center (x, y) found at each step, before its correction
        target_xy: target (x, y) of the cross
        xt_values: xt at the start and after each step's correction
        final_center_xy: cross center (x, y) of the final check at the last C3 offset
        applied: whether each step's correction was applied (default all)
    """
    global CrossSteps_C3, CrossSteps_Centers, CrossSteps_Target, CrossSteps_XT, CrossSteps_Final, CrossSteps_Applied
    CrossSteps_C3.append(np.asarray(c3_offsets, dtype=float))
    CrossSteps_Centers.append(np.asarray(centers_xy, dtype=float))
    CrossSteps_Target.append(np.asarray(target_xy, dtype=float))
    CrossSteps_XT.append(np.asarray(xt_values, dtype=float))
    CrossSteps_Final.append(np.asarray(final_center_xy, dtype=float))
    CrossSteps_Applied.append(np.ones(len(c3_offsets), dtype=bool) if applied is None else np.asarray(applied, dtype=bool))

def resetCrossSteps():
    global CrossSteps_C3, CrossSteps_Centers, CrossSteps_Target, CrossSteps_XT, CrossSteps_Final, CrossSteps_Applied
    CrossSteps_C3 = []
    CrossSteps_Centers = []
    CrossSteps_Target = []
    CrossSteps_XT = []
    CrossSteps_Final = []
    CrossSteps_Applied = []

def showCrossSteps(panel_width=4.5, title='', savefig=False):
    """
    One row per on-cross correction run:
    left: cross offset from the target (x, y) at each step before its
    correction (hollow: correction skipped as unreliable), and the final
    check (star) at the last C3 offset.
    middle: xt (x, y) at the start and after each step's correction.
    right: path of the cross center in the image with the target.
    """
    import matplotlib.pyplot as plt
    number_of_runs = len(CrossSteps_C3)
    fig, ax = plt.subplots(number_of_runs, 3, figsize=(3 * panel_width, number_of_runs * panel_width * 0.8), squeeze=False)
    fig.suptitle(title or "On-cross correction steps", fontweight="bold")
    for n in range(number_of_runs):
        c3 = CrossSteps_C3[n]
        centers = CrossSteps_Centers[n]
        target = CrossSteps_Target[n]
        xt = CrossSteps_XT[n]
        final = CrossSteps_Final[n]
        applied = CrossSteps_Applied[n]
        step = np.arange(1, len(c3) + 1)
        offsets = centers - target
        final_offset = final - target
        a = ax[n][0]
        for k, (label, color) in enumerate((('x', 'tab:blue'), ('y', 'tab:orange'))):
            a.plot(step, offsets[:, k], '-', color=color, label=f'{label}')
            a.plot(step[applied], offsets[applied, k], 'o', color=color)
            a.plot(step[~applied], offsets[~applied, k], 'o', markerfacecolor='none', color=color)
            a.plot(len(c3) + 1, final_offset[k], '*', color=color, markersize=14, markeredgecolor='black',
                   label=f'{label} final {final_offset[k]:.1f}')
        a.axhline(0, color='gray', linewidth=0.8)
        a.set_xticks(list(step) + [len(c3) + 1])
        a.set_xticklabels([f'{v:.0f}' for v in c3] + [f'final\n{c3[-1]:.0f}'], fontsize=7)
        a.set_xlabel('C3 offset of step')
        a.set_ylabel('cross offset from target (px)')
        a.legend(fontsize=7)
        b = ax[n][1]
        xt_step = np.arange(len(xt))
        for k, (label, color) in enumerate((('xt x', 'tab:blue'), ('xt y', 'tab:orange'))):
            b.plot(xt_step, xt[:, k], 'o-', color=color, label=label)
        b.set_xticks(xt_step)
        b.set_xticklabels(['start'] + [f'{v:.0f}' for v in c3], fontsize=7)
        b.set_xlabel('after correction at C3 offset')
        b.set_ylabel('xt (rad)')
        b.legend(fontsize=7)
        c = ax[n][2]
        path = np.vstack([centers, final])
        c.plot(path[:, 0], path[:, 1], '-', color='gray', linewidth=0.8)
        sc = c.scatter(centers[:, 0], centers[:, 1], c=c3, cmap='plasma', s=25)
        c.plot(final[0], final[1], '*', color='green', markersize=14, markeredgecolor='black', label='final')
        c.plot(target[0], target[1], '+', color='red', markersize=14, markeredgewidth=2.5, label='target')
        fig.colorbar(sc, ax=c, fraction=0.046, pad=0.04, label='C3 offset')
        c.invert_yaxis()
        c.set_aspect('equal', adjustable='datalim')
        c.set_xlabel('x (px)')
        c.set_ylabel('y (px)')
        c.legend(fontsize=7)
    plt.tight_layout()
    if savefig:
        plt.savefig((title or 'cross_steps') + '.png', dpi=300, bbox_inches="tight")
    plt.show()
