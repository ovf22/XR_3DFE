import matplotlib.pyplot as plt
import numpy as np
import sys


def tomo_crop_center(data_crop, width=100, height=100, depth=100):
    """ Central cropping function for 3D data

    Parameters
    ----------
    data_crop: 3D X-ray data for cropping  [3D Array] \n
    width: cropping width  [int] \n
    height: cropping height  [int] \n
    depth: cropping depth [int] \n

    Returns
    -------
    data: 3D X-ray CT data after cropping [3D Array] \n
    """
    data_shape = data_crop.shape

    max_model_width = data_shape[0]
    if width > max_model_width:
        width = max_model_width
        print('The chosen model width was to large, '
              'and was set to the maximum allowed height of ', max_model_width)
        print('The new aspect ratio width/height = ', width/height)

    max_model_depth = data_shape[2]
    if depth > max_model_depth:
        depth = max_model_depth
        print('The chosen model depth was to large, '
              'and was set to the depth of the available data: ',
              max_model_depth)
        print('The new aspect ratio width/depth = ', width/depth)

    # Write tomography data into array
    cut_arr = np.array([width, height, depth])

    # If the cropping dimensions are bigger than the tomography data set,
    # then kill the script
    if np.any(cut_arr > data_shape):
        print('Too much data was removed. Adjust cut parameters')
        print('Data_shape = ' + str(data_shape))
        print('Crop array shape = ' + str(cut_arr.shape()))
        sys.exit()

    # If any cutting parameters are chosen, then crop the volume accordingly
    # Else return the full data set
    xc, yc, zc = np.array(data_shape[:])/2
    if np.any(cut_arr > 0):
        print('Cropping data')
        x_slice = slice(int(xc-width/2), int(xc+width/2))
        y_slice = slice(int(yc-height/2), int(yc+height/2))
        z_slice = slice(int(zc-depth/2), int(zc+depth/2))
        data = data_crop[x_slice, y_slice, z_slice]
    else:
        data = data_crop
        print('No cropping performed')
    return data


def tomo_crop(data_crop, xcut=0, ycut=0, zcut=0):
    """ Crop the tomography data according to the cutting parameters

    Parameters
    ----------
    data_crop : Tomography data file [Array] \n
    xcut : Amount of voxels removed from each side of the volume in the
    x-direction [int]. Default is 0.\n
    ycut : Amount of voxels removed from each side of the volume in the
    y-direction [int]. Default is 0.\n
    zcut : Amount of voxels removed from each side of the volume in the
    x-direction [int]. Default is 0.

    Returns
    -------
    data : Tomography data from ROI [Array]

    """
    # Write tomography data into array
    no_slicing = slice(None)
    data_shape = data_crop.shape
    cut_arr = 2 * np.array([xcut, ycut, zcut])

    # If the cropping dimensions are bigger than the tomography data set,
    # then kill the script
    if np.any(cut_arr > data_shape):
        print('Too much data was removed. Adjust cut parameters')
        sys.exit()

    # If any cutting parameters are chosen, then crop the volume accordingly
    # Else return the full data set
    if np.any(cut_arr > 0):
        print('Cropping data')
        if xcut != 0:
            x_slice = slice(xcut, -xcut)
        else:
            x_slice = no_slicing
        if ycut != 0:
            y_slice = slice(ycut, -ycut)
        else:
            y_slice = no_slicing
        if zcut != 0:
            z_slice = slice(zcut, -zcut)
        else:
            z_slice = no_slicing

        data = data_crop[x_slice, y_slice, z_slice]
    else:
        data = data_crop
        print('No cropping performed')
    return data


def tomo_plot_3(data, x_slice=0, y_slice=0, z_slice=0, title='',
                fig_name='Tomo_fig', fig_path='', figsize=(10, 10),
                FrameOn=True, LineOn=True, VOXEL_SIZE=None):
    """ PLotting three mutually orthogonal planes of X-ray CT data

    Parameters
    ----------
    data: 3D X-ray CT data  [3D array] \n
    x_slice: x slice index [int] \n
    y_slice: y slice index [int] \n
    z_slice: z slice index [int] \n
    title: figure title [str] \n
    figsize: figure size [1x2 tuple] \n
    FrameOn: Include colored frame for different planes \n
    LineOn: Include colored lines at slice locations  \n
    """

    fig3 = plt.figure(constrained_layout=True, figsize=figsize)
    gs = fig3.add_gridspec(2, 2, width_ratios=([1.5, 1]),
                           height_ratios=(1, data.shape[2]/data.shape[1]),
                           wspace=-0.1, hspace=0.05, left=0.1, right=0.9,
                           bottom=0.1, top=0.9)
    f3_ax1 = fig3.add_subplot(gs[0, 0])
    f3_ax2 = fig3.add_subplot(gs[1, 0], sharex=f3_ax1)
    f3_ax3 = fig3.add_subplot(gs[:, 1])
    L, W, T = data.shape

    liney = np.array([[0, L], [y_slice, y_slice]])
    linex = np.array([[x_slice, x_slice], [0, W]])
    frame2 = np.array([[1, L-1, L-1, 1, 1], [1, 1, W-1, W-1, 1]])
    if LineOn:
        f3_ax1.plot(liney[0, :], liney[1, :], 'g', linewidth=2)
        f3_ax1.plot(linex[0, :], linex[1, :], 'r', linewidth=2)
    if FrameOn:
        f3_ax1.plot(frame2[0, :], frame2[1, :], 'b', linewidth=3)
    img = data[:, :, z_slice].T
    f3_ax1.imshow(img.squeeze(), cmap='gray')
    f3_ax1.set_ylim(f3_ax1.get_ylim()[::-1])
    f3_ax1.set(xlabel='x-axis  [Voxels]',
               ylabel='y-axis  [Voxels]')

    # Plot xz-slice
    linex = np.array([[x_slice, x_slice], [0, T]])
    linez = np.array([[0, L], [z_slice, z_slice]])
    frame1 = np.array([[2, L-2, L-2, 2, 2], [2, 2, T-2, T-2, 2]])
    img = data[:, y_slice, :].T
    if LineOn:
        f3_ax2.plot(linex[0, :], linex[1, :], 'r', linewidth=2)
        f3_ax2.plot(linez[0, :], linez[1, :], 'b', linewidth=2)
    if FrameOn:
        f3_ax2.plot(frame1[0, :], frame1[1, :], 'g', linewidth=4)
    f3_ax2.imshow(img.squeeze(), cmap='gray')
    f3_ax2.set_ylim(f3_ax2.get_ylim()[::-1])
    f3_ax2.set(xlabel='x-axis  [Voxels]',
               ylabel='z-axis  [Voxels]')

    # Plot yz-slice
    liney = np.array([[0, T], [y_slice, y_slice]])
    linez = np.array([[z_slice, z_slice], [0, W]])
    frame0 = np.array([[1, T-1, T-1, 1, 1], [1, 1, W-1, W-1, 1]])
    if LineOn:
        f3_ax3.plot(liney[1, :], liney[0, :], 'g', linewidth=2)
        f3_ax3.plot(linez[1, :], linez[0, :], 'b', linewidth=2)
    if FrameOn:
        f3_ax3.plot(frame0[1, :], frame0[0, :], 'r', linewidth=3)
    img = data[x_slice, :, :].T
    f3_ax3.imshow(img.squeeze(), cmap='gray')
    f3_ax3.set_ylim(f3_ax3.get_ylim()[::-1])
    f3_ax3.set(xlabel='y-axis  [Voxels]',
               ylabel='z-axis  [Voxels]')

    if not VOXEL_SIZE:
        pass
    else:
        Nticks = []
        for idx, dim in enumerate(data.shape):
            if dim < 500:
                Nticks.append(100)
            elif dim > 1200:
                Nticks.append(5000)
            else:
                Nticks.append(200)

        xticks = np.arange(0, data.shape[0], Nticks[0]/VOXEL_SIZE)
        xlabels = ["{:.0f}".format(round(x*VOXEL_SIZE, -1)) for x in xticks]

        yticks = np.arange(0, data.shape[1], Nticks[1]/VOXEL_SIZE)
        ylabels = ["{:.0f}".format(round(y*VOXEL_SIZE, -1)) for y in yticks]

        zticks = np.arange(0, data.shape[2], Nticks[2]/VOXEL_SIZE)
        zlabels = ["{:.0f}".format(round(z*VOXEL_SIZE, -1)) for z in zticks]

        f3_ax1.set_xticks(xticks, labels=xlabels)
        f3_ax1.set_yticks(yticks, labels=ylabels)
        f3_ax1.set(xlabel='x-axis  [$\\mu m$]',
                   ylabel='y-axis  [$\\mu m$]')

        f3_ax2.set_yticks(zticks, labels=zlabels)
        f3_ax2.set(xlabel='x-axis  [$\\mu m$]',
                   ylabel='z-axis  [$\\mu m$]')

        f3_ax3.set_xticks(yticks, labels=ylabels)
        f3_ax3.set_yticks(zticks, labels=zlabels)
        f3_ax3.set(xlabel='y-axis  [$\\mu m$]',
                   ylabel='z-axis  [$\\mu m$]')

    plt.savefig(fig_path + fig_name + '.png')


def scaling_data_fade(data, PADDING, VOXEL_SIZE, FVF_measured=None):
    """ Function used for scaling field variables at the start and end along
        the x-axis. This is used to avoid stress concentrations at boundary
        conditions.


    Parameters
    ----------
    data : Dataset being scaled at start and end of x-axis
    PADDING : Length of padding at start and end [int]
    VOXEL_SIZE : X-rat CT voxel size [float] \n
    FVF_measured : Measured Fiber Volume Fraction [float] \n
        Default value is None. This value is only used if the field variable
        being scaled is the FVF. \n

    Returns
    -------
    data_new : Updated dataset with with scaling along the
        start and end of the x-axis

    """
    data_copy = data.copy()
    PADDING_idx = int(PADDING/VOXEL_SIZE)

    scale1d = np.linspace(1, 0, PADDING_idx)
    beta = 2
    scale1dS = 1 / (1 + (scale1d[1: -1] / (1 - scale1d[1: -1]))**(-beta))
    scale1dS = np.insert(scale1dS, [0, scale1dS.size], [1, 0])

    scale2d = np.repeat(scale1dS[:, np.newaxis], data_copy.shape[1], axis=1)
    scale3d = np.repeat(scale2d[:, :, np.newaxis], data_copy.shape[2], axis=2)

    if FVF_measured is not None:
        scale1dSrev = np.flip(scale1dS)

        scale2drev = np.repeat(scale1dSrev[:, np.newaxis],
                               data_copy.shape[1], axis=1)
        scale3drev = np.repeat(scale2drev[:, :, np.newaxis],
                               data_copy.shape[2], axis=2)

        data_right = data_copy[-PADDING_idx:, :, :] * scale3d + \
            scale3drev * FVF_measured
        data_left = data_copy[:PADDING_idx, :, :] * np.flip(scale3d, 0) + \
            np.flip(scale3drev, 0) * FVF_measured

    else:
        data_right = data_copy[-PADDING_idx:, :, :] * scale3d
        data_left = data_copy[:PADDING_idx, :, :] * np.flip(scale3d, 0)

    data_new = np.concatenate((data_left,
                               data[PADDING_idx:-PADDING_idx, :, :],
                               data_right), axis=0)
    return data_new


def tomo_plot_3_overlay(data, field1, field2, y_slice=0, z_slice=0,
                        alpha=0.5, vmin=-5, vmax=5, variable='$\\phi_{ij}$',
                        unit='[$^\\circ$]', cmap='coolwarm', title='',
                        fig_name='Tomo_fig', fig_path='', figsize=(10, 10)):
    """ Figure for plotting field variable on top of X-ray CT data

    Parameters
    ----------
    data: 3D X-ray dataset  [3D array] \n
    field1: Field variable data 1 [3D array] \n
    field2: Field variable data 2 [3D array] \n
    y_slice: y slice  [int] \n
    z_slice: z-slice  [int] \n
    alpha: Transparency factor  [float between 0 and 1] \n
    vmin: lower limit for field variable  [float] \n
    vmax: upper limit for field variable  [float] \n
    variable: Field variable string  [str] \n
    unit: Field variable unit  [str] \n
    cmap: cooler map  [str] \n
    title: figure title  [str] \n
    figsize: figure sixe  [tuple 1x2] \n
    """

    fig3 = plt.figure(constrained_layout=True, figsize=figsize)
    gs = fig3.add_gridspec(2, 1,
                           height_ratios=(1, data.shape[2]/data.shape[1]),
                           wspace=-0.1, hspace=0.05, left=0.1, right=0.9,
                           bottom=0.1, top=0.9)
    f3_ax1 = fig3.add_subplot(gs[0, 0])
    f3_ax2 = fig3.add_subplot(gs[1, 0], sharex=f3_ax1)
    L, W, T = data.shape

    img = data[:, :, z_slice].T
    f3_ax1.imshow(img.squeeze(), cmap='gray')
    f3_ax1.set_ylim(f3_ax1.get_ylim()[::-1])
    f3_ax1.set(xlabel='x-axis  [Voxel index]',
               ylabel='y-axis  [Voxel index]')
    f3_ax1.title.set_text('Out-of-plane')
    im = f3_ax1.imshow(field1[:, :, z_slice].T, alpha=alpha, cmap=cmap,
                       vmin=vmin, vmax=vmax)

    img = data[:, y_slice, :].T
    f3_ax2.imshow(img.squeeze(), cmap='gray')
    f3_ax2.set_ylim(f3_ax2.get_ylim()[::-1])
    f3_ax2.title.set_text('In-plane')
    f3_ax2.set(xlabel='x-axis  [Voxel index]',
               ylabel='z-axis  [Voxel index]')
    im = f3_ax2.imshow(field2[:, y_slice, :].T, alpha=alpha, cmap=cmap,
                       vmin=vmin, vmax=vmax)

    cb_ax = fig3.add_axes([1.01, 0.32, 0.02, 0.5])
    cbar = fig3.colorbar(im, cax=cb_ax)

    # Set the colorbar ticks and tick labels
    cbar.set_ticks([vmin, vmin/2, 0, vmax/2, vmax])
    cbar.set_ticklabels(['$<$'+str(vmin), str(vmin/2), '$0$',
                         str(vmax/2), '$>$'+str(vmax)])
    cbar.set_label(variable + unit, rotation=0, y=1.08, labelpad=-21)
    plt.savefig(fig_path + fig_name + '.png', dpi=70)
