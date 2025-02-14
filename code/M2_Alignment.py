from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd


def align_with_x(vec, bc=False, edge_crop=100):
    """ Rotate vectors to be aligned with the global x-axis using
    Rodrigues rotation formula
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

    Parameters
    ----------
    vec: Vector field  [4D array] \n

    Returns
    -------
    vec_new: Rotated vector [4D array] \n
    angle: Angle of rotation [float]

    """

    if bc is False:
        vec_mean = vec.mean(axis=(1, 2, 3))
        vec_mean *= 1/np.sqrt((vec_mean**2).sum())
        z, y = vec_mean[1:]/np.sqrt((vec_mean[1:] ** 2).sum())
        angle = np.arccos(vec_mean[0])
    elif bc is True:
        vec_crop = vec[:, :,
                       slice(edge_crop, -edge_crop),
                       slice(edge_crop, -edge_crop)]
        vec_mean = vec_crop.mean(axis=(1, 2, 3))
        vec_mean *= 1/np.sqrt((vec_mean**2).sum())
        z, y = vec_mean[1:]/np.sqrt((vec_mean[1:] ** 2).sum())
        angle = np.arccos(vec_mean[0])

    K = np.array([[0, z, y], [-z, 0, 0], [-y, 0, 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K

    vec_inter = (R @ vec.reshape(3, -1))
    vec_norm = np.linalg.norm(vec_inter, axis=0, keepdims=True)
    vec_new = np.true_divide(vec_inter, vec_norm,
                             where=(vec_inter != 0) | (vec_norm != 0))

    vec_new = vec_new.reshape(vec.shape)
    angle *= 180 / np.pi

    return vec_new, angle


def angles_from_vec_atan2(vec):
    """ Calculate orientations from vector field


    Parameters
    ----------
    vec: Vector field  [4D array] \n

    Returns
    -------
    theta: Azimuth angle [float - degrees] \n
    phi: Elevation angle [float - degrees] \n
    phi_xy: Orientations projected onto the xy-plane [float - degrees] \n
    phi_xz: Orientations projected onto the xz-plane [float - degrees]

"""
    phi_xy = np.arcsin(vec[1]) * 180 / np.pi
    phi_xz = np.arcsin(vec[2]) * 180 / np.pi
    theta = np.arctan2(vec[2], vec[1]) * 180 / np.pi
    phi = np.arccos(vec[0]) * 180 / np.pi
    return theta, phi, phi_xy, phi_xz


def fig_with_colorbar(data, misalignment, title='', alpha=0.5, cmap=None,
                      vmin=None, vmax=None, variable='$\\phi$ [$^\\circ$]',
                      x_str='x', y_str='y', VOXEL_SIZE=None,
                      fig_name='Misalignment_overlay', fig_path=''):
    """ Creates a figure with data, fiber misalignment overlay and color bar.

    Parameters
    ----------
    data : Tomography data 2D slice [Array]\n
    misalignment : Fiber misalignment 2D overlay [Array]\n
    title : Figure title [str]. Default is empty string\n
    alpha : overlay plot bleeding parameter [float]. The default is 0.5\n
    cmap : Color map for fiber misalignment overlay plot [str]\n
    vmin : Minimum limit for fiber misalignment overlay plot.
    The default is None.\n
    vmax : Minimum limit for fiber misalignment overlay plot.
    The default is None.\n
    fig_name : Name of saved figure without extension [str].
    Default is Tomo_fig\n
    fig_path : Directory for saving figure [str].
    Default is working directory

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    # cax.tick_params(labelsize=16)

    # plot gray scale tomogram image
    ax.imshow(data, cmap='gray')
    ax.set_ylim(ax.get_ylim()[::-1])

    # plot reg overlay of fibermisalignment
    im = ax.imshow(misalignment, alpha=alpha, cmap=cmap, vmin=vmin, vmax=vmax)

    # add colorbar for misalignment range
    clb = fig.colorbar(im, cax=cax, orientation='vertical',
                       ticks=[vmin, vmin/2, 0, vmax/2, vmax])
    clb.ax.set_title(variable)
    # vertically oriented colorbar
    clb.ax.set_yticklabels(['$<$'+str(vmin), str(vmin/2), '$0$',
                            str(vmax/2), '$>$'+str(vmax)])
    # add plot formatting
    ax.set_xlabel(x_str+'-axis [Volxels]')
    ax.set_ylabel(y_str+'-axis [Volxels]')
    plt.savefig(fig_path + fig_name + '.png')
    plt.show()


def FVF_limits(data, FVF_measured=0.67,
               fig_path='', fig_name='distribution_hist'):
    """ Function for estimating the lower and upper threshold for estimating
        the Fiber Volume Fraction of a dataset with a weak contrast.


    Parameters
    ----------
    data : X-ray CT data [3D array] \n
    FVF_measured : Measured Fiber Volume Fraction [float] \n
    fig_path : Figure directory for saving [str] \n
    fig_name : Figure name used for saving [str] \n

    Returns
    -------
    data : Updated X-ray CT data [3D array] \n
    LT : Lower threshold used for assigning matrix properties \n
    UT : Upper threshold used for assigning fiber properties \n
    FVF_model : Fiber Volume Fraction used for the model \n

    """
    figure, ax = plt.subplots(figsize=(6, 6))
    y, x, _ = plt.hist(data.ravel(), bins=255,
                       alpha=.3, label='data', density=True)
    x = (x[1:] + x[:-1]) / 2

    def gauss(x, mu, sigma, A):
        """ Gaussian dsitribution function for normal distribution


        Parameters
        ----------
        x: Bin edges of the histogram  [Array] \n
        mu: Mean value of distribution  [float] \n
        sigma: Standard deviation of distribution  [float] \n
        A: Constant multiplied with the normal distribution  [float] \n

        Returns
        -------
        Gaussian distribution function

        """
        return A*np.exp(-(x-mu)**2/2/sigma**2)

    def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
        """ Sum of Gaussian distribution functions


        Parameters
        ----------
        x: Bin edges of the histogram  [Array] \n
        mu1: Mean value of distribution 1 [float] \n
        sigma1: Standard deviation of distribution 1 [float] \n
        A1: Constant multiplied with the normal distribution 1  [float] \n
        mu2: Mean value of distribution 2  [float] \n
        sigma2: Standard deviation of distribution 2 [float] \n
        A2: Constant multiplied with the normal distribution 2 [float] \n

        Returns
        -------
        Sum of Gaussian distribution functions
        """
        return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

    expected = (75, 20, 20000, 170, 50, 94000)
    params, cov = curve_fit(bimodal, x, y, expected)
    sigma = np.sqrt(np.diag(cov))
    x_fit = np.linspace(x.min(), x.max(), 500)
    ax.plot(x_fit, bimodal(x_fit, *params), color='red',
            lw=3, label='Fitted model')
    ax.plot(x_fit, gauss(x_fit, *params[:3]), color='red',
            lw=1, ls="--", label='Matrix')
    ax.plot(x_fit, gauss(x_fit, *params[3:]), color='red',
            lw=1, ls=":", label='Fiber')

    plt.legend()
    print(pd.DataFrame(data={'params': params, 'sigma': sigma},
                       index=bimodal.__code__.co_varnames[1:]))
    plt.savefig(fig_path + fig_name + '.png', dpi=70)

    thres_lower = params[3] - 3 * params[4]
    data[data > 254] = thres_lower

    bin_div = 5
    fx_init = np.linspace(0, 254, 254*bin_div + 1) 
    hist, bin_edges = np.histogram(data.ravel(),
                                   bins=fx_init.size, density=True)

    Vf = 1.
    dVf = 1 - FVF_measured
    LT = thres_lower
    UT = LT + 50
    break_while = False
    while Vf > FVF_measured - 0.2:
        fxm = np.zeros(fx_init[fx_init < LT].size)
        fxc = np.linspace(0, 1,
                          fx_init[(fx_init >= LT) & (fx_init <= UT)].size)
        fxf = np.ones(fx_init[fx_init > UT].size)
        fx = np.concatenate((fxm, fxc, fxf))
        if break_while is True:
            print('Fiber Volume Fraction used for modelling = ' + str(Vf))
            break

        if abs(np.sum(hist*fx)/bin_div - FVF_measured) > dVf:
            UT -= 1/bin_div
            break_while = True
        else:
            dVf = abs(np.sum(hist*fx)/bin_div - FVF_measured)
            Vf = np.sum(hist*fx)/bin_div
            UT += 1/bin_div

    FVF_model = Vf
    return data, LT, UT, FVF_model
