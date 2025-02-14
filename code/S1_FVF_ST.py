# This script is executed from a shell script and contains constants
# that are controlled by it. If you want to use this script independently
# of the shell script, simply change the constants "model_height" and
# "kernel_multiplier" to a desired value and run the script.

# %% Load modules for runnig the script
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import os
from structure_tensor import eig_special_3d, structure_tensor_3d
from scipy.signal import fftconvolve

import M1_TomoHandling as M1TH
import M2_Alignment as M2A

# %%  File names and directory
sample_name = 'shell_sample_name'
data_name = 'S50crop'
file_name = data_name + '.nii'

data_path = '../data'
result_path = '../results/'+sample_name+'_files/figures/'
data_file_path = os.path.join(data_path, file_name)

# Assign sample material coordinate system:
# [x, y, z] = [length, thickness, width]
sample_coor_axis = [1, 2, 0]

# %% Load NIfTI data.
nii_file = nib.load(data_file_path)
data_full = np.array(nii_file.dataobj[slice(None), slice(None), slice(None)])


# Change the tomogram coordinate system to the material coordinate system
data_full = np.moveaxis(data_full, [0, 1, 2], sample_coor_axis)

# Define model dimensions used for cropping the original dataset
model_height  = shell_crop

model_aspect_ratio = 1.75  # Aspect ratio between length and height
specimen_aspect_ratio = 5 / 6  # Aspect ratio between height and width

model_width = int(model_height * model_aspect_ratio)
model_depth = int(model_height / specimen_aspect_ratio)

data = M1TH.tomo_crop_center(data_full,
                             width=model_width,
                             height=model_height,
                             depth=model_depth)

# Read meta data.
data_shape = data.shape
data_type = nii_file.get_data_dtype()
VOXEL_SIZE = abs(nii_file.affine[0, 0])

# Set known fiber diameter in micro meters.
FIBER_DIAMETER = 7

# Normalize intensity to grayscale range
I_min, I_max = np.quantile(data, [0.0005, 0.9995])

data[data < I_min] = I_min
data = (data - I_min) * (255) / (I_max - I_min)
data[data > 255] = 255

M1TH.tomo_plot_3(data,
                 x_slice=int(data_shape[0]/2),
                 y_slice=int(data_shape[1]/2),
                 z_slice=int(data_shape[2]/2),
                 figsize=(12, 10),
                 fig_name='fig1_' + sample_name+'_Tomo_fig',
                 fig_path=result_path)

# %% Structure tensor analysis
r = FIBER_DIAMETER / 2 / VOXEL_SIZE
sigma = round(np.sqrt(r**2 / 2), 2)
rho = 4 * sigma

# Copy volume and cast it to 64-bit floating point.
data_f = data.astype(np.float64)

# Calculate S.
S = structure_tensor_3d(data_f, sigma, rho)

truncate = 4
kernel_radius = int(max(sigma, rho) * truncate + 0.5)
print('kernel_radius:', kernel_radius)

S = S[:, kernel_radius:-kernel_radius,
      kernel_radius:-kernel_radius,
      kernel_radius:-kernel_radius]
S.shape

data_s = data_f[kernel_radius:-kernel_radius,
                kernel_radius:-kernel_radius,
                kernel_radius:-kernel_radius]

# %% Eigendecomposition
val, vec = eig_special_3d(S, full=False)

del S
val = val.astype(np.float16)
vec = vec.astype(np.float32)

# Smallest eigenvalue corresponds to dominant direction i.e., fiber direction
# Eigenvectors are returned as vec=[z,y,x], this is flipped back to vec=[x,y,z]
val = np.flip(val, axis=0)
vec = np.flip(vec, axis=0)


# Eigenvectors desribe the dominant material orientation. This is not a unique
# direction as an opposite vector share the same orientation as the eigenvector
# This will cause a noisy image in the visualization of material orientations.
# All eigenvectors are thus multiplied with the sign of the x-component to
# align the vectors in the positive x-direction.
vec *= np.sign(vec[0])[np.newaxis]

# The global fiber orientation may not be aligned with the global material
# coordinate system.
# All eigenvectors may thus be rotated to a new reference axis.

vec_new, rot_angle = M2A.align_with_x(vec)


# %% Calculate orientations
# Calculate fiber mialignment (Azimuth and elevation angles)
theta, phi, phi_xy, phi_xz = M2A.angles_from_vec_atan2(vec_new)

# %% Overlay plot of fiber misalignment
Sz = data_s.shape[2]//2
Sy = data_s.shape[1]//2
M1TH.tomo_plot_3_overlay(data_s, phi_xy, phi_xz, y_slice=Sy, z_slice=Sz,
                         alpha=0.5, vmin=-5, vmax=5, variable='$\\phi_{ij}$',
                         unit='[$^\\circ$]', cmap='coolwarm', title='',
                         fig_name='fig2_' + sample_name + '_Fiber_misalignment_overlay',
                         fig_path=result_path, figsize=(10, 10))

# Calculate the spherical constant from the eigenvectors and define a threshold
# for masking out high density particals
cs = np.true_divide(val[0], val[2], where=(val[0] != 0) | (val[2] != 0))
cs_thres = np.quantile(cs, 0.95)

# %% FVF estimation
FVF_measured = 0.67  # Measured using a burn-off test

# Define the upper and lower thresholds for estimating the FVF distribution
data_s, LT, UT, FVF_model = M2A.FVF_limits(data_s, FVF_measured=FVF_measured,
                                           fig_path=result_path, 
                                           fig_name='fig3_' + sample_name + '_data_intensity_hist')

# Copy the gray scale dataset and normalize according the ramping function
data_bin = np.copy(data_s).astype(float)
data_bin -= LT
data_bin[data_bin <= 0.] = 0.
data_bin[(data_bin > 0.) & (data_bin <= (UT-LT))] /= (UT-LT)
data_bin[data_bin >= (UT-LT)] = 1.

# Plot the normalized dataset (localized FVF distribution)
M1TH.tomo_plot_3(data_bin,
                 x_slice=int(data_bin.shape[0]/2),
                 y_slice=int(data_bin.shape[1]/2),
                 z_slice=int(data_bin.shape[2]/2),
                 figsize=(12, 10),
                 fig_name='fig4_' + sample_name + '_data_bin_FVF',
                 fig_path='',
                 FrameOn=False)

# Plot the localized FVF ontop of the original CT data
M2A.fig_with_colorbar(
    data_s[:, 50, :].T,
    data_bin[:, 50, :].T,
    'Fiber misalignment relative to global fiber direction',
    cmap='coolwarm',
    alpha=0.3,
    vmin=0.0,
    vmax=1.,
    variable='$V_f$ [-]',
    x_str='z',
    y_str='y',
    fig_name='fig5_' + sample_name + '_bin_FVF_overlay',
    fig_path='')

# Define the kernel window used for convolution
kernel_multiplier = shell_kernel
if kernel_multiplier > 99:
    k_param = int(FIBER_DIAMETER / VOXEL_SIZE) * 16
else:
    k_param = int(FIBER_DIAMETER / VOXEL_SIZE) * kernel_multiplier
k_dim = 2*k_param+1
kern = np.ones([k_dim, k_dim, k_dim]) / (k_dim)**3

# Add padding to the dataset to avoid boundary effects from the mean filter
pad_dim = k_dim//2
data_pad = np.pad(data_bin, ((pad_dim, pad_dim),
                             (pad_dim, pad_dim),
                             (pad_dim, pad_dim)),
                  'reflect')

# Convolve the mean filter with the normalized dataset
# to estimate the distributed FVF
FVF = fftconvolve(data_pad.astype('float32'), kern, mode='same')
FVF = FVF[pad_dim:-pad_dim, pad_dim:-pad_dim, pad_dim:-pad_dim]

FVF_mean = FVF.mean()
print('Average FVF after masking', str(round(FVF_mean, 3)))

if kernel_multiplier > 99:
    FVF = np.ones(data_s.shape) * FVF_mean

figure, ax = plt.subplots(figsize=(6, 6))
ax.hist(FVF.ravel(), bins=256, density=True, label='$K_{size}$ = ' +
        str(round(k_dim / FIBER_DIAMETER, 1)) + '[Number of fibers]')
ax.set_xlabel('FVF [-]')
ax.set_ylabel('Density [-]')
plt.savefig(result_path + 'fig6_' + sample_name + '_FVF_hist.png')
plt.show()

# %% Mask field variables and plot
# Mask the elevation angles based on the spherical threshold
phi[cs > cs_thres] = 0.

M1TH.tomo_plot_3_overlay(data_s, FVF, FVF, y_slice=Sy, z_slice=Sz,
                         alpha=0.5, vmin=0, vmax=1, variable='$V_f$',
                         unit='[$-$]', cmap='PRGn', title='',
                         fig_name='fig7_' + sample_name + '_FVF_overlay',
                         fig_path=result_path, figsize=(10, 10))

# %% Scaling padded regions
# Define padding length
PADDING = 2 * 70

# Pad the elevation angles along at the start and end of the dataset (x-axis)
phi_pad = M1TH.scaling_data_fade(phi,
                                 PADDING,
                                 VOXEL_SIZE)
# Pad the FVF along at the start and end of the dataset (x-axis)
FVF_pad = M1TH.scaling_data_fade(FVF,
                                 PADDING,
                                 VOXEL_SIZE,
                                 FVF_measured=FVF_mean)

# %% FE-Model dimensions
# Define the FE-model dimensions based on the data geometries
L_MODEL = len(data_s[:, 0, 0]) * VOXEL_SIZE
T_MODEL = len(data_s[0, :, 0]) * VOXEL_SIZE
W_MODEL = len(data_s[0, 0, :]) * VOXEL_SIZE

# Save model dimensions to file. The file is used in S2_Cube.py
MODEL_DIM = [L_MODEL, T_MODEL, W_MODEL]
MODEL_DIM_FILE = sample_name+'_TomoDim.txt'
np.savetxt(MODEL_DIM_FILE, MODEL_DIM, delimiter=';')

# %% Save variables and constants for mapping orientations to integration
np.savez(sample_name+'_MAP_VAR.npz', VOXEL_SIZE=VOXEL_SIZE,
         MODEL_DIM=MODEL_DIM, phi=phi_pad, theta=theta, FVF=FVF_pad)
