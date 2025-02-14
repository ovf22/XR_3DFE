import matplotlib.pyplot as plt
import numpy as np

import M4_IntegrationPoints as M4IP

sample_name = 'shell_sample_name'
result_path = '../results/'+sample_name+'_files/figures/'

# %% Define function for reading Abaqus data


def load_AbqData(file_name, start_line, end_line, start_offset, end_offset):
    """
    Parameters
    ----------
    file_name : Abaqus file name [str]
    start_line : Text in line from where to start loading data [str]
    end_line : Text in line from where to end loading data [str]
    start_offset : Offset from start line from where to load data [int]
    end_offset : Offset from end line from where to load data [int]

    Returns
    -------
    ip_data : Integration point data [array]

    """

    with open(file_name, 'r', encoding='cp1257') as f:
        lines = f.readlines()

    # Get start of integration point (IP) list.
    # Get number of first line that contains...
    start_ip = next(i for i, line in enumerate(lines) if start_line in line)
    start_ip += start_offset  # Add offset.
    # Find stop of IPs.
    # Get the first line searching backwards that contains...
    end_ip = len(lines) - next(i for i, line in enumerate(lines[::-1])
                               if end_line in line)
    end_ip -= end_offset  # Subtract offset.

    # Get IP lines.
    ip_lines = lines[start_ip:end_ip]
    # Create array with IPs.
    ip_data = np.array([line.split() for line in ip_lines], dtype=float)
    # ip_indices = ip_coords[:, :2].astype(int)
    ip_data = ip_data[:, 2:]
    return ip_data


# %% Load integration point coordinates
ip_coords = load_AbqData(file_name='Out-'+sample_name+'_INTCOOR.dat',
                         start_line='THE FOLLOWING TABLE IS PRINTED AT '
                         'THE INTEGRATION POINTS FOR ELEMENT TYPE',
                         end_line='THE ANALYSIS HAS BEEN COMPLETED',
                         start_offset=6,
                         end_offset=3)

# %% Load stress components at integration points
S_local = load_AbqData(file_name='Out-'+sample_name+'_S_local.out',
                       start_line='Field Output reported at integration '
                       'points for part',
                       end_line='Minimum',
                       start_offset=5,
                       end_offset=3)

E_local = load_AbqData(file_name='Out-'+sample_name+'_E_local.out',
                       start_line='Field Output reported at '
                       'integration points for part',
                       end_line='Minimum',
                       start_offset=5,
                       end_offset=3)

S_global = load_AbqData(file_name='Out-'+sample_name+'_S_global.out',
                        start_line='Field Output reported at integration '
                        'points for part',
                        end_line='Minimum',
                        start_offset=5,
                        end_offset=3)


# %% Plot stresses in local and global coordinate systems
stresses = ['S11', 'S22', 'S12', 'S13']

SLmax = [max(idx) for idx in zip(*S_local)]
SLmin = [min(idx) for idx in zip(*S_local)]
M4IP.Int_point_plotting(ip_coords, S_local,
                        S_local.shape[1],
                        vmin=SLmin, vmax=SLmax,
                        unit='MPa', variable='',
                        fig_name='fig10_' + sample_name+'_local_S',
                        fig_path=result_path,
                        fig_title='Stress in local coordinate system',
                        figsize=(25, 6),
                        sp_title=stresses)

SGmax = [max(idx) for idx in zip(*S_global)]
SGmin = [min(idx) for idx in zip(*S_global)]
M4IP.Int_point_plotting(ip_coords, S_global,
                        S_global.shape[1],
                        vmin=SGmin, vmax=SGmax,
                        unit='MPa', variable='',
                        fig_name='fig11_' + sample_name+'_global_S',
                        fig_path=result_path,
                        fig_title='Stress in global coordinate system',
                        figsize=(25, 6),
                        sp_title=stresses)

E_local *= 100
strains = ['E11', 'E22', 'E12', 'E13', 'E_max_principal']
Emax = [max(idx) for idx in zip(*E_local)]
Emin = [min(idx) for idx in zip(*E_local)]

M4IP.Int_point_plotting(ip_coords, E_local,
                        E_local.shape[1],
                        vmin=Emin, vmax=Emax,
                        unit='%', variable='$\\epsilon_{ij}$',
                        fig_name='fig12_' + sample_name + '_local_E',
                        fig_path=result_path,
                        fig_title='Strain in local coordinate system',
                        figsize=(25, 6),
                        sp_title=strains)

# %% Plot stress-strain


def FuncRange(x, xRange):
    """ Save data from list within a defined range.

    Parameters
    ----------
    x: Data for extracting range of data [Array] \n
    xRange: Defined range to extract data [list] \n

    Returns
    -------
    ia: Data point from start of range [float] \n
    ib: Data point from end of range [float]
    """
    ia = 0
    ib = len(x)
    for i in range(len(x)):
        if x[i] > xRange[1]:
            ib = i - 1
            break
    for i in range(ib):
        if x[i] < xRange[0]:
            ia = i
    ib = ib + 1
    return ia, ib


LD_file = 'Out-' + sample_name + '_load-disp.out'

disp, load = np.loadtxt(open(LD_file, 'r'), delimiter=' ', usecols=(0, 1),
                        dtype=float, unpack=True)

length, thickness, width = np.loadtxt(sample_name + '_TomoDim.txt') * 1e-3
A = thickness * width

strain = -disp / length * 100
stress = -load / A

stress_max = - np.min(load) / A
strain_max = strain[np.argmin(load)]

StrainRange = [0.05, 0.25]
(ia, ib) = FuncRange(strain, StrainRange)
pEmod = np.polyfit(strain[ia:ib], stress[ia:ib], 1)
Emod = pEmod[0] / 10   # Emod in GPa
Emod_str = r'E-modulus =' + str(format(Emod, '1.2f')) + '$GPa$'

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax = plt.plot(strain, stress,
              'r--', linewidth=2, label=Emod_str)
ax = plt.plot(strain_max, stress_max, 'rx', markersize=12)

plt.xlim([0, np.max(strain) * 1.05])
plt.ylim([0, np.max(stress) * 1.05])
plt.grid()
plt.legend()
plt.xlabel('Strain [$\\%$]')
plt.ylabel('Stress [MPa]')

plt.savefig(result_path + 'fig13_' + sample_name + '_Stress_strain.png')
