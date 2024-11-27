from __future__ import division

import sys
print(sys.path)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import readgadget
import readfof
import MAS_library as MASL
import Pk_library as PKL
import redshift_space_library as RSL



kmin = float(input('kmin = '))
kmax = float(input('kmax = '))
grid_input = int(input('grid ='))

df_path = '/home/gremlin/Msc Project/Bk_tables/redshift_space_new/equil/0.01_0.1'

#----------------------------SAVING TO TABLE------------------------------------------


def save_to_csv(data, file_name_prefix, df_path):
    """
    Saves the data to a CSV file. If the file already exists, adds the data as a new column.

    :param data: DataFrame containing the data to save.
    :param file_name_prefix: Prefix for the CSV file name.
    :param df_path: Directory path where the CSV file will be saved.
    """
    # Construct the file name and path
    file_name = f"{file_name_prefix}.csv"
    file_path = os.path.join(df_path, file_name)

    # Convert the data to a DataFrame if it isn't already
    data_df = pd.DataFrame(data)

    # Flatten any list/array entries to comma-separated strings
    for col in data_df.columns:
        data_df[col] = data_df[col].apply(lambda x: ','.join(map(str, x)) if isinstance(x, (list, np.ndarray)) else x)

    if os.path.exists(file_path):
        # Load the existing CSV into a DataFrame
        existing_df = pd.read_csv(file_path)

        # Check if the column count matches
        if existing_df.shape[0] == data_df.shape[0]:
            # Add the data as a new column
            new_column_name = f'{file_name_prefix}_{len(existing_df.columns) + 1}'
            existing_df[new_column_name] = data_df.iloc[:, 0]
        else:
            raise ValueError("The existing file and new data have different numbers of rows. Cannot append column.")

        # Save the updated DataFrame back to CSV
        existing_df.to_csv(file_path, index=False)
        print(f"New column '{new_column_name}' added to the existing CSV file '{file_name}'.")
    else:
        # Save as a new file if it doesn't exist
        data_df.to_csv(file_path, index=False)
        print(f"CSV file '{file_name}' created with the data.")

    # Print the saved DataFrame
    saved_df = pd.read_csv(file_path)
    print(f"\nSaved DataFrame in {file_name_prefix}:\n", saved_df)






#-----------------------------------------------------------------------
#------------------------------------------------------------------------
# Function to process a single snapshot
def process_snapshot(snapshot):
    #---------------------density field parameters------------------------
    grid    = grid_input  # the density field will have grid^3 voxels
    MAS     = 'CIC'  # Mass-assignment scheme: 'NGP', 'CIC', 'TSC', 'PCS'
    verbose = True   # whether to print information about the progress

    # power spectrum parameters
    axis = 0    # axis along which redshift-space distortions have been placed
    threads = 1 # number of openmp threads to compute the power spectrum

    # -------------------Read Particle positions-----------------------------
    # read the redshift, boxsize, cosmology...etc in the header
    header   = readgadget.header(snapshot)
    BoxSize  = header.boxsize/1e3  # Mpc/h
    Nall     = header.nall         # Total number of particles
    Masses   = header.massarr*1e10 # Masses of the particles in Msun/h
    Omega_m  = header.omega_m      # value of Omega_m
    Omega_l  = header.omega_l      # value of Omega_l
    h        = header.hubble       # value of h
    redshift = header.redshift     # redshift of the snapshot
    Hubble   = 100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)# Value of H(z) in km/s/(Mpc/h)

    print('Processing snapshot: %s' % snapshot)
    print('BoxSize = %.3f Mpc/h'%BoxSize)
    print('Number of particles in the snapshot:', Nall)
    print('Masses of the particles in Msun/h =', Masses)
    print('Omega_m = %.3f'%Omega_m)
    print('Omega_l = %.3f'%Omega_l)
    print('h = %.3f'%h)
    print('redshift = %.1f'%redshift)
    print('\nk1 = ', kmin)
    print('k2 = ', kmax)
    print('grid = ', grid)

    # ------------------Compute 3D density field--------------------------------
    
    # particle positions in 3D
    pos = readgadget.read_block(snapshot, "POS ", [1])/1e3#[1]=DM, [2]=neutrinos
    print('\nDone pos')

    #particle velocities
    vel   = readgadget.read_block(snapshot, "VEL ", [1])
    print('\nDone vel')

    # putting particles into redshift space
    # pos_original = np.copy(pos)
    RSL.pos_redshift_space(pos, vel, BoxSize, Hubble, redshift, axis)#moves particles to redshift-space
    print('\nDone real-space -> redshift-space')

    # Define 3D density field for-----RSD--------------------
    #fill a box with zeros
    delta = np.zeros((grid, grid, grid), dtype=np.float32)#3D array filled with zeros

    # Construct 3D density field
    MASL.MA(pos, delta, BoxSize, MAS, verbose=verbose)

    delta /= np.mean(delta, dtype=np.float64)
    delta -= 1.0

    # ----------------Compute bispectrum-------------------------------------------
    #comment out equilateral/folded/squeezed config as needed
    
    k_range = np.linspace(kmin, kmax, 50) #range for k in h/Mpc

    data = []

    # Iterate over k values
    for k in k_range:

        #equilateral
        k1 = k
        k2 = k
        k3 = k
        
        # #folded
        # k1 = k
        # k2 = k
        # k3 = k/2
        
        # # Squeezed
        # epsilon = 0.2
        # k1 = k
        # k2 = k
        # k3 = epsilon * k  # epsilon is a small number, e.g., 0.01

        #theta = np.array([np.pi/3])#keep theta as 60 deg for equilateral
        theta = np.array([np.arccos((k3**2 - k1**2 - k2**2)/(2*k1*k2))])


        # Compute bispectrum
        BBk = PKL.Bk(delta, BoxSize, k1, k2, theta, MAS, threads)
        Bk = BBk.B  # bispectrum
        Qk = BBk.Q  # reduced bispectrum
        k_Bk = BBk.k  # k-bins for power spectrum
        Pk = BBk.Pk  # power spectrum
        
        # Append results to data list
        data.append({
            'k': k,
            'Bk': Bk,
            'Qk': Qk,
            'k_Bk': k_Bk,
            'Pk': Pk
            })

    # Convert list of dictionaries to DataFrame
    df_results = pd.DataFrame(data)


    #---------------------------------------------------------------


    #save results to csv
    save_to_csv(df_results['Bk'], 'Bk', df_path)
    save_to_csv(df_results['Qk'], 'Qk', df_path)
    #save_to_csv(df_results['Power_Spectrum'], 'Pk', df_path)
    #save_k_bk_to_csv(df_results['k_Bins'], df_path)

    


# Main function to iterate over snapshots from my dir
if __name__ == "__main__":
    base_path = '/home/gremlin/Msc Project/Quijote_simulations2/Snapshots/fiducial'
    for i in range(1000, 1050):
        snapshot_path = os.path.join(base_path, str(i), 'snapdir_004', 'snap_004')
        process_snapshot(snapshot_path)
