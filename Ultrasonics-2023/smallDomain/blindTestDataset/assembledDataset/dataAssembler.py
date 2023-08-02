# Import the required libraries
import numpy as np
import mat4py
import scipy

from scipy.io import savemat
import matplotlib.pyplot as plt
import os

# Define necessary function
# Define function to clear terminal window
def clear_screen():
    os.system('clear')

# Get data from different voids
def get_displacement_data(file_directory, num_samples):

    num_sensors = 41
    num_timesteps = 701
    displacement_data = np.zeros((num_samples, num_sensors, num_timesteps))
    void_data = np.zeros((num_samples, 5000))

    for n_sample in range(num_samples):

        print("Receiving data from sample %d..." %(n_sample+1))

        uy_data = scipy.io.loadmat(file_directory + "/uy_ansys_blind_test_v" + str(n_sample+1) + ".mat")
        uy_data = uy_data["uy_ansys"]

        displacement_data[n_sample, :, :] = uy_data

        indiv_void_data = scipy.io.loadmat(file_directory + "/Blind_test_Ex_" + str(n_sample+1) + ".mat")
        indiv_void_data = indiv_void_data["Crack_data_yes_no_1"]

        void_data[n_sample, :] = indiv_void_data

    return displacement_data, void_data
    # return displacement_data

def displacement_data_realign(displacement_data):

    print("\nRealigning data...")

    # Specify how many sensors the displacement data has
    num_sensors = np.shape(displacement_data)[1]

    # Time stamp count
    time_step_count = 3 # Skip every steps
    displacement_data = displacement_data[:, :, 1::time_step_count]

    # To change the channel and timestep positions
    displacement_data_reshaped = np.zeros((displacement_data.shape[0], displacement_data.shape[2], displacement_data.shape[1]))

    # Make the timestep values the second dimension and the sensor data as the channels
    print("Changing sensor and timestamp positions...")
    for i in range(displacement_data.shape[0]):
        for j in range(displacement_data.shape[1]):
            for k in range(displacement_data.shape[2]):
                displacement_data_reshaped[i, k, j] = displacement_data[i, j, k]

    return displacement_data_reshaped

# Define function to save the dataset
def save_dataset(displacement_data, displacement_file_name, file_directory_to_save_data):

    savemat(file_directory_to_save_data + displacement_file_name, {"data": displacement_data})
    print("\nSaved " + displacement_file_name)

# Define main function
def main():

    # Clear screen before startup
    clear_screen()
    print("Running dataRandomizer.py\n")

    # Assign the file directory as a string to a variable for re-usability
    main_file_directory = '/media/swimlab/8e0a5339-75ae-4b57-aaae-375e5bb09ac3/ML_Projects/anisotropic/iteration4/blindTestDataset'
    raw_data_file_directory = main_file_directory + '/rawDataset'
    file_directory_to_save_data = main_file_directory + '/assembledDataset'

    # Get displacement and void data 
    displacement_data, void_data = get_displacement_data(raw_data_file_directory, 6)
    # displacement_data = get_displacement_data(raw_data_file_directory, 6)
    
    # Realign the displacement data
    displacement_data = displacement_data_realign(displacement_data)

    # Save the divided dataset
    print("Saving datasets...")
    save_dataset(displacement_data, "/blindDisplacementData.mat", file_directory_to_save_data)
    save_dataset(void_data, "/blindVoidData.mat", file_directory_to_save_data)


if __name__ == '__main__':
    main()
