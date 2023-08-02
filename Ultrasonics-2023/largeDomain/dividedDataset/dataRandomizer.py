# Import the required libraries
import numpy as np
import mat4py
import scipy
import matplotlib.pyplot as plt
import os

from scipy.io import savemat

# Define necessary function
# Define function to clear terminal window
def clear_screen():
    os.system('clear')

# Get data from different voids
def get_void_data(file_directory, void_number):

    print("\nReceiving data from void %c..." %(void_number))

    print("Retrieving Displacement data...")
    displacement_data = scipy.io.loadmat(file_directory + '/Void' + void_number + '/signal_data.mat')
    displacement_data = displacement_data["signalData"]

    print("Retrieving Void data...")
    void_data = scipy.io.loadmat(file_directory + '/Void' + void_number + '/void_data.mat')
    void_data = void_data["voidData"]

    return displacement_data, void_data

def displacement_data_realign(displacement_data, void_number):

    print("\nRealigning data for void %c..." %(void_number))

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

# Define function to divide dataset into training, validation, and test dataset
def dataset_divider(min_num_data, displacement_data, void_data):

    start_size = 0
    training_index = int((80/100)*min_num_data)
    val_index = training_index + int((10/100)*min_num_data)
    test_index = min_num_data # val_index + int((10/100)*min_num_data)

    displacement_training_data = displacement_data[start_size:training_index, :, :]
    void_training_data = void_data[start_size:training_index, :]

    displacement_val_data = displacement_data[training_index:val_index, :, :]
    void_val_data = void_data[training_index:val_index, :]

    displacement_test_data = displacement_data[val_index:test_index, :, :]
    void_test_data = void_data[val_index:test_index, :]

    return displacement_training_data, void_training_data, displacement_val_data, void_val_data, displacement_test_data, void_test_data

# Define function to stack the dataset vertically
def stack_dataset(void_0_displacement_data, void_1_displacement_data, void_2_displacement_data, void_0_void_data, void_1_void_data, void_2_void_data):

    return np.row_stack((void_0_displacement_data, void_1_displacement_data, void_2_displacement_data)), np.row_stack((void_0_void_data, void_1_void_data, void_2_void_data))

# Define function to randomize dataset
def randomize_dataset(displacement_data, void_data):

    indices_displacement_data = np.arange(displacement_data.shape[0])
    np.random.shuffle(indices_displacement_data)

    return displacement_data[indices_displacement_data], void_data[indices_displacement_data]

# Define function to save the dataset
def save_dataset(displacement_data, void_data, displacement_file_name, void_file_name, file_directory_to_save_data):

    savemat(file_directory_to_save_data + displacement_file_name, {"data": displacement_data})
    print("\nSaved " + displacement_file_name)

    savemat(file_directory_to_save_data + void_file_name, {"data": void_data})
    print("Saved " + void_file_name)

# Define main function
def main():

    # Clear screen before startup
    clear_screen()
    print("Running dataRandomizer.py")

    # Assign the file directory as a string to a variable for re-usability
    main_file_directory = '/media/swimlab/8e0a5339-75ae-4b57-aaae-375e5bb09ac3/ML_Projects/anisotropic/iteration3/dataset'
    raw_data_file_directory = main_file_directory + '/rawDataset'
    file_directory_to_save_data = main_file_directory + '/randomizedDataset'

    # Get displacement and void data 
    displacement_data_void_0, void_data_void_0 = get_void_data(raw_data_file_directory, '0')
    displacement_data_void_1, void_data_void_1 = get_void_data(raw_data_file_directory, '1')
    displacement_data_void_2, void_data_void_2 = get_void_data(raw_data_file_directory, '2')

    # Realign the displacement data
    displacement_data_void_0 = displacement_data_realign(displacement_data_void_0, '0')
    displacement_data_void_1 = displacement_data_realign(displacement_data_void_1, '1')
    displacement_data_void_2 = displacement_data_realign(displacement_data_void_2, '2')

    # Find the minimum data available
    num_data_void_0 = np.shape(displacement_data_void_0)[0]
    num_data_void_1 = np.shape(displacement_data_void_1)[0]
    num_data_void_2 = np.shape(displacement_data_void_2)[0]

    min_num_data = min(num_data_void_0, num_data_void_1, num_data_void_2)
    
    # Divide the data into training, validation, and test dataset
    void_0_displacement_train, void_0_void_train, void_0_displacement_val, void_0_void_val, void_0_displacement_test, void_0_void_test = \
        dataset_divider(min_num_data, displacement_data_void_0, void_data_void_0)
    void_1_displacement_train, void_1_void_train, void_1_displacement_val, void_1_void_val, void_1_displacement_test, void_1_void_test = \
        dataset_divider(min_num_data, displacement_data_void_1, void_data_void_1)
    void_2_displacement_train, void_2_void_train, void_2_displacement_val, void_2_void_val, void_2_displacement_test, void_2_void_test = \
        dataset_divider(min_num_data, displacement_data_void_2, void_data_void_2)
    
    # Assemble the dataset
    displacement_train, void_train = stack_dataset(void_0_displacement_train, void_1_displacement_train, void_2_displacement_train, void_0_void_train, void_1_void_train, void_2_void_train)
    displacement_val, void_val = stack_dataset(void_0_displacement_val, void_1_displacement_val, void_2_displacement_val, void_0_void_val, void_1_void_val, void_2_void_val)
    displacement_test, void_test = stack_dataset(void_0_displacement_test, void_1_displacement_test, void_2_displacement_test, void_0_void_test, void_1_void_test, void_2_void_test)

    # Randomize the displacement and void dataset uniformly
    displacement_train, void_train = randomize_dataset(displacement_train, void_train)
    displacement_val, void_val = randomize_dataset(displacement_val, void_val)
    displacement_test, void_test = randomize_dataset(displacement_test, void_test)

    # Save the divided dataset
    print("Saving datasets...")
    save_dataset(displacement_train, void_train, "/displacementTrainingData.mat", "/voidTrainingData.mat", file_directory_to_save_data)
    save_dataset(displacement_val, void_val, "/displacementValidationData.mat", "/voidValidationData.mat", file_directory_to_save_data)
    save_dataset(displacement_test, void_test, "/displacementTestData.mat", "/voidTestData.mat", file_directory_to_save_data)


if __name__ == '__main__':
    main()
