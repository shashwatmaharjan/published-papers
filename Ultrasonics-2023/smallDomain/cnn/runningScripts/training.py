# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from tensorflow.keras import optimizers, models, layers, metrics, losses
# import tensorflow_addons as tfa
import time
import os

# Define Fonts for graphs
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
FONT_SIZE = 15
FONT_SIZE2 = 12

# Define necessary function
# Define function to clear terminal window
def clear_screen():
    os.system('clear')

# Define function to get displacement and void data
def get_data(file_directory):

    # Get the training displacement and void data
    displacement_data_training = loadmat(file_directory + '/displacementTrainingData.mat')
    displacement_data_training = np.array(tuple(displacement_data_training['data']))

    void_data_training = loadmat(file_directory + '/voidTrainingData.mat')
    void_data_training = np.array(tuple(void_data_training['data']))

    # Get the validation displacement and void data
    displacement_data_validation = loadmat(file_directory + '/displacementValidationData.mat')
    displacement_data_validation = np.array(tuple(displacement_data_validation['data']))

    void_data_validation = loadmat(file_directory + '/voidValidationData.mat')
    void_data_validation = np.array(tuple(void_data_validation['data']))

    # Get the test displacement and void data
    displacement_data_test = loadmat(file_directory + '/displacementTestData.mat')
    displacement_data_test = np.array(tuple(displacement_data_test['data']))

    void_data_test = loadmat(file_directory + '/voidTestData.mat')
    void_data_test = np.array(tuple(void_data_test['data']))

    print('\nOrganized Data Received.\n')

    return displacement_data_training, void_data_training, displacement_data_validation, void_data_validation, displacement_data_test, void_data_test

# Define function to normalize the input dataset
def normalize_displacement_data(displacement_data_used_to_normalize, data_to_be_normalized, type_of_data_subset):

    displacement_data_mean = np.mean(displacement_data_used_to_normalize)
    displacement_data_range = np.max(displacement_data_used_to_normalize) - np.min(displacement_data_used_to_normalize)

    data_normalized = (data_to_be_normalized - displacement_data_mean) / displacement_data_range
    normalizing_parameters = np.array(list([displacement_data_mean, displacement_data_range]))

    print('Normalized %s dataset.' %(type_of_data_subset))

    return data_normalized, normalizing_parameters

def save_mat_file(values, file_name, file_directory):

    # Save the file as a .mat file
    savemat(file_directory + '/' + file_name, {'data': values})

    print('Saved data to %s' %(file_directory))

# Define function to set up a CNN network
def cnn_model(train_x, train_y, val_x, val_y, test_x, test_y):

    # Record time that it takes to train
    start_train = time.time()

    model = models.Sequential([

        # Input Layer
        layers.Input(shape=train_x.shape[1:]),

        # Convolution Layer
        layers.Conv1D(filters=150, kernel_size=10, padding='same',
                                kernel_initializer='he_uniform'),

        layers.LeakyReLU(),

        # Add a MaxPooling layer
        layers.MaxPooling1D(pool_size=7, strides=None, padding='same'),

        # Add BatchNormalization
        layers.BatchNormalization(),

        layers.Dropout(0.5),

        # Add a Flatten layer
        layers.Flatten(),

        # Add hidden layers
        layers.Dense(4250),

        layers.LeakyReLU(),

        # Add BatchNormalization
        layers.BatchNormalization(),

        layers.Dense(len(train_y[1]), activation='sigmoid')])

    model.compile(optimizer=optimizers.Nadam(learning_rate=0.01),
                    loss=losses.BinaryCrossentropy(name='binary_crossentropy'),
                        metrics=metrics.Precision(name='precision'))

    history_cnn = model.fit(train_x, train_y,
                                epochs=20, batch_size=32,
                                    validation_data=(val_x, val_y))

    # End recording training time
    end_train = time.time()

    # Evaluate performance of CNN on test dataset
    model.evaluate(test_x, test_y)

    # Observe summary of the model
    model.summary()

    training_time = end_train - start_train
    print('CNN trains for %f seconds' %(training_time))

    return model, history_cnn

# Define function to plot loss function and evaluation metric for training and validation dataset
def plot_loss(history_cnn, file_directory_to_save_data):

    loss_name = list(history_cnn.history.keys())[0]

    loss = history_cnn.history[loss_name]
    val_loss = history_cnn.history['val_' + loss_name]

    loss_plot = plt.figure()
    epochs = range(1, len(loss) + 1)
    # plt.yscale('log')
    plt.plot(epochs, loss, 'bo--', label='Training Data', markersize=2)
    plt.plot(epochs, val_loss, 'go--', label='Validation Data', markersize=2)
    # plt.title('Training and Validation Binary Crossentropy')
    plt.xlabel('Epochs', fontsize=FONT_SIZE)
    plt.ylabel('Binary Crossentropy', fontsize=FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
    plt.legend(['Training', 'Validation'], fontsize = FONT_SIZE)
    plt.savefig(file_directory_to_save_data + '/Loss.pdf', bbox_inches='tight')

    return loss_plot

def plot_evaluation_metric(history_cnn, file_directory_to_save_data):

    evaluation_metric_name = list(history_cnn.history.keys())[1]

    evaluation_metric_plot = plt.figure()
    metric = history_cnn.history[evaluation_metric_name]
    val_metric = history_cnn.history['val_' + evaluation_metric_name]
    epochs = range(1, len(metric) + 1)
    # plt.yscale('log')
    plt.plot(epochs, metric, 'bo--', label='Training Data', markersize=2)
    plt.plot(epochs, val_metric, 'go--', label='Validation Data', markersize=2)
    # plt.title('Training and Validation Precision')
    plt.xlabel('Epochs', fontsize=FONT_SIZE)
    plt.ylabel('Precision', fontsize=FONT_SIZE)
    plt.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
    plt.legend(['Training', 'Validation'], fontsize = FONT_SIZE)
    plt.savefig(file_directory_to_save_data + '/Precision.pdf', bbox_inches='tight')

    return evaluation_metric_plot

# Define a main function
def main():

    # Clear screen before startup
    clear_screen()
    print('Running training.py')

    # Assign the file directory as a string to a variable for re-usability
    main_file_directory = 'YOUR_FILE_DIRECTORY_HERE'
    randomized_data_file_directory = main_file_directory + '/dataset/divideddDataset'
    file_directory_to_save_results = main_file_directory + '/cnn/trainingResults'

    # Get displacement and void data 
    displacement_data_training, void_data_training, displacement_data_validation, void_data_validation, \
        displacement_data_test, void_data_test = get_data(randomized_data_file_directory)
    
    # Normalize the displacement datasets
    displacement_training_data_normalized, normalizing_parameters = \
        normalize_displacement_data(displacement_data_used_to_normalize = displacement_data_training, \
                                        data_to_be_normalized = displacement_data_training, type_of_data_subset = 'training')
    displacement_val_data_normalized = \
        normalize_displacement_data(displacement_data_used_to_normalize = displacement_data_training, \
                                        data_to_be_normalized = displacement_data_validation, \
                                            type_of_data_subset = 'validation')[0]
    displacement_test_data_normalized = \
        normalize_displacement_data(displacement_data_used_to_normalize = displacement_data_training, \
                                        data_to_be_normalized = displacement_data_test, \
                                            type_of_data_subset = 'test')[0]

    # Train the CNN model
    model, history_cnn = cnn_model(displacement_training_data_normalized, void_data_training, displacement_val_data_normalized, void_data_validation,
        displacement_data_test, void_data_test)
    
    # Plot loss and evaluation metrics
    loss_plot = plot_loss(history_cnn, file_directory_to_save_results)
    evaluation_metric_plot = plot_evaluation_metric(history_cnn, file_directory_to_save_results)

    # Save normalizing parameters
    save_mat_file(normalizing_parameters, 'normalizing_parameters.mat', file_directory_to_save_results)

    # Save the trained model
    model.save(file_directory_to_save_results + '/trainedCNNModel.h5')
    print('Saved trainedCNNModel.h5 to disk.')

    # Show all existing graphs
    plt.show()

if __name__ == '__main__':
    main()
