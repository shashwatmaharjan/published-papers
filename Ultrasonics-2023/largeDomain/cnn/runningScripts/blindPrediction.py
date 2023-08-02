# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.io import loadmat
from tensorflow.keras import optimizers, models, layers
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

# Define function to get the saved trained model
def get_model(file_directory):

    # return models.load_model(file_directory + '/trainedCNNModel.h5')
    return models.load_model(file_directory + '/bestTrainedCNNModel.h5')

# Define function to get data for specific voids
def get_data(file_directory):

    displacement_data = loadmat(file_directory + '/blindDisplacementData.mat')
    displacement_data = np.array(tuple(displacement_data['data']))

    void_data = loadmat(file_directory + '/blindVoidData.mat')
    void_data = np.array(tuple(void_data['data']))

    return displacement_data, void_data

# Define function to normalize the input dataset
def normalize_displacement_data(data_to_be_normalized, normalizing_parameters):

    displacement_data_mean = normalizing_parameters[0, 0]
    displacement_data_range = normalizing_parameters[0, 1]

    data_normalized = (data_to_be_normalized - displacement_data_mean) / displacement_data_range

    return data_normalized

# Define a function to predict test data using trained model
def predict_from_data(model, normalized_displacement_data, threshold):

    # Make a prediction on the specific void dataset
    probabilistic_predicted_data = model.predict(normalized_displacement_data)
    exact_predicted_data = np.where(probabilistic_predicted_data > threshold, 1, 0)

    return probabilistic_predicted_data, exact_predicted_data

# Define function to calculate the classification metrics
def metrics_calculator(target_data, predicted_data):

    # True positive
    tp = 0
    # True negative
    tn = 0
    # False positive
    fp = 0
    # False positive
    fn = 0

    # Accuracy, Precision, Recall, and F1-score
    accuracy = []
    precision = []
    recall = []

    for i in range(target_data.shape[0]):
        for j in range(target_data.shape[1]):
            if target_data[i, j] == predicted_data[i, j]:
                if target_data[i, j] == 0 and predicted_data[i, j] == 0:
                    tn = tn + 1

                elif target_data[i, j] == 1 and predicted_data[i, j] == 1:
                    tp = tp + 1

            else:
                if target_data[i, j] == 1 and predicted_data[i, j] == 0:
                    fp = fp + 1

                elif target_data[i, j] == 0 and predicted_data[i, j] == 1:
                    fn = fn + 1

        accuracy.append((tp + tn) / (tp + tn + fp + fn))
        precision.append((tp)/(tp+fp))
        recall.append((tp) / (tp + fn))

    accuracy = np.array(accuracy)
    precision = np.array(precision)
    recall = np.array(recall)
    f1_score = 2 * precision * recall / (precision + recall)

    return accuracy, precision, recall, f1_score

# Plot prediction using a particular metric
def plot_prediction(sample_number, predicted_data, file_directory):

    plot = plt.figure()
    ax = plt.axes()
    im = plt.imshow(predicted_data[sample_number-1, :, :], cmap='binary', origin='lower', vmin=0, vmax=1, extent=[0, 1, 0, 0.5])
    plt.axhline(y=0.15, color='k')
    plt.axhline(y=0.35, color='k')
    plt.xlabel('x [m]', fontsize=FONT_SIZE)
    plt.ylabel('y [m]', fontsize=FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
    cax = plot.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax)
    plt.savefig(file_directory + '/' + 'pred_blind_sample_' + str(sample_number) + '.pdf', bbox_inches='tight')

    return plot

# Define a main function
def main():

    # Clear screen before startup
    clear_screen()
    print('Running blindPrediction.py')

    # Assign the file directory as a string to a variable for re-usability
    main_file_directory = '/media/swimlab/8e0a5339-75ae-4b57-aaae-375e5bb09ac3/ML_Projects/anisotropic/iteration2'
    organized_data_file_directory = main_file_directory + '/dataset/organizedDataset'
    blind_test_data_file_directory = main_file_directory + '/blindTestDataset/assembledDataset'
    trained_results_file_directory = main_file_directory + '/CNN/trainingResults'
    file_directory_to_save_data = main_file_directory + '/CNN/predictionResults/blindTestData'

    # Get the trained CNN model
    model = get_model(trained_results_file_directory)
    print(model.summary())

    # Get displacement and void data
    displacement_data, void_data = get_data(blind_test_data_file_directory)

    # Get the normalizing parameters
    normalizing_parameters = loadmat(trained_results_file_directory + '/normalizing_parameters.mat')
    normalizing_parameters = np.array(tuple(normalizing_parameters['data']))

    # Normalize the displacement data
    displacement_test_data_normalized = normalize_displacement_data(displacement_data, normalizing_parameters)

    # Predict using the trained model
    probabilistic_predicted_data, exact_predicted_data = predict_from_data(model, displacement_test_data_normalized, threshold = 0.4)
    
    # Get the four classficationmetrics after prediction
    accuracy, precision, recall, f1_score = metrics_calculator(target_data = void_data, predicted_data = exact_predicted_data)

    # Reshape the flattened domain to the original 2D domain
    void_data = np.reshape(void_data, (np.shape(void_data)[0], 50, 100))
    probabilistic_predicted_data = np.reshape(probabilistic_predicted_data, (np.shape(probabilistic_predicted_data)[0], 50, 100))

    for n_sample in range(np.shape(displacement_data)[0]):
    
        sample_number = n_sample+1
        plot = plot_prediction(sample_number, probabilistic_predicted_data, file_directory_to_save_data)

        print('Accuracy = %0.2f, Precision = %0.2f, Recall = %0.2f, F1-score = %0.2f' %(accuracy[n_sample]*100, precision[n_sample]*100, recall[n_sample]*100, f1_score[n_sample]*100))

    # plt.show()

if __name__ == '__main__':
    main()
