# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.io import loadmat
from tensorflow.keras import optimizers, models, layers
import os
from tensorflow.keras.utils import plot_model

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
def get_data(file_directory, void_number):
    void_number = str(void_number)

    displacement_data = loadmat(file_directory + '/void_' + void_number + '_displacement_test.mat')
    displacement_data = np.array(tuple(displacement_data['data']))

    void_data = loadmat(file_directory + '/void_' + void_number + '_void_test.mat')
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

# Define a function to just plot void 0 prediction
def plot_void_0(target_data, predicted_data, file_directory):

    sample_num = 10

    actual_plot = plt.figure()
    ax = plt.axes()
    im1 = plt.imshow(target_data[sample_num, :, :], cmap='binary', origin='lower', vmin=0, vmax=1, extent=[0, 0.1, 0, 0.05])
    plt.axhline(y=0.015, color='k')
    plt.axhline(y=0.035, color='k')
    plt.xlabel('x [m]', fontsize=FONT_SIZE)
    plt.ylabel('y [m]', fontsize=FONT_SIZE)
    cax = actual_plot.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.title(type + ' Actual')
    cb = plt.colorbar(im1, cax=cax)
    cb.ax.tick_params(labelsize=FONT_SIZE)
    plt.savefig(file_directory + '/void_0_actual.pdf', bbox_inches='tight')

    # See how the predicted crack looks in the samples
    predicted_plot = plt.figure()
    ax = plt.axes()
    im2 = plt.imshow(predicted_data[sample_num, :], cmap='binary', origin='lower', vmin=0, vmax=1, extent=[0, 0.1, 0, 0.05])
    plt.axhline(y=0.015, color='k')
    plt.axhline(y=0.035, color='k')
    plt.xlabel('x [m]', fontsize=FONT_SIZE)
    plt.ylabel('y [m]', fontsize=FONT_SIZE)
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
    cax = predicted_plot.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # plt.title(type + ' Prediction')
    cb = plt.colorbar(im2, cax=cax)
    cb.ax.tick_params(labelsize=FONT_SIZE)
    plt.savefig(file_directory + '/void_0_predicted.pdf', bbox_inches='tight')

    return actual_plot, predicted_plot

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
def plot_prediction(void_number, target_data, predicted_data, accuracy, precision, recall, f1_score, file_directory, type):

    chosen_metric = f1_score

    if type == 'Best':
        comparison_metric = np.max(chosen_metric)
    
    elif type == 'Worst':
        comparison_metric = np.min(chosen_metric)
    
    print("%s prediction:" %(type))

    for i in range(chosen_metric.shape[0]):

        if chosen_metric[i] == comparison_metric:

            actual_plot = plt.figure()
            ax = plt.axes()
            im1 = plt.imshow(target_data[i, :, :], cmap='binary', origin='lower', vmin=0, vmax=1, extent=[0, 0.1, 0, 0.05])
            plt.axhline(y=0.015, color='k')
            plt.axhline(y=0.035, color='k')
            plt.xlabel('x [m]', fontsize=FONT_SIZE)
            plt.ylabel('y [m]', fontsize=FONT_SIZE)
            cax = actual_plot.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            # plt.title(type + ' Actual')
            cb = plt.colorbar(im1, cax=cax)
            cb.ax.tick_params(labelsize=FONT_SIZE)
            ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
            plt.savefig(file_directory + '/void_' + void_number + '_' + type.lower() + '_actual.pdf', bbox_inches='tight')

            # See how the predicted crack looks in the samples
            predicted_plot = plt.figure()
            ax = plt.axes()
            im2 = plt.imshow(predicted_data[i, :], cmap='binary', origin='lower', vmin=0, vmax=1, extent=[0, 0.1, 0, 0.05])
            plt.axhline(y=0.015, color='k')
            plt.axhline(y=0.035, color='k')
            plt.xlabel('x [m]', fontsize=FONT_SIZE)
            plt.ylabel('y [m]', fontsize=FONT_SIZE)
            ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE2)
            cax = predicted_plot.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            # plt.title(type + ' Prediction')
            cb = plt.colorbar(im2, cax=cax)
            cb.ax.tick_params(labelsize=FONT_SIZE)
            plt.savefig(file_directory + '/void_' + void_number + '_' + type.lower() + '_predicted.pdf', bbox_inches='tight')
            print("Accuracy: %.2f, Precision: %.2f, Recall: %.2f, F1-score: %.2f" %(accuracy[i]*100, precision[i]*100, recall[i]*100, f1_score[i]*100))

            return actual_plot, predicted_plot
        
            # worst_plot, (ax1, ax2) = plt.subplots(2, 1, constrained_layout = True)
            # axlist = [ax1, ax2]
            # subfig1 = ax1.imshow(target_data[i, :, :], cmap='binary', origin='lower', vmin=0, vmax=1)
            # ax1.set_title(type + " Actual")
            # subfig2 = ax2.imshow(predicted_data[i, :], cmap='binary', origin='lower', vmin=0, vmax=1)
            # ax2.set_title(type + " Prediction")
            # # worst_plot = plt.colorbar()
            # worst_plot.colorbar(subfig1, ax=axlist)
            # # worst_plot.ax.tick_params(labelsize=FONT_SIZE)
            # plt.savefig(file_directory + '/' + type.lower() + '_void_' + void_number + '.pdf', bbox_inches='tight')
            # print("Accuracy: %.2f, Precision: %.2f, Recall: %.2f, F1-score: %.2f" %(accuracy[i]*100, precision[i]*100, recall[i]*100, f1_score[i]*100))

            # return worst_plot


# Define a main function
def main():

    # Clear screen before startup
    clear_screen()
    print("Running prediction.py")

    # Assign the file directory as a string to a variable for re-usability
    main_file_directory = "/media/swimlab/8e0a5339-75ae-4b57-aaae-375e5bb09ac3/ML_Projects/anisotropic/iteration4"
    organized_data_file_directory = main_file_directory + '/dataset/organizedDataset'
    trained_results_file_directory = main_file_directory + '/CNN/trainingResults'
    file_directory_to_save_data = main_file_directory + '/CNN/predictionResults/testData'

    # Get the trained CNN model
    model = get_model(trained_results_file_directory)

    # Get displacement and void data
    void_number = '0'
    displacement_data, void_data = get_data(organized_data_file_directory, void_number)

    # Get the normalizing parameters
    normalizing_parameters = loadmat(trained_results_file_directory + '/normalizing_parameters.mat')
    normalizing_parameters = np.array(tuple(normalizing_parameters['data']))

    # Normalize the displacement data
    displacement_test_data_normalized = normalize_displacement_data(displacement_data, normalizing_parameters)


    # Predict using the trained model
    probabilistic_predicted_data, exact_predicted_data = predict_from_data(model, displacement_test_data_normalized, threshold = 0.4)
    
    if void_number != "0":

        # Get the four classficationmetrics after prediction
        accuracy, precision, recall, f1_score = metrics_calculator(target_data = void_data, predicted_data = exact_predicted_data)
    
    # Reshape the flattened domain to the original 2D domain
    void_data = np.reshape(void_data, (np.shape(void_data)[0], 50, 100))
    probabilistic_predicted_data = np.reshape(probabilistic_predicted_data, (np.shape(probabilistic_predicted_data)[0], 50, 100))

    if void_number == "0":
        actual_plot, predicted_plot = plot_void_0(void_data, probabilistic_predicted_data, file_directory_to_save_data)
    
    else:

        # Plot the best prediction
        actual_best_plot, predicted_best_plot = \
            plot_prediction(void_number, void_data, probabilistic_predicted_data, accuracy, precision, recall, f1_score, file_directory_to_save_data, 'Best')
        actual_worst_plot, predicted_worst_plot = \
            plot_prediction(void_number, void_data, probabilistic_predicted_data, accuracy, precision, recall, f1_score, file_directory_to_save_data, 'Worst')
        
        # best_plot = plot_prediction(void_number, void_data, probabilistic_predicted_data, accuracy, precision, recall, f1_score, file_directory_to_save_data, 'Best')
        # worst_plot = plot_prediction(void_number, void_data, probabilistic_predicted_data, accuracy, precision, recall, f1_score, file_directory_to_save_data, 'Worst')

    plt.show()

if __name__ == '__main__':
    main()
