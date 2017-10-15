#from sklearn import datasets
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import Feature_generation as FG

def load_data(file_directory):

    target_dictionary, protein_list = load_target("/Users/lucasminghu/Desktop/Pharmacogenomics/all.csv")

    protein_similarity_dictionary = {}
    number_file = 0
    for file in os.listdir(file_directory):
        if file.endswith(".csv"):
            with open(file_directory + "/" + file) as fp:
                next(fp)
                for line in fp:
                    line = line.rstrip()
                    items = line.split(",")
                    protein1 = items[0]
                    protein2 = items[1]
                    score = items[2]
                    edge = "\t".join(sorted([protein1, protein2]))
                    if edge not in protein_similarity_dictionary:
                        protein_similarity_dictionary[edge] = []
                        if number_file == 0:
                            protein_similarity_dictionary[edge].append(score)
                        else:
                            protein_similarity_dictionary[edge] += number_file * [0]
                            protein_similarity_dictionary[edge].append(score)
                    else:
                        no_items_before = len(protein_similarity_dictionary[edge])
                        protein_similarity_dictionary[edge] += [0] * (number_file - no_items_before)
                        protein_similarity_dictionary[edge].append(score)
                number_file += 1

    for keys in protein_similarity_dictionary:
        if len(protein_similarity_dictionary[keys]) < number_file:
            protein_similarity_dictionary[keys] += [0] * (number_file - len(protein_similarity_dictionary[keys].split("\t")))

    for keys in protein_similarity_dictionary:
        if len(protein_similarity_dictionary[keys]) != number_file:
            print len(protein_similarity_dictionary[keys])

    #get all the drugs
    all_drug_list = []
    for key in protein_similarity_dictionary:
        keys = key.split("\t")
        if keys[0] not in all_drug_list:
            all_drug_list.append(keys[0])
        if keys[1] not in all_drug_list:
            all_drug_list.append(keys[1])

    #create a matrix to store all similarities for all drugs with all drugs...
    similarity_matrix = np.zeros((len(all_drug_list), len(all_drug_list) * number_file))

    #create target matrix
    target_matrix = np.zeros((len(all_drug_list), len(protein_list)))

    #loop through the list to generate the matrix
    for i in range(len(all_drug_list)):
        drug = all_drug_list[i]

        # the similarity matrix part
        for keys in protein_similarity_dictionary:
            if drug in keys:
                the_two_drugs = keys.split("\t")
                the_two_drugs.remove(drug)
                the_paired_drug = the_two_drugs[0]
                index = all_drug_list.index(the_paired_drug)
                similarity_matrix[i,index * number_file : (index + 1) * number_file] = protein_similarity_dictionary[keys]
                #print protein_similarity_dictionary[keys]

        # the target martix part
        for keys in target_dictionary:
            if drug in target_dictionary[keys]:
                target_matrix[i, protein_list.index(keys)] = 1

        print ("druge " + drug + " is done!")

    np.savetxt('test.txt', target_matrix,delimiter="\t")

    print target_matrix

    #boston = datasets.load_boston()
    #X = boston.data
    #y = boston.target
    #features = boston.feature_names
    #return X, y, features

# load target file, return a target dictionary and a protein list
def load_target(file_dirctory):
    target_dictionary = {}
    with open(file_dirctory) as fp:
        next(fp)
        for line in fp:
            line = line.rstrip()
            items = line.split(",")
            target = items[0]
            drugs = (items[-1]).replace(" ", "") #some drugbankID has white space in the front
            drugs_list = drugs.split(";")

            if target not in target_dictionary:
                target_dictionary[target] = drugs_list

    protein_list = []
    for keys in target_dictionary:
        if keys not in protein_list:
            protein_list.append(keys)

    return target_dictionary, protein_list

def generate_target_matrix(file_dirctory, drugID_list):
    target_dictionary, protein_list = load_target(file_dirctory)

    # set a np.array has the size as No_Drugs * No_Proteins
    target_matrix = np.zeros((len(drugID_list), len(protein_list)))

    for i in range(len(drugID_list)):
        for keys in target_dictionary:
            if drugID_list[i] in target_dictionary[keys]:
                index = protein_list.index(keys)
                target_matrix[i, index] = 1

    return target_matrix, protein_list




def visualize(X, y, features):
    plt.figure(figsize=(20, 5))
    feature_count = X.shape[1]

    # i: index
    for i in range(feature_count):
        plt.subplot(3, 5, i + 1)
        # TODO: Plot feature i against y
        plt.scatter(X[:, i], y)

        # Add x-axis and y-axis labels.
        plt.xlabel(features[i])
        plt.ylabel("house price")

    plt.tight_layout()
    plt.show()


def fit_regression(X, Y):
    # Remember to use np.linalg.solve instead of inverting!

    # Add bias term to the design matrix
    new_X = add_bias(X)

    # Use Moore-Penrose pseudo-inverse of the design matrix to calculate optimal weight vector

    # First, calculate the inverse of transpose(design_matrix) * design_matrix
    inverse_matrix = np.linalg.solve(np.matmul(np.transpose(new_X), new_X), np.eye(np.shape(new_X)[1]))

    # Second, calculate the Moore-Penrose pseudo-inverse of the design matrix
    MPP_inverse = np.matmul(inverse_matrix, np.transpose(new_X))

    # Third, calculate the optimal w vector
    w = np.matmul(MPP_inverse, Y)

    # Return weight vector
    return w


# Build a function to to add bias to the design matrix
def add_bias(X):
    added_term = np.ones(np.shape(X)[0])
    new_X = np.insert(X, 0, added_term, axis=1)
    return new_X


# Build a function to calculate error measurement metric: MSE
def compute_MSE(predicted_y, target_y):
    sum = 0
    for i in range(np.size(target_y)):
        sum = sum + math.pow((predicted_y[i] - target_y[i]), 2)
    MSE = sum / np.size(target_y)
    return MSE


# Build a function to calculate error measurement metric: root_mean_square_error (RMS)
def compute_RMS(predicted_y, target_y):
    RMS = math.sqrt(compute_MSE(predicted_y, target_y))
    return RMS


# Build a function to calculate error measurement metric: Mean absolute error (MAE)
def compute_MAE(predicted_y, target_y):
    sum = 0
    for i in range(np.size(target_y)):
        sum = sum + abs(predicted_y[i] - target_y[i])
    MAE = sum / np.size(target_y)
    return MAE

# Build a function for 5_fold cross_validation
def five_fold_cross_validation(X, y):
    MSE = 0
    RMS = 0
    MAE = 0
    for i in range(5):
        # Generate the indices for training and testing sets based on the 80% vs. 20% splitting using np.random.choice function
        training_index = np.random.choice(np.arange(np.shape(X)[0]), int(round(np.shape(X)[0] * 0.8)),
                                      replace=False)
        testing_index = np.setdiff1d(np.arange(np.shape(X)[0]), training_index)

        # Generate training and testing sets using the previously generated indices
        training_data = X[training_index, :]
        training_targets = y[training_index]

        test_data = X[testing_index, :]
        test_targets = y[testing_index]

        # Fit regression model by training data
        new_testing_data = add_bias(test_data)  # add bias term
        w2 = fit_regression(training_data, training_targets)
        print ("the weight vector generated from training data: ")
        print w2

        predicted_y = np.matmul(new_testing_data, w2)
        print predicted_y

        # calculate MSE
        MSE += compute_MSE(predicted_y, test_targets)
        #print ("Mean Square Error is: " + str(MSE))

        # calculate RMS
        RMS += compute_RMS(predicted_y, test_targets)
        #print ("Root Mean Square Error is: " + str(RMS))

        # calculate MAE
        MAE += compute_MAE(predicted_y, test_targets)
        #print ("Mean Absolute Error is: " + str(MAE))

    average_MSE = MSE / 5
    average_RMS = RMS / 5
    average_MAE = MAE / 5

    print ("Mean Square Error is: " + str(average_MSE))
    print ("Root Mean Square Error is: " + str(average_RMS))
    print ("Mean Absolute Error is: " + str(average_MAE))

    return average_MSE, average_RMS, average_MAE


def main():
    # Load the data
    #X, y, features = load_data("/Users/lucasminghu/Desktop/Pharmacogenomics/similarity_scores_for_fingerprints")
    #print("Features: {}".format(features))

    combined_matrix, drugbank_ID = FG.main()

    target_matrix, protein_list = generate_target_matrix("/Users/lucasminghu/Desktop/Pharmacogenomics/all.csv", drugbank_ID)

    X = combined_matrix
    y = target_matrix[:,0] #just test the first protein


    # Visualize the features
    #visualize(X, y, features)


    # Fit regression model
    #w = fit_regression(X, y)
    #print ("the weight vector generated from all data: ")
    #print w
    text_file = open("MAE_for_each_protein.txt", "w")
    for i in range(len(protein_list)):
        y = target_matrix[:,i]
        average_MSE, average_RMS, average_MAE = five_fold_cross_validation(X,y)
        text_file.write(protein_list[i] + "\t" + str(np.sum(y)) + "\t" + str(average_MSE))
        text_file.write("\n")


    sys.exit()

    # create a file to store the feature selection result.
    text_file = open("feature_selection.txt", "w")
    text_file.write("the feature ommited is:" + "\t" + "resulting MSE")
    text_file.write("\n")

    # Test each feature individually to find out the most important features
    # in each loop, omit one feature, and re-do the linear regression and compute the MSE
    # if the the higher the MSE is, the more important the omitted feature is.
    for i in range(np.size(feature_list)):
        column_number = i + 1
        selected_training = np.concatenate(
            (training_data[:, 0:i], training_data[:, i + 1:np.shape(training_data)[1]]), axis=1)
        selected_testing = np.concatenate((new_testing_data[:, 0:column_number], new_testing_data[:,
                                                                                 column_number + 1:
                                                                                 np.shape(new_testing_data)[
                                                                                     1]]), axis=1)
        print ("the omitted feature: ")
        print ("feature" + str(feature_list[i]))
        w3 = fit_regression(selected_training, training_targets)
        new_predicted_y = np.matmul(selected_testing, w3)
        new_MSE = compute_MSE(new_predicted_y, test_targets)
        print (new_MSE)
        text_file.write("feature"+str(feature_list[i]))
        text_file.write("\t")
        text_file.write(str(new_MSE))
        text_file.write("\n")

if __name__ == "__main__":
    main()

