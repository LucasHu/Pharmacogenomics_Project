from __future__ import print_function

import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import RMSprop
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import sys
import numpy as np



def deep_learning(x_train, x_test, y_train, y_test):

    batch_size = 128
    num_classes = 2
    epochs = 10

    # the data, split between train and test sets

    #x_train = x_train.reshape(60000, 784)
    #x_test = x_test.reshape(10000, 784)
    x_train = x_train.astype('float32')
    x_test = x_test.astype('float32')

    print(x_train.shape[0], 'train samples')
    print(x_test.shape[0], 'test samples')

    # y_train = keras.utils.to_categorical(y_train, num_classes)
    # y_test = keras.utils.to_categorical(y_test, num_classes)
    #
    # model = Sequential()
    # model.add(Dense(512, activation='relu', input_shape=(3594,)))
    # model.add(Dropout(0.2))
    # model.add(Dense(512, activation='relu'))
    # model.add(Dropout(0.2))
    # model.add(Dense(num_classes, activation='softmax'))
    #
    # model.summary()
    #
    # model.compile(loss = 'categorical_crossentropy',
    #             optimizer = RMSprop(),
    #             metrics=['accuracy'])
    #
    # history = model.fit(x_train, y_train,
    #                     batch_size=batch_size,
    #                     epochs=epochs,
    #                     verbose = 1,
    #                     validation_data=(x_test,y_test))
    #
    #score=model.evaluate(x_test, y_test, verbose=0)
    #print('Test loss', score[0])
    #print('Test accuracy:', score[1])

    #########
    # Generate dummy data
    #x_train = np.random.random((1000, 20))
    #y_train = np.random.randint(2, size=(1000, 1))
    #x_test = np.random.random((100, 20))
    #y_test = np.random.randint(2, size=(100, 1))

    model = Sequential()
    model.add(Dense(64, input_dim=x_train.shape[1], activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))

    model.compile(loss='binary_crossentropy',
                  optimizer='rmsprop',
                  metrics=['accuracy'])

    print (x_train)
    print (y_train)
    print (type(x_train))
    print (type(y_train))
    print (np.shape(x_train))
    print (np.shape(y_train))

    model.fit(x_train, y_train,
              epochs=20,
              batch_size=128)

    score = model.evaluate(x_test, y_test, batch_size=128)
    print('Test loss', score[0])
    print('Test accuracy:', score[1])

def random_forest(x_train, x_test, y_train, y_test):
    rf = RandomForestClassifier(n_estimators=100, oob_score=True)
    rf.fit(x_train, y_train)
    predicted = rf.predict(x_test)
    accuracy = accuracy_score(y_test, predicted)
    print('Mean accuracy score: '+ str(accuracy) )

def random_forest_prediction(x, y, exp_data, input_molecule_list, protein_list, out_file):
    rf = RandomForestClassifier(n_estimators=100, oob_score=True)
    rf.fit(x, y)
    predicted_score = rf.predict_proba(exp_data)

    text_file = open(out_file, 'w')

    len_protein = len(protein_list)

    for j in range(len(input_molecule_list)):
        for i in range(len(protein_list)):
            proba = predicted_score[j * len_protein + i]
            molecule = input_molecule_list[j]
            protein = protein_list[i]
            text_file.write(str(molecule)+ "\t" + str(protein) + "\t" + str(proba[1]))
            text_file.write("\n")


    return predicted_score



