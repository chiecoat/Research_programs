# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:54:37 2022

@author: chiec
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:19:44 2022

@author: chiec


This file will read in data for bubble shapes, sizes and shape parameters 
"""
#basic needs for math and reading text, respectively
import numpy as np
import pandas as pd
#this is needed to prepare the data for the neural net
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
#TensorFlow performs the ML 
import tensorflow as tf
#needed to plot outputs of ML data
import matplotlib.pyplot as plt
from matplotlib import rcParams

#needed to test the success of the algorithm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, precision_score, recall_score

#read in data from csv file
bub_data_train_in=pd.read_csv('E:/Chieco/NASA foam/Images foam/ttab_25%/machine learning lists/train_foam_ttab25_bubble_data_all_ages.csv')

#there are some values that are Infinite. They are replaced with NA and then 
#removed
bub_data_train_in.replace([np.inf, -np.inf], np.nan, inplace=True)
# Drop rows with NaN
bub_data_train_in.dropna(inplace=True)
bub_data_train_in = bub_data_train_in.dropna()
bub_data_train=bub_data_train_in.drop(['frame','x_cens','y_cens'],axis=1)

#this elimiates the bubble_yes_no column from the data
X_train = bub_data_train.drop('bubble_yes_no', axis=1)
#the only result we train for is whether a region is a bubble or not
Y_train = bub_data_train['bubble_yes_no']

"""
this will output [x,y] columns ot training data and test data where the split
is 80% training data and 20% test data. The random state variable shuffles the
data randomly but by picking an ineger we can return to the same random state 
for reproducability

X_train, X_test, y_train, y_test = train_test_split(
    X, y, 
    test_size=0.2, random_state=42
)
"""

#here we rescale the data so the values are on the same order of magnitude
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)

#Now we perform the learning
#randomizes the initial learning seed
tf.random.set_seed(42)

"""
there are three layers of the neural net. The first passes in the inital data 
into a 128 neuron, then thwo dense layers with 256 neurons each. The final 
output is a one layer sigmoid to separate data in
"""
model = tf.keras.Sequential([
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dense(256, activation='relu'),
    tf.keras.layers.Dense(256, activation='relu'),
    tf.keras.layers.Dense(1, activation='sigmoid')
])
"""
These are the outputs from the model. Since we have a yes/no question we are
asking for binary_crossentropy. 'Adam' is the optimize function and a learning
rate of 0.03. The last three are metrics we can measure for success.
"""
model.compile(
    loss=tf.keras.losses.binary_crossentropy,
    optimizer=tf.keras.optimizers.Adam(lr=0.03),
    metrics=[
        tf.keras.metrics.BinaryAccuracy(name='accuracy'),
        tf.keras.metrics.Precision(name='precision'),
        tf.keras.metrics.Recall(name='recall')
    ]
)

n_epochs=100
history = model.fit(X_train_scaled, Y_train, epochs=n_epochs)


#now we plot the results from training
#this just makes the figure a aprticular size
rcParams['figure.figsize'] = (18, 8)

#if loss plateau's in the plot we have run enough epochs
plt.plot(
    np.arange(1, n_epochs+1), 
    history.history['loss'], label='Loss'
)

#these are the metrics for our success. We'd like them as close to 1 as 
#possible
plt.plot(
    np.arange(1, n_epochs+1), 
    history.history['accuracy'], label='Accuracy'
)
plt.plot(
    np.arange(1, n_epochs+1), 
    history.history['precision'], label='Precision'
)
plt.plot(
    np.arange(1, n_epochs+1), 
    history.history['recall'], label='Recall'
)
#This is perfct accuracy
plt.plot(np.arange(1, n_epochs+1),1+0*np.arange(1, n_epochs+1))

plt.title('Evaluation metrics', size=20)
plt.xlabel('Epoch', size=14)
plt.legend();

"""
now we test how the model does when given new data. This is from our test set
data which is a different data set than the training data
"""
#read in data from csv file
bub_data_test_in=pd.read_csv('E:/Chieco/NASA foam/Images foam/ttab_25%/machine learning lists/test_foam_ttab25_bubble_data_all_ages.csv')

#there are some values that are Infinite. They are replaced with NA and then 
#removed
bub_data_test_in.replace([np.inf, -np.inf], np.nan, inplace=True)
# Drop rows with NaN
bub_data_test_in.dropna(inplace=True)
bub_data_test_in = bub_data_test_in.dropna()
bub_data_test=bub_data_test_in.drop(['frame','x_cens','y_cens'],axis=1)

#this elimiates the bubble_yes_no column from the data
X_test = bub_data_test.drop('bubble_yes_no', axis=1)
#the only result we train for is whether a region is a bubble or not
Y_test = bub_data_test['bubble_yes_no']

#here we rescale the data so the values are on the same order of magnitude
X_test_scaled = scaler.transform(X_test)

#predictions outputs the probability a region is or is not a bubble. We need
#to switch the value to binary 1 for is a bubble and 0 for is not a bubble
predictions = model.predict(X_test_scaled)

prediction_classes = [
    1 if prob > 0.5 else 0 for prob in np.ravel(predictions)
]

"""
Now we test the model success. To see if data is overfit
"""
print(confusion_matrix(Y_test, prediction_classes))
print(f'Accuracy: {accuracy_score(Y_test, prediction_classes):.2f}')
print(f'Precision: {precision_score(Y_test, prediction_classes):.2f}')
print(f'Recall: {recall_score(Y_test, prediction_classes):.2f}')

predictions_arr = [float(i) for i in prediction_classes]
predictions_arr=np.atleast_2d(predictions_arr).T

bub_data_test_out=np.append(bub_data_test_in, predictions_arr, 1)

filename_out='E:/Chieco/NASA foam/Images foam/ttab_25%/machine learning lists/ML_out_test_foam_ttab25_bubble_data_all_ages.csv'
pd.DataFrame(bub_data_test_out).to_csv(filename_out, index_label='bub_id', header  = ['frame','x_cens','y_cens','areas','equiv_radius','perimeter','circularity','eccentricity','major_axis','minor_axis','bubble_yes_no_test','bubb_yes_no_ML'])