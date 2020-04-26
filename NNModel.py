
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

from keras.optimizers import SGD
from keras.models import Sequential
from keras.layers import Dense, Activation,Dropout
from keras.models import model_from_json

data=pd.read_csv('/Users/siobhanmcloughlin/RubyProjects/BioHackCPH/Pinguestigators/PancreasDataSio.csv')
        model.add(Dense(1,  kernel_initializer='uniform', activation='softmax'))
        return model

net=makeModel()
# net.compile(loss='categorical_crossentropy', optimizer='Adam')
net.compile(loss='mean_squared_error', optimizer='sgd')
Xt, Xv, Yt, Yv= train_test_split(X_train,Y_train, test_size=.2)
history = net.fit(Xt, Yt, epochs=300, verbose=1, validation_data=(Xv,Yv) )
scores = net.evaluate(X_test, Y_test, verbose=1)
print(scores)



import matplotlib.pyplot as plt

from IPython import embed
embed()

plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
# plt.show()

