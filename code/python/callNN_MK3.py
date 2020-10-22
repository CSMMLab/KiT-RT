# imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras.backend as K

plt.style.use("kitish")  # or "kitishnotex" to avoid text rendering with TeX.

# in-project imports
#import sphericalquadpy as sqp




def initialize_network():
    # initalizes the MK3 network

    # Specify the Problem specific parameters
    quadOrder = 10  # QuadOrder for training
    BasisDegree = 2  # Determines the number of moments used

    # Load model
    model = tf.keras.models.load_model('neural_network_model/MK3_nnM_1/model')  # ,custom_objects={'loss': cLoss_FONC_varD(quadOrder, BasisDegree)})

    # Check its architecture
    model.summary()
    return model


# make the network a gobal variable here
model = initialize_network()


def main():
    # 1) Read in test data
    (u_test, alpha_test) = read_data("testData_Solver_M2.csv")
    print("Read test data")
    # 2) Print test data for overview
    plot_data(u_test, alpha_test)

    # 3)  Evaluate Network
    print("Evaluate Network")

    # initialize_network()
    print(call_network(u_test[0, :]))
    print("Single input tested")

    # a) Evaluate model with custom loss
    print(model.evaluate(u_test, u_test, verbose=1))
    print("Custom Loss evaluated")

    # b) Evaluate model with RMSE loss
    modelRMSE = tf.keras.models.load_model('nn_MK3_M2/nn_model', custom_objects={'loss': tf.keras.losses.mse})
    print(model.evaluate(u_test, alpha_test, verbose=1))
    print("RMSE Loss evaluated")

    print("Network testing finished!")
    return 0


def call_network(input):
    # Input: input.shape = (,nMaxMoment), nMaxMoment = 9 in case of MK3
    #         ==> check this!

    inputNP = np.asarray([input])
    predictions = model.predict(inputNP)
    return predictions[0]


def call_networkBatchwise(input):
    # Input: input.shape = (batchSize,nMaxMoment), nMaxMoment = 9 in case of MK3
    #         ==> check this!

    #print(input)
    #print("printed raw input")
    inputNP = np.asarray(input)
    #print(inputNP.shape)
    #print("printed input shape")
    #print(inputNP)
    #print("printed nparray input")

    predictions = model.predict(inputNP)

    # print(predictions)

    size = predictions.shape[0] * predictions.shape[1]
    test = np.zeros(size)
    for i in range(0, size):
        test[i] = predictions.flatten(order='C')[i]
    #print("print result")
    #print(test)
    #print(test.shape)
    return test

    # return predictions.flatten(order='C')[0]


def read_data(filename):
    # Note! This is essentially the same code as preprocess_data(filename) of trainNN_MK1.py. Needs to be refactored!

    ### ------- Code Starts Here ----- ####

    # reading csv file
    dataFrameInput = pd.read_csv(filename)  # outputs a dataframe object

    idx = 0
    xDataList = list()
    yDataList = list()
    data = dataFrameInput.values  # numpy array

    ## Somehow t = 0 has an odd number of elements, just cut it out

    for row in data:
        if (row[2] > 0):
            if (idx % 2 == 0):
                xDataList.append(row)
            else:
                yDataList.append(row)

            idx = idx + 1

    # merge the lists
    DataList = list()
    for rowX, rowY in zip(xDataList, yDataList):
        DataList.append([rowX, rowY])

    # Shuffle data
    # random.shuffle(DataList)

    DataArray = np.asarray(DataList)

    # Strip off header information, i.e. the first 3 cols
    DataArraySlim = DataArray[:, :, 3:]

    xData = DataArraySlim[:, 0, :]
    yData = DataArraySlim[:, 1, :]

    return (xData, yData)


def plot_data(xDataTrain, yDataTrain):
    fig, axs = plt.subplots(4)
    axs[0].plot(xDataTrain[:, 0])
    axs[0].set_title('$m_0^0$')
    axs[0].set_ylabel('value')
    axs[0].set_xlabel('sample id')

    axs[1].plot(xDataTrain[:, 1])
    axs[1].set_title('$m_1^-1$')
    axs[1].set_ylabel('value')
    axs[1].set_xlabel('sample id')

    axs[2].plot(xDataTrain[:, 2])
    axs[2].set_title('$m_1^0$')
    axs[2].set_ylabel('value')
    axs[2].set_xlabel('sample id')

    axs[3].plot(xDataTrain[:, 3])
    axs[3].set_title('$m_1^1$')
    axs[3].set_ylabel('value')
    axs[3].set_xlabel('sample id')

    plt.show()

    return 0

# Custom Loss
'''
def cLoss_FONC_varD(qOrder, bDegree):
    # 1) Compute basis at quadrature points of GaussLegendreQuadrature
    def computeBasis(quadorder, basisDegree):
        Q = sqp.gausslegendre.GaussLegendre(order=quadorder)

        # a) compute quadrature points.
        quadPts = Q.computequadpoints(quadorder)
        quadWeights = Q.computequadweights(quadorder)
        # b) Compute spherical harmonics basis at quadPts
        mBasisList = list()
        for l in range(0, basisDegree + 1):
            for k in range(-l, l + 1):
                temp = sqp.tools.sphericalharmonics.ylm(k, l, quadPts[:, 0], quadPts[:, 1], quadPts[:, 2])
                if (k < 0):
                    mBasisList.append(temp.imag)
                elif (k == 0):
                    mBasisList.append(temp.real)
                else:  # (k > 0):
                    mBasisList.append(temp.real)

        temp = np.array(mBasisList)

        ### Manual computation of spherical harmonics up to order 2###
        # the first nine harmonics:
        def Y0_0(mu, phi):
            return np.sqrt(1 / (4 * np.pi))

        def Y1_m1(mu, phi):
            return -np.sqrt(3 / (4 * np.pi)) * np.sqrt(1 - mu * mu) * np.sin(phi)

        def Y1_0(mu, phi):
            return np.sqrt(3 / (4 * np.pi)) * mu

        def Y1_1(mu, phi):
            return -np.sqrt(3 / (4 * np.pi)) * np.sqrt(1 - mu * mu) * np.cos(phi)

        def Y2_m2(mu, phi):
            return np.sqrt(15 / (16 * np.pi)) * (1 - mu * mu) * np.sin(2 * phi)

        def Y2_m1(mu, phi):
            return -1 * np.sqrt(15 / (4 * np.pi)) * mu * np.sqrt(1 - mu * mu) * np.sin(phi)

        def Y2_0(mu, phi):
            return np.sqrt(5 / (16 * np.pi)) * (3 * mu * mu - 1)

        def Y2_1(mu, phi):
            return -1 * np.sqrt(15 / (4 * np.pi)) * mu * np.sqrt(1 - mu * mu) * np.cos(phi)

        def Y2_2(mu, phi):
            return np.sqrt(15 / (16 * np.pi)) * (1 - mu * mu) * np.cos(2 * phi)

        # Transform karth coordinates to shperical coordinates:
        thetaMu = sqp.tools.xyz2thetaphi(quadPts)  # mu in [0,pi]
        phi = thetaMu[:, 0]
        mu = np.cos(thetaMu[:, 1])
        nPts = quadPts.shape[0]
        basisManual = np.zeros(temp.shape)

        for i in range(0, nPts):  # manual computation...

            basisManual[0, i] = Y0_0(mu[i], phi[i])
            basisManual[1, i] = Y1_m1(mu[i], phi[i])
            basisManual[2, i] = Y1_0(mu[i], phi[i])
            basisManual[3, i] = Y1_1(mu[i], phi[i])
            basisManual[4, i] = Y2_m2(mu[i], phi[i])
            basisManual[5, i] = Y2_m1(mu[i], phi[i])
            basisManual[6, i] = Y2_0(mu[i], phi[i])
            basisManual[7, i] = Y2_1(mu[i], phi[i])
            basisManual[8, i] = Y2_2(mu[i], phi[i])
        #### End Manual compuation ####

        return (basisManual, quadPts, quadWeights)

    # 2) Helper functions
    def innerProd(alphaX, mBasisY):
        # alphaX.shape = (batchSize,NmaxMoment)
        # mBasisY.shape = (NmaxMoment,1)
        # output.shape = (batchSize,1)
        return K.dot(alphaX, mBasisY)

    def entropyDualPrimeMaxwellBoltzmann(alpha, mBasisPt):
        # alpha.shape = (batchSize,NmaxMoment)
        # mBasisPt.shape = (NmaxMoment,1)
        # output.shape = (batchSize,NmaxMoment)

        mBasisPtTensor = tf.constant(mBasisPt, dtype=tf.float32, shape=(mBasisPt.shape[0], 1),
                                     name='Const')  # shape = (N,1)

        temp = K.exp(innerProd(alpha, mBasisPtTensor))
        return K.dot(temp, K.transpose(mBasisPtTensor))

    def integrate(mBasis, alpha, quadWeight):
        # input: mBasis.shape = (NmaxMoment, numQuadPts)
        #         alpha.shape = (batchSize, NmaxMoment)
        #         quadWeight.shape = (numQuadPts)
        # output: integral.shape = (batchSize, NmaxMoment)

        (NmaxMoment, numPts) = mBasis.shape
        (batchSize, NmaxMoment_redundant) = alpha.shape

        integral = tf.constant(np.zeros((batchSize, NmaxMoment)), dtype=tf.float32, name='Const')

        for i in range(0, numPts):
            integral += quadWeight[i] * entropyDualPrimeMaxwellBoltzmann(alpha, mBasis[:, i])

        return integral

    (mbasis, quadPts, quadWeights) = computeBasis(qOrder, bDegree)

    # 3) Compute the Loss
    def loss(u_input, alpha_pred):
        # input: u_input.shape = (batchSize, NmaxMoment)
        #         alpha_pred.shape = (batchSize, NmaxMoment)
        # output: out.shape = (batchSize, 1)

        temp = integrate(mbasis, alpha_pred, quadWeights)
        return K.sum(K.square(temp - u_input), axis=1)

    return loss
'''

if __name__ == '__main__':
    main()
