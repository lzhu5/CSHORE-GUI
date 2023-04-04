import scipy.io
import numpy as np

def trained_model(input, model_parameters ):
    # input = [f1, f2, f3]: should be a list
    # model_parameters: a dictionary contains parameters of ANN model

    input_nparray = np.array([input]).transpose()

#     ## load model data
    netLW1 = model_parameters['netLW1']
    netLW2 = model_parameters['netLW2']
    netLW3 = model_parameters['netLW3']
    netLW4 = model_parameters['netLW4']
    netLW5 = model_parameters['netLW5']
    netb1 = model_parameters['netb1']
    netb2 = model_parameters['netb2']
    netb3 = model_parameters['netb3']
    netb4 = model_parameters['netb4']
    netb5 = model_parameters['netb5']

    lb_x = model_parameters['lb_x']
    ub_x = model_parameters['ub_x']
    lb_y = model_parameters['lb_y']
    ub_y = model_parameters['ub_y']


    # mapminmax
    input_norm = (input_nparray - lb_x)/(ub_x - lb_x)*2.0 - 1.0

    # construct
    # L1
    tmp = np.matmul(netLW1, input_norm) + netb1
    tmpL1 = 2.0 / (1+np.exp(-2*tmp)) - 1.0

    # L2
    tmp = np.matmul(netLW2, tmpL1) + netb2
    tmpL2 = 2.0 / (1+np.exp(-2*tmp)) - 1.0

    # L3
    tmp = np.matmul(netLW3, tmpL2) + netb3
    tmpL3 = 2.0 / (1+np.exp(-2*tmp)) - 1.0

    # L4
    tmp = np.matmul(netLW4, tmpL3) + netb4
    tmpL4 = 2.0 / (1+np.exp(-2*tmp)) - 1.0

    # output layer
    output = np.matmul(netLW5, tmpL4) + netb5

    # reverse mapminmax
    output_reverse_norm = (output + 1.0) / 2.0*(ub_y - lb_y) + lb_y # reverse-normalization

    return output_reverse_norm[0][0]



def load_ANN_model_parameters():
    # the neural network model was trained using MATLAB fitnet.

    model = scipy.io.loadmat('./Auxiliary/net_params.mat')
    netLW1 = model['netLW1']
    netLW2 = model['netLW2']
    netLW3 = model['netLW3']
    netLW4 = model['netLW4']
    netLW5 = model['netLW5']
    netb1 = model['netb1']
    netb2 = model['netb2']
    netb3 = model['netb3']
    netb4 = model['netb4']
    netb5 = model['netb5']

    lb_x = model['lb_x'].transpose()
    ub_x = model['ub_x'].transpose()
    lb_y = model['lb_y'].transpose()
    ub_y = model['ub_y'].transpose()

    model_parameters = {'netLW1': netLW1,
                        "netLW2": netLW2,
                        "netLW3": netLW3,
                        "netLW4": netLW4,
                        "netLW5": netLW5,
                        "netb1": netb1,
                        "netb2": netb2,
                        "netb3": netb3,
                        "netb4": netb4,
                        "netb5": netb5,
                        "lb_x": lb_x,
                        "ub_x": ub_x,
                        "lb_y": lb_y,
                        "ub_y": ub_y}

    return model_parameters
