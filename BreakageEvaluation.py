import numpy as np


def muconv(m, v):
    # need to convert mean and std. see explanation here:
    # https://www.mathworks.com/matlabcentral/answers/270716-generate-a-random-variable-that-is-log-normal-distributed-and-has-a-correlation-coefficent
    # v is variance, v = std^2
    return np.log(m/np.sqrt(1.0+v/m**2.0))


def sigmacov(m, v):
    # need to convert mean and std. see explanation here:
    # https://www.mathworks.com/matlabcentral/answers/270716-generate-a-random-variable-that-is-log-normal-distributed-and-has-a-correlation-coefficent
    # v is variance, v = std^2
    return np.sqrt(np.log(1.0+v/m**2.0))



def MvLogNRand(Mu , Sigma , Simulations , CorrMat):

    # MVLOGNRAND MultiVariant Lognormal random numbers with correlation
    #
    #   Mu: The Lognormal parameter Mu  (can be column or row vector)
    #
    #   Sigma: The Lognormal parameter Sigma (can be column or row vector)
    #
    #   Simulations:  The Number of simulations to run (scalar)
    #
    #   CorrMat:  OPTIONAL A square matrix with the number of rows and columns
    #   equal to the number of elements in Mu/Sigma.  Each element on the
    #   diagonal is equal to one, with the off diagonal cells equal to the
    #   correlation of the marginal Lognormal distributions. If not specified,
    #   then assume zero correlation.
    #
    #   To check the simulation run corrcoef(Y) and that should be the same as
    #   your CorrMat.
    #
    #   For more information see: Aggregration of Correlated Risk Portfolios:
    #   Models and Algorithms; Shaun S. Wang, Phd.  Casualty Actuarial Society
    #   Proceedings Volume LXXXV www.casact.org
    #
    #   Author: Stephen Lienhard

    # Calculate the covariance structure
    sigma_down = np.array([Sigma,]*len(Sigma))
    sigma_acrs = np.array([Sigma,]*len(Sigma)).transpose()

    covv = np.log( CorrMat * np.sqrt(np.exp(sigma_down**2.0)-1.0) * np.sqrt(np.exp(sigma_acrs**2.0)-1.0) + 1.0 );

    # The Simulation
    y = np.exp( np.random.multivariate_normal( Mu , covv , Simulations ))
    random_hv = y[:, 0]
    random_bv = y[:, 1]
    random_flex = y[:, 2]

    return random_hv, random_bv, random_flex
