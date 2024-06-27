import numpy as np
from scipy import linalg
from scipy.optimize import minimize

def MLE_fit(y, X, D, cov_model):
    eta_ini = np.array([16.03, 0.24])
    lo_bound = np.array([1e-10, 1e-10])
    up_bound = np.array([200,5])
    bound = []
    nug = True
    for i in range(eta_ini.size):
        bound.append((lo_bound[i], up_bound[i]))   
    eta = minimize(innerf, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="Nelder-Mead").x
    cmat = (1-eta[eta.size-1])*np.exp(-D/eta[0])
    np.fill_diagonal(cmat, 1+1e-8)
    white_X = np.linalg.solve(np.linalg.cholesky(cmat), X) 
    white_y = np.linalg.solve(np.linalg.cholesky(cmat), y)
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym') 
    sill = (white_y - white_X @ beta_est).T @ (white_y - white_X @ beta_est)/y.size
    theta_est = np.copy(eta) 
    nug_effect = sill*eta[eta_ini.size-1]
    psill = sill*(1-eta[eta_ini.size-1])
    theta_est = np.append(theta_est[0:(eta_ini.size-1)], [psill, nug_effect])
    return theta_est, beta_est

def innerf(eta, y, X, D, cov_model, nug):
    cmat = (1-eta[eta.size-1])*np.exp(-D/eta[0])
    np.fill_diagonal(cmat, 1+1e-8) 
    white_X = np.linalg.solve(np.linalg.cholesky(cmat), X)
    white_y = np.linalg.solve(np.linalg.cholesky(cmat), y)
    beta = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym')
    return (y.size*np.log(2*np.pi) + y.size + 2*np.sum(np.log(np.diag(np.linalg.cholesky(cmat)))) + y.size*np.log((white_y - white_X @ beta).T @ (white_y - white_X @ beta)/y.size))/y.size


