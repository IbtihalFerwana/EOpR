import pandas as pd
import numpy as np

def get_ubar(AA,u,xx, reg):
    '''a function to calculate the recovered signal u_bar'''
    yy = u[xx]
    I = np.eye(len(AA[0]),len(AA[0]))
    R = (AA.T@AA + I*reg)
    
    Q = np.linalg.pinv(R,hermitian=False)
    PHI = R[:,xx].T@Q@R[:,xx]
    w = np.linalg.inv(PHI)@yy
    ubar = R[:,xx]@w
    return xx,ubar,PHI, R, Q

def error_bars(xx,ubar,S, PHI, R, Q):
    PHI_inv = np.linalg.inv(PHI)
    T_size = len(S.T)
    # calculate cbar
    cbar = np.zeros((len(xx),T_size))
    for ii in range(0, T_size):
        cbar[:,ii] = PHI_inv@R[:,xx].T@Q@R[:,ii] 
    
    # calculate ybar
    eps = 0.001
    ybar = np.zeros((T_size,T_size))
    for ii in range(0,T_size):
        ybar[:,ii] = R[:,np.append(ii,xx)]@np.append(1,-cbar[:,ii])
        if sum(ybar[:,ii]) > eps:
            ybar[:,ii] = ybar[:,ii]/np.sqrt(ybar[:,ii].T@Q@ybar[:,ii])
    
    # error bars
    a = np.abs(np.random.randn(len(S),1))
    scale = np.sqrt(abs(a.T@a - ubar.T@Q@ubar))
    
    uworst1 = np.zeros((T_size,T_size))
    uworst2 = np.zeros((T_size,T_size))
    normvect = np.zeros(T_size)
    for ii in range(0,T_size):
        uworst1[:,ii] = ubar.T + scale*ybar[:,ii]
        uworst2[:,ii] = ubar.T - scale*ybar[:,ii]
        normvect[ii] = np.sqrt(uworst1[:,ii].T@Q@uworst1[:,ii])
    max_error = scale*abs(np.diag(ybar))
    return max_error

def apply_op(noisy_M0, noisy_m0, TrainingEnd):
    # optimize for \lambda
    vals = []
    op_errors = []
    regs = np.arange(0.01,1, 0.01)
    xx = np.array([i for i in range(TrainingEnd)])
    for reg in regs:
        xx,yy,u, ubar,w, Q = get_ubar(noisy_M0,noisy_m0,xx,reg)
        vals.append(np.sqrt(np.mean((ubar[:TrainingEnd] - noisy_m0[:TrainingEnd])**2)))

    mm = np.argmin(vals)
    min_reg = regs[mm]
    print(min_reg)
    xx,yy,u, ubar,w, Q = get_ubar(noisy_M0,noisy_m0,xx,min_reg)
    
    pre_err  = np.sqrt(np.mean((ubar[:TrainingEnd] - noisy_m0[:TrainingEnd])**2))
    post_err  = np.sqrt(np.mean((ubar[TrainingEnd:] - noisy_m0[TrainingEnd:])**2))
    
    return pre_err, post_err

 