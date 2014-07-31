# -*- coding: utf-8 -*-
from numpy import random
from numpy.random import normal
from numpy import concatenate,array,zeros,r_,ones,log,sqrt,pi,exp,abs
import sys
import math
import logging
_norm_pdf_C = math.sqrt(2*pi)
def normpdf(x, mu=0.0, sigma=1.0):
    u = (x-mu)/abs(sigma)
    y = exp(-u*u/2.0)/(_norm_pdf_C*abs(sigma))
    return y

def model_str(p):
    return "mu_1:{u1:.3f},sig_1:{s1:.3f},pi_1:{p1:.3f}\tmu_2:{u2:.3f},sig_2:{s2:.3f},pi_2:{p2:.3f}".format(u1=p[0],s1=p[1],p1=p[4],u2=p[2],s2=p[3],p2=1-p[4])
def pdf_model(x, p):
    mu1, sig1, mu2, sig2, pi_1 = p
    return pi_1*normpdf(x, mu1, sig1) + (1.0-pi_1)*normpdf(x, mu2, sig2)
def bayes_p1(x,p): # bayes prob on modle 1 ,  
    mu1, sig1, mu2, sig2, pi_1 = p
    p1=pi_1*normpdf(x,mu1,sig1)
    p2=(1.0-pi_1)*normpdf(x,mu2,sig2)
    return p1/(p1+p2)
def bayes_p2(x,p):
    return 1.0-bayes_p1(x,p)
def log_likelihood_two_1d_gauss(p, sample):
    return -log(pdf_model(sample, p)).sum()
def sim_two_gauss_mix(p, N=1000): 
    s1 =  normal(p[0], p[1], size=N*p[4])
    s2 = normal(p[2], p[3], size=N*(1-p[4]))
    return concatenate([s1,s2])
def fit_two_peaks_EM(sample, sigma=None, weights=False, p0=array([0.2,0.2,0.7,0.2,0.5]), 
        max_iter=300, tollerance=1e-5):
    if not weights: w = ones(sample.size)
    else: w = 1./(sigma**2)
    w *= 1.*w.size/w.sum() # renormalization so they sum to N
    # Initial guess of parameters and initializations
    mu = array([p0[0], p0[2]])
    sig = array([p0[1], p0[3]])
    pi_ = array([p0[4], 1-p0[4]])
    gamma, N_ = zeros((2, sample.size)), zeros(2)
    p_new = array(p0)
    N = sample.size
    # EM loop
    counter = 0
    converged, stop_iteration = False, False
    while not stop_iteration:
        p_old = p_new
        # Compute the responsibility func. and new parameters
        for k in [0,1]:
            gamma[k,:] = w*pi_[k]*normpdf(sample, mu[k], sig[k])/pdf_model(sample, p_new) # SCHEME1
            #gamma[k,:] = pi_[k]*normpdf(sample, mu[k], sig[k])/pdf_model(sample, p_new) # SCHEME2
            N_[k] = gamma[k,:].sum()
            mu[k] = sum(gamma[k]*sample)/N_[k] # SCHEME1
            #mu[k] = sum(w*gamma[k]*sample)/sum(w*gamma[k]) # SCHEME2
            sig[k] = sqrt( sum(gamma[k]*(sample-mu[k])**2)/N_[k] )
            pi_[k] = 1.*N_[k]/N
        p_new = array([mu[0], sig[0], mu[1], sig[1], pi_[0]])
        assert abs(N_.sum() - N)/float(N) < 1e-6 
        assert abs(pi_.sum() - 1) < 1e-6
        # Convergence check
        counter += 1
        max_variation = max((p_new-p_old)/p_old)
        converged = True if max_variation < tollerance else False
        stop_iteration = converged or (counter >= max_iter)
    #print "Iterations:", counter
    if not converged: logging.warning("WARNING: Not converged")
    return p_new



if __name__=="__main__":
    import logging
    from types import *
    logging.basicConfig(level=logging.WARNING)
    p=[0.2,0.3,0.7,0.2,0.3]
    s = sim_two_gauss_mix(N=1000, p=p)
    print "sample:",s
    print "sample:",type(s)
    print "data:",model_str(p)
    model=fit_two_peaks_EM(s)
    print "em  :",model_str(model)
    r=normal(p[0],p[1],size=10)
    print r,"\n",pdf_model(r,model)
    print bayes_p1(r,model)
    print bayes_p2(r,model)

