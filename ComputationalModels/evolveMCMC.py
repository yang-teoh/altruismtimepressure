#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: teohyiyang
"""
#Packages.

import os
import scipy as sp
import numpy as np
from scipy import io as input
from scipy.special import erf as erf
from scipy.special import erfinv as erfinv
import sys
import scipy.stats as sps
import time
import concurrent.futures
import numba as nb
from numba import jit
import pickle
import copy 

#########START OF MCMC##################

########## COMPILED FUNCTION FOR CREATING PROPOSALS FOR NEXT STEP IN MCMC ##################
@jit(nb.types.Tuple((nb.float64[:],nb.float64[:,:]))(nb.int64, nb.float64[:,:],nb.float64[:],nb.float64[:],nb.float64, nb.float64,nb.float64[:], nb.int64[:]), nopython = True, nogil = True, cache = True)
def createpropsevolving(nchains, param_old, LLold, priorold, MCMCgamma, noise_size,probabilitydist, freeparams):
    #NCHAINS: number of chains in the MCMC
    #param_old: current parameter combination
    #ll_old: the loglikelihood of that parameter combination fitting your data
    #prior_old: the prior probability of that combination
    #MCMCgamma is a constant factor based on number of parameters that determines step size for the MCMC in parameter search
    #noise_size is random factor that introduces stochasticity to the parameter search
    #probabilitydist: the vector of probability distribution for the priors 
    proposalprob = np.zeros(nchains) #vector of prior probabilities for the proposed parameter combinations
    shaped = param_old.shape #shape of the parameters vector
    proposals8 = np.zeros(shaped) #empty matrix of possible parameters
    for nc in np.arange(nchains): #loop through each chain
        proposals8[:,nc] = param_old[:,nc] #get current position
        chains = np.arange(nchains) #find number of chains
        chains = chains[np.where(np.not_equal(chains,nc))] #identify all other chains
        chains = np.random.permutation(chains) #randomly permutate the chains
        d1 = chains[0] #select 2
        d2 = chains[1]
        proposals8[:,nc] += (param_old[:,d1]-param_old[:,d2])*MCMCgamma #the new proposal for chain nc is the current position + the scaled difference between two other random chains
        proposals8[:,nc] *= freeparams
        runif = np.zeros(shaped[0]) #create empty vector for noise
        for i in np.arange(shaped[0]): #numba uses for loops faster than numpy random generators
            runif[i] = np.random.uniform(0,1) - .5
        proposals8[:,nc] +=  (noise_size*runif) #add noise

        paramset = proposals8[:,nc]  
        getnorm = paramset/.001 #find position of paramset in the probability distribution
        if np.any(np.logical_or(np.greater_equal(getnorm,  probabilitydist.size/2),np.less_equal(getnorm, - probabilitydist.size/2))):
            probab = 0 # if not between -10 to 10, assume 0
        else:
            getnorm1 = np.zeros(getnorm.shape)
            probab = 1
            for j in np.arange(getnorm.shape[0]):
                getnorm1[j] = int(getnorm[j] +  probabilitydist.size/2)
                probab *= probabilitydist[int(getnorm1[j])]    
        proposalprob[nc] = probab #prior probability of the parameter
            
    return(proposalprob, proposals8) #proposalprob (vector nchains long with the prior probability of each combination), #proposals8 (matrix nchains*nparams with a parameter combination for each chain)
        


#compiled prior probability retrieval of a combination given a parameter set and the probability distribution            
@jit(nb.float64(nb.float64[:],nb.float64[:],nb.int64[:]), nopython = True, nogil=True,cache =True)
def priorpjit(paramset,probabilitydist,freeparams):
    paramset = paramset[np.where(np.greater(freeparams,0))[0]]
    getnorm = paramset/.001
    if np.any(np.logical_or(np.greater_equal(getnorm,  probabilitydist.size/2),np.less_equal(getnorm, - probabilitydist.size/2))):
        probab = 0
    else:
        getnorm1 = np.zeros(getnorm.shape)
        probab = 1
        for j in np.arange(getnorm.shape[0]):
            getnorm1[j] = int(getnorm[j]+  probabilitydist.size/2)
            probab *= probabilitydist[int(getnorm1[j])]    
    prob = probab
    return(prob)
