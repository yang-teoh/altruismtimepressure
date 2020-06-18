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

################## MODEL FOR EYE TRACKING ADDM ##################

#Collect a matrix of eyetracking data, and create a vector of drifts that are dependent on the aoi given some parameter combination  
def eyetodriftvec(eye, finalsim, nfix,timevec, percept, wself, selfamt, wother, otheramt, wfair, theta):
    #eye: tuple of eyetrackingdata [col 0:subject, 1: fixation number, 2: start,3: stop, 4: duration of trial, 5: AOI]
    #finalsim is length of trial in ms
    #timevec is a vector the length of finalsim
    #percept is minimum time for visual transformation of information
    #theta is the attentional discount
    if not np.equal(eye.size,0):
        eye2 = eye[0]
        for i in range(1,nfix):
            eye2 = np.vstack([eye2, eye[i]]) #stacking the tuple into an array
        prev = nfix -1  #last recorded fixation 
        fillfix = nfix + 1 #last recorded fixation may not extend into the length of time vector (round up to the nearest 100ms)
        if np.equal(nfix,  1):
            nextrow = np.asarray([eye2[0],fillfix, eye2[3] + 1, finalsim, eye2[4], eye2[5]]) #assume that continuation of previous (single row array (1D)) 
        else:
            nextrow = np.asarray([eye2[prev,0],fillfix, eye2[prev,3] + 1, finalsim, eye2[prev,4], eye2[prev,5]]) #fill in information until the nearest 
        
        eye2 = np.vstack([eye2,nextrow]) #stack data with the simulated data
        
        for i in range(fillfix): 
            p = int(eye2[i,2]) #start of the fixation
            if p == 0: #if fixation in AOI starts the moment trial begins
                p = 1 
            q = int(eye2[i,3]) + 1 #get fixation end (+1 for indexing i python)
            timevec[p:q] = eye2[i,5] #fill in information for that particular period
        selfproc = np.where(np.equal(eye2[:,5] , 1)) # find whe participant is looking at self
    
        if np.greater(np.size(selfproc[0]) , 1): #if look at all
            selfproc = selfproc[0][0] #find the first fixation on self
            selfprocblindstart = np.int(eye2[selfproc, 2]) #start of that fixation
            selfprocblindend = np.int(selfprocblindstart + (percept*1000)) #percept ms after end 
            if np.greater(selfprocblindend , finalsim): #check if it overlaps with end of the simulated trials
                selfprocblindend = finalsim + 1
            timevec[selfprocblindstart:selfprocblindend] = np.nan #no processing because perceptual transform
    
        otherproc = np.where(eye2[:,5] == 2) #do same for the first fixation on other
    
    
        if np.greater(np.size(otherproc[0]) , 1):
            otherproc = otherproc[0][0]
            otherprocblindstart = np.int(eye2[otherproc, 2])
            otherprocblindend = np.int(otherprocblindstart + (percept*1000))
            if np.greater(otherprocblindend , finalsim):
                otherprocblindend = finalsim + 1
            timevec[otherprocblindstart:otherprocblindend] = np.nan
    
    
        driftvec = np.zeros(finalsim) #create a vector of drift rates same length at timevec that vary from ms to ms depending on the AOI indicated by time vec
    
        
        p = np.where(np.isfinite(timevec)) #check if there is any fixation at all
        if np.size(p[0]) == 0:  #if never look at any information
            dstart = 0 #dstart = 0
        else: 
            dstart = p[0][1] #first time you look at information
            AOIs  = np.arange(1,3) #AOI - 1 is self, 2 is other
            #different  drifts for differet contingencies
            selfonly = wself*selfamt 
            otheronly = wother*otheramt
            notfocusself = (1-theta)*wself*selfamt
            notfocusother = (1-theta)*wother*otheramt
            notfocusboth = (1-theta)*wself*selfamt + (1-theta)*wother*otheramt \
                           - wfair*(np.abs(selfamt-otheramt))
            selfdom = wself*selfamt + (1-theta)*wother*otheramt - wfair*(np.abs(selfamt-otheramt))
            otherdom = (1-theta)*wself*selfamt + wother*otheramt - wfair*(np.abs(selfamt-otheramt))
            
            driftvec = filldriftvecjit(driftvec,dstart,AOIs,timevec,selfdom,otherdom,notfocusboth,notfocusself,notfocusother,otheronly,selfonly)
       
    else:
        driftvec = np.zeros(finalsim)
        dstart = 0
          
    return(driftvec,dstart)
    
#compiled commands for filling in driftvector given the time vector   
@jit(nb.float64[:](nb.float64[:], nb.int64, nb.int64[:], nb.float64[:], nb.float64, nb.float64, nb.float64, nb.float64, nb.float64, nb.float64, nb.float64), nopython = True,nogil=True,cache =True) 
def filldriftvecjit(driftvec,dstart,AOIs,timevec,selfdom,otherdom,notfocusboth,notfocusself,notfocusother,otheronly,selfonly):
    #driftvec is a vector of zeros the length of trial in milliseconds
    #dstart is the time at which the first piece of information is fixated or 0
    #AOIs is the list of possible AOIs
    #timevec is a vector the length of the trial in milliseconds with information of which AOI is being fixated upon at that specific time
    #selfdom is the drift rate when fixating on $Self after having acquired $Self & $Other
    #otherdom is the drift rate when fixating on $Other after having acquired $Self & $Other
    #notfocusboth is the drift rate when not fixating on $Self or $Other after having acquired both $Self & $Other
    #notfocusself is the drift rate when not fixating on $Self or $Other after only having acquired $Self
    #notfocusother is the drift rate when not fixating on $Self or $Other after only having acquired $Other
    #otheronly is the drift rate when fixating on $Other after having only acquired $Other
    #selfonly  is the drift rate when fixating on $Self after having only acquired $Self
    if np.isfinite(dstart):
        info2 = AOIs[np.where(np.not_equal(AOIs, timevec[dstart]))[0][0]]
        infot2 = np.where(np.equal(timevec, info2))[0]
    if np.equal(infot2.size ,0):
        aoiswitch = np.nan
        first = driftvec
        firsttime = timevec
    else:
        aoiswitch = infot2[0]
        first = driftvec[0:aoiswitch]
        firsttime = timevec[0:aoiswitch]
        bothtime = timevec[aoiswitch:]
        both = driftvec[aoiswitch:]
        for j in np.arange(both.size):
            if np.isfinite(bothtime[j]):
                if np.equal(bothtime[j],1):
                    both[j] = selfdom
                elif np.equal(bothtime[j],2):
                    both[j] = otherdom
            else:
                both[j] = notfocusboth

    if np.equal(timevec[dstart] ,1):
        for j in np.arange(first.size):
            if np.isfinite(firsttime[j]):
                first[j] = selfonly
            else:
                first[j] = notfocusself
    elif np.equal(timevec[dstart] ,2):
        for j in np.arange(first.size):
            if np.isfinite(firsttime[j]):
                first[j] = otheronly
            else:
                first[j] = notfocusother

    if np.greater(dstart , 1):
        first[0:dstart] = 0

    if np.isfinite(aoiswitch):
        driftvec[0:aoiswitch] = first
        driftvec[aoiswitch:] = both
    else:
        driftvec = first 
        
    return driftvec 


#main looping mechanism for simulations 
@jit(nb.types.UniTuple(nb.float64[:],2)(nb.int64, nb.float64[:], nb.int64, nb.float64[:], nb.float64,nb.int64,nb.float64,nb.float64[:]), nopython = True,nogil=True,cache =True)    
def ddmloopjit(nsims, boundaries, dstart, driftvec,stbias, finalsim,precision, noisematrix):
    #nsims = number of simulations
    #boundaries is a vector the length of the simulation in milliseconds filled with the decision threshold values at the given time
    #dstart is the point at which accumulation begins
    #driftvec is the vector of drift rates as long as the lenght of the simulation
    #stbias is the starting point bias towards accepting the proposal on this trial
    # finalsim is the length of the simulations in milliseconds
    #precision is the stepsize of each time-step in seconds (.001 is millisecond)
    #noisematrix is a pregenerated matrix filled with gaussian noise
    ddms = np.zeros(nsims) + stbias
    t = dstart-1
    withinbounds = np.ones(nsims)
    resp = np.full(nsims,np.nan)
    rt = np.full(nsims,np.nan)

    while (np.logical_and(np.less(t,finalsim-1), np.any(withinbounds))):
        t += 1
        stillgoing = np.where(np.equal(withinbounds,1))[0]
        dy = np.full(stillgoing.size,driftvec[t]*precision) + np.random.choice(noisematrix,stillgoing.size,replace=True)
        ddms[stillgoing] += dy
        bcoll = boundaries[t]
        boundcompare = np.where(np.greater_equal(np.abs(ddms), bcoll))[0]
        withinbounds[boundcompare] = 0
        justdone = np.where(np.logical_and(np.equal(withinbounds,0),np.isnan(resp)))[0]
        resp[justdone] = np.sign(ddms[justdone])
        rt[justdone] = t + 1
    return(resp,rt)
   
def maxmax(finish,maxRT, motor):
    #this function determines the maximum amount of time to simulate the choice 
    #usually 100ms after choiceRT - assumed motor response time ('nerve conduction latency')
    #if finish is within 100ms of the time limit, simulate up to time limit - assumed motor response time.

    #finish is the length of the trial (choice RT)
    #maxRT is the time limit on the choice
    #motor is the assumed motor response time : 80ms
    if not(np.isnan(finish)):
        out = finish + .1 - motor
        if np.greater(out, maxRT-motor):
            out = maxRT - motor
    else:
        out = maxRT - motor
    return(out)

#initialise a vector of boundaries given a collapse rate and starting bounnds
@jit(nb.float64[:](nb.float64, nb.float64, nb.int64),nopython = True,nogil=True,cache=True)
def collapseboundjit(bound,collapse,finalsim):
    #this function produces a vector the length of the simulation filled with the value of the boundary threshold at each time point given exponential decay
    #bound is the initial height of the threshold
    #collapse is the decay rate per second
    #finalsim is the length of the simulation

    out = np.full(finalsim, bound)
    for j in np.arange(finalsim):
        out[j] *= np.exp(-1* collapse * .001 *(j+1))
    return out

#parent command that takes into account parameters and eyetracking data  to return simulations 
def simul_giADDM_jit(params,Sims,noisematrix, precision): #Sims : number of simulations
    #s1 = time.perf_counter()
    nsims = int(Sims)

    selfamt =  params['selfamt']
    otheramt =  params['otheramt']
    eye = params['eye']
    respoptions = (0,1)

    wself = params['self']
    wother = params['other']
    wfair = params['fair']
    bound = params['bound']
    collapse = params['collapse']
    stbias = params['modbias']
    percept = params['percept']
    theta = params['theta']

    nfix = len(eye)    
    finalsim = int(params['simmaxrt']*1000)
    timevec = np.full(finalsim, np.nan)

    driftvec,dstart = eyetodriftvec(eye, finalsim, nfix,timevec, percept, wself, selfamt, wother, otheramt, wfair, theta)
    
    stbias *= bound
    boundaries = collapseboundjit(bound,collapse,finalsim)
    
    resp,rt = ddmloopjit(nsims, boundaries, dstart, driftvec,stbias, finalsim,precision,noisematrix)

    
    rt = rt*precision
    resp[np.where(resp == -1)] = respoptions[0]
    resp[np.where(resp == 1)] = respoptions[1]
#    
    outcomes = {'resp' : resp, 'rt' : rt}
    return(outcomes)
    
    

#parent command that takes into account parameters and eyetracking data  to return simulations 
def simul_multiDDM_jit(params,Sims,noisematrix, precision): #Sims : number of simulations
    #s1 = time.perf_counter()
    nsims = int(Sims)

    selfamt =  params['selfamt']
    otheramt =  params['otheramt']
    respoptions = (0,1)

    wself = params['self']
    wother = params['other']
    wfair = params['fair']
    bound = params['bound']
    collapse = params['collapse']
    stbias = params['modbias']
    dstart = int(params['ndt']*1000)

    finalsim = int(params['simmaxrt']*1000)

    drift  = wself*selfamt + wother*otheramt - wfair*np.abs(selfamt-otheramt)
    
    driftvec = np.full(finalsim, drift)
    
    stbias *= bound
    boundaries = collapseboundjit(bound,collapse,finalsim)
    
    resp,rt = ddmloopjit(nsims, boundaries, dstart, driftvec,stbias, finalsim,precision,noisematrix)

    
    rt = rt*precision
    resp[np.where(resp == -1)] = respoptions[0]
    resp[np.where(resp == 1)] = respoptions[1]
#    
    outcomes = {'resp' : resp, 'rt' : rt}
    return(outcomes)
    

#obtain NLL for the given parameter set in fitting the data using the KDE method with (silverman's factor )/2 : recommended 10k sims
def get_giADDM_NLLs(i, splitparams, Sims, timelimit, tempselfamt,tempotheramt, tempfinish, tempfix, tempchoice,noisematrix,precision):
    simParams = splitparams
    modbias = np.zeros(len(tempfinish), dtype = float) #modbias is the composite starting point bias on a trial given both (stbias & genbias)
    modbias[np.where(np.greater(tempselfamt,  tempotheramt))] = simParams[5] - simParams[6] #if acceptance is a selfish choice, a genbias is bias towards rejection (-)
    modbias[np.where(np.less(tempselfamt, tempotheramt))] = simParams[5] + simParams[6] #if acceptance is a generous choice, a genbias is bias towards acceptance (-)
    params = {'self' : simParams[0], 'other' : simParams[1],'fair' : simParams[2],\
              'bound' : simParams[3],'collapse' : simParams[4],\
              'ndt': simParams[7]}
    maxRT = timelimit

    params['selfamt'] = (tempselfamt[i]-50)/5   #scaling self amount
    params['otheramt'] = (tempotheramt[i]-50)/5     #scaling other amount
    params['modbias'] = modbias[i]
    precision = .001
    
    params['simmaxrt'] = np.min(tempfinish[i]+.1, maxRT)  
    
    sims = simul_multiDDM_jit(params,Sims,noisematrix, precision) 
    
    nsims = len(sims['resp'])
    
    decisiontime = tempfinish[i]
    ###CORRECT FOR MOTOR COMMAND TIMING "params['motor']", assuming choice is made well before response is recorded due to this timing.
    
    if tempchoice[i]==0: #IF REJECT
        no = np.asarray(np.where(sims['resp'] == 0))[0]
        nno = no.size
        if nno < 4:  # at least four simulated responses
            prob = (nno/nsims)/(params['simmaxrt']*1000)
        else:
            hsmooth = (nno*(3)/4)**(-1/5)/2
            nort = sims['rt'][no]
            try:
                nodist = sps.gaussian_kde(nort,hsmooth)
                prob = nodist.evaluate(decisiontime)*precision
                if prob > 1:
                    prob = 1
                elif prob < 0:
                    prob = 0
            except: #if only one simulated RT 
                if np.equal(decisiontime,np.unique(nort)): 
                    prob = 1
                else:
                    prob = 0
            prob *= (nno/nsims)
    elif tempchoice[i]==1: #IF ACCEPT
        yes = np.asarray(np.where(sims['resp'] == 1))[0]
        nyes = yes.size
        if nyes < 4:
            prob = (nyes/nsims)/(params['simmaxrt']*1000)
        else:
            hsmooth = (nyes*(3)/4)**(-1/5)/2
            yesrt = sims['rt'][yes]
            try:
                yesdist = sps.gaussian_kde(yesrt,hsmooth)
                prob = yesdist.evaluate(decisiontime)*precision
                if prob > 1:
                        prob = 1
                elif prob < 0:
                        prob = 0
            except:
                if np.equal(decisiontime, np.unique(yesrt)):
                    prob = 1
                else:
                    prob = 0
            prob*=(nyes/nsims)
    else: #IF NO RESPONSE GATHER SIMULATED DATA PAST THE TIME LIMITS TO CONSTRUCT PROBABILITY
        fail = np.where(np.isnan(sims['resp'])) 
        nfail = fail[0].size
        prob = nfail/nsims
        
    out = np.log(prob + sys.float_info.min)     #LOG TRANSFORM
    return(out)



#obtain NLL for the given parameter set in fitting the data using the KDE method with (silverman's factor )/2 : recommended 10k sims
def get_giADDM_NLLs(i, splitparams, Sims, timelimit, tempselfamt,tempotheramt, tempfinish, tempfix, tempchoice,noisematrix,precision):
    simParams = splitparams
    modbias = np.zeros(len(tempfinish), dtype = float) #modbias is the composite starting point bias on a trial given both (stbias & genbias)
    modbias[np.where(np.greater(tempselfamt,  tempotheramt))] = simParams[5] - simParams[6] #if acceptance is a selfish choice, a genbias is bias towards rejection (-)
    modbias[np.where(np.less(tempselfamt, tempotheramt))] = simParams[5] + simParams[6] #if acceptance is a generous choice, a genbias is bias towards acceptance (-)
    params = {'self' : simParams[0], 'other' : simParams[1],'fair' : simParams[2],\
              'bound' : simParams[3],'collapse' : simParams[4],\
              'theta': simParams[7],'percept' : .08, 'motor': .08}
    maxRT = timelimit

    params['selfamt'] = (tempselfamt[i]-50)/5   #scaling self amount
    params['otheramt'] = (tempotheramt[i]-50)/5     #scaling other amount
    params['modbias'] = modbias[i]
    params['eye'] = tempfix[i]
    precision = .001
    
    params['simmaxrt'] = maxmax(tempfinish[i], maxRT, params['motor'])  
    
    sims = simul_giADDM_jit(params,Sims,noisematrix, precision) 
    
    nsims = len(sims['resp'])
    
    decisiontime = tempfinish[i]-params['motor']
    ###CORRECT FOR MOTOR COMMAND TIMING "params['motor']", assuming choice is made well before response is recorded due to this timing.
    
    if tempchoice[i]==0: #IF REJECT
        no = np.asarray(np.where(sims['resp'] == 0))[0]
        nno = no.size
        if nno < 4:  # at least four simulated responses
            prob = (nno/nsims)/(params['simmaxrt']*1000)
        else:
            hsmooth = (nno*(3)/4)**(-1/5)/2
            nort = sims['rt'][no]
            try:
                nodist = sps.gaussian_kde(nort,hsmooth)
                prob = nodist.evaluate(decisiontime)*precision
                if prob > 1:
                    prob = 1
                elif prob < 0:
                    prob = 0
            except: #if only one simulated RT 
                if np.equal(decisiontime,np.unique(nort)): 
                    prob = 1
                else:
                    prob = 0
            prob *= (nno/nsims)
    elif tempchoice[i]==1: #IF ACCEPT
        yes = np.asarray(np.where(sims['resp'] == 1))[0]
        nyes = yes.size
        if nyes < 4:
            prob = (nyes/nsims)/(params['simmaxrt']*1000)
        else:
            hsmooth = (nyes*(3)/4)**(-1/5)/2
            yesrt = sims['rt'][yes]
            try:
                yesdist = sps.gaussian_kde(yesrt,hsmooth)
                prob = yesdist.evaluate(decisiontime)*precision
                if prob > 1:
                        prob = 1
                elif prob < 0:
                        prob = 0
            except:
                if np.equal(decisiontime, np.unique(yesrt)):
                    prob = 1
                else:
                    prob = 0
            prob*=(nyes/nsims)
    else: #IF NO RESPONSE GATHER SIMULATED DATA PAST THE TIME LIMITS TO CONSTRUCT PROBABILITY
        fail = np.where(np.isnan(sims['resp'])) 
        nfail = fail[0].size
        prob = nfail/nsims
        
    out = np.log(prob + sys.float_info.min)     #LOG TRANSFORM
    return(out)


    
#transform parameters from normal to uniform between some bound            
def transform(splitparams, fixedparamsvals, fixedparams, paramnames):
    pars11 = splitparams[0:8] - splitparams[8:]  
    pars21 = splitparams[0:8] + splitparams[8:]
    pars1 = .5*(1 + erf(pars11/np.sqrt(2)))    #probit transformations to get pars1 : the low time pressure condition
    pars2= .5*(1 + erf(pars21/np.sqrt(2)))    #probit transformations to get pars2 : the high time pressure condition
    pars1[:3] *= .2  #first three parameters are weights (wself, wother, wfair) with a possible range of (-.1, .1), this step creates the range of .2
    pars2[:3] *= .2
    pars1[:3] -= .1 #centers this range around 0
    pars2[:3] -= .1
    pars1[3] *= .5 #range for boundary (0, .5)
    pars2[3] *= .5
    pars1[4] *= 2.5 #range for collapse (0, 2.5)
    pars2[4] *= 2.5
    pars1[5] -= .5  #range for stbias (-.5, .5)
    pars2[5] -= .5
    pars1[6] -= .5 #range for genbias (-.5, .5)
    pars2[6] -= .5

    fixpos = np.where(np.equal(fixedparams, 1))[0]

    for i in np.arange(fixpos.size):
        pars1[fixpos[i]] = fixedparamsvals[paramnames[fixpos[i]]]
        pars2[fixpos[i]] = fixedparamsvals[paramnames[fixpos[i]]]

    #pars1[7] & pars1[8] : parameters that range (0, 1) : default range of parameters after transformation
    return(pars1, pars2) #pars1: low time, pars2: high time

#return MCMCparams given giADDM params in each condition (pars11: low time, pars21: high time )
def invtransform(pars11, pars21):
    pars11[:3] += .1
    pars21[:3] += .1
    pars11[:3] /= .2
    pars21[:3] /= .2
    pars11[3] /= .5 #range for boundary (0, .5)
    pars21[3] /= .5
    pars11[4] /= 2.5 #range for collapse (0, 2.5)
    pars21[4] /= 2.5
    pars11[5] += .5  #range for stbias (-.5, .5)
    pars21[5] += .5
    pars11[6] += .5 #range for genbias (-.5, .5)
    pars21[6] += .5

    pars1 = np.sqrt(2)*erfinv(2*pars11-1)     #probit transformations to get pars1 : the low time pressure condition

    pars2 = np.sqrt(2)*erfinv(2*pars21-1)     #probit transformations to get pars2 : the high time pressure condition

    splitparams = np.concatenate(((pars1+pars2)/2,
    								(pars2-pars1)/2))

    return(splitparams)

