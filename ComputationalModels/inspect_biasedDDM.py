#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: teohyiyang
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
import copy
from scipy import io as input
from scipy.special import erf as erf
from scipy.special import erf as erfinv
import numba as nb
from numba import jit
import scipy.stats as sps
import sys
import concurrent.futures
import time
import giADDM
import evolveMCMC
from pathlib import Path

mainfolder = str(Path(__file__).resolve().parent.parent)
                           
modelname = 'biasedDDM'
kfolds = 1
nk = 1
ntrials = 160
ncores = 4
if kfolds == 1:
    kk = 1
    modifier = 'fullfit'
else:
    kk = (kfolds-1)/kfolds
    modifier = str(nk) + '_of_' + str(kfolds) + '_folds'


direct = 'Analysis' + '/migrateMCMC/'+ modelname +'/' + modifier +'/'

paramstofit = {'wself':1,'wother':1,'wfair':1,'bound':1,'collapse':1,'stbias':1, 'genbias':1, 'ndt':1, 
    'dwself':1,'dwother':1,'dwfair':1,'dbound':1,'dcollapse':1,'dstbias':1, 'dgenbias':1, 'dndt':1}

paramnames = ('wself','wother','wfair','bound','collapse','stbias', 'genbias', 'ndt', 
    'dwself','dwother','dwfair','dbound','dcollapse','dstbias', 'dgenbias', 'dndt')

fixedparamsvals = {}

fixedparams = np.zeros((8,8))

for j in np.arange(len(list(fixedparamsvals.keys()))):
    for i in np.arange(8):
        fixedparams[j,i] = paramnames[i]==list(fixedparamsvals.keys())[j]

fixedparams = np.sum(fixedparams,axis = 0)

freeparams = np.asarray(list(paramstofit.values()))

sub = [2,\
             3,\
             4,5,6,8,
             10,11,12,13,14,16,\
             #17,\
             18,\
             19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,\
             #34,
             35,\
             #37,\
             39,40, \
            41,42,43,44,45,46,47,48,49,50,51,53,54,55,56,57,\
            58,\
            59 ,60]

firstpd = 0

ss = time.perf_counter()
dffull = np.zeros([len(sub), 43+24])

subjectfitpars = np.full([len(sub),21],np.nan)

converge1 = np.zeros((len(sub), 5))

cond = {0:'long', 1:'short'}
maxtime = {0: 10, 1: 1.5}
nsimspfit = 10000 #number of simulations for kde

#INITIALISE PRIOR PROBABILITY DISTRIBUTIONS (NORMAL DISTRIBUTION VAR = .5, MEAN = 0)  : FITTING PARAMETERS AS A MEAN +/- (CHANGE) -> PROBIT TRANSFORM -> UNIFORM 
### SO THE RESULTING FITTED PARAMETERS WOULD BE A UNIFORM DISTRIBUTION SINCE A PROBIT TRANSFORM OF A NORMAL(0,1) DISTRIBUTION -> FULLY UNIFORM BETWEEN (0,1)
### FITTING PARAMETERS BETWEEN CONDITIONS AS MEAN (1 PARAMETER) +/- CHANGE (1 PARAMETER) EACH WITH VARIANCE (0.5) RESULTS IN A DISTRIBUTION OF (0,1) FOR EACH PARAMETER 
### TRANSFORMED


#FOR EFFICIENCY, STORE PROBABILITIES OF A NORMAL DISTRIBUTION (0,0.5) AS A VECTOR WITH A PRECISION OF .001, FROM -10 TO 10 
vector1 = np.arange(-10000,10001)  
vector1 = vector1/1000
vector2 = sps.norm(0,np.sqrt(.5)).cdf(vector1)
probabilitydist = vector2[1:] - vector2[:-1]
probabilitydist[10000:] = np.fliplr(probabilitydist[np.newaxis,:10000])

precision = float(.001)
s = float(.1)#ddm noise
noisematrix = np.random.normal(0,1,100000)*s*np.sqrt(precision)


subjectchangepars = np.zeros([len(sub), 8])
subjectavgpars = np.zeros([len(sub), 8])

    
burn = 500
steps= 1000

for sub1 in np.arange(len(sub)):
    subject = sub[sub1]
    subjectfitpars[sub1,0] = sub[sub1]
    dffull[sub1,0] = subject
    substr = str(subject)
    folder = direct + str(subject) 
    converge1[sub1,0] = subject
    chainfile = folder + '/fit_' + str(subject)+'_bothcond'
    with open(chainfile, 'rb') as f:
        out=pickle.load(f)
    
    
    chainmem = out['chainmem']
    chainll = out['chainLL']
    nreset = out['nreset']
    
    chainll = chainll[burn:steps,:]
    chainmem = chainmem[burn:steps,:,:]
    
    wself = chainmem[:,0,:]
    wother = chainmem[:,1,:]
    wfair = chainmem[:,2,:]
    bound = chainmem[:,3,:]
    collapse = chainmem[:,4,:]
    stbias = chainmem[:,5,:]
    genbias = chainmem[:,6,:]
    ndt = chainmem[:,7,:]
    dwself = chainmem[:,8,:]
    dwother = chainmem[:,9,:]
    dwfair = chainmem[:,10,:]
    dbound = chainmem[:,11,:]
    dcollapse = chainmem[:,12,:]
    dstbias = chainmem[:,13,:]
    dgenbias = chainmem[:,14,:]
    dndt = chainmem[:,15,:]
    
    ##Calculating the gelman-rubin r-hat
    Rhat = np.zeros(chainmem.shape[1])
    
    chainlength = steps - burn
    
    x = np.zeros(chainmem.shape[1])
    y = np.zeros(chainmem.shape[1])
    for j in np.arange(chainmem.shape[1]):
        for i in np.arange(chainmem.shape[2]):
            x[j] += (np.nanmean(chainmem[:,j,i]) - np.nanmean(chainmem[:,j,:]))**2
            y[j] += i/(chainlength - 1)*np.sum((chainmem[:,j,i] - np.nanmean(chainmem[:,j,i]))**2)
            
    B = chainlength*x/(chainmem.shape[2]-1)
    W = y/chainmem.shape[2] 
    
    varndtdata = (chainlength-1)*W/chainlength + B/chainlength
    Rhat = np.sqrt(varndtdata/W)
    
    Rhat < 1.1
    
    converge1[sub1,1] = np.sum(Rhat > 1.1)
    converge1[sub1,2] = np.max(Rhat)
    
    converge1[sub1,3] = steps 
#    
    converge1[sub1,4] = nreset
        
    
    temp = input.loadmat(mainfolder + '/SubjectData/' + substr + '/Data.' + substr + '.wfix.choice.mat',\
                         struct_as_record = False)
    x = temp['Data']
     
    finish = np.round(x[0,0].ChoiceRT[0].astype(float),3)
    finish = tuple(finish)
    choice = x[0,0].Resp[0]
    choice[np.where(np.asarray(choice == 'NULL'))] = np.asarray([[np.nan]])
    choice[np.where(np.equal(choice,1))] = 0
    choice[np.where(np.equal(choice,2))] = 1
    choice = tuple(choice.astype(float))
    SelfProposal = tuple(x[0,0].SelfProposal[0].astype(int))
    OtherProposal = tuple(x[0,0].OtherProposal[0].astype(int))
    TimeLimit = tuple(x[0,0].TimeLimit[0].astype(float))
    Fixations = tuple(x[0,0].Fix[0])
    
    
    # plt.figure()
    # plt.plot(np.arange(len(chainll)),chainll, alpha = 1)
    # plt.suptitle('ChainNLL subject #' + substr)
    # plt.savefig(folder + '/chainll.png')
    
    # fig, axs = plt.subplots(4,4,figsize=(10,8))
    # plt.subplots_adjust(wspace = .3,  # the amount of width reserved for space between subplots,
    #           hspace = .5)# expressed as a fraction of the average axis width
    # axs[0,0].plot(np.arange(len(wself)),wself)
    # axs[0,1].plot(np.arange(len(wother)),wother)
    # axs[0,2].plot(np.arange(len(wfair)),wfair)
    # axs[0,3].plot(np.arange(len(bound)),bound)
    # axs[1,0].plot(np.arange(len(collapse)),collapse)
    # axs[1,1].plot(np.arange(len(stbias)),stbias)
    # axs[1,2].plot(np.arange(len(genbias)),genbias)
    # axs[1,3].plot(np.arange(len(ndt)),ndt)
    # axs[2,0].plot(np.arange(len(dwself)),dwself)
    # axs[2,1].plot(np.arange(len(dwother)),dwother)
    # axs[2,2].plot(np.arange(len(dwfair)),dwfair)
    # axs[2,3].plot(np.arange(len(dbound)),dbound)
    # axs[3,0].plot(np.arange(len(dcollapse)),dcollapse)
    # axs[3,1].plot(np.arange(len(dstbias)),dstbias)
    # axs[3,2].plot(np.arange(len(dgenbias)),dgenbias)
    # axs[3,3].plot(np.arange(len(dndt)),dndt)
    # for j in np.arange(16):
    # 	axs[int(np.floor(j/4)), int(np.mod(j,4))].set_title(paramnames[j])
    # fig.suptitle('Chain Traces Subject #' + substr)
    # plt.savefig(folder + '/chaintraces.png')
    
    wself = np.concatenate(wself, axis = 0)
    wother = np.concatenate(wother, axis = 0)
    wfair = np.concatenate(wfair, axis = 0)
    bound = np.concatenate(bound, axis = 0)
    collapse = np.concatenate(collapse, axis = 0)
    stbias = np.concatenate(stbias, axis = 0)
    genbias = np.concatenate(genbias, axis = 0)
    ndt = np.concatenate(ndt, axis = 0)
    dwself = np.concatenate(dwself, axis = 0)
    dwother = np.concatenate(dwother, axis = 0)
    dwfair = np.concatenate(dwfair, axis = 0)
    dbound = np.concatenate(dbound, axis = 0)
    dcollapse = np.concatenate(dcollapse, axis = 0)
    dstbias = np.concatenate(dstbias, axis = 0)
    dgenbias = np.concatenate(dgenbias, axis = 0)
    dndt = np.concatenate(dndt, axis = 0)

    
    # plt.figure(figsize=(10,8))
    # plt.subplot(4,4,1)
    # plt.hist(wself,100)
    # plt.subplot(4,4,2)
    # plt.hist(wother,100)
    # plt.subplot(4,4,3)
    # plt.hist(wfair,100)
    # plt.subplot(4,4,4)
    # plt.hist(bound,100)
    # plt.subplot(4,4,5)
    # plt.hist(collapse,100)
    # plt.subplot(4,4,6)
    # plt.hist(stbias,100)
    # plt.subplot(4,4,7)
    # plt.hist(genbias,100)
    # plt.subplot(4,4,8)
    # plt.hist(ndt,100)
    # plt.subplot(4,4,9)
    # plt.hist(dwself,100)
    # plt.subplot(4,4,10)
    # plt.hist(dwother,100)
    # plt.subplot(4,4,11)
    # plt.hist(dwfair,100)
    # plt.subplot(4,4,12)
    # plt.hist(dbound,100)
    # plt.subplot(4,4,13)
    # plt.hist(dcollapse,100)
    # plt.subplot(4,4,14)
    # plt.hist(dstbias,100)
    # plt.subplot(4,4,15)
    # plt.hist(dgenbias,100)
    # plt.subplot(4,4,16)
    # plt.hist(dndt,100)
    
    mwself = np.nanmean(wself)
    mwother = np.nanmean(wother)
    mwfair = np.nanmean(wfair)
    mbound = np.nanmean(bound)
    mcollapse=np.nanmean(collapse)
    mstbias = np.nanmean(stbias)
    mgenbias = np.nanmean(genbias)
    mndt = np.nanmean(ndt)
    mdwself = np.nanmean(dwself)
    mdwother=np.nanmean(dwother)
    mdwfair = np.nanmean(dwfair)
    mdbound=np.nanmean(dbound)
    mdcollapse=np.nanmean(dcollapse)
    mdstbias=np.nanmean(dstbias)
    mdgenbias=np.nanmean(dgenbias)
    mdndt = np.nanmean(dndt)
    
    
    
    #plt.hist(wself,100)
    sdwself = np.nanstd(wself)
    #plt.hist(wother,100)
    sdwother = np.nanstd(wother)
    #plt.hist(wfair,100)
    sdwfair = np.nanstd(wfair)
    #plt.hist(bound,100)
    sdbound = np.nanstd(bound)
    #plt.hist(collapse,100)
    sdcollapse=np.nanstd(collapse)
    #plt.hist(stbias,100)
    sdstbias = np.nanstd(stbias)
    #plt.hist(genbias,100)
    sdgenbias = np.nanstd(genbias)
    #plt.hist(ndt,100)
    sdndt = np.nanstd(ndt)
    #plt.hist(dwself,100)
    sddwself = np.nanstd(dwself)
    #plt.hist(dwother,100)
    sddwother=np.nanstd(dwother)
    #plt.hist(dwfair,100)
    sddwfair = np.nanstd(dwfair)
    #plt.hist(dbound,100)
    sddbound=np.nanstd(dbound)
    #plt.hist(dcollapse,100)
    sddcollapse=np.nanstd(dcollapse)
    #plt.hist(dstbias,100)
    sddstbias=np.nanstd(dstbias)
    #plt.hist(dgenbias,100)
    sddgenbias=np.nanstd(dgenbias)
    #plt.hist(dndt,100)
    sddndt = np.nanstd(dndt)
    
    
    
    fitpars = [mwself, mwother, mwfair,mbound,mcollapse,mstbias,mgenbias,mndt, mdwself, mdwother, mdwfair, mdbound,mdcollapse,mdstbias,mdgenbias,mdndt]
    fitparssd = [sdwself, sdwother, sdwfair,sdbound,sdcollapse,sdstbias,sdgenbias,sdndt, sddwself, sddwother, sddwfair, sddbound,sddcollapse,sddstbias,sddgenbias,sddndt]
    ci = 1.95*np.asarray(fitparssd)
    
    
    fitpars = np.asarray(fitpars)
    fitparssd = np.asarray(fitparssd)
    
    fittedparams = fitpars
       
    subjectfitpars[sub1,1:17] = copy.deepcopy(fitpars)
    acceptance = out['acceptance']
    acceptance = acceptance[burn:steps,:]
    acceptancerate = np.sum(acceptance)/48/(steps-burn)
    subjectfitpars[sub1,-1] = acceptancerate


    #@jit
    pars1, pars2 = giADDM.transform(fittedparams, fixedparamsvals, fixedparams, paramnames)


    iteration mechanism for each step given some parameter for parallelisation
    def iterate1(x): # x is a vector of iterables
        i = int(x[0]) #trial order
        splitparams = x[1:] #parameters untransformed
        pars1 , pars2 = giADDM.transform(splitparams, fixedparamsvals, fixedparams, paramnames)

        selfamt1 = np.asarray(tempselfamt)[selected1] #selected1 is the defined outside this function in a wrapper where it is the index of trials that is in low time pressure
        otheramt1 = np.asarray(tempotheramt)[selected1]
        finish1 = np.asarray(tempfinish)[selected1]
        fix1 = np.asarray(tempfix)[selected1]
        choice1 = np.asarray(tempchoice)[selected1]
        first = giADDM.get_multiDDM_NLLs(i, pars1, nsimspfit, maxtime[0],selfamt1, otheramt1,finish1,fix1,choice1,noisematrix,precision)

        selfamt2 = np.asarray(tempselfamt)[selected2] #selected2 is the defined outside this function in a wrapper where it is the index of trials that is in high time pressure
        otheramt2 = np.asarray(tempotheramt)[selected2]
        finish2 = np.asarray(tempfinish)[selected2]
        fix2 = np.asarray(tempfix)[selected2]
        choice2 = np.asarray(tempchoice)[selected2]
        out2 = giADDM.get_multiDDM_NLLs(i, pars2, nsimspfit, maxtime[1],selfamt2, otheramt2,finish2,fix2,choice2,noisematrix,precision)
        out = first + out2 #summed likelihood for the two trials
        return(out)

    
    
    cond = {0:'long', 1:'short'}
    maxtime = {0: 10, 1: 1.5}
    nsimspfit = 10000
    
    tempchoice = tuple(np.asarray(choice))
    tempfinish = tuple(np.asarray(finish))
    tempselfamt = tuple(np.asarray(SelfProposal))
    tempotheramt = tuple(np.asarray(OtherProposal))
    tempfix = tuple(np.asarray(Fixations))
    
    tempfinish = np.asarray(tempfinish)
    tempchoice = np.asarray(tempchoice)
    templimit = np.asarray(TimeLimit)

    if kfolds == 1:
        selected1 = np.where(np.equal(TimeLimit,maxtime[0]))[0] #low time pressure trials
        selected2 = np.where(np.equal(TimeLimit,maxtime[1]))[0] #high time pressure tirals
    else:
        selected1 = np.where(np.logical_and(np.equal(TimeLimit,maxtime[0]), np.not_equal(np.mod(np.arange(len(TimeLimit)),kfolds), kfolds-nk)))[0] #low time pressure trials
        selected2 = np.where(np.logical_and(np.equal(TimeLimit,maxtime[1]), np.not_equal(np.mod(np.arange(len(TimeLimit)),kfolds), kfolds-nk)))[0] #high time pressure tirals

    paramset = fittedparams
    
    it =  np.arange(int(ntrials*(kk/2)))
    x = np.hstack([it[:,np.newaxis], np.tile(paramset[np.newaxis,:], [int(ntrials*(kk/2)),1])])
    with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
        future = executor.map(iterate1, x ,chunksize = 1)#int(ncombo/100))
    future = tuple(future)
    LL = sum(future)
    
    Dndtbar = -2*LL[0]
    Dbar = -2*np.nanmean(chainll[:,:])
    DIC = 2*Dbar - Dndtbar
    
    
    
    subjectfitpars[sub1,17] =Dndtbar
    subjectfitpars[sub1,18] =Dbar
    subjectfitpars[sub1,19] =DIC
    
    def iterate2(i):
        simulated = np.zeros(3)
        pars = {0: pars1, 1: pars2}
        selfamt = np.asarray(tempselfamt)[i]
        otheramt = np.asarray(tempotheramt)[i]
        fix = np.asarray(tempfix)[i]
        if (np.asarray(TimeLimit)[i] > 5):
            simParams = pars[0]
        else:
            simParams = pars[1]
        if (selfamt > otheramt):
            modbias = simParams[5] - simParams[6]
        else:
            modbias = simParams[5] + simParams[6]
        params = {'self' : simParams[0], 'other' : simParams[1],'fair' : simParams[2],\
                  'bound' : simParams[3],'collapse' : simParams[4],\
                  'ndt': simParams[7],'percept' : .08, 'motor': .08}
    
        params['selfamt'] = (selfamt-50)/5
        params['otheramt'] = (otheramt-50)/5
        params['modbias'] = modbias
        params['eye'] = fix
        params['simmaxrt'] = np.asarray(TimeLimit)[i] - params['motor']
        sims = giADDM.simulEADDMjit(params,nsimspfit,noisematrix,precision) 
        sims['rt'] += params['motor']
        simulated[0] = np.nanmean(sims['resp'])
        simulated[1] = np.nanmean(sims['rt'][np.where(np.equal(sims['resp'],1))])
        simulated[2] = np.nanmean(sims['rt'][np.where(np.equal(sims['resp'],0))])
        return(simulated)
    it2 = np.arange(ntrials)     
    with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
        simulated = executor.map(iterate2, it2)#int(ncombo/100))
    simulated = tuple(simulated)
    simulated = np.asarray(simulated)
        
        
        
    simulated[np.where(np.equal(simulated[:,0],1)),2] = 0
    simulated[np.where(np.equal(simulated[:,0],0)),1] = 0
    avert =  (simulated[:,0]*simulated[:,1]+ (1-simulated[:,0])*simulated[:,2])
    simulated = np.hstack([simulated,avert[:,np.newaxis]])    
    
    import pandas as pd
    
    simuldata = pd.DataFrame(data = simulated, columns = ['accept', 'yesrt','nort','avert'])
    
    simuldata['r.accept'] = tempchoice
    simuldata['r.rt'] = tempfinish
    simuldata['time'] = templimit
    simuldata['self'] = tempselfamt
    simuldata['other'] = tempotheramt
    simuldata['genoffer'] = (simuldata['self'] < simuldata['other']).astype(int) 
    simuldata['genoffer'][np.where(simuldata['genoffer'] < .5)[0]] = -1
    simuldata['gen'] = simuldata['accept']*simuldata['genoffer']
    simuldata['gen'][np.where(simuldata['gen'] < 0)[0]] +=1
    simuldata['r.gen'] = tempchoice *  2- 1
    simuldata['r.gen'] *= simuldata['genoffer']
    simuldata['r.gen'][np.where(simuldata['r.gen'] < 0)[0]]  = 0
    simuldata['comp.rt'] = np.nan
    simuldata['comp.rt'][np.where(simuldata['r.accept']==1)[0]] = simuldata['yesrt'][np.where(simuldata['r.accept']==1)[0]]
    simuldata['comp.rt'][np.where(simuldata['r.accept']==0)[0]] = simuldata['nort'][np.where(simuldata['r.accept']==0)[0]]


    short1 = simuldata['time'] < 5
    long1 = simuldata['time'] > 5
    simuldata['subjval'] = np.nan
    simuldata['subjval'][short1] = pars2[0] * simuldata['self'][short1] +pars2[1] * simuldata['other'][short1] - pars2[2]*np.abs(simuldata['self'][short1] - simuldata['other'][short1])
    simuldata['subjval'][long1] = pars1[0] * simuldata['self'][long1] +pars1[1] * simuldata['other'][long1] - pars1[2]*np.abs(simuldata['self'][long1] - simuldata['other'][long1])
    
    simuldata['subj'] = subject 

    if firstpd == 0:
        allsimulate = simuldata
        firstpd = 1
    else:
        allsimulate = pd.concat([allsimulate, simuldata])
        
    
    # for e in np.arange(chainlength*chainmem.shape[2]):
        
    #       pars1, pars2 = giADDM.transform( np.asarray([wself[e],wother[e],
    #                                                   wfair[e],bound[e],
    #                                                   collapse[e],stbias[e],
    #                                                   genbias[e],ndt[e],
    #                                                   dwself[e], dwother[e],
    #                                                   dwfair[e],dbound[e],
    #                                                   dcollapse[e],dstbias[e],
    #                                                   dgenbias[e],dndt[e]]), 
    #                                       fixedparamsvals, fixedparams, paramnames)
    #       if e == 0:
    #         params1 = (pars1+pars2)/2
    #         params2 = pars2-pars1
    #       else:
    #         params1 = np.vstack((params1, (pars1+pars2)/2))
    #         params2 = np.vstack((params2, pars2-pars1))
    
    
    
    
    p1, p2 = giADDM.transform(np.asarray(subjectfitpars)[sub1,1:17], fixedparamsvals, fixedparams, paramnames)
    subjectchangepars[sub1,] = p2-p1
    subjectavgpars[sub1,] = (p1+p2)/2
    
    # fig, axs = plt.subplots(8,2,figsize=(8,16))
    # plt.subplots_adjust(wspace = .3,  # the amount of width reserved for space between subplots,
    #           hspace = .5)# expressed as a fraction of the average axis width
    # plt.suptitle('Parameters Subject #' + substr)
    
    # for plot in np.arange(8):
    #     bins=50
    #     axs[int(np.floor(plot)), 0].hist(params1[:, plot], bins,alpha = .5, color = 'darkcyan')
    #     axs[int(np.floor(plot)), 1].hist(params2[:, plot], bins,alpha = .4, color = 'crimson')
    #     axs[int(np.floor(plot)), 0].axvline(x = subjectavgpars[sub1,plot], color = 'darkcyan', alpha = .9, linestyle = ':')
    #     axs[int(np.floor(plot)), 1].axvline(x = subjectchangepars[sub1,plot], color = 'crimson',alpha = .9, linestyle = ':')
    #     axs[int(np.floor(plot)), 0].set_title('mean' + paramnames[plot])
    #     axs[int(np.floor(plot)), 1].set_title('change in' + paramnames[plot])
    
    # plt.savefig(folder + '/posteriors.png')
    

import pandas as pd



summarypath = direct + '/summary'
if not os.path.exists(summarypath):
    os.makedirs(summarypath)


allsimulate.to_csv(direct + '/summary/allsimulate.csv')



subjectfitpars = pd.DataFrame(subjectfitpars)
subjectfitpars.columns = ['subj','wself','wother','wfair','bound','collapse','stbias','genbias','ndt','dwself','dwother','dwfair','dbound','dcollapse','dstbias','dgenbias','dndt','Dndtbar','Dbar','DIC', 'acceptance']


subjectfitpars.to_csv(direct+ '/summary'+'/fittedmcmcparams.csv')

    
dftranspars = pd.DataFrame(np.hstack([subjectfitpars.subj[:,np.newaxis],subjectavgpars, subjectchangepars]))

dftranspars.columns = ['subj', 'wself','wother','wfair','bound','collapse','stbias','genbias','ndt','dwself','dwother','dwfair','dbound','dcollapse','dstbias','dgenbias','dndt']

dftranspars.to_csv(direct+ '/summary'+'/transformedparams.csv')



convergence = pd.DataFrame(converge1, columns = ['subj', 'N_nonconverge','max_rhat','nsteps','nreset'])

convergence.to_csv(direct+ '/summary'+'/convergencestat.csv')
    