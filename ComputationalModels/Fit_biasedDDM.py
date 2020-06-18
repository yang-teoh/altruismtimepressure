#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: teohyiyang
"""
#Packages.
import os
import numpy as np
from scipy import io as input
import sys
import scipy.stats as sps
import time
import concurrent.futures
import pickle
import copy 
import giADDM
import evolveMCMC

#############INITIALIZATION#############

mainfolder = str(Path(__file__).resolve().parent.parent)

###SUBJECT INFORMATION###
subject = int(sys.argv[1])

##number of processing cores
ncores =80 
##Current split-half, fit on odd trials
kfolds  = 1 #KFOLD? 
nk = 1 #which of the kfolds
ntrials = 160 
if kfolds == 1:
    kk = 1
else:
    kk = (kfolds-1)/kfolds

##working directory
path1 = os.path.dirname(os.path.realpath(__file__))
path1 = os.getcwd()

##Data Location
path = '/SubjectData/'
substr = str(subject)

#burstbuffer directory on supercomputer allowing rapid i/o (save each step)
# burstpath = '/bb/c/chutcher/teohyi2/giADDM/' 

modelname = 'giADDM'

if kfolds == 1:
    modifier = 'fullfit'
else:
    modifier = str(nk) + '_of_' + str(kfolds) + '_folds'


#create directory for analysis if does not exist
direct2 = burstpath + '/Analysis/' + 'migrateMCMC/' + modelname + '/' + modifier + '/' + substr
direct = path1 + '/Analysis/'  + 'migrateMCMC/'+ modelname +  '/' + modifier + '/' + substr

if not os.path.exists(direct):
    os.makedirs(direct)
if not os.path.exists(direct2):
    os.makedirs(direct2)


#loading savefile is exists
savefilenamescratch = direct +  '/fit_' + substr + '_' + 'bothcond'
savefilenameburst = direct2 +  '/fit_' + substr + '_' + 'bothcond'
if os.path.isfile(savefilenameburst):
    savefilename = savefilenameburst
else:
    savefilename  = savefilenamescratch


### giADDM model parameters of interest ###

###Control MCMC params: 
### toggle both average and change valeus to 0 if you are fixing the values across both conditions
# paramstofit = {'wself':1,'wother':1,'wfair':1,'bound':1,'collapse':1,'stbias':1, 'genbias':0, 'ndt':1, 
#     'dwself':1,'dwother':1,'dwfair':1,'dbound':1,'dcollapse':1,'dstbias':1, 'dgenbias':0, 'dndt':1}


###only toggle change params here to determine if only estimating 1 parameter value across both conditions
### e.g. if 'ndt': 1 and 'dndt':0, model will estimate one ndt parameter for high and low time pressure.
# paramstofit = {'wself':1,'wother':1,'wfair':1,'bound':1,'collapse':1,'stbias':1, 'genbias':1, 'ndt':1, 
#     'dwself':1,'dwother':1,'dwfair':1,'dbound':1,'dcollapse':1,'dstbias':1, 'dgenbias':1, 'dndt':0}

paramstofit = {'wself':1,'wother':1,'wfair':1,'bound':1,'collapse':1,'stbias':1, 'genbias':1, 'ndt':1, 
    'dwself':1,'dwother':1,'dwfair':1,'dbound':1,'dcollapse':1,'dstbias':1, 'dgenbias':1, 'dndt':1}

paramnames = ('wself','wother','wfair','bound','collapse','stbias', 'genbias', 'ndt', 
    'dwself','dwother','dwfair','dbound','dcollapse','dstbias', 'dgenbias', 'dndt')

##Fix value of model params across both conditions: wself, wother, wfair, bound, collapse, stbias, genbias, ndt, 
###if you want to fix specific parameters to be one value across both conditions
#fixedparamsvals = {'genbias': 0} #no generosity bias at all.
fixedparamsvals = {}

fixedparams = np.zeros((8,8))

for j in np.arange(len(list(fixedparamsvals.keys()))):
    for i in np.arange(8):
        fixedparams[j,i] = paramnames[i]==list(fixedparamsvals.keys())[j]
fixedparams = np.sum(fixedparams,axis = 0)

freeparams = np.asarray(list(paramstofit.values()))

##Maximum number of free parameters 
nparams = np.sum(freeparams) 
##Number of free parameters to fit
nfreeparams = np.sum(freeparams) 

##DRIFT DIFFUSION SETTINGS
precision = float(.001)
s = float(.1)#ddm noise
noisematrix = np.random.normal(0,1,100000)*s*np.sqrt(precision)
migrateprob = .1 #probability of a migration step


###deMCMC Settings
MCMCgamma = 2.38/np.sqrt(2*nparams) #step size parameter formula in MCMC sampling.
niter = 2500 #number of max iterations (size of savefile)
nburn = 500 #burn-in
nchains = 3*nfreeparams # 3x number of free params to fit
noise_size = .001 #mcmcnoise
clen = nburn + niter #data structure size (up to #clen#  samples can add more manually)
samplelength = nburn + 500 #chain length 


############ FITTING GI-ADDM MODEL WITH DE-MCMC#############
###PARALLELISATION IN THIS CODE HAPPENS AT THE LEVEL OF SIMULATING TRIALS AND EVALUATING THE KERNEL DENSITY ESTIMATION. 80 CORES, 160 (80 pairs) TRIALS. MOST "OPTIMAL WAY"
###DEPENDING ON NUMBER OF CORES AVAILABLE, OPTIMISATING MAY HAPPEN AT THE LEVEL OF CHAINS OR SOMEWHERE ELSE.


#load data from .mat file
temp = input.loadmat(mainfolder + path + substr + '/Data.' + substr + '.wfix.choice.mat',\
                     struct_as_record = False)
x = temp['Data']

##conversion of data into tuples
finish = np.round(x[0,0].ChoiceRT[0].astype(float),3) #reaction times
finish = tuple(finish)
choice = x[0,0].Resp[0] #choice responses
choice[np.where(np.asarray(choice == 'NULL'))] = np.asarray([[np.nan]])
choice[np.where(np.equal(choice,1))] = 0
choice[np.where(np.equal(choice,2))] = 1
choice = tuple(choice.astype(float))
SelfProposal = tuple(x[0,0].SelfProposal[0].astype(int))  #self amount
OtherProposal = tuple(x[0,0].OtherProposal[0].astype(int)) #other amount
TimeLimit = tuple(x[0,0].TimeLimit[0].astype(float)) #time pressure (response limit)
Fixations = tuple(x[0,0].Fix[0]) #eye gaze data

cond = {0:'long', 1:'short'} #TIME PRESSURE CONDITIONS
maxtime = {0: 10, 1: 1.5} #TIME LIMIT
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

###Identifying training data based on kfold split
if kfolds == 1:
    selected1 = np.where(np.equal(TimeLimit,maxtime[0]))[0] #low time pressure trials
    selected2 = np.where(np.equal(TimeLimit,maxtime[1]))[0] #high time pressure tirals
else:
    selected1 = np.where(np.logical_and(np.equal(TimeLimit,maxtime[0]), np.not_equal(np.mod(np.arange(len(TimeLimit)),kfolds), kfolds-nk)))[0] #low time pressure trials
    selected2 = np.where(np.logical_and(np.equal(TimeLimit,maxtime[1]), np.not_equal(np.mod(np.arange(len(TimeLimit)),kfolds), kfolds-nk)))[0] #high time pressure tirals

tempchoice = tuple(np.asarray(choice))
tempfinish = tuple(np.asarray(finish))
tempselfamt = tuple(np.asarray(SelfProposal))
tempotheramt = tuple(np.asarray(OtherProposal))
tempfix = tuple(np.asarray(Fixations))


#iteration mechanism for each step given some parameter for parallelisation: gets the summed loglikelihood of the xth high time pressure trial & xth low time pressure trial
def iterate1(x): # x is a vector of iterables
    i = int(x[0]) #trial order
    splitparams = x[1:] #parameters untransformed
    pars1 , pars2 = giADDM.transform(splitparams, fixedparamsvals, fixedparams, paramnames)

    selfamt1 = np.asarray(tempselfamt)[selected1] #selected1 is the defined outside this function in a wrapper where it is the index of trials that is in low time pressure
    otheramt1 = np.asarray(tempotheramt)[selected1]
    finish1 = np.asarray(tempfinish)[selected1]
    fix1 = np.asarray(tempfix)[selected1]
    choice1 = np.asarray(tempchoice)[selected1]
    first = giADDM.getNLLs(i, pars1, nsimspfit, maxtime[0],selfamt1, otheramt1,finish1,fix1,choice1,noisematrix,precision)

    selfamt2 = np.asarray(tempselfamt)[selected2] #selected2 is the defined outside this function in a wrapper where it is the index of trials that is in high time pressure
    otheramt2 = np.asarray(tempotheramt)[selected2]
    finish2 = np.asarray(tempfinish)[selected2]
    fix2 = np.asarray(tempfix)[selected2]
    choice2 = np.asarray(tempchoice)[selected2]
    out2 = giADDM.getNLLs(i, pars2, nsimspfit, maxtime[1],selfamt2, otheramt2,finish2,fix2,choice2,noisematrix,precision)

    out = first + out2 #summed likelihood for the two trials
    return(out)


#allow to resume, check if savefile of analysis already exists
if os.path.isfile(savefilename):
    with open(savefilename, 'rb') as f:
        saved = pickle.load(f)
    laststep = saved['laststep'] #chain position 
    chainLL = saved['chainLL'] #history of likelihoods
    chainmem = saved['chainmem'] #history of parameter combinations
    acceptance = saved['acceptance'] #acceptance rate of proposals after burn-in
    priorold = saved['priorold'] #prior of the current combinations
    rejection = saved['rejection'] #how many sequential rejections per chain
    migratevevolve = saved['migratevevolve']
    nmigrations = saved['nmigrations']
    migratingchains = saved['migratingchains']
    nreset = saved['nreset'] #number of resets after burn-in if already passed
else:
    ###initalize MCMC chains###
    migratingchains = [np.nan]*clen
    nmigrations = np.full(clen, np.nan)
    migratevevolve = np.full(clen, np.nan)
    acceptance = np.full([clen,nchains], np.nan)
    chainmem = np.full([clen, nparams, nchains], np.nan) #initialise matrix to record parameters
    chainLL = np.zeros([clen,nchains]) #matrix to record likelihoods
    chainoldprior = np.zeros(nchains) #matrix to record current prior probabilities
    chainoldparams= np.zeros([nparams,nchains]) #matrix to record current parameter values
    iterative = np.zeros(int(1+nparams))
    iterative = iterative[np.newaxis,:] 
    migratevevolve[0] = np.nan
    acceptance[0,:] = 1 
    nreset = 0

    ###INITIALISING THE POSITION OF CHAINS### 
    for nc in np.arange(nchains):
        ###STARTING POINT OF CHAINS WITH RANDOM NOISE####
        chainoldparams[0,nc] = .1 + np.random.uniform(0,1)*.2
        chainoldparams[1,nc] = .1 + np.random.uniform(0,1)*.2
        chainoldparams[2,nc] = .1 + np.random.uniform(0,1)*.1
        chainoldparams[3,nc] = -.5 + np.random.uniform(-1,1)*.1
        chainoldparams[4,nc] = -.5 + np.random.uniform(-1,1)*.1
        chainoldparams[5,nc] = np.random.uniform(-1,1)*.15
        chainoldparams[6,nc] = np.random.uniform(-1,1)*.15
        chainoldparams[7,nc] = -.7 + np.random.uniform(-1,1)*.1

        ##THESE SUBSEQUENT PARAMETERS GOVERN CHANGES ACROSS THE CONDITION (CENTERED AROUND 0: NULL HYPOTHESIS NO CHANGE)
        chainoldparams[8,nc] = np.random.uniform(-1,1)*.025
        chainoldparams[9,nc] = np.random.uniform(-1,1)*.025
        chainoldparams[10,nc] = np.random.uniform(-1,1)*.025
        chainoldparams[11,nc] =  np.random.uniform(-1,1)*.0125
        chainoldparams[12,nc] =  np.random.uniform(-1,1)*.05
        chainoldparams[13,nc] =  np.random.uniform(-1,1)*.025
        chainoldparams[14,nc] =  np.random.uniform(-1,1)*.025
        chainoldparams[15,nc] =  np.random.uniform(-1,1)*.025
        
        paramset = chainoldparams[:,nc] #CURRENT PARAMSET
        paramset *= freeparams

        it = np.arange(int(ntrials*(kk/2))) #HOW MANY ITERATIONS

        x = np.hstack([it[:,np.newaxis], np.tile(paramset[np.newaxis,:], [int(ntrials*(kk/2)),1])]) #CREATE ITERABLE

        iterative = np.vstack([iterative, x])

        chainoldprior[nc] = evolveMCMC.priorpjit(paramset,probabilitydist,freeparams) #GET PRIOR PROBABILITY


    iterative = iterative[1:,:]
    ss = time.perf_counter()
    #PARALLELISE, EACH ITERATION COMPUTES THE LIKELIHOOD OF ONE HIGH TIME PRESSURE TRIAL & ONE LOW TIMEPRESSURE TRIAL 
    with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
        future = executor.map(iterate1, iterative,  chunksize = int(np.ceil((nchains*kk)))) #int(ncombo/100)) 
    future = tuple(future) #ALL THE LIKELIHOODS COMPUTED FOR THE GIVEN DATA FOR THE CURRENT PARAMETER COMBINATION
    ee = time.perf_counter()
    ee-ss
    print(ee-ss)

    #compute summed liklihood for each parameter combination
    for nc in np.arange(nchains):
        chainLL[0,nc] = sum(future[(nc*int(ntrials*(kk/2))):((nc+1)*int(ntrials*(kk/2)))]) #STORE IN MEMORY
        
    chainmem[0,:,:] = chainoldparams #STORE MEMORY
    
    rejection = np.zeros(nchains) #INITIALISE matrix to keep track of how many rejections 

    priorold = copy.deepcopy(chainoldprior) #(prob of params given prior of previous step for each of chains)
    laststep = 0
    nreset = 0

LLold = copy.deepcopy(chainLL[laststep]) #(loglikelihood of previous step for each of chains)
param_old = copy.deepcopy(chainmem[laststep]) 
   
 #stepping through the MCMC
for tstep in np.arange(laststep+1,samplelength):
    ss = time.perf_counter()
    ####determine id it is a migrating step.
    if np.random.uniform(0,1) > migrateprob: #not migrating
        migratevevolve[tstep] = 0
        #create proposals
        proposalprob, proposals8 = evolveMCMC.createpropsevolving(nchains, param_old, LLold, priorold, MCMCgamma, noise_size,probabilitydist,freeparams) #create proposals and get their priors
        #initialize acceptance as 0
        acceptance[tstep,:] = 0    
        #assume rejection
        rejection += 1 #add count 
        #identify possible proposals
        notimpossible = np.where(np.greater(proposalprob, 0))[0]
        if notimpossible.size > 0:
            iterative1 = np.zeros(int(1+nparams))
            iterative1 = iterative1[np.newaxis,:]   
            for nc in notimpossible: #if prior is not 0, evaluate the LL
                paramset = proposals8[:,nc] 
                it = np.arange(int(ntrials*(kk/2)))
                x = np.hstack([it[:,np.newaxis], np.tile(paramset[np.newaxis,:], [int(ntrials*(kk/2)),1])])


                iterative1 = np.vstack([iterative1, x])

            iterative1 = iterative1[1:,:]

            with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
                future = executor.map(iterate1, iterative1, chunksize = int(np.ceil(notimpossible.size*kk)))#int(ncombo/100))
            future = tuple(future)

            for n1 in np.arange(notimpossible.size):
                LL = sum(future[(n1*int(ntrials*(kk/2))):((n1+1)*int(ntrials*(kk/2)))])
                nc = notimpossible[n1]
                log_accept = LL- LLold[nc] + np.log(proposalprob[nc]) - np.log(priorold[nc]) #the acceptance probability given by bayes' rule
                acceptno = np.log(np.random.uniform(0,1)) 
                if np.less(acceptno, log_accept): #if accept
                    acceptance[tstep,nc] = 1 #add to acceptance
                    param_old[:,nc] = proposals8[:,nc]  #replace current parameter combination
                    priorold[nc] = proposalprob[nc] #update the prior
                    LLold[nc] = LL #update likelihood
                    rejection[nc] = 0 #reset rejection counter
    else:  #migrating step
        migratevevolve[tstep] = 1
        chains = np.arange(nchains)
        #randomly select number of chains to migrate
        nmigrate = np.random.permutation(chains)[0]+1
        #randomly identify which chains to migrate
        selectedchains = np.random.permutation(chains)[:nmigrate]


        nmigrations[tstep] = nmigrate #record number of migrations
        migratingchains[tstep] = tuple(selectedchains) #record chaings to migrate
        acceptance[tstep,:] = 0
        rejection += 1 #add count 
        
        #copying existing state of the chains
        proposals8 = np.zeros(param_old.shape)
        proposalprob = np.zeros(priorold.shape)
        proposals8 += param_old
        proposalprob += priorold
        llvec = np.zeros(LLold.shape)
        llvec += LLold
        
        #migrating chains
        for n2 in np.arange(nmigrate):
            #inital position of chain
            nc = selectedchains[n2]
            n3 = n2 + 1
            if n3 > nmigrate-1:
                n3 = 0
            #target chain
            tc = selectedchains[n3]
            LL = LLold[tc]
            
            #proability of migration
            log_accept = LL - LLold[nc] + np.log(priorold[tc]) - np.log(priorold[nc])
            acceptno = np.log(np.random.uniform(0,1))
            
            # if migration is accepted
            if np.less(acceptno, log_accept):
                acceptance[tstep,nc] = 1
                proposals8[:,nc] = param_old[:,tc]
                proposalprob[nc] = priorold[tc]
                llvec[nc] = LLold[tc]
                rejection[nc] = 0
                
        param_old = proposals8
        priorold = proposalprob
        LLold = llvec
   
    
    #to prevent sticking of chains from a poor estimation of the kde of the true likelihoods every 3 rejections, resimulate 
    rejected = np.where(np.greater_equal(rejection, 3))[0]
    if rejected.size > 0:
        iterative2 = np.zeros(int(1+nparams))
        iterative2 = iterative2[np.newaxis,:]
        for nc in rejected: 
            paramset = param_old[:,nc]
            it = np.arange(int(ntrials*(kk/2)))
            x = np.hstack([it[:,np.newaxis], np.tile(paramset[np.newaxis,:], [int(ntrials*(kk/2)),1])])
            iterative2 = np.vstack([iterative2, x])

        iterative2 = iterative2[1:,:]
        with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
            future = executor.map(iterate1, iterative2, chunksize = int(np.ceil(rejected.size*kk)))#int(ncombo/100))
        future = tuple(future)

        for r1 in np.arange(rejected.size):
            LLold[rejected[r1]] = sum(future[(r1*int(ntrials*(kk/2))):((r1+1)*int(ntrials*(kk/2)))])
            rejection[rejected[r1]] = 0
           

       # resetting any chains with mean more than 1.96 * (sd of all other samples of all chains except those) to the mean of the samples outside those chains
    if np.equal(tstep, nburn/2-1): #at the midpoint of burn-in, reset chains that are stuck and outliers
        p = np.zeros(16)
        for j in np.arange(nchains):
            subset = copy.deepcopy(chainmem)
            subset[:(tstep+1),:, j] = np.nan
            meanpars = np.concatenate(subset, axis = 1)
            meanpars = np.nanmean(meanpars, axis = 1 )
            sdpars = np.concatenate(subset, axis = 1)
            sdpars = np.nanstd(sdpars, axis = 1 )
            chainposition = np.nanmean(chainmem[:(tstep+1),:,j],0)
            p = np.vstack((p, np.greater(np.abs(chainposition-meanpars),1.96*sdpars)))
        p = p[1:,:]    

        resetchain = np.where(np.greater(np.sum(p,1),0))[0]
        noreset =  np.where(~np.greater(np.sum(p,1),0))[0]

        subset = copy.deepcopy(chainmem)
        subset = subset[:(tstep+1),:,noreset]
        meanpars = np.concatenate(subset, axis = 1)
        meanpars = np.nanmean(meanpars, axis = 1 )

        paramset = meanpars
        it = np.arange(int(ntrials*(kk/2)))
        x = np.hstack([it[:,np.newaxis], np.tile(paramset[np.newaxis,:], [int(ntrials*(kk/2)),1])])
        ss = time.perf_counter()
        with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
            future = executor.map(iterate1, x)#int(ncombo/100))
        future = tuple(future)
        ee = time.perf_counter()
        ee-ss
        print(ee-ss)
        ll = sum(future)

        for nc in resetchain:
            param_old[:,nc] = meanpars

        LLold[resetchain] = ll
        rejection[resetchain] = 0

        #end of burn-in
        nreset = noreset.size

    #update the likelihood and the parameter combination in memory
    chainLL[tstep,:] = LLold 
    chainmem[tstep,:] = param_old
   
   
    #saving every step in burst
    saved = {'laststep': tstep,  'chainLL' :chainLL, 'chainmem':chainmem,'acceptance':acceptance, 'priorold':priorold, 'paramold':param_old, 'llold':LLold, 'rejection':rejection, 'nreset' : nreset, 'migratevevolve' :migratevevolve, 'nmigrations':nmigrations, 'migratingchains':migratingchains}
    with open(savefilenameburst, 'wb') as f:
        pickle.dump(saved, f)
    ee = time.perf_counter()
    print(ee-ss)
    
    #saving every 100 steps in scratch
    if np.equal((tstep+1)%100, 0):
        with open(savefilenamescratch, 'wb') as f:
            pickle.dump(saved, f)







