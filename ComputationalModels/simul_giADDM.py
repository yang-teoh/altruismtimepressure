#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: teohyiyang
"""
#Packages.
import os
from scipy import io
import numpy as np
import pandas as pd
import giADDM
from scipy import io as input
import sklearn
from sklearn.linear_model import LogisticRegression


#Simulates 20 subjects; 160 trials (80 trials high tp, 80 trials low tp)

for sub1 in np.arange(30):
    ##For each subject
    sub = str(sub1)
    
    path1 = os.getcwd()
    
    precision = float(.001)
    s = float(.1)#ddm noise
    noisematrix = np.random.normal(0,1,100000)*s*np.sqrt(precision)
    
    resp_options = [1,2] #1 is reject, 2 is accept
    
    
    wself_init = [.01, .01]
    wother_init = [.005, .01]
    wfair_init = [.005, .005]
    bound_init = [.32, .1]
    collapse_init = [.7, .1]
    stbias_init = [.025, .1]
    genbias_init = [.025, .1]
    theta_init = [.45, .25]
    
    
    #Params to simulate
    wself = np.repeat(np.max([0.005, np.random.normal(wself_init[0], wself_init[1],1)]),2)
    wself += np.random.uniform(0,.01,2)
    wother = np.repeat(np.max([0.003,np.random.normal(wother_init[0], wother_init[1],1)]),2)
    wother += np.random.uniform(0,.01,2)
    wfair =  np.repeat(np.max([0.003,np.random.normal(wfair_init[0], wfair_init[1],1)]),2)
    wfair += np.random.uniform(0,.01,2)
    bound = np.repeat(np.random.normal(bound_init[0], bound_init[1],1),2)
    bound += np.random.uniform(-.05,.05,2)
    collapse = np.repeat(np.random.normal(collapse_init[0], collapse_init[1],1),2)
    collapse = np.asarray([collapse[0] + np.random.uniform(.2,1), collapse[1] - np.random.uniform(.1,.5)])
    stbias = np.repeat(np.random.normal(stbias_init[0], stbias_init[1],1),2)
    stbias  += np.random.uniform(-.05,.05,2)
    genbias = np.repeat(np.random.normal(genbias_init[0], genbias_init[1],1),2)
    genbias += np.random.uniform(-.05,.05,2)
    theta = np.repeat(np.min([np.max([0.25,np.random.normal(theta_init[0], theta_init[1],1)]),.9]),2)
    theta += np.random.uniform(-.1,.1,2)
    
    motor = .08
    percept = .08
    
    ffixbias = [np.nanmean(wself)/(np.nanmean(wself)+np.nanmean(wother)), \
                np.max([0,np.min([np.random.normal(0.5, .05), 1])])]
    
    
    # probabilistic selfish fixation bias
    fixdurm = np.log(200) #mean log(length in ms) of fixations
    fixdursd = .2 # sd of log(fixation length)
    sacdur = [50, 150] #uniform distribution of saccade length in ms
    
    timelimit = [1.5,10] #time limit for trial
    fix = list()
    fixdf = list()

    outcomes = list() #simulated choice & rt
    lastj = -1
    lastf = -1
    ntrialpblock = 10
    trialn = 0
    for block in np.arange(8):
        if block < 2:
            selfamt = np.asarray([0,5,3,25,28,30,48,48,48,50,50,53,53,53,80,85,75,98,100,100])  #self vals 
            otheramt = np.asarray([68,83,100,55,85,100,60,80,100,100,43,40,23,0,45,15,0,38,15,0]) # other vals
            initamts = np.hstack((selfamt[:,np.newaxis], otheramt[:,np.newaxis]))
            if block == 0:
                amts = np.tile(initamts[:ntrialpblock,:], (2,1))
            else:
                amts = np.vstack((amts, np.tile(initamts[ntrialpblock:,:], (2,1))))
        else:
            #create list of fixation data 
            gen = np.logical_or(np.logical_and(np.equal(np.asarray(outcomes)[:,1],resp_options[0]), np.greater(amts[:,0], amts[:,1])),\
                                np.logical_and(np.equal(np.asarray(outcomes)[:,1],resp_options[1]), np.less(amts[:,0], amts[:,1])))
            gen = gen.astype(int)
            
            if np.nanmean(gen) >= .25:
                trialsahd = np.greater(amts[:,0], amts[:,1])
                selfahd = amts[:,0]* trialsahd.astype(int)
                selfbeh = amts[:,0] * np.logical_not(trialsahd).astype(int)
                otherahd = amts[:,1]* np.logical_not(trialsahd).astype(int)
                otherbeh = amts[:,1] * trialsahd.astype(int)
                model = LogisticRegression(solver='liblinear', random_state=0).fit(np.hstack((selfahd[:,np.newaxis],\
                                                                                              selfbeh[:,np.newaxis],\
                                                                                              otherahd[:,np.newaxis],\
                                                                                              otherbeh[:,np.newaxis])),gen)
                
                posselfamt = np.arange(0,101,1) #range of self vals + 50
                posotheramt = np.arange(0,101,1) # range of other vals + 50
                
                posamts = np.array(np.meshgrid(posselfamt, posotheramt)).T.reshape(-1,2) #all unique combinations
                
                posamts = posamts[np.where(np.logical_or(np.logical_and(np.greater(posamts[:,0],50),np.less(posamts[:,1],50)),
                                                   np.logical_and(np.less(posamts[:,0],50),np.greater(posamts[:,1],50))))] #only retrieve combinations where conflict between self and other
            
                postrialsahd = np.greater(posamts[:,0], posamts[:,1])
                
                posselfahd = posamts[:,0]* postrialsahd.astype(int)
                posselfbeh = posamts[:,0] * np.logical_not(postrialsahd).astype(int)
                posotherahd = posamts[:,1]* np.logical_not(postrialsahd).astype(int)
                posotherbeh = posamts[:,1] * postrialsahd.astype(int)
                
                probs = model.predict_proba(np.hstack((posselfahd[:,np.newaxis],
                                               posselfbeh[:,np.newaxis],
                                               posotherahd[:,np.newaxis],
                                               posotherbeh[:,np.newaxis])))
                
                numcan = 3
                
                potentialg = np.where(probs[:,1] > .9)[0]
                if len(potentialg) < numcan:
                    cang =  np.argsort(probs[:,1])[-3:]
                else:
                    cang = np.random.choice(potentialg,numcan)  
                potentials = np.where(probs[:,1] < .1)[0]
                if len(potentials) < numcan:
                    cans =  np.argsort(probs[:,1])[:3]
                else:
                    cans = np.random.choice(potentials,numcan)  
                rest = np.delete(np.arange(probs.shape[0]),np.concatenate((cang, cans)))
                restcan = np.random.choice(rest,10-(numcan*2))
                canamts = posamts[np.random.permutation(np.hstack((cang, cans, restcan))),:]
            else:
                posselfamt = np.arange(0,101,1) #range of self vals + 50
                posotheramt = np.arange(0,101,1) # range of other vals + 50
                
                posamts = np.array(np.meshgrid(posselfamt, posotheramt)).T.reshape(-1,2) #all unique combinations
                
                posamts = posamts[np.where(np.logical_or(np.logical_and(np.greater(posamts[:,0],50),np.less(posamts[:,1],50)),
                                                   np.logical_and(np.less(posamts[:,0],50),np.greater(posamts[:,1],50))))] #only retrieve combinations where conflict between self and other
                
                potentialg = np.where(np.logical_and(np.logical_and(posamts[:,0] > 50 , posamts[:,0] <55),
                                                     np.logical_and(posamts[:,1] > 0 , posamts[:,1] <10)))[0]
                cang = np.random.choice(potentialg,numcan)
                rest = np.delete(np.arange(probs.shape[0]),cang)
                restcan = np.random.choice(rest,10-(numcan))
                
                canamts = posamts[np.random.permutation(np.hstack((cang, restcan))),:]
                
            amts = np.vstack((amts, canamts, canamts))
            
        for cond in np.arange(len(timelimit)):
            for f in np.arange(ntrialpblock):
                
                rstfix = np.random.uniform(0,1)
                stfix = rstfix < ffixbias[cond] #get first fixation
                if stfix:
                    fixnow = 1 # 1 is self
                else:
                    fixnow = 2 # 2 is other
                nfix = 1
                #pre-fixation entry is always np.nan 
                fixmat = np.asarray([cond*amts.shape[0]+f+1,nfix,0, np.round(np.random.uniform(sacdur[0], sacdur[1])), timelimit[cond]*1000,np.nan])
                                   
                nfix +=1 
                
                #next fixation (first in AOI)
                ed = fixmat[3] +1 + np.round(np.exp(np.random.normal(fixdurm, fixdursd)))
                st = fixmat[3] + 1
                                    
                #simulate subsequent fixations (alternating between self & other) with a lagtime between fixations for saccades
                while ed < (timelimit[cond]*1000 - 1):
                    
                    fixmat1 = np.asarray([cond*amts.shape[0]+f+1,nfix,st, ed, timelimit[cond]*1000, fixnow])
                    fixmat = np.vstack((fixmat,fixmat1))
                    nfix += 1
                    
                    if fixnow == 1:
                        fixnow = 2
                    else:
                        fixnow = 1
                    
                    st = ed + 1+np.round(np.random.uniform(sacdur[0], sacdur[1]))
                    ed = st + 1 +np.round(np.exp(np.random.normal(fixdurm, fixdursd)))
                    
                fix.append(fixmat)
                fixmat = pd.DataFrame(fixmat, columns = ['trial', 'nfix','start','end','triallength','fixpos'])
                fixdf.append(fixmat) #each trial has a dataframe of fixations.
                
                if amts[trialn,0] > amts[trialn, 1]:
                    modbias = stbias[cond] - genbias[cond]
                else:
                    modbias = stbias[cond] + genbias[cond]
                    
                simmaxrt = timelimit[cond] - motor
                
                params = {'selfamt' : (amts[trialn,0]-50)/5 ,\
                          'otheramt' : (amts[trialn,1]-50)/5 ,\
                          'eye' : np.asarray(fix[trialn]) ,\
                          'self' : wself[cond] ,\
                          'other' : wother[cond] ,\
                          'fair' : wfair[cond] ,\
                          'bound' : bound[cond] ,\
                          'collapse' : collapse[cond] ,\
                          'modbias' : modbias,\
                          'theta' : theta[cond] ,\
                          'percept' : percept ,\
                          'simmaxrt' : simmaxrt}
                
                outcome1 = giADDM.simulEADDMjit(params, 1, noisematrix, precision)
                
                if outcome1['resp'] == 0:
                    outcome1['resp'] = [resp_options[0]]
                elif outcome1['resp'] == 1:
                    outcome1['resp'] = [resp_options[1]]
                    
                
                outcomes.append(np.asarray([trialn, outcome1['resp'][0], outcome1['rt'][0]+motor, timelimit[cond]])) #adding motor response time to rts
      
                trialn += 1
        
    fix = tuple(fix) #tuple of dataframes for all trials
    fixdf = tuple(fixdf)
    outcomes = np.asarray(outcomes)
    
    
    
    simuldata = pd.DataFrame(data = np.hstack((outcomes,amts)),
                 columns = ['trial', 'resp','rt','time','selfamt','otheramt'])
    
    
    
    Data = {'ChoiceRT' :  np.asarray(simuldata['rt']), \
            'Resp' : np.asarray(simuldata['resp']),\
            'SelfProposal' :  np.asarray(simuldata['selfamt']),\
            'OtherProposal' :  np.asarray(simuldata['otheramt']),\
            'TimeLimit' : np.asarray( simuldata['time']),\
            'Fix' : np.asarray(fix),\
            'Generating_Params': [wself, wother, wfair, bound, collapse, stbias, genbias, theta],\
            'Fixbias': ffixbias}
        
    
    Data = {'Data':Data}    
    if not os.path.exists(path1+'/SimulatedData/'):
        os.mkdir(path1+'/SimulatedData/')
    if not os.path.exists(path1+'/SimulatedData/'+sub):
        os.mkdir(path1+'/SimulatedData/'+sub)
    
    filename = path1+'/SimulatedData/'+sub + '/Data.' + sub + '.wfix.choice.mat'
    io.savemat(filename, Data)


