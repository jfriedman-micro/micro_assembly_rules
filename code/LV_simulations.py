#!/usr/bin/env python

'''
@author: jonathanfriedman

'''


import pandas as pd
import scipy.integrate as spi
import scipy.stats as stats
import copy

from itertools import chain, combinations
from numpy import ones, vstack, full, diag, arange, nan, array
from assembly_rule import get_survival_options
from scipy.special import binom


########################
#### LV simulations ####
########################
def LV_fun(x, t, r, k, A):
    dx = r*x*(1-A.dot(x)/k)
    return dx

def LV_simulate(y0, t, r, k, A):
    args = (r,k,A)
    y = spi.odeint(LV_fun, y0, t, args)
    return y


######################################
#### generate random interactions ####
######################################
def get_A_pair(A_dist=stats.norm(.6,.46), allow_bis=False):
    if allow_bis:
        return A_dist.rvs(2)
    
    found = False
    while not found:
        a = A_dist.rvs(2)
        found = not all(a>1)
    return a
        
def get_A(n, A_dist=stats.norm(.6,.46), allow_bis=False):
    sps = arange(n)
    A = pd.DataFrame(1, index=sps, columns=sps)
    for s1,s2 in combinations(sps, 2):
        a = get_A_pair(A_dist, allow_bis)
        A.loc[s1,s2] = a[0]
        A.loc[s2,s1] = a[1]
    return A

def A2otucomes(A):
    outcomes = pd.DataFrame(index=A.index, columns=A.columns)
    for s1,s2 in combinations(A.index,2):
        a12, a21 = A.loc[s1,s2], A.loc[s2,s1]
        if a12<1 and a21<1:
            o = 'COX'
        elif a12>1 and a21<1:
            o = s2
        elif a12<1 and a21>1:
            o = s1
        else:
            o = 'BIS'
        outcomes.loc[s1,s2] = o
        outcomes.loc[s2,s1] = o
    return outcomes


#############################
#### infer modifications ####
#############################
def outcomes_given_survival(survival, original_outcomes):
    '''
    Get a simple pairwise outcome matrix that's consistent with the observed survival.
    Currently only works for 3 species
    '''
    sps = survival.index
    if len(sps)!=3:
        raise NotImplementedError, '%d species communities are not yet supported'%len(sps)
    outcomes = original_outcomes.copy()
    sps1 = survival.index[survival==1]
    sps0 = survival.index[survival==0]
    for s1,s2 in combinations(sps1,2):
        outcomes.loc[s1,s2] = 'COX'
        outcomes.loc[s2,s1] = 'COX'
    for s0 in sps0:
        for s1 in sps1:
#             print s0,s1
            if outcomes.loc[s0, s1] == s1:
#                 print 'b1'
                break
            elif outcomes.loc[s0, s1] == s0:
                outcomes.loc[s0,s1] = s1
                outcomes.loc[s1,s0] = s1
#                 print 'b2'
                break
        outcomes.loc[s0,s1] = s1
        outcomes.loc[s1,s0] = s1

    return outcomes

def update_modifications(survival, pair_outcomes, modifications):
    '''
    Update trio modification
    
    Modifications are stored as a dictionary where:
    keys   = the modifying species and 
    values = list of 3-tuple storing the interacting pair and the new outcome 
    '''
    new_out = outcomes_given_survival(survival, pair_outcomes)
    sps = new_out.index
    for s1,s2 in combinations(sps,2):
        s3 = [s for s in sps if s not in [s1,s2]][0]
        if s3 not in modifications:
            modifications[s3] = []
        o = pair_outcomes.loc[s1,s2]
        new_o = new_out.loc[s1,s2]
        if o!=new_o:
            modifications[s3].append((s1,s2, new_o))
    return modifications
        

#############################
#### prediction accuracy ####
#############################
def get_survival(y, sps, th=1e-5):
    return pd.Series((y[-1]>th).astype(int), index=sps)

def get_prediction_accuracy(pred, obs):
    FN = ((pred==0) &  (obs==1)).sum()
    FP = ((pred==1) &  (obs==0)).sum()
    return pd.Series(array([FP,FN]), index=['FP', 'FN'])

def get_pair_prediction_error(preds, A, err_type='mean'):    
    if isinstance(preds, str):
        if preds=='majority':
            make_preds = True
            n_preds = 1
        else:
            raise ValueError, 'Unsupported preds argument encountered: %s' %preds
    else:
        preds = pd.DataFrame(preds)

        if preds.shape[1]==0:
            return pd.Series(nan, index=['FP', 'FN']), nan
    
        n_preds = preds.shape[1]
        make_preds = False

    n = A.shape[0]
    r = ones(n)
    k = ones(n)    
    t = arange(500)
    y0s = vstack( (full(n, 1./n), .05*ones((n,n)) + diag([1 - .05*n]*n)) )
    n_y0s = y0s.shape[0]
    
    err = pd.Series(0, index=['FP', 'FN'])
    for y0 in y0s:
        y = LV_simulate(y0, t, r, k, A.values)
        survival = get_survival(y, A.index)
        if make_preds:
            f_survival = 1.*survival.sum()/n
            if f_survival >= .5:
                preds = pd.DataFrame(1, columns=[0], index=A.index)
            else:
                preds = pd.DataFrame(0, columns=[0], index=A.index)
        
        es = pd.DataFrame(nan, index=range(n_preds), columns=['FP','FN'])
        for i,p in enumerate(preds):
            e = get_prediction_accuracy(preds[p], survival)
#             print pd.concat((preds[p], survival), axis=1), '\n'
            es.iloc[i] = e
        if err_type is 'min':
            err += es.min()
        elif err_type is 'mean':
            err += es.mean()
        elif err_type is 'median':
            err += es.median()
    return err/(n+1.), survival

def get_group_error(n, sps, outcomes, A, modifications={}, infer_modifications=False, err_type='mean'):
    '''
    Ge the average prediction error for all groups of size n from given sps.
    
    n : int
        Group size.
    sps : list
        All the species. 
    outcome : DataFrame
        pairwise outcomes. 
        All sps should be in index/column labels.
    A : DataFrame
        LV interaction coeeficients.
        All sps should be in index/column labels.
    '''
    if infer_modifications and n!=3:
        raise ValueError, 'Modifications can only be inferred when group size is 3. Given group size is %d'%n
    modifications_new = copy.deepcopy(modifications)
    err = pd.Series(0, index=['FP', 'FN'])
    err_all_survive = err.copy()
    for grp in combinations(sps,n):
        out = outcomes.reindex(grp,grp)
        
        preds = get_survival_options(grp, outcomes, modifications)
        if not preds.shape[1]: 
            preds = pd.DataFrame(0, columns=[0], index=grp)
        e, obs = get_pair_prediction_error(preds, A.reindex(grp,grp), err_type=err_type)
        err += e
        
        pred_all_survive = obs.copy()
        pred_all_survive[:] = 1
        err_all_survive += get_prediction_accuracy(pred_all_survive, obs)
        
        if infer_modifications:
            modifications_new = update_modifications(obs, out, modifications_new)
           
    n_grp = binom(len(sps), n)
    err_avg = err/n_grp/n
    err_all_survive_avg = err_all_survive/n_grp/n
    return err_avg, err_all_survive_avg, modifications_new

###############
#### demos ####
###############
def accuracy_demo(n_tot=8, n_group=7, A_dist=stats.norm(.6,.46)):
	'''
	Get accuracy of prediction made use pairwise outcomes only, or including trio outcomes as well.
	Accuracy is avergaed over all groups of size n_group in a set of n_tot species.

	Species interacte according to the gLV model, with random interaction strengths drawn from A_dist. 
	'''
    A = get_A(n_tot, A_dist)
    sps = A.index
    outcomes = A2otucomes(A)

    ## infer modifications
    e_tmp, _, modifications = get_group_error(3, sps, outcomes, A, infer_modifications=True)
    
    ## get prediction error for all 7 species sets
    e2, e_all_survive, _ = get_group_error(n_group, sps, outcomes, A) #pair only prediction
    e3, _, _ = get_group_error(n_group, sps, outcomes, A, modifications=modifications) #pair + trio prediction
    
    print '\n'
    print 'Average prediction accuray for all %d species competitions from an %d species set:'%(n_group, n_tot)
    print 'pair only accuracy   = %.2f'%(1-e2.sum())
    print 'pair + trio accuracy = %.2f'%(1-e3.sum())

if __name__ == '__main__':
	accuracy_demo(5,4)