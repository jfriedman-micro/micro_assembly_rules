#!/usr/bin/env python

'''
@author: jonathanfriedman

'''


import pandas as pd
from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def modify_outcomes(outcomes, modifiers):
    new_outcomes = outcomes.copy()
    nodes = outcomes.index
    losers = []
    for n_mod, mod in modifiers.iteritems():
        if n_mod not in nodes:
            continue
        for m in mod:
            n1,n2,o = m
            if (n1 not in nodes) or (n2 not in nodes):
                continue
            if o =='COX':
                new_outcomes.loc[n1,n2] = o
                new_outcomes.loc[n2,n1] = o
            else:
                winner1 = n_mod
                winner2 = n1 if n1==o else n2
                loser = n2 if n1==o else n1
                losers.append(loser)
    return new_outcomes, losers

def is_excluded(s, sps, outcomes, losers):
    '''
    Is species s excluded by anybody in a community composed of species set sps.
    Outcomes already include coexistence modifications.
    Losers are species that are excluded by a pair, rather than by a single species.
    '''
    if s in losers:
        return True
    
    others = list(sps)
    if s in others:
        others.remove(s)
    try:
        os = set(outcomes.loc[s,others])
    except:
        print s, others
        print outcomes, '\n'
        raise ValueError
#     print os
    for other in others:
        if other in os:
            return True
    
    return False

def can_invade(invader, residents, outcomes, losers):
    '''
    Can species invader invade species residents.
    A species can invade if it is not excluded by any of the residents or residents pairs.
    
    Outcomes is the modified outcomes in the presence of both residents and invader.
    '''
    return not is_excluded(invader, residents, outcomes, losers)

def can_any_invade(all_species, residents, outcomes, modifications):
    '''
    Can any of the non-resident species invade?
    
    Outcomes should be the full set of unmofified outcomes for both residents and invaders
    '''
    invaders = set(all_species) - set(residents)
    for invader in invaders:
        ss_inv = list(residents) + [invader]
        os = outcomes.reindex(ss_inv,ss_inv)
        new_outcomes, losers = modify_outcomes(os, modifications)
        if can_invade(invader, residents, new_outcomes, losers):
#             print invader
            return True
    return False

def is_any_excluded(sps, outcomes, losers):
    return any([is_excluded(s, sps, outcomes, losers) for s in sps])

def get_survival_options(sps, outcomes, modifications):
    '''
    Find all potential outcomes of competition between a set of species with speciefied pairwise outcomes and trio modifications.

    Parameters
    ----------
    sps : iterable
    	names of competing species.
    outcomes: DataFrame
    	Pairwise outcomes of competing speices.
    	Each outcome should be either the name of the surviving species, or 'COX' if both species survive.
    	Bistability is not currently supported.
    modifications: mappable
    	Effective modifications of pairwise outcomes. These are derived from trio competitions.
    	Key =  name of the modifying species.
    	Values = list of modification caused by the key species. Each modification is a 3-tuple: (species1, species2, modified outcome).
    
    	Example: 
    	- Species A and B coexist in the presence of species C.
    	- Species A excludes D in the presence of species C.
    	- Species D and E coexist in the presence of species B.
    	 modifications = {'C':[('A', 'B', 'COX'), ('A', 'D', 'A'),], 'B':[('D', 'E', 'COX')]}

    Returns
    -------
    survival : DataFrame
    	Predicted survival of all species. 
    	Values denote extinction (0), or survival (1). 
    	rows = species
    	columns = predicted outcomes
    '''
    survival = pd.DataFrame(index=sps)
    i = 0
    for ss in powerset(sps):
        os = outcomes.reindex(ss,ss)
        ##losers are species that are excluded by a pair, rather than by a single species
        new_outcomes, losers = modify_outcomes(os, modifications) 
        if not is_any_excluded(ss, new_outcomes, losers) and not can_any_invade(sps, ss, outcomes, modifications):
            survival[i] = 0
            survival.loc[list(ss), i] = 1
            i += 1
#             print ss
    return survival


def demo_unique_prediction():
	sps = ('A', 'B', 'C', 'D', 'E')
	outcomes = pd.DataFrame('COX', index=sps, columns=sps)
	outcomes.loc['A','B'] = 'A'
	outcomes.loc['B','A'] = 'A'
	outcomes.loc['D','E'] = 'D'
	outcomes.loc['E','D'] = 'D'

	modifications = {'C':[('A', 'B', 'COX'), ('A', 'D', 'A'),], 'B':[('D', 'E', 'COX')]}  
	pred = pd.DataFrame(0, index=sps, columns=['pair only', 'pair + trio'])
	pred['pair only'] = get_survival_options(sps, outcomes, {})
	pred['pair + trio'] = get_survival_options(sps, outcomes, modifications)
	print pred

def demo_multiple_predictions():
	sps = ('A', 'B', 'C', 'D')
	outcomes = pd.DataFrame('COX', index=sps, columns=sps)

	modifications = {'B':[('D', 'C', 'C')], 'A':[('D', 'C', 'D')]}  
	pred_temp = pd.DataFrame(0, index=['type', 'pred #'] + list(sps), columns=(0,1,2))
	pred_temp.loc['type'] = ['pair only', 'pair + trio', 'pair + trio']
	pred_temp.loc['pred #'] = [1,1,2]
	pred = pred_temp.T.set_index(['type', 'pred #']).T
	pred.iloc[:,0] = get_survival_options(sps, outcomes, {})
	pred.iloc[:,1:] = get_survival_options(sps, outcomes, modifications).values
	print pred


if __name__ == '__main__':
	demo_multiple_predictions()