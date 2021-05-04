'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.3
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  modpilgrim         |
| Sub-module :  kmc                |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''


import numpy as np
import random
import time
import common.physcons as pc

from common.Exceptions import NoReacMol
from common.criteria import EPS_FSTRICT



#=====================================================#
def calculate_propensities(dxvec, processes):
    '''
    Calculates the propensitie of each reaction
    '''
    propensities = []
    total_propensity = 0.0
    for Rs,Ps,ks in processes:
        prop = ks
        for Ri in Rs: prop *= dxvec[Ri]
        propensities.append( prop )
        total_propensity += prop
    return propensities, total_propensity
#-----------------------------------------------------#
def generate_random(eps=EPS_FSTRICT):
    '''
    Generates a random number, excluding those smaller than eps
    '''
    while True:
      randx=random.random() 
      if randx >= eps: return randx
#=====================================================#




#=====================================================#
def kmc(ipops, processes, excess_species=None, volume=1.0/pc.ML, nstpdata=1000):
    '''
    ipops   : dictionary with initial populations (only needed those != 0.0)
    processes: a list with the reactants, products and the rate constant
               processes[idx] = (Rs,Ps,k)

    INPUT UNITS: atomic units
                 concentrations: molecules/bohr**3
                 volume        : bohr**3
                 rate constants: in au
    '''

    #----------------------------#
    # set initial concentrations #
    # and reactant molecules     #
    #----------------------------#
    reactants = []
    for Rs,Ps,k in processes:
        for species in Rs+Ps:
            if species not in ipops.keys(): ipops[species] = 0.0
            bool1 = species in excess_species
            bool2 = species in reactants
            if (not bool1) and (not bool2): reactants.append(species)

    #------------------------#
    # rate constants to s^-1 #
    #------------------------#
    for idx,(Rs,Ps,k) in enumerate(processes):
        nR = len(Rs)
        k /= volume**(nR-1)
        processes[idx] = (Rs,Ps,k)

    #-----------------------------------------#
    # Get dict of xvec and limiting molecules #
    #-----------------------------------------#
    dxvec = ipops.copy()
    try   : N0 = min([pop for pop in dxvec.values() if pop != 0.0])
    except: N0 = 0.0
    if N0 == 0.0: raise NoReacMol(Exception)
        
    # --------------------------#
    # START KINETIC MONTE-CARLO #
    # --------------------------#
    # initialize variables
    tau, tx = 0.0, 0.0
    jcount  = 0

    # data
    xvalues = [tx]
    yvalues = {key:[val] for key,val in dxvec.items()}

    Nj = N0
    xi = np.array([ dxvec[species] for species in sorted(dxvec.keys())])
    while Nj > 0.0:
       # compare each 1000 steps
       if jcount % 1000 == 0 and jcount > 0:
          xj = np.array([ dxvec[species] for species in sorted(dxvec.keys())])
          diff = np.linalg.norm(xj-xi)/1000
          if diff < 1e-4: break
          xi = xj
       # calculate propensities
       propensities, tot_propensity = calculate_propensities(dxvec,processes)
       #if tot_propensity < 1.e-15: break
       # generate two random numbers 
       r1=generate_random()    # 1) Tau (time)
       r2=generate_random()    # 2) Changes in population
       # Select process (stacking)
       value = tot_propensity * r2
       sum_props = 0.0
       for target,prop in enumerate(propensities):
           sum_props += prop
           if sum_props >= value: break
       # Modify populations (except those in excess)
       Rs,Ps,ks = processes[target]
       for Ri in Rs:
           if Ri not in excess_species: dxvec[Ri] -= 1.0
       for Pi in Ps:
           if Pi not in excess_species: dxvec[Pi] += 1.0
       # Time step
       if tot_propensity == 0.0: break
       tau = np.log(1./r1)/tot_propensity
       tx += tau
       # Keep data
       jcount += 1
       if jcount%nstpdata==0:
         xvalues.append(tx)
         for specie in dxvec.keys():
             yvalues[specie].append( dxvec[specie] )
    
       # Calculate current number of reactant molecules
       Nj = sum([dxvec[specie] for specie in reactants])

    return xvalues, yvalues
#=====================================================#


