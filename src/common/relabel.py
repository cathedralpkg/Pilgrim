'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.4
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
| Module     :  common             |
| Sub-module :  relabel            |
| Last Update:  2021/05/27 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''
#==================================================#
import os
import numpy           as     np
#--------------------------------------------------#
import common.fncs     as     fncs
import common.internal as     intl
from   common.fncs     import center_and_orient
from   common.Ugraph   import UGRAPH
#==================================================#


#===================================================#
def generate_enantio(xcc):
    return [-xcc[idx] if idx%3==0 else xcc[idx] for idx in range(len(xcc))]
#===================================================#

#===================================================#
def get_assigned(dcorr):
    # first, get assigned nodes
    assigned = [node_0 for node_t,node_0 in dcorr.items() if type(node_0) == int]
    # assert no repetitions in assigned
    num_assigned = len(assigned)
    assigned     = set(assigned)
    if len(assigned) != num_assigned: raise Exception
    return assigned,num_assigned
#---------------------------------------------------#
def get_num_of_not_assigned(dcorr):
    return sum([1 if type(nodes_0) != int else 0 for node_t,nodes_0 in dcorr.items()])
#---------------------------------------------------#
def assign_check(dcorr):
    # get assigned
    assigned,numassig1 = get_assigned(dcorr)
    # now, clean them up
    for node_t,node_0 in dcorr.items():
        # node already assigned
        if type(node_0) == int: continue
        # remove already assigned
        node_0 = node_0.difference(assigned)
        dcorr[node_t] = node_0
        # check length of node_0
        if   len(node_0 ) == 0: raise Exception
        elif len(node_0 ) == 1: dcorr[node_t] = list(node_0)[0]
    # any new assignation?
    assigned,numassig2 = get_assigned(dcorr)
    if numassig1 != numassig2: dcorr = assign_check(dcorr)[0]
    # return data
    return dcorr,get_num_of_not_assigned(dcorr) 
#---------------------------------------------------#
def assign_concatenate(dcorr,graph_0,graph_t):
    '''
    It tries to assign neighbors of an assigned node using connectivity.
    This process is executed until no new assignations.
    '''
    dcorr,na = assign_check(dcorr)
    while True:
        na_i = get_num_of_not_assigned(dcorr)
        dcorr2 = {k:v for k,v in dcorr.items()}
        for node_t,node_0 in dcorr.items():
            # node_t not assigned
            if type(node_0) != int: continue
            # node_t assigned --> check not assigned neighbors
            neighbors_t = [neigh for neigh in graph_t.neighbors(node_t) if type(dcorr[neigh]) != int]
            # correspoding neighbors in graph_0
            neighbors_0 = graph_0.neighbors(node_0)
            for neigh2 in neighbors_t:
                dcorr2[neigh2] = dcorr2[neigh2].intersection(neighbors_0)
        # check
        dcorr,na_j = assign_check(dcorr2)
        # return
        if na_i == na_j or na_j == 0: return dcorr,na_j
#===================================================#



#===================================================#
def assign_initialize(symbols_t,symbols_0):
    dcorr = {idx2:set([idx1 for idx1,s1 in enumerate(symbols_0) if s1 == s2]) \
                            for idx2,s2 in enumerate(symbols_t)}
    return assign_check(dcorr)
#---------------------------------------------------#
def assign_neighbors(dcorr,graph_0,graph_t,symbols_0,symbols_t):
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # Compare symbols for neighbors
    for node_t,nodes_0 in list(dcorr.items()):
        # already assigned?
        if type(nodes_0) == int: continue
        # get symbol of neighbors
        neighbors_t = sorted([symbols_t[idx] for idx in graph_t.neighbors(node_t)])
        # compare with graph_0
        dcorr[node_t] = set([])
        for node_0 in nodes_0:
            # get symbol of neighbors
            neighbors_0 = sorted([symbols_0[idx] for idx in graph_0.neighbors(node_0)])
            # compare neighbors
            if neighbors_0 == neighbors_t: dcorr[node_t].add(node_0)
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_layers(dcorr,graph_0,graph_t,symbols_0,symbols_t):
    '''
    if this is executed, assign_neighbors is not required,
    as the second layer (layer[1]) corresponds to neighors
    '''
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # Compare symbols for each layer
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        # layers node 2
        layers_t =  [sorted([symbols_t[at] for at in l]) for l in graph_t.get_layers(node_t)]
        # comparison graph_0 <--> graph_t
        dcorr[node_t] = set([])
        for node_0 in nodes_0:
            # layers node 1
            layers_0 = [sorted([symbols_0[at] for at in l]) for l in graph_0.get_layers(node_0)]
            if layers_0 == layers_t: dcorr[node_t].add(node_0)
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_pyramid(dcorr,graph_0,graph_t,symbols_0,symbols_t,xcc_0,xcc_t):
    '''
     A1  B1    |    B2  A2    X1 --> X1 or X2?
      \ /            \ /
       X1      |      X2     Ai!=Bi!=Ci!=Di
      / \            / \     
     C1  D1    |    D2  C2   
    '''
    for X_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        if  len(nodes_0) !=  2 : continue
        neighbors_t = graph_t.neighbors(X_t)
        # number of neighbors
        nn_t = len(neighbors_t)
        if nn_t not in [3,4]: continue
        # check all neghbors are different
        nsymbols_t = len(set([symbols_t[idx] for idx in neighbors_t]))
        if nsymbols_t != nn_t: continue
        # sort neighbors by symbol
        neighbors_t = sorted(neighbors_t, key=lambda idx:symbols_t[idx])
        # dihedral A-B-C-D or X-A-B-C?
        if nn_t == 4: atoms_t = neighbors_t
        else        : atoms_t =  [X_t]+list(neighbors_t)
        # calculate dihedral
        xs_t     = (xcc_t[3*at:3*at+3] for at in atoms_t)
        angle_t  = np.rad2deg(fncs.dihedral(*xs_t))%360
        confi_t  = angle_t > 180
        # Check same distribution
        toremove = []
        for X_0 in nodes_0:
            neighbors_0 = graph_0.neighbors(X_0)
            nn_0 = len(neighbors_0)
            if nn_0 != nn_t: continue
            # assert they are different
            nsymbols_0 = len(set([symbols_0[idx] for idx in neighbors_0]))
            if nsymbols_0 != nn_0: continue
            # sort neighbors by symbol
            neighbors_0 = sorted(neighbors_0, key=lambda idx:symbols_0[idx])
            # dihedral A-B-C-D or X-A-B-C?
            if nn_0 == 4: atoms_0 = neighbors_0
            else        : atoms_0 =  [X_0]+list(neighbors_0)
            # calculate dihedral
            xs_0    = (xcc_0[3*at:3*at+3] for at in atoms_0)
            angle_0 = np.rad2deg(fncs.dihedral(*xs_0))%360
            confi_0 = angle_0 > 180
            if confi_0 != confi_t: toremove.append(X_0)
        # remove
        for X_0 in toremove: nodes_0.remove(X_0)
        dcorr[X_t] = nodes_0
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_via_improper_case1(dcorr,graph_0,graph_t,xcc_0,xcc_t):
    '''
     C   D           correlation between A1, A2, ..., An
      \ /
  A1 - B - An        in this situation, B, C and D
      / \            are assigned so Ai-B-C-D improper
    A2  ...          torsion can be used
    '''
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # Compare spatial organization (improper torsions)
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        # check neighbors
        atB2,atC2,atD2 = None,None,None
        for neigh in graph_t.neighbors(node_t):
            if type(dcorr[neigh]) != int: continue
            atB2,atC2,atD2 = neigh,None,None
            for nneigh in graph_t.neighbors(atB2):
               #print(node_t+1,neigh+1,nneigh+1)
                if nneigh == node_t: continue
                if type(dcorr[nneigh]) != int: continue
                if atC2 is None: atC2 = nneigh; continue
                if atD2 is None: atD2 = nneigh; break
            if atD2 is not None: break
        if atD2 is None: continue
        # now check spatial configuration
        atoms_t  = [node_t,atB2,atC2,atD2]
        xs_t     = (xcc_t[3*at:3*at+3] for at in atoms_t)
        angle_t  = np.rad2deg(fncs.dihedral(*xs_t))%360
        dcorr[node_t] = set([])
        # correlate to closest one
        mindist, minnode = float("inf"),None
        for node_0 in nodes_0:
            atoms_0  = [node_0,dcorr[atB2],dcorr[atC2],dcorr[atD2]]
            xs_0     = (xcc_0[3*at:3*at+3] for at in atoms_0)
            angle_0  = np.rad2deg(fncs.dihedral(*xs_0))%360
            dist     = fncs.angular_dist(angle_0,angle_t,u="deg")
            if dist < mindist:
                mindist = dist
                minnode = node_0
        dcorr[node_t] = minnode
        # break to avoid assign same node
        break
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # do again?
    if na_i != na_j and na_j != 0:
       dcorr = assign_via_improper_case1(dcorr,graph_0,graph_t,xcc_0,xcc_t)[0]
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_via_improper_case2(dcorr,graph_0,graph_t,xcc_0,xcc_t):
    '''
     C   D           correlation between A1 & A2
      \ /
       B             in this situation, B, C are
      / \            assigned, but D is not!
    A1  A2           Assignation is done via A1-B-C-A2
    '''
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # Compare spatial organization (improper torsions)
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        if  len(nodes_0) !=  2 : continue
        BC = None
        # check neighbors
        atB2,atC2,atD2 = None,None,None
        for neigh in graph_t.neighbors(node_t):
            if type(dcorr[neigh]) != int: continue
            atB2,atC2,atD2 = neigh,None,None
            for nneigh in graph_t.neighbors(atB2):
               #print(node_t+1,neigh+1,nneigh+1)
                if nneigh == node_t: continue
                if type(dcorr[nneigh]) != int: continue
                if atC2 is None: atC2 = nneigh; continue
                if atD2 is None: atD2 = nneigh; break
            if atD2 is not None: break
            if atC2 is not None: BC = (atB2,atC2)
        # assert D is not assigned
        if atD2 is not None: continue
        # assert B and C are assigned
        if BC   is     None: continue
        # Unpack
        atB2,atC2 = BC
        # look for the other node with same assigments in graph_0
        node_ta = node_t
        node_tb = None
        for node_t_,nodes_0_ in dcorr.items():
            if node_t_ == node_ta : continue
            if nodes_0 != nodes_0_: continue
            node_tb = node_t_
            break
        if node_tb is None: continue
        # assert both nodes are connected to the same atoms!!
        if graph_t.neighbors(node_ta) != graph_t.neighbors(node_tb): continue
        # config in graph_t
        atoms2  = [node_ta,node_tb,atB2,atC2]
        xs2     = (xcc_t[3*at:3*at+3] for at in atoms2)
        angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
        config2 = angle2 < 180
        # config in graph_0
        node_0a,node_0b = nodes_0
        atoms1 = [node_0a,node_0b,dcorr[atB2],dcorr[atC2]]
        xs1     = (xcc_0[3*at:3*at+3] for at in atoms1)
        angle1  = np.rad2deg(fncs.dihedral(*xs1))%360
        config1 = angle1 < 180
        if config1 == config2:
           dcorr[node_ta] = node_0a
           dcorr[node_tb] = node_0b
        else:
           dcorr[node_ta] = node_0b
           dcorr[node_tb] = node_0a
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
##  def assign_spatial(dcorr,graph_0,graph_t,xcc_0,xcc_t):
##      assign_print(dcorr)
##      # number of non-assigned (initial)
##      na_i = get_num_of_not_assigned(dcorr)
##      # Compare spatial organization (improper torsions)
##      for node_t,nodes_0 in dcorr.items():
##          if type(nodes_0) == int: continue
##          if len(nodes_0) != 2: continue
##          print("*",node_t+1,[ii+1 for ii in nodes_0])
##          BC = None
##          # check neighbors
##          atB2,atC2,atD2 = None,None,None
##          for neigh in graph_t.neighbors(node_t):
##              if type(dcorr[neigh]) != int: continue
##              atB2,atC2,atD2 = neigh,None,None
##              for nneigh in graph_t.neighbors(atB2):
##                 #print(node_t+1,neigh+1,nneigh+1)
##                  if nneigh == node_t: continue
##                  if type(dcorr[nneigh]) != int: continue
##                  if atC2 is None: atC2 = nneigh; continue
##                  if atD2 is None: atD2 = nneigh; break
##              if atD2 is not None: break
##              if atC2 is not None: BC = (atB2,atC2)
##          # now check spatial configuration
##          if atD2 is not None:
##             atoms2  = [node_t,atB2,atC2,atD2]
##             xs2     = (xcc_t[3*at:3*at+3] for at in atoms2)
##             angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
##             config2 = angle2 < 180
##             dcorr[node_t] = set([])
##            #print("*",node_t+1,[i+1 for i in atoms2],angle2)
##             for node_0 in nodes_0:
##                 atoms1  = [node_0,dcorr[atB2],dcorr[atC2],dcorr[atD2]]
##                 xs1     = (xcc_0[3*at:3*at+3] for at in atoms1)
##                 angle1  = np.rad2deg(fncs.dihedral(*xs1))%360
##                 config1 = angle1 < 180
##                #print("*",node_0+1,[i+1 for i in atoms1],angle1)
##                 if config1 == config2:
##                    dcorr[node_t].add(node_0)
##            #print()
##          elif BC is not None and len(nodes_0) == 2:
##             atB2,atC2 = BC
##             # look for the other node with same assigments in graph_0
##             node_ta = node_t
##             node_tb = None
##             for node_t_,nodes_0_ in dcorr.items():
##                 if node_t_ == node_ta : continue
##                 if nodes_0 != nodes_0_: continue
##                 node_tb = node_t_
##                 break
##             if node_tb is None: continue
##             # assert both nodes are connected to the same atoms!!
##             if graph_t.neighbors(node_ta) != graph_t.neighbors(node_tb): continue
##             # config in graph_t
##             atoms2  = [node_ta,node_tb,atB2,atC2]
##             xs2     = (xcc_t[3*at:3*at+3] for at in atoms2)
##             angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
##             config2 = angle2 < 180
##             # config in graph_0
##             node_0a,node_0b = nodes_0
##             atoms1 = [node_0a,node_0b,dcorr[atB2],dcorr[atC2]]
##             xs1     = (xcc_0[3*at:3*at+3] for at in atoms1)
##             angle1  = np.rad2deg(fncs.dihedral(*xs1))%360
##             config1 = angle1 < 180
##             if config1 == config2:
##                dcorr[node_ta] = node_0a
##                dcorr[node_tb] = node_0b
##             else:
##                dcorr[node_ta] = node_0b
##                dcorr[node_tb] = node_0a
##      # number of non-assigned(final)
##      dcorr,na_j = assign_check(dcorr)
##      # return data
##      return assign_concatenate(dcorr,graph_0,graph_t)
##      if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph_0,graph_t)
##      else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_torsions(dcorr,graph_0,graph_t,xcc_0,xcc_t):
    '''
    only one at a time
    '''
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # assignment based in proper torsions
    minad = float("inf")
    minnode_0 = None
    minnode_t = None
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        BC = None
        # check neighbors
        atA2,atB2,atC2,atD2 = node_t,None,None,None
        for atB2 in graph_t.neighbors(atA2):
            if type(dcorr[atB2]) != int: continue
            for atC2 in graph_t.neighbors(atB2):
                if type(dcorr[atC2]) != int: continue
                for neigh in graph_t.neighbors(atC2):
                    if type(dcorr[neigh]) != int: continue
                    if len(set([atA2,atB2,atC2,neigh])) != 4: continue
                    atD2 = neigh
                    break
                if atD2 is not None: break
            if atD2 is not None: break
        # torsion located?
        if atD2 is None: continue
        # Now, check A-B-C-D
        atoms2  = [atA2,atB2,atC2,atD2]
        xs2     = (xcc_t[3*at:3*at+3] for at in atoms2)
        angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
        config2 = angle2 < 180
        # compare torsions
        for node_0 in nodes_0:
            atoms1 = [node_0,dcorr[atB2],dcorr[atC2],dcorr[atD2]]
            xs1    = (xcc_0[3*at:3*at+3] for at in atoms1)
            angle1 = np.rad2deg(fncs.dihedral(*xs1))
            # difference
            ad = fncs.angular_dist(-angle2,angle1,u="deg")
            if ad < minad:
                minad    = ad
                minnode_0 = node_0
                minnode_t = node_t
    # Any new assignation to do??
    if minnode_0 is None: return dcorr, na_i
    # assign
    dcorr[minnode_t] = minnode_0
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph_0,graph_t)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_cx3(dcorr,graph_0,graph_t,xcc_0,xcc_t,symbols_0,symbols_t):
    # number of non-assigned (initial)
    na_i = get_num_of_not_assigned(dcorr)
    # Assign Xs in CX3 groups (CH3, CF3, CCl3, et cetera)
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        if  len(nodes_0) !=   3: continue
        neighbors_t = graph_t.neighbors(node_t)
        atC = neighbors_t.pop()
        # check we have a C atom
        if len(neighbors_t) != 0: continue
        if symbols_t[atC] != "C": continue
        # determine torsion X-C-B-A
        atA = None
       #print("*",node_t+1,atC+1)
        for atB in graph_t.neighbors(atC):
            if type(dcorr[atB]) != int: continue
            for neigh in graph_t.neighbors(atB):
                if type(dcorr[neigh]) != int: continue
                if neigh == node_t: continue
                if neigh == atC  : continue
                atA = neigh
               #print(" ",atB+1,atA+1)
                break
            if atA is not None: break
        if atA is None: continue
        # calculate dihedral
        atoms2  = [node_t,atC,atB,atA]
       #print([i+1 for i in atoms2])
        xs2     = (xcc_t[3*at:3*at+3] for at in atoms2)
        angle2  = np.rad2deg(fncs.dihedral(*xs2))
        # compare with the others
        minad = float("inf")
        minnode = None
        for node_0 in nodes_0:
            atoms1 = [node_0,dcorr[atC],dcorr[atB],dcorr[atA]]
            xs1    = (xcc_0[3*at:3*at+3] for at in atoms1)
            angle1 = np.rad2deg(fncs.dihedral(*xs1))
            # difference
            ad = fncs.angular_dist(-angle2,angle1,u="deg")
            if ad < minad:
                minad   = ad
                minnode = node_0
        # apply
        dcorr[node_t] = minnode
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_force(dcorr,graph_0,graph_t):
    '''
    force assignment
    tries to assign a non-terminal node
    '''
    target2  = None
    target1  = None
    mindoubt = float("inf")
    ttarget2 = None # terminals
    ttarget1 = None # terminals
    tmindoubt = float("inf")
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int: continue
        nneighbors = len(graph_t.neighbors(node_t))
        # look for minimum doubt
        if nneighbors == 1 and len(nodes_0) < tmindoubt:
           tmindoubt = len(nodes_0)
           ttarget1,ttarget2 = sorted(list(nodes_0))[0], node_t
        if nneighbors != 1 and len(nodes_0) < mindoubt:
           mindoubt = len(nodes_0)
           target1,target2 = sorted(list(nodes_0))[0], node_t
    # apply target
    if tmindoubt <= mindoubt:
       if   ttarget2 is not None: dcorr[ttarget2] = ttarget1
       elif  target2 is not None: dcorr[target2 ] =  target1
    else:
       if    target2 is not None: dcorr[target2 ] =  target1
       elif ttarget2 is not None: dcorr[ttarget2] = ttarget1
    # return
    return assign_concatenate(dcorr,graph_0,graph_t)
#---------------------------------------------------#
def assign_print(dcorr):
    for node_t,nodes_0 in dcorr.items():
        if type(nodes_0) == int:
            print("%2i : %2i"%(node_t+1,nodes_0+1))
        else:
            print("%2i :"%(node_t+1),[node_0+1 for node_0 in nodes_0])
    print()
#===================================================#



#===================================================#
def correlate_amatrices(amatrix_t,amatrix_0,symbols_t,symbols_0,\
                        xcc_t=None,xcc_0=None,onlyconn=False):

    if not onlyconn and (xcc_t is None or xcc_0 is None): raise Exception
    # Graphs
    graph_0 = UGRAPH(); graph_0.set_from_amatrix(amatrix_0); graph_0.gen_laplacian()
    graph_t = UGRAPH(); graph_t.set_from_amatrix(amatrix_t); graph_t.gen_laplacian()

    
    # Correlations
    dcorr,na_0 = assign_initialize(symbols_t,symbols_0)
    nforce = 0
    for ii in range(10): # it should be a while True but... just in case
        # use symbols of neighbors
        dcorr,na_1  = assign_neighbors(dcorr,graph_0,graph_t,symbols_0,symbols_t)
        if na_1 == 0: break
        # use symbols of layers
        dcorr,na_2  = assign_layers(dcorr,graph_0,graph_t,symbols_0,symbols_t)
        if na_2 == 0: break
        # only asked for connectivity
        if onlyconn: break
        # use pyramid with 4 different vertices
        dcorr,na_j = assign_pyramid(dcorr,graph_0,graph_t,symbols_0,symbols_t,xcc_0,xcc_t)
        if na_j == 0: break
        # use improper torsions (case 1)
        dcorr,na_3 = assign_via_improper_case1(dcorr,graph_0,graph_t,xcc_0,xcc_t)
        if na_3 == 0: break
        # use improper torsions (case 2)
        dcorr,na_4 = assign_via_improper_case2(dcorr,graph_0,graph_t,xcc_0,xcc_t)
        if na_4 == 0: break
        # assign CX3 groups
        dcorr,na_5  = assign_cx3(dcorr,graph_0,graph_t,xcc_0,xcc_t,symbols_0,symbols_t)
        if na_5 == 0: break
        # assign based on torsions
        dcorr,na_6  = assign_torsions(dcorr,graph_0,graph_t,xcc_0,xcc_t)
        if na_6 == 0: break
        # force assignation
        if na_0 == na_6:
           nforce += 1
           dcorr,na_6 = assign_force(dcorr,graph_0,graph_t)
        # update na_0
        na_0 = na_6
    # return dcorr
    return dcorr
#---------------------------------------------------#
def equivalent_atoms_by_connectivity(xcc,symbols,fconnect=1.3):
    amatrix = np.matrix(intl.get_adjmatrix(xcc,symbols,fconnect,"int")[0])
    amatrix = np.matrix(intl.link_fragments(xcc,amatrix.tolist(),1)[0])
    dcorr   = correlate_amatrices(amatrix,amatrix,symbols,symbols,onlyconn=True)
    return dcorr
#---------------------------------------------------#
def correlate_xccs(xcc_t,symbols_t,xcc_0,symbols_0,fconnect=1.3,pp=False):
    '''
    xcc_0 --> geometry
    xcc_t --> enantiomer
    '''

    if len(xcc_0) != len(xcc_t)      : raise Exception
    if len(xcc_0) != len(symbols_0)*3: raise Exception
    
    # Adjacency matrices
    amatrix_0 = np.matrix(intl.get_adjmatrix(xcc_0,symbols_0,fconnect,"int")[0])
    amatrix_t = np.matrix(intl.get_adjmatrix(xcc_t,symbols_t,fconnect,"int")[0])

    # autoconnect!
    amatrix_0 = np.matrix(intl.link_fragments(xcc_0,amatrix_0.tolist(),1)[0])
    amatrix_t = np.matrix(intl.link_fragments(xcc_t,amatrix_t.tolist(),1)[0])

    # correlate amatrices
    dcorr = correlate_amatrices(amatrix_t,amatrix_0,symbols_t,symbols_0,xcc_t,xcc_0)

    if pp:
       print("forced: %i"%nforce)
       assign_print(dcorr)

    return dcorr
#---------------------------------------------------#
def reorder_xcc(symbols_old,xcc_old,dcorr):
    xcc     = [0.0  for ii in range(len(xcc_old    ))]
    symbols = [None for ii in range(len(symbols_old))]
    for node_t,node_0 in dcorr.items():
        xcc[3*node_0:3*node_0+3] = xcc_old[3*node_t:3*node_t+3]
        symbols[node_0]          = symbols_old[node_t]
    return symbols,xcc
#===================================================#



#===================================================#
def correlate(symbols_t,xcc_t,symbols_0,xcc_0,fconnect=1.3):
    '''
    reference geometry: symbols_0,xcc_0
    target    geometry: symbols_t,xcc_t
    '''
    # prepare reference geometry
    masses_0 = fncs.symbols2masses(symbols_0)
    xcc_0    = center_and_orient(xcc_0,None,None,masses_0)[0]
    # prepare target geometry
    masses_t = fncs.symbols2masses(symbols_t)
    xcc_t    = center_and_orient(xcc_t,None,None,masses_t)[0]
    # correlate enantiomer
    dcorr = correlate_xccs(xcc_t,symbols_t,xcc_0,symbols_0,fconnect=fconnect,pp=False)
    # non-assigned?
    if get_num_of_not_assigned(dcorr) != 0: return None
    # Re-order xcc_t
    symbols,xcc = reorder_xcc(symbols_t,xcc_t,dcorr)
    return symbols,xcc
#===================================================#




