'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.5
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
| Sub-module :  enantor            |
| Last Update:  2021/05/20 (Y/M/D) |
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
def get_not_assigned(dcorr):
    return sum([1 if type(nodes1) != int else 0 for node2,nodes1 in dcorr.items()])
#---------------------------------------------------#
def assign_check(dcorr):
    assigned = []
    for node2,node1 in dcorr.items():
        # node already assigned
        if   type(node1) == int:
            if node1 in assigned: raise Exception
            assigned.append(node1)
            continue
        # remove already assigned
        node1 = node1.difference(assigned)
        dcorr[node2] = node1
        # check length of node1
        if   len(node1 ) == 0  : raise Exception
        elif len(node1 ) == 1  : dcorr[node2] = list(node1)[0]
    return dcorr,get_not_assigned(dcorr) 
#---------------------------------------------------#
def assign_concatenate(dcorr,graph1,graph2):
    '''
    It tries to assign neighbors of an assigned node using connectivity.
    This process is executed until no new assignations.
    '''
    dcorr,na = assign_check(dcorr)
    while True:
        na_i = get_not_assigned(dcorr)
        dcorr2 = {k:v for k,v in dcorr.items()}
        for node2,node1 in dcorr.items():
            # node2 not assigned
            if type(node1) != int: continue
            # node2 assigned --> check not assigned neighbors
            neighbors2 = [neigh for neigh in graph2.neighbors(node2) if type(dcorr[neigh]) != int]
            # correspoding neighbors in graph1
            neighbors1 = graph1.neighbors(node1)
            for neigh2 in neighbors2:
                dcorr2[neigh2] = dcorr2[neigh2].intersection(neighbors1)
        # check
        dcorr,na_j = assign_check(dcorr2)
        # return
        if na_i == na_j or na_j == 0: return dcorr,na_j
#===================================================#



#===================================================#
def assign_initialize(symbols):
    dcorr = {idx2:set([idx1 for idx1,s1 in enumerate(symbols) if s1 == s2]) \
                            for idx2,s2 in enumerate(symbols)}
    return assign_check(dcorr)
#---------------------------------------------------#
def assign_neighbors(dcorr,graph1,graph2,symbols):
    # number of non-assigned (initial)
    na_i = get_not_assigned(dcorr)
    # Compare symbols for neighbors
    for node2,nodes1 in list(dcorr.items()):
        # already assigned?
        if type(nodes1) == int: continue
        # get symbol of neighbors
        neighbors2 = sorted([symbols[idx] for idx in graph2.neighbors(node2)])
        # compare with graph1
        dcorr[node2] = set([])
        for node1 in nodes1:
            # get symbol of neighbors
            neighbors1 = sorted([symbols[idx] for idx in graph1.neighbors(node1)])
            # compare neighbors
            if neighbors1 == neighbors2: dcorr[node2].add(node1)
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph1,graph2)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_layers(dcorr,graph1,graph2,symbols):
    '''
    if this is executed, assign_neighbors is not required,
    as the second layer (layer[1]) corresponds to neighors
    '''
    # number of non-assigned (initial)
    na_i = get_not_assigned(dcorr)
    # Compare symbols for each layer
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int: continue
        # layers node 2
        layers2 =  [sorted([symbols[at] for at in l]) for l in graph2.get_layers(node2)]
        # comparison graph1 <--> graph2
        dcorr[node2] = set([])
        for node1 in nodes1:
            # layers node 1
            layers1 = [sorted([symbols[at] for at in l]) for l in graph1.get_layers(node1)]
            if layers1 == layers2: dcorr[node2].add(node1)
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph1,graph2)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_spatial(dcorr,graph1,graph2,xcc1,xcc2,symbols):
    # number of non-assigned (initial)
    na_i = get_not_assigned(dcorr)
    # Compare spatial organization (improper torsions)
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int: continue
        BC = None
        # check neighbors
        atB2,atC2,atD2 = None,None,None
        for neigh in graph2.neighbors(node2):
            if type(dcorr[neigh]) != int: continue
            atB2 = neigh
            atC2,atD2 = None,None
            for nneigh in graph2.neighbors(atB2):
               #print(node2+1,neigh+1,nneigh+1)
                if nneigh == node2: continue
                if type(dcorr[nneigh]) != int: continue
                if atC2 is None: atC2 = nneigh; continue
                if atD2 is None: atD2 = nneigh; break
            if atD2 is not None: break
            if atC2 is not None: BC = (atB2,atC2)
        # now check spatial configuration
        if atD2 is not None:
           atoms2  = [node2,atB2,atC2,atD2]
           xs2     = (xcc2[3*at:3*at+3] for at in atoms2)
           angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
           config2 = angle2 < 180
           dcorr[node2] = set([])
          #print("*",node2+1,[i+1 for i in atoms2],angle2)
           for node1 in nodes1:
               atoms1  = [node1,dcorr[atB2],dcorr[atC2],dcorr[atD2]]
               xs1     = (xcc1[3*at:3*at+3] for at in atoms1)
               angle1  = np.rad2deg(fncs.dihedral(*xs1))%360
               config1 = angle1 < 180
              #print("*",node1+1,[i+1 for i in atoms1],angle1)
               if config1 == config2:
                  dcorr[node2].add(node1)
          #print()
        elif BC is not None and len(nodes1) == 2:
           atB2,atC2 = BC
           # look for the other node with same assigments in graph1
           node2a = node2
           node2b = None
           for node2_,nodes1_ in dcorr.items():
               if node2_ == node2a : continue
               if nodes1 != nodes1_: continue
               node2b = node2_
               break
           if node2b is None: continue
           # assert both nodes are connected to the same atoms!!
           if graph2.neighbors(node2a) != graph2.neighbors(node2b): continue
           # config in graph2
           atoms2  = [node2a,node2b,atB2,atC2]
           xs2     = (xcc2[3*at:3*at+3] for at in atoms2)
           angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
           config2 = angle2 < 180
           # config in graph1
           node1a,node1b = nodes1
           atoms1 = [node1a,node1b,dcorr[atB2],dcorr[atC2]]
           xs1     = (xcc1[3*at:3*at+3] for at in atoms1)
           angle1  = np.rad2deg(fncs.dihedral(*xs1))%360
           config1 = angle1 < 180
           if config1 == config2:
              dcorr[node2a] = node1a
              dcorr[node2b] = node1b
           else:
              dcorr[node2a] = node1b
              dcorr[node2b] = node1a
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph1,graph2)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_torsions(dcorr,graph1,graph2,xcc1,xcc2,symbols):
    '''
    only one at a time
    '''
    # number of non-assigned (initial)
    na_i = get_not_assigned(dcorr)
    # assignment based in proper torsions
    minad = float("inf")
    minnode1 = None
    minnode2 = None
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int: continue
        BC = None
        # check neighbors
        atA2,atB2,atC2,atD2 = node2,None,None,None
        for atB2 in graph2.neighbors(atA2):
            if type(dcorr[atB2]) != int: continue
            for atC2 in graph2.neighbors(atB2):
                if type(dcorr[atC2]) != int: continue
                for neigh in graph2.neighbors(atC2):
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
        xs2     = (xcc2[3*at:3*at+3] for at in atoms2)
        angle2  = np.rad2deg(fncs.dihedral(*xs2))%360
        config2 = angle2 < 180
        # compare torsions
        for node1 in nodes1:
            atoms1 = [node1,dcorr[atB2],dcorr[atC2],dcorr[atD2]]
            xs1    = (xcc1[3*at:3*at+3] for at in atoms1)
            angle1 = np.rad2deg(fncs.dihedral(*xs1))
            # difference
            ad = fncs.angular_dist(-angle2,angle1,u="deg")
            if ad < minad:
                minad    = ad
                minnode1 = node1
                minnode2 = node2
    # Any new assignation to do??
    if minnode1 is None: return dcorr, na_i
    # assign
    dcorr[minnode2] = minnode1
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph1,graph2)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_cx3(dcorr,graph1,graph2,xcc1,xcc2,symbols):
    # number of non-assigned (initial)
    na_i = get_not_assigned(dcorr)
    # Assign in CX3 groups (CH3, CF3, CCl3, et cetera)
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int: continue
        if  len(nodes1) !=   3: continue
        neighbors2 = graph2.neighbors(node2)
        atC = neighbors2.pop()
        # check we have a C atom
        if len(neighbors2) != 0: continue
        if symbols[atC] != "C": continue
        # determine torsion X-C-B-A
        atA = None
       #print("*",node2+1,atC+1)
        for atB in graph2.neighbors(atC):
            if type(dcorr[atB]) != int: continue
            for neigh in graph2.neighbors(atB):
                if type(dcorr[neigh]) != int: continue
                if neigh == node2: continue
                if neigh == atC  : continue
                atA = neigh
               #print(" ",atB+1,atA+1)
                break
            if atA is not None: break
        if atA is None: continue
        # calculate dihedral
        atoms2  = [node2,atC,atB,atA]
       #print([i+1 for i in atoms2])
        xs2     = (xcc2[3*at:3*at+3] for at in atoms2)
        angle2  = np.rad2deg(fncs.dihedral(*xs2))
        # compare with the others
        minad = float("inf")
        minnode = None
        for node1 in nodes1:
            atoms1 = [node1,dcorr[atC],dcorr[atB],dcorr[atA]]
            xs1    = (xcc1[3*at:3*at+3] for at in atoms1)
            angle1 = np.rad2deg(fncs.dihedral(*xs1))
            # difference
            ad = fncs.angular_dist(-angle2,angle1,u="deg")
            if ad < minad:
                minad   = ad
                minnode = node1
        # apply
        dcorr[node2] = minnode
    # number of non-assigned(final)
    dcorr,na_j = assign_check(dcorr)
    # return data
    if na_j != 0 and na_j != na_i: return assign_concatenate(dcorr,graph1,graph2)
    else                         : return assign_check(dcorr)
#---------------------------------------------------#
def assign_force(dcorr,graph1,graph2):
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
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int: continue
        nneighbors = len(graph2.neighbors(node2))
        # look for minimum doubt
        if   nneighbors == 1 and len(nodes1) < tmindoubt:
            tmindoubt = len(nodes1)
            ttarget1,ttarget2 = sorted(list(nodes1))[0], node2
        elif nneighbors != 1 and len(nodes1) < mindoubt:
            mindoubt = len(nodes1)
            target1,target2 = sorted(list(nodes1))[0], node2
    # apply target
    if mindoubt <= tmindoubt:
       if   target2  is not None: dcorr[target2 ] = target1
       elif ttarget2 is not None: dcorr[ttarget2] = ttarget1
    else:
       if   ttarget2 is not None: dcorr[ttarget2] = ttarget1
       elif target2  is not None: dcorr[target2 ] = target1
    # return
    return assign_concatenate(dcorr,graph1,graph2)
#---------------------------------------------------#
def assign_print(dcorr):
    for node2,nodes1 in dcorr.items():
        if type(nodes1) == int:
            print("%2i : %2i"%(node2+1,nodes1+1))
        else:
            print("%2i :"%(node2+1),[node1+1 for node1 in nodes1])
    print()
#---------------------------------------------------#
def equivalent_atoms_in_graphs(amatrix1,amatrix2,symbols):
    # Graphs
    graph1 = UGRAPH(); graph1.set_from_amatrix(amatrix1)
    graph2 = UGRAPH(); graph2.set_from_amatrix(amatrix2)
    
    # Correlations
    dcorr,na_0 = assign_initialize(symbols)
    for ii in range(10): # it should be a while True but... just in case
        # use symbols of neighbors
        dcorr,na_1  = assign_neighbors(dcorr,graph1,graph2,symbols)
        if na_1 == 0: break
        # use symbols of layers
        dcorr,na_2  = assign_layers(dcorr,graph1,graph2,symbols)
        if na_2 == 0: break
        # update na_0
        na_0 = na_2

    return dcorr
#---------------------------------------------------#
def equivalent_atoms(xcc,symbols,fconnect=1.3):
    xcc1 = [xi for xi in xcc]
    xcc2 = [xi for xi in xcc]

    # Adjacency matrices
    amatrix1 = np.matrix(intl.get_adjmatrix(xcc1,symbols,fconnect,"int")[0])
    amatrix2 = np.matrix(intl.get_adjmatrix(xcc2,symbols,fconnect,"int")[0])

    # autoconnect!
    amatrix1 = np.matrix(intl.link_fragments(xcc1,amatrix1.tolist(),1)[0])
    amatrix2 = np.matrix(intl.link_fragments(xcc2,amatrix2.tolist(),1)[0])

    # Correlate
    dcorr = equivalent_atoms_in_graphs(amatrix1,amatrix2,symbols)
    return dcorr
#===================================================#



#===================================================#
def generate_enantio(xcc):
    return [-xcc[idx] if idx%3==0 else xcc[idx] for idx in range(len(xcc))]
#---------------------------------------------------#
def correlate_enantio(xcc1,xcc2,symbols,fconnect=1.3,pp=False):
    '''
    xcc1 --> geometry
    xcc2 --> enantiomer
    '''

    if len(xcc1) != len(xcc2)     : raise Exception
    if len(xcc1) != len(symbols)*3: raise Exception
    
    # Adjacency matrices
    amatrix1 = np.matrix(intl.get_adjmatrix(xcc1,symbols,fconnect,"int")[0])
    amatrix2 = np.matrix(intl.get_adjmatrix(xcc2,symbols,fconnect,"int")[0])

    # autoconnect!
    amatrix1 = np.matrix(intl.link_fragments(xcc1,amatrix1.tolist(),1)[0])
    amatrix2 = np.matrix(intl.link_fragments(xcc2,amatrix2.tolist(),1)[0])

    # Graphs
    graph1 = UGRAPH(); graph1.set_from_amatrix(amatrix1)
    graph2 = UGRAPH(); graph2.set_from_amatrix(amatrix2)

    
    # Correlations
    dcorr,na_0 = assign_initialize(symbols)
    nforce = 0
    for ii in range(10): # it should be a while True but... just in case
        # use symbols of neighbors
        dcorr,na_1  = assign_neighbors(dcorr,graph1,graph2,symbols)
        if na_1 == 0: break
        # use symbols of layers
        dcorr,na_2  = assign_layers(dcorr,graph1,graph2,symbols)
        if na_2 == 0: break
        # use spatial disposition
        dcorr,na_3  = assign_spatial(dcorr,graph1,graph2,xcc1,xcc2,symbols)
        if na_3 == 0: break
        # assign CX3 groups
        dcorr,na_4  = assign_cx3(dcorr,graph1,graph2,xcc1,xcc2,symbols)
        if na_4 == 0: break
        # assign based on torsions
        dcorr,na_5  = assign_torsions(dcorr,graph1,graph2,xcc1,xcc2,symbols)
        if na_5 == 0: break
        # force assignation
        if na_0 == na_5:
           nforce += 1
           dcorr,na_5 = assign_force(dcorr,graph1,graph2)
        # update na_0
        na_0 = na_5

    if pp:
       print("forced: %i"%nforce)
       assign_print(dcorr)

    return dcorr
#---------------------------------------------------#
def reorder_enantio(xcc_enantio,dcorr):
    xcc = [0.0 for ii in range(len(xcc_enantio))]
    for node2,node1 in dcorr.items():
        xcc[3*node1:3*node1+3] = xcc_enantio[3*node2:3*node2+3]
    return xcc
#---------------------------------------------------#
def gen_enantio_and_correlate(xcc,symbols,masses=None,fconnect=1.3,pp=False):
    # generate list of masses
    if masses is None: masses = fncs.symbols2masses(symbols)
    # xcc1: center and reorient
    xcc1 = center_and_orient(xcc,None,None,masses)[0]
    # xcc2: xcc for enantiomer
    xcc2 = generate_enantio(xcc1)
    xcc2 = center_and_orient(xcc2,None,None,masses)[0]
    # correlate enantiomer
    dcorr = correlate_enantio(xcc1,xcc2,symbols,fconnect=fconnect,pp=pp)
    # non-assigned?
    if get_not_assigned(dcorr) != 0: return None
    # Generate xcc for enantiomer
    try:
       xcc_enantio = [0.0 for ii in range(len(xcc))]
       for node2,node1 in dcorr.items():
           xcc_enantio[3*node1:3*node1+3] = xcc2[3*node2:3*node2+3]
       return xcc_enantio
    except: return None
#===================================================#



#===================================================#
if __name__ == '__main__':

   from   common.files    import read_xyz
   from   common.files    import read_xyz_zmat
   from   common.files    import write_xyz
   from   common.pgs      import get_pgs

   files = [fname for fname in os.listdir(".") if fname.endswith(".xyz")]
   print(files)
   for xyz in files:
       #if "xyz_prueba2_en.xyz" not in xyz: continue
       print(xyz)
       # read file
       try:
          xcc, symbols, masses = read_xyz(xyz)
       except:
          (lzmat,zmatvals,zmatatoms), symbols, masses = read_xyz_zmat(xyz)
          xcc = intl.zmat2xcc(lzmat,zmatvals)
       # point group
       pg, rotsigma = get_pgs(symbols,masses,xcc)
       print("   -",pg)
       # Name for molden file
       xyz2 = xyz.replace(".xyz","_en.xyz")
       # Get enantiomer
       try   : xcc_enantio = gen_enantio_and_correlate(xcc,symbols,masses,pp=True)
       except: xcc_enantio = None
       if   xcc_enantio is None : print("unable to correlate!\n")
       elif os.path.exists(xyz2): print("%s already exists!\n"%xyz2)
       else                     : write_xyz(xyz2,xcc_enantio,symbols)
#===================================================#



