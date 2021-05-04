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
| Module     :  common             |
| Sub-module :  MolGraph           |
| Last Update:  2020/05/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the MolGraph class
'''

#=============================================#
import numpy as np
#---------------------------------------------#
import common.fncs     as     fncs
import common.internal as     intl
from   common.Ugraph   import UGRAPH
#=============================================#


class MolGraph(UGRAPH):

      def __init__(self,xcc,symbols,cscal=1.3,nfrags=1,epslin=4.5):
          # initialize data for UGRAPH
          self._ugdict = {}
          self._nnodes = 0
          self._nedges = 0
          self._cnumber= 0
          # calculate connectivity matrix
          cmatrix = intl.get_adjmatrix(xcc,symbols,cscal,mode="int")[0]
          # autoconnect
          cmatrix = intl.link_fragments(xcc,cmatrix,nfrags)[0]
          # save cmatrix
          cmatrix = np.matrix(cmatrix)
          # set graph
          self.set_from_amatrix(cmatrix)
          # detect atoms in cycles
          nodes_cycles = self.nodes_in_cycles()
          # save data
          self._xcc     = [xi for xi in xcc]
          self._symbols = symbols
          self._incycle = nodes_cycles
          self._cmatrix = cmatrix
          self._epslin  = epslin
          self._cscal   = cscal
          self._angles  = {}

          # see if dummy is required by checking all atoms with two connections (B-A-C)
          dummies = []
          nodeX   = int(self._nnodes) - 1
          for nodeA,neighbors in self._ugdict.items():
              if len(neighbors) != 2: continue
              nodeB,nodeC = neighbors
              # add dummy atom??
              if self.islinear( (nodeB,nodeA,nodeC) ):
                 xB = self._xcc[3*nodeB:3*nodeB+3]
                 xA = self._xcc[3*nodeA:3*nodeA+3]
                 # coordinates of dummy atom
                 xX = intl.zmat_thirdatom(xB,xA,dist=1.0,angle=np.pi/2)
                 # save data
                 nodeX += 1
                 dummies.append( (nodeA,nodeX,xX) )

          #---------------------#
          # Include dummy atoms #
          #---------------------#
          nn = int(self._nnodes)
          nd = len(dummies)
          if nd > 0:
             # increase cmatrix
             zerocols      = np.zeros( (nn,nd   ) , dtype=int)
             zerorows      = np.zeros( (nd,nn+nd) , dtype=int)
             self._cmatrix = np.hstack( (self._cmatrix,zerocols) )
             self._cmatrix = np.vstack( (self._cmatrix,zerorows) )
             # Add dummy atoms!!
             for nodeA,nodeX,xX in dummies:
                 # add connection in graph
                 self.add_node(nodeX)
                 self.add_edge(nodeX,nodeA)
                 # add connection in cmatrix
                 self._cmatrix[nodeA,nodeX] = 1
                 self._cmatrix[nodeX,nodeA] = 1
                 # add dummy to symbols and xcc
                 self._symbols.append("X")
                 self._xcc += [xi for xi in xX]

      def __str__(self):
          return "(n,e)=(%i,%i)"%(self._nnodes,self._nedges)

      def islinear(self,triad):
          if triad not in self._angles:
             nodeA,nodeB,nodeC = triad
             xA = self._xcc[3*nodeA:3*nodeA+3]
             xB = self._xcc[3*nodeB:3*nodeB+3]
             xC = self._xcc[3*nodeC:3*nodeC+3]
             angABC = fncs.angle(xA,xB,xC)
             # save angles
             self._angles[(nodeA,nodeB,nodeC)] = angABC
             self._angles[(nodeC,nodeB,nodeA)] = angABC
          else: angABC = self._angles[triad]
          # return boolean (True if linear)
          return abs(180-np.rad2deg(angABC)) < self._epslin

      def nonlinearpath(self,start_idx,visited=[],ppoint=None):
          '''
          Get longest path without three connected atoms in a line
          To do so, it always choose the dummy atom (if not visited yet)
          '''
          # Get neighbors, excluding previously visited ones
          neighbors = [node for node in self._ugdict[start_idx] if node not in visited]

          if len(neighbors) == 0: return [start_idx]

          # see if any neighbor is X (dummy atom), it must be the end of path
          symbol_neighbors = [self._symbols[node] for node in neighbors]
          if symbol_neighbors.count("X") != 0:
             nodeX = neighbors[symbol_neighbors.index("X")]
             return [start_idx,nodeX]


          # Get longest from non-visited neighbors
          length = - float("inf")
          for neighbor in neighbors:
              # assert path does not contain linear angles
              if ppoint is not None and self.islinear( (ppoint,start_idx,neighbor) ): continue
              # now we now it is not 180, we can continue
              visited_i = visited + [start_idx,neighbor]
              path_i    = self.nonlinearpath(neighbor,visited=visited_i,ppoint=start_idx)
              symbol_path = [self._symbols[node] for node in path_i]
              ppath = [start_idx]+path_i
              if len(path_i) > length:
                 length = len(path_i)
                 the_path = path_i
          return [start_idx] + the_path

      def construct_zmatrix(self,seed=None,visited=set([]),previous=None):
          '''
          Returns $zmatrix and $visited
             * $zmatrix --> a list of tuples
             * $visited --> a list of visited nodes

          Tuples in $zmatrix contains the following elements
                 (at1,at2,at3,at4,boolean1,boolean2)
          where:
             * at1 to at4 are node indices
             * boolean1 is True is the tuple corresponds to a PROPER torsion

          '''
          #---------------------------------------------------------------------------------#
          # For some reason, Python messes when the following line is not included          #
          # If ommited and construct_zmatrix is inside a loop to analize different graphs,  #
          # the $visited variable DO NOT initialize when calling it as:                     #
          #          construct_zmatrix(graph)                                               #
          # I guess this is a Python bug or something....                                   #
          #---------------------------------------------------------------------------------#
          visited = set(visited) # for some reason, this is required
          #---------------------------------------------------------------------------------#

          zmatrix = []
          # find longest path
          if seed is None:
             for node in range(self._nnodes):
                 if len(self.neighbors(node)) == 1: break
             # DFS to find one end point of longest path
             path = self.nonlinearpath(node)
             # DFS to find the actual longest path
             path = self.nonlinearpath(path[-1])
          else:
             path = self.nonlinearpath(seed,list(visited))

    ###########################
    #     for idx in range(1,len(path)-1):
    #         triad = (idx-1,idx,idx+1)
    #         if self.islinear( triad ): print("*linear:",triad)
    ###########################

          tovisit = set([])
          links = {}
          zmatrix2 = []
          added_as_itor = set([])
          for idx,atom in enumerate(path):
              if   idx == 0: atoms = ( atom,   -1       ,  -1        ,  -1         ,False)
              elif idx == 1: atoms = ( atom, path[0]    ,  -1        ,  -1         ,False)
              elif idx == 2: atoms = ( atom, path[1]    , path[0]    ,  -1         ,False)
              else         : atoms = ( atom, path[idx-1], path[idx-2], path[idx-3] ,True )
              zmatrix.append(atoms)
              visited.add(atom)
              # link rest of bonded atoms to path
              for vecino in self.neighbors(atom):
                  if vecino in path         : continue
                  if vecino in added_as_itor: continue
                  if vecino in visited      : continue
                  if idx == 0:
                     # select node in previous such as angle != 180
                     for prev_i in previous:
                         if not self.islinear( (vecino,atom,prev_i) ): break
                     atoms = (vecino,atom,prev_i     ,path[idx+1],False)
                  else:
                     if self.islinear( (vecino,atom,path[idx-1]) ):
                        atoms = (vecino,atom,path[idx+1],path[idx-1],False)
                     else:
                        atoms = (vecino,atom,path[idx-1],path[idx+1],False)
                  zmatrix2.append( atoms )
                  tovisit.add(vecino)
                  added_as_itor.add(vecino)
                  # save link (if X, save X)
                  link1 = path[idx-1]
                  link2 = path[idx-1]
                  if idx+1 != len(path): link2 = path[idx+1]
                  if self._symbols[link2] == "X": links[vecino] = (atom,link2,link1)
                  else                          : links[vecino] = (atom,link1,link2)

          zmatrix += zmatrix2
          # now work with neighbors
          while len(tovisit) != 0:
              node3 = tovisit.pop()
              # update visited
              visited.add(node3)
              # recurrence!!
              zmatrix_i, visited_i = self.construct_zmatrix(node3,visited,links[node3])
              # update visited and tovisit
              visited = visited.union(visited_i)
              tovisit = tovisit.difference(visited)
              # update zmatrix with torsions in ramifications
              for torsion in zmatrix_i:
                  numnegs = torsion.count(-1)
                  if   numnegs == 3: continue
                  # requires connection to main path
                  elif numnegs == 2:
                       atom,link1,link2 = links[node3]
                       triad1 = (node3,atom,link1)
                       triad2 = (node3,atom,link2)
                       if self.islinear( triad1 ):
                          torsion = [torsion[0],node3,atom,link2,True]
                       else:
                          torsion = [torsion[0],node3,atom,link1,True]
                  # requires connection to main path
                  elif numnegs == 1: torsion = torsion[0:3] + links[node3][0:1]+(True,)
                  zmatrix.append(torsion)
          return zmatrix, visited


