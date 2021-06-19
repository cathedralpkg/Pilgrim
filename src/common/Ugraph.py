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
| Sub-module :  Ugraph             |
| Last Update:  2020/05/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the Ugraph class
'''

#=============================================#
import numpy as np
#=============================================#


#>>>>>>>>>>>>>>>>>>*
# CLASS: Queue     *
#>>>>>>>>>>>>>>>>>>*
class Queue:
    """
    A simple implementation of a FIFO queue.
    """
    def __init__(self):
        self._items = []
    def __len__(self):
        return len(self._items)
    def __iter__(self):
        for item in self._items:
            yield item
    def __str__(self):
        return str(self._items)
    def enqueue(self, item):
        self._items.append(item)
    def dequeue(self):
        return self._items.pop(0)
    def clear(self):
        self._items = []


#>>>>>>>>>>>>>>>>>>*
# CLASS: Stack     *
#>>>>>>>>>>>>>>>>>>*
class Stack:
    """
    A simple implementation of a LIFO stack
    """
    def __init__(self):
        self._items = []
    def __len__(self):
        return len(self._items)
    def __iter__(self):
        for item in self._items:
            yield item
    def __str__(self):
        return str(self._items)
    def push(self, item):
        self._items = [item] + self._items
    def pop(self):
        return self._items.pop(0)
    def clear(self):
        self._items = []


#>>>>>>>>>>>>>>>>>>*
# CLASS: UGRAPH    *
#>>>>>>>>>>>>>>>>>>*
class UGRAPH:
      """
      A simple implementation of a undirected graph
      """

      def __init__(self):
          self._ugdict = {}
          self._nnodes = 0
          self._nedges = 0
          self._cnumber= 0 # cycle number
          self._lapla  = None

      def __str__(self):
          return "(n,e)=(%i,%i)"%(self._nnodes,self._nedges)

      #-----------------#
      # Add/Remove node #
      #-----------------#
      def add_node(self,node):
          if node not in self._ugdict.keys():
             self._nnodes += 1
             self._ugdict[node] = set([])

      def remove_node(self,node1):
          # Remove node
          self._ugdict.pop(node1)
          self._nnodes -= 1
          # Remove edges with that node
          for node2 in self._ugdict.keys():
              self._ugdict[node2].discard(node1)
      #-----------------#
          
      #-----------------#
      # Add/Remove edge #
      #-----------------#
      def add_edge(self,node1,node2):
          self.add_node(node1)
          self.add_node(node2)
          if node2 not in self._ugdict[node1]:
             self._ugdict[node1].add(node2)
             self._ugdict[node2].add(node1)
             self._nedges += 1

      def remove_edge(self,node1,node2):
          self._ugdict[node1].discard(node2)
          self._ugdict[node2].discard(node1)
          self._nedges -= 1
      #-----------------#

      def set_from_amatrix(self,amatrix):
          '''
          set graph from adjacency matrix
          '''
          nn = len(amatrix)
          for node1 in range(nn):
              self.add_node(node1)
              for node2 in range(node1+1,nn):
                  if amatrix[node1,node2] in [True,1]:
                     self.add_edge(node1,node2)

      def get_amatrix(self):
          # create adj matrix
          amatrix = np.zeros((self._nnodes,self._nnodes),dtype=int)
          # complete it
          for node1,neighbors in self._ugdict.items():
              for node2 in neighbors:
                  amatrix[node1,node2] = 1
                  amatrix[node2,node1] = 1
          return amatrix

      #-------------------------#
      # get different variables #
      #-------------------------#
      def get_nnodes(self):
          '''
          Returns number of nodes in the ugraph
          '''
          return self._nnodes

      def get_nedges(self):
          '''
          Returns number of edges in the ugraph
          '''
          return self._nedges

      def get_nodes(self):
          '''
          Returns the nodes in the ugraph
          '''
          return list(self._ugdict.keys())

      def get_edges(self):
          '''
          Returns the edges in the ugraph
          '''
          edges = set([])
          for node1 in self._ugdict.keys():
              for node2 in self._ugdict[node1]:
                  edge = tuple(sorted((node1,node2)))
                  edges.add(edge)
          return edges
      #-------------------------#


      def neighbors(self,node):
          return self._ugdict[node].copy()


      #-------------#
      # BFS and DFS #
      #-------------#
      def bfsearch(self,start_idx):
          '''
          Breadth First Search for undirected graph
          Input:
            * graph_dict: a dict of the graph representing the
                          adjacency list
                          - key  : integer
                          - value: list of integers
            * start_idx : the index where to start the BFS
          '''
          # Initialize queue
          queue   = Queue()
          visited = [start_idx]
          queue.enqueue(start_idx)
          # Start BFS
          while len(queue) != 0:
               # Take node out of queue
               target_idx = queue.dequeue()
               # Get neighbors
               neighbors  = self._ugdict[target_idx]
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   queue.enqueue(neighbor)
          return visited

      def dfsearch(self,start_idx):
          '''
          Depth First Search
          Breadth First Search for undirected graph
          Input:
            * graph_dict: a dict of the graph representing the
                          adjacency list
                          - key  : integer
                          - value: list of integers
            * start_idx : the index where to start the BFS
          '''
          # Initialize queue
          stack   = Stack()
          visited = [start_idx]
          stack.push(start_idx)
          # Start BFS
          while len(stack) != 0:
               # Take node out of queue
               target_idx = stack.pop()
               # Get neighbors
               neighbors  = self._ugdict[target_idx]
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   stack.push(neighbor)
          return visited
      #-------------#

      def get_fragments(self):
          fragments     = []
          nodes         = list(self._ugdict.keys())
          visited_nodes = set([])
          for node in nodes:
              if node in visited_nodes: continue
              # explore graph with BFS from node
              fragment      = self.bfsearch(node)
              # update visited nodes
              visited_nodes = visited_nodes.union(fragment)
              # store fragment
              fragments.append(fragment)
          return fragments

      def bfsearch1d(self,idx1,idx2):
          '''
          Using a BFS algorithm, goes through
          the graph.
          However, it does it in the idx1-->idx2
          directions.
          '''
          # Initialize queue
          queue      = Queue()
          neighbors1 = self._ugdict[idx1]
          old2       = None

          # idx1 and idx2 are not bonded, there is a node in the middle (idx1--idxJ--idx2)
          if idx2 not in neighbors1:
             neighbors2 = self._ugdict[idx2]
             idxJ = list(neighbors1.intersection(neighbors2))
             if idxJ == []:
                return None
             else:
                old2 = idx2
                idx2 = idxJ[0]

          visited = [idx2]
          queue.enqueue(idx2)
          # Start BFS
          while len(queue) != 0:
               # Take node out of queue
               target_idx = queue.dequeue()
               # Get neighbors
               neighbors  = list(self._ugdict[target_idx])
               if target_idx == idx2:
                  neighbors.remove(idx1)
               # Visit neighbors
               for neighbor in neighbors:
                   if neighbor in visited: continue
                   visited.append(neighbor)
                   queue.enqueue(neighbor)
          visited.remove(idx2)
          if old2 is not None: visited.remove(old2)
          return visited

      def dfs_cycle(self,node2,node1,color,mark,par):
          # node2 was completely visited
          if color[node2] == 2: return

          # node2 was visited once. Now, backtrack based on parents to find complete cycle
          if color[node2] == 1:
              self._cnumber += 1
              cur = node1
              mark[cur] = self._cnumber

              # backtrack the vertex which are
              # in the current cycle thats found
              while cur != node2:
                  cur = par[cur]
                  mark[cur] = self._cnumber

              return

          par[node2] = node1

          # partially visited.
          color[node2] = 1

          # simple dfs on graph
          for node3 in self.neighbors(node2):

              # if it has not been visited previously
              if node3 == par[node2]:
                  continue
              self.dfs_cycle(node3, node2, color, mark, par)

          # completely visited.
          color[node2] = 2

      def nodes_in_cycles(self):
          '''
          graph coloring method --> uses dfs_cycle
          Good when the graph is not very dense in edges; in such case, use v2
          For molecules, v1 is the best option
          '''
          if self._nnodes - self._nedges == 1: return []

          color = [0] * self._nnodes
          par   = [0] * self._nnodes # parent of node
          mark  = [0] * self._nnodes
          self.dfs_cycle(0,-1,color=color,mark=mark,par=par)
          nodes_in_cycles = [node_i for node_i,mark_i in enumerate(mark) if mark_i != 0]
          return nodes_in_cycles

      def remove_external(self):
          '''
          remove external nodes iteratively
          If there are no cycles, it should returns an empty list
          '''
          cmatrix = self.get_amatrix()
          # remove atoms with connectivity smaller than 2 until convergence
          while True:
              # check which nodes have only 1 neighbor (or none)
              toremove = []
              for node in range(self._nnodes):
                  if sum(cmatrix[node,:]) == 1: toremove.append(node)
              # nothing to remove --> so no cycle
              if len(toremove) == 0: break
              # remove them
              for node in toremove:
                  cmatrix[node,:] = 0
                  cmatrix[:,node] = 0
          # see nodes
          remaining_nodes = [node for node in range(self._nnodes) if sum(cmatrix[node,:] != 0)]
          return remaining_nodes

      def longest_path_from_node(self,start_idx,visited=[]):
          '''
          Naive algorithm to explore the graph, starting at start_idx,
          and return the longest path
          '''
          # Get neighbors, excluding previously visited ones
          neighbors = [node for node in self._ugdict[start_idx] if node not in visited]

          if len(neighbors) == 0: return [start_idx]

          # Get longest from non-visited neighbors
          length = - float("inf")
          for neighbor in neighbors:
              visited_i = visited + [start_idx,neighbor]
              path_i    = self.longest_path_from_node(neighbor,visited=visited_i)
              if len(path_i) > length:
                 length = len(path_i)
                 the_path = path_i
          return [start_idx] + the_path

      def longest_path(self):
          # DFS to find one end point of longest path
          lnode = self.longest_path_from_node(0)[-1]
          # DFS to find the actual longest path
          path = self.longest_path_from_node(lnode)
          return path

      def get_layers(self,center):
          '''
           returns a list of layers for the node center
              * 1st layer: neighbors of node center
              * 2nd layer: neighbors of neighbors of center
                           (excluding repetitions of previous layers)
          '''
          layers  = [set([center])]
          current = [center]
          visited = set([center])
          nnodes  = len(self._ugdict.keys())

          while len(visited) != nnodes:
                layer = []
                for node in current:
                    neighbors = self._ugdict[node]
                    layer     = layer + list(neighbors)
                layer = set(layer).difference(visited)
                visited = visited.union(layer)
                layers.append(layer)
                current = list(layer)
          return layers
                
      #----------------------------#
      # Get matrix representations #
      #----------------------------#
      def gen_laplacian(self):
          self._lapla = np.zeros((self._nnodes,self._nnodes))
          for node in self._ugdict.keys():
              neighbors = self._ugdict[node]
              for neighbor in neighbors:
                  self._lapla[node,node] = self._lapla[node,node] + 1
                  self._lapla[node,neighbor] = -1

          # Eigenvalues
          vals, vecs = np.linalg.eigh(self._lapla)

          # Degenerancies?
          degs = [0]*len(vals)
          for i in range(len(vals)):
              val_i = vals[i]
              for j in range(len(vals)):
                 val_j = vals[j]
                 if abs(val_i-val_j) < 1e-3: degs[i] += 1

          # Data for each node
          dict_vecs = {}
          for node in self._ugdict.keys():
              vector = [ abs(float(vecs[node,idx])) for idx in range(len(vals)) if degs[idx] == 1]
              dict_vecs[node] = vector
          return dict_vecs


if __name__ == "__main__":
    graph = UGRAPH()
    graph.add_edge(1-1 ,  2-1)
    graph.add_edge(2-1 ,  3-1)
    graph.add_edge(3-1 ,  4-1)
    graph.add_edge(4-1 ,  6-1)
    graph.add_edge(4-1 ,  7-1)
    graph.add_edge(3-1 ,  5-1)
    graph.add_edge(7-1 ,  8-1)
    graph.add_edge(6-1 , 10-1)
    graph.add_edge(9-1 , 10-1)
    graph.add_edge(5-1 ,  9-1)
    graph.add_edge(11-1, 12-1)
    graph.add_edge(11-1, 13-1)
    graph.add_edge(5-1 , 11-1)
    graph.add_edge(6-1 , 13-1)
    print([at+1 for at in graph.nodes_in_cycles()])
