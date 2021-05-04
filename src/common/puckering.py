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
| Sub-module :  puckering          |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different functions
related to puckering coordinates
'''

#=============================================#
import numpy        as np
#---------------------------------------------#
import common.fncs  as fncs
import common.dicts as dd
#=============================================#

#===============================================#
# Puckering coordinates                         #
#===============================================#
def puckering_rjs(xvec,atoms_in_ring):
    ''' JACS 1975, 97:6, 1354-1358; eq (4) '''
    gc = np.array(fncs.get_centroid(xvec,atoms_in_ring))
    list_Rj = []
    for at in atoms_in_ring:
        rj = np.array(fncs.xyz(xvec,at))
        Rj = rj-gc
        list_Rj.append(Rj)
    return list_Rj, gc
#-----------------------------------------------#
def puckering_normal(list_Rj):
    ''' JACS 1975, 97:6, 1354-1358; eq (10) '''
    N   = len(list_Rj)
    Rp  = np.array([0.0,0.0,0.0])
    Rpp = np.array([0.0,0.0,0.0])
    for j in range(N):
        Rj   = list_Rj[j]
        sin  = np.sin(2*np.pi*j/N)
        cos  = np.cos(2*np.pi*j/N)
        Rp  += sin * Rj
        Rpp += cos * Rj
    normal = np.cross(Rp,Rpp)
    return fncs.normalize_vec(normal)
#-----------------------------------------------#
def puckering_zjs(list_Rj,normal):
    ''' JACS 1975, 97:6, 1354-1358; eq (11) '''
    list_zj = []
    for Rj in list_Rj:
        zj = np.dot(Rj,normal)
        list_zj.append(zj)
    return list_zj
#-----------------------------------------------#
def puckering_coords(xvec,atoms_in_ring=None):
    '''
    Generalized Ring-Puckering Coordinates
    JACS 1975, 97:6, 1354-1358; eqs (12)-(14)
    '''
    nat = fncs.howmanyatoms(xvec)
    if atoms_in_ring is None: atoms_in_ring = range(nat)
    list_Rj, gc = puckering_rjs(xvec,atoms_in_ring)
    normal      = puckering_normal(list_Rj)
    list_zj     = puckering_zjs(list_Rj,normal)
    N = len(list_zj)
    if N%2 != 0: mmax = (N-1)//2
    else       : mmax = N//2 -1
    list_qm   = []
    list_phim = []
    for m in range(2,mmax+1,1):
       qcos = 0.0
       qsin = 0.0
       for j in range(N):
           zj = list_zj[j]
           qcos += zj*np.cos(2*np.pi*m*j/N)
           qsin += zj*np.sin(2*np.pi*m*j/N)
       qcos *= +np.sqrt(2.0/N)
       qsin *= -np.sqrt(2.0/N)
       # Get q_m and phi_m
       phi_m = fncs.sincos2angle(qsin,qcos)
       q_m   = qcos / np.cos(phi_m)
       list_qm.append(q_m)
       list_phim.append(phi_m)
    if N%2 == 0:
       qNhalf = 0.0
       for j in range(N):
           zj      = list_zj[j]
           qNhalf += ((-1)**j) * zj
       qNhalf *= np.sqrt(1.0/N)
       list_qm.append(qNhalf)
    # total puckering amplitude
    Q = np.sqrt(sum([zj**2 for zj in list_zj]))
    return list_qm, list_phim, Q
#-----------------------------------------------#
def puckering_5ring(xvec,atoms_in_ring=None):
    [q2], [phi2], Q = puckering_coords(xvec,atoms_in_ring)
    return q2, phi2
#-----------------------------------------------#
def puckering_6ring(xvec,atoms_in_ring=None):
    [q2,q3], [phi2], Q = puckering_coords(xvec,atoms_in_ring)
    # Special case for 6 member-ring
    sintheta = q2/Q
    costheta = q3/Q
    theta    = fncs.sincos2angle(sintheta,costheta)
    return Q, theta, phi2
#-----------------------------------------------#
def classify_6ring(theta,phi):
    '''
    Determines conformation of 6-member ring (Boeyens)
    '''
    # reduces dictionary
    #dicCONF6 = {k:v for k,v in dd.dicCONF6.items() if k in "25B,6H1,B25,1H6,E3"}


    p1 = (theta,phi)
    mindist = float("inf")
    for conf, p0 in dd.dicCONF6.items():
   #for conf, p0 in dicCONF6.items():
   #for conf, p0 in dd.dicCONF6_v2.items():
        distance = fncs.angdist_sphere( p0,p1 )
        if distance < mindist:
           mindist      = distance
           conformation = conf
    return conformation.strip()
#===============================================#

