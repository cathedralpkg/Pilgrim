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
| Sub-module :  fncs               |
| Last Update:  2021/05/19 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains some functions that may help to reduce
the number of lines in other modules or programs.
''' 


#==============================================#
import os
import sys
import numpy    as np
#==============================================#
import common.physcons as     pc
from   common.dicts    import dpt_z2s, dpt_z2m
from   common.dicts    import dpt_s2z, dpt_s2m
from   common.criteria import EPS_CCIC
from   common.criteria import EPS_FLOAT
from   common.criteria import EPS_SCX
from   common.criteria import EPS_GEOM
from   common.criteria import EPS_INERTIA
from   common.criteria import EPS_NORM
#==============================================#

PROJECT_TRA = True
PROJECT_ROT = True

#==============================================#
alphUC   = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"        #
alphLC   = "abcdefghijklmnopqrstuvwxyz"        #
numbers  = "0123456789"                        #
#==============================================#


#===============================================#
def exp128(arg):
    try   :
        if arg < -700 or arg > 700: return np.exp(np.float128(arg))
        return np.exp(arg)
    except: return np.exp(arg)
#---------------------------------------------#
def eformat(value,ndec):
    try:
      # value inside normal limits
      if 1E-300 < abs(value) < 1E+300 or value == 0.0:
         string = "%%.%iE"%ndec%value
      # value outside of normal limits
      else:
         string = str(value).replace("e","E")
         coef, exponent = string.split("E")
         coef = ("%%.%if"%ndec)%float(coef)
         string = "%sE%s"%(coef,exponent)
      return string
    except: return str(value)
#===============================================#

#===============================================#
# Some useful things to make better software    #
#===============================================#
def classify_args(user_args):
    '''
    classify the arguments
    '''
    dargs = {None:[]}
    # in case -h is used instead of --help
    user_args = [arg if arg != "-h" else "--help" for arg in user_args]
    # the first arguments until one with --
    for ii,arg in enumerate(user_args):
        if "--" in arg: break
        dargs[None].append(arg)
    # the rest of arguments
    for jj,arg in enumerate(user_args):
        if jj < ii: continue
        if arg.startswith("--"):
           current = arg[2:]
           dargs[current] = []
           continue
        try   : dargs[current].append(arg)
        except: pass
    return dargs
#-----------------------------------------------#
def do_parallel(parallel=False):
    yes = [True,"yes","Yes","YES","y","Y",1]
    if parallel in yes:
       try   : import multiprocessing
       except: return False
       return True
    else: return False
#-----------------------------------------------#
def set_parallel(parallel):
    '''
    use it with these two line:
        >> global PARALLEL, delayed, multiprocessing, Parallel
        >> PARALLEL, delayed, multiprocessing, Parallel = fncs.set_parallel(parallel)
    '''
    global PARALLEL
    global delayed
    global multiprocessing
    global Parallel
    PARALLEL        = False
    delayed         = None
    multiprocessing = None
    Parallel        = None
    # Parallization required
    if parallel in [True,"yes","Yes","YES","y","Y",1]:
        try:
            from   joblib import delayed, Parallel
            import multiprocessing
            PARALLEL = True
        except: PARALLEL = False
    # return
    return PARALLEL, delayed, multiprocessing, Parallel
#===============================================#


#=============================================#
# Functions related to (lists of) strings     #
#=============================================#
def add_iblank(string,nib):
    return "".join([" "*nib+line+"\n" for line in string.split("\n")])
#---------------------------------------------#
def fill_string(string,length):
    if len(string) >= length: return string
    while len(string) < length: string = " "+string+" "
    if len(string) > length: string = string[:-1]
    return string
#---------------------------------------------#
def clean_line(line,cs="#",strip=False):
    ''' cs : comment symbol '''
    line = line.split(cs)[0]
    if strip: line = line.strip()
    if not line.endswith("\n"): line += "\n"
    return line
#---------------------------------------------#
def clean_lines(lines,cs="#",strip=False):
    '''
    from a list of lines, it removes the comment part
    - cs stands for comment symbol
    '''
    return [clean_line(line,cs,strip) for line in lines]
#---------------------------------------------#
def extract_lines(lines,key1,key2,fmatch=False,cleanline=True,ignorecase=False):
    '''
    Returns the lines between the one starting with key1 and the one starting with key2
    if fmatch is True, returns the first match
    '''
    if ignorecase:
       key1 = key1.lower()
       key2 = key2.lower()
    record   = False
    selected = []
    for line in lines:
        if cleanline : line = clean_line(line,cs="#",strip=True)
        if ignorecase: line = line.lower()
        if line.startswith(key2):
           record = False
           if fmatch: break
        if record: selected.append(line)
        if line.startswith(key1):
           record   = True
           selected = []
    return selected
#---------------------------------------------#
def extract_string(lines,key1,key2,accumulate=False):
    '''
    Returns a string with lines between key1 and key2 (both included)
    '''
    strings  = []
    # Initialize variables
    record,selected  = False, []
    # Loop along lines
    for line in lines:
        # start-key found!
        if key1 in line: record = True
        # save line?
        if record: selected.append(line)
        # end-key found!
        if key2 in line and record:
           # Save data
           strings.append( "".join(selected) )
           # initialize
           record, selected = False, []
    # return data
    if accumulate: return strings
    else         : return strings[-1]
#---------------------------------------------#
def is_string_valid(string,allowed=alphUC+alphLC+numbers,extra=""):
    allowed = allowed+extra
    for character in string:
        if character not in allowed: return False
    return True
#=============================================#

#=============================================#
# Dealing with float numbers                  #
#=============================================#
def is_smaller(num1,num2,eps):
    ''' check if num1 < num2 (for floats) '''
    return num1-num2 < -abs(eps)
#---------------------------------------------#
def is_greater(num1,num2,eps):
    ''' check if num1 > num2 (for floats) '''
    return num1-num2 > +abs(eps)
#---------------------------------------------#
def is_smallereq(num1,num2,eps):
    ''' check if num1 < num2 (for floats) '''
    return num1-num2 < +abs(eps)
#---------------------------------------------#
def is_greatereq(num1,num2,eps):
    ''' check if num1 >= num2 (for floats) '''
    return num1-num2 > -abs(eps)
#=============================================#


#===============================================#
# Some basic MATHematics                        #
#===============================================#
def sign(number):
    if number >= 0.0: return +1
    else:             return -1
#-----------------------------------------------#
def delta_ij(i,j):
   if i==j: return 1.0
   else:    return 0.0
#-----------------------------------------------#
def float_in_domain(value,domain):
    '''
    domain is a list of tuples, each one with initial
    and end points
    '''
    value = float(value)
    for ival,fval in domain:
        if ival <= value <= fval: return True
    return False
#-----------------------------------------------#
def distance(x1,x2):
    ''' returns distance between two points'''
    diff = [x1_i-x2_i for x1_i,x2_i in zip(x1,x2)]
    return norm(diff)
#-----------------------------------------------#
def angle_vecs(u,v):
    ''' returns angle between two vectors'''
    cos = np.dot(u,v) / norm(u) / norm(v)
    if abs(cos-1.0) <= EPS_SCX: cos = +1.0
    if abs(cos+1.0) <= EPS_SCX: cos = -1.0
    return np.arccos(cos)
#-----------------------------------------------#
def angle(x1,x2,x3):
    '''returns angle between three points (1-2-3)'''
    u = np.array(x1)-np.array(x2)
    v = np.array(x3)-np.array(x2)
    return angle_vecs(u,v)
#-----------------------------------------------#
def dihedral(x1,x2,x3,x4):
    '''returns dihedral angle between 4 points (1-2-3-4)'''
    vec_23 = np.array(x3) - np.array(x2)
    vec_23 = normalize_vec(vec_23)
    # Compute plane vectors
    n1 = get_normal(x1,x2,x3)
    n2 = get_normal(x2,x3,x4)
    # Vector perpendicular to (n1,vec23)
    m1 = np.cross(n1,vec_23)
    # Coordinates of n2 in this frame (n1,vec23,m1)
    x = np.dot(n1,n2)
    y = np.dot(m1,n2)
    # Angle
    return -np.arctan2(y,x)
#-----------------------------------------------#
def dihedral_of_torsion(xcc,torsion_atoms):
    at1,at2,at3,at4 = torsion_atoms
    x1 = xcc[3*at1:3*at1+3]
    x2 = xcc[3*at2:3*at2+3]
    x3 = xcc[3*at3:3*at3+3]
    x4 = xcc[3*at4:3*at4+3]
    return dihedral(x1,x2,x3,x4)
#-----------------------------------------------#
def get_zmatvals_from_xcc(xcc,zmatatoms):
    zmatvals = {}
    for key,atoms in zmatatoms.items():
        # ask for a distance
        if len(atoms) == 2:
            x1 = xcc[3*atoms[0] : 3*atoms[0]+3]
            x2 = xcc[3*atoms[1] : 3*atoms[1]+3]
            value = distance(x1,x2) * pc.ANGSTROM
        # ask for an angle
        elif len(atoms) == 3:
            x1 = xcc[3*atoms[0] : 3*atoms[0]+3]
            x2 = xcc[3*atoms[1] : 3*atoms[1]+3]
            x3 = xcc[3*atoms[2] : 3*atoms[2]+3]
            value = np.rad2deg(angle(x1,x2,x3))
        # ask for a dihedral angle
        elif len(atoms) == 4:
            x1 = xcc[3*atoms[0] : 3*atoms[0]+3]
            x2 = xcc[3*atoms[1] : 3*atoms[1]+3]
            x3 = xcc[3*atoms[2] : 3*atoms[2]+3]
            x4 = xcc[3*atoms[3] : 3*atoms[3]+3]
            value = np.rad2deg(dihedral(x1,x2,x3,x4))
        # none of the previous
        else: value = None
        zmatvals[key] = value
    return zmatvals
#-----------------------------------------------#
def norm(vec):
    return np.linalg.norm(vec)
#-----------------------------------------------#
def normalize_vec(vec):
    norm_vec = norm(vec)
    if norm_vec == 0.0: return vec
    return np.array(vec)/norm_vec
#-----------------------------------------------#
def get_normal(p1,p2,p3):
    '''returns the normal vector to the plane defined by three points'''
    v12 = np.array(p2)-np.array(p1)
    v23 = np.array(p3)-np.array(p2)
    v12 = normalize_vec(v12)
    v23 = normalize_vec(v23)
    normal = np.cross(v12,v23)
    return normalize_vec(normal)
#-----------------------------------------------#
def angle_in_interval(ang,case="0,360"):
    ''' cases = '-180,180', '0,360', '-pi,pi', '0,2pi' '''
    if "pi" in case: ang = ang%pc.TWOPI
    else           : ang = ang%360.0
    if case == "-pi,pi"   and ang > pc.PI  : ang = ang-pc.TWOPI
    if case == "-180,180" and ang > 180.: ang = ang-360.0
    return ang
#-----------------------------------------------#
def angular_dist(ang1,ang2,u="rad",limit=None):
    '''
    if 'limit' is defined, 'u' lacks of sense
    '''
    if limit is None:
       if   u == "rad": limit = pc.TWOPI
       elif u == "deg": limit = 360.0
       else           : return None
    # Get distances
    diff1 = (ang1-ang2)%limit
    diff2 = (ang2-ang1)%limit
    # Return distance
    return min(diff1,diff2)
#-----------------------------------------------#
def angular_dist_with_sign(angle1,angle2):
    '''
    angle1, angle2 in [0,2pi]
    '''
    d1 = angle2-2.0*np.pi-angle1
    d2 = angle2          -angle1
    d3 = angle2+2.0*np.pi-angle1
    if abs(d1) <= abs(d2) and abs(d1) <= abs(d3): return d1
    if abs(d2) <= abs(d1) and abs(d2) <= abs(d3): return d2
    if abs(d3) <= abs(d1) and abs(d3) <= abs(d2): return d3
#-----------------------------------------------#
def sincos2angle(sin,cos):
    # Get angle in first quadrant
    if   cos == 0.0: angle = np.pi/2.0
    elif sin == 0.0: angle = 0.0
    else           : angle = np.arctan( abs(sin/cos) )
    # quadrant 1 and 2
    if sin >= 0.0:
       if cos >= 0.0: angle = angle
       else         : angle = np.pi - angle
    # quadrant 3 and 4
    else:
       if cos >= 0.0: angle = 2*np.pi - angle
       else         : angle = angle + np.pi
    return angle
#-----------------------------------------------#
def angdist_sphere(p1,p2,units="rad"):
  '''
  Angular distance between points p1 and p2 on a sphere
  Each point is a tuple (theta,phi) with
     theta = th in [0, pi]
     phi   = ph in [0,2pi]
  Information:
      * spherical coordinate system (theta,phi definition)
  if units = "rad": input and output angles in radians
  if units = "deg": input and output angles in degrees
  '''
  # Points in sphere
  th1, ph1 = p1
  th2, ph2 = p2
  if units == "deg":
     th1 = np.deg2rad(th1)
     th2 = np.deg2rad(th2)
     ph1 = np.deg2rad(ph1)
     ph2 = np.deg2rad(ph2)
  sth1, cth1 = np.sin(th1), np.cos(th1)
  sth2, cth2 = np.sin(th2), np.cos(th2)
  sph1, cph1 = np.sin(ph1), np.cos(ph1)
  sph2, cph2 = np.sin(ph2), np.cos(ph2)
  v1 = [sth1*cph1,sth1*sph1,cth1]
  v2 = [sth2*cph2,sth2*sph2,cth2]
  angle = angle_vecs( v1,v2 )
  if units == "deg": angle = np.rad2deg(angle)
  return angle
#-----------------------------------------------#
def getPerim_circle(radius):
    return pc.TWOPI*radius
#-----------------------------------------------#
def getArea_triangle(pA,pB,pC):
    '''
    area of triangle by shoelace algorithm
    pA represents the vertex coordinates (xA,yA)
    pB represents the vertex coordinates (xB,yB)
    pC represents the vertex coordinates (xC,yC)
    '''
    Area = + (pB[0]*pC[1] - pC[0]*pB[1]) \
           - (pA[0]*pC[1] - pC[0]*pA[1]) \
           + (pA[0]*pB[1] - pB[0]*pA[1])
    Area = 0.5 * abs(Area)
    return Area
#-----------------------------------------------#
def getArea_circle(radius) :
    return pc.PI*radius**2
#-----------------------------------------------#
def getArea_sphere(radius) :
    return 4.*pc.PI*radius**2
#-----------------------------------------------#
def getVol_sphere(radius)  :
    return 4./3.*pc.PI*radius**3
#===============================================#


#===============================================#
# Symmetry operations                           #
#===============================================#
def getmatrix_inversion():
    return np.diag([-1.,-1.,-1.])
#-----------------------------------------------#
def getmatrix_Cn(n,u):
    ''' n: order; u: vector'''
    t    = pc.TWOPI/n
    cosu = np.cos(t)
    sinu = np.sin(t)
    ux,uy,uz = u
    Cngen = np.zeros( (3,3) )
    Cngen[0,0] = cosu + (ux**2)*(1.-cosu)
    Cngen[1,1] = cosu + (uy**2)*(1.-cosu)
    Cngen[2,2] = cosu + (uz**2)*(1.-cosu)
    Cngen[0,1] = ux*uy*(1.-cosu)-uz*sinu
    Cngen[0,2] = ux*uz*(1.-cosu)+uy*sinu
    Cngen[1,2] = uy*uz*(1.-cosu)-ux*sinu
    Cngen[1,0] = ux*uy*(1.-cosu)+uz*sinu
    Cngen[2,0] = ux*uz*(1.-cosu)-uy*sinu
    Cngen[2,1] = uy*uz*(1.-cosu)+ux*sinu
    return Cngen
#===============================================#


#===============================================#
# Integrating functions                         #
#===============================================#
def intg_trap(function,x0,xf,args=None,n=100,dx=None):
    '''
    function: the function to integrate
    args    : arguments of the function; function(x,args)
    '''
    if dx is None:
       dx = 1.0 * (xf-x0)/(n-1)
    else:
       n = int(round((xf-x0)/dx))+1
    xvalues = [x0+ii*dx for ii in range(n)]
    if args is not None: integrand = [function(x,*args) for x in xvalues]
    else               : integrand = [function(x)       for x in xvalues]
    integral  = sum([y*dx for y in integrand])
    return integral
#-----------------------------------------------#
def intg_gau(function,x0,xf,args=None,n=80):
    points, weights = np.polynomial.legendre.leggauss(n)
    suma     = (xf+x0)/2.0
    resta    = (xf-x0)/2.0
    # Points to evaluate and weights
    points   = [resta*xi+suma for xi in points]
    weights  = [wi*resta for wi in weights]
    # calculate integral
    if args is not None: integral = sum([w*function(x,*args) for x,w in zip(points,weights)])
    else               : integral = sum([w*function(x)       for x,w in zip(points,weights)])
    del points
    del weights
    return integral
#===============================================#





#===============================================#
# Diverse set of useful functions               #
#===============================================#
def frange(start,end,dx,include_end=True):
    '''
    returns a list of floats
    '''
    start = float(start)
    end   = float(end  )
    dx    = float(dx   )
    nsteps = int( round( (end-start)/dx ) )
    if include_end: nsteps = nsteps + 1
    return [ start+i*dx for i in range(nsteps)]
#-----------------------------------------------#
def prod_list(lists):
    '''
    lists is a tuple with lists of floats
    '''
    prod = [1.0 for ii in lists[0]]
    for list_i in lists:
        prod = [v1*v2 for v1,v2 in zip(prod,list_i)]
    return prod
#-----------------------------------------------#
def same_lfloats(list1,list2,eps=EPS_FLOAT):
    if len(list1) != len(list2): return False
    for f1,f2 in zip(list1,list2):
        diff = abs(f1-f2)
        if diff > eps: return False
    return True
#-----------------------------------------------#
def flatten_llist(list_of_lists):
    ''' converts list of lists to list'''
    flattened_list = [y for x in list_of_lists for y in x]
    return flattened_list
#-----------------------------------------------#
def uniquify_flist(list_of_floats,eps=EPS_FLOAT):
    ''' removes duplicates in a list of float numbers'''
    if len(list_of_floats) < 2: return list_of_floats
    # Copy list and sort it
    list2 = list(list_of_floats)
    list2.sort()
    # Get position of duplicates
    repeated = []
    if len(list2) > 1:
       for idx in range(1,len(list2)):
           previous_float = list2[idx-1]
           current_float  = list2[idx]
           if abs(current_float - previous_float) < eps: repeated.append(idx)
    # Change values of positions to Nones
    for idx in repeated: list2[idx] = None
    # Remove Nones from list
    removals = list2.count(None)
    for removal in range(removals): list2.remove(None)
    return list2
#-----------------------------------------------#
def ll2matrix(ll,varcomp=None,pos=0):
    '''
    Converts a list of lists (of different sizes)
    to a list of list with equal sizes.
    Begin N the maximum length:
          N = max([len(alist) for alist in ll])
    all the lists in ll are modified to present
    the same size. To do so, the variable 'varcomp'
    is added to each list until len(alist) = N.
    The variable is added at the beggining (pos=0)
    or at the end (pos=-1) of the list
    '''
    if pos not in [0,-1]: raise Exception("Wrong 'pos' variable given in ll2matrix")
    # empty list?
    if len(ll) == 0: return None
    # maximum size
    N = max([len(alist) for alist in ll])
    if N == 0: return None
    # Initialize matrix
    shape = (len(ll),N)
    matrix = [ [0.0 for col in range(shape[1])] for row in range(shape[0])]
    # Now, complete
    for row in range(shape[0]):
        arow = ll[row]
        while len(arow) < N:
              if pos ==  0: arow = [varcomp]+arow
              if pos == -1: arow = arow+[varcomp]
        for col in range(shape[1]):
            matrix[row][col] = arow[col]
    # return list of list
    return matrix
#-----------------------------------------------#
def remove_float(thefloat,thelist,eps=EPS_FLOAT):
    return [ii for ii in thelist if abs(ii-thefloat)>eps]
#-----------------------------------------------#
def uppt2matrix(upptriangle):
    if upptriangle is None: return None
    l = len(upptriangle)
    N = int( (-1 + (1+8*l)**0.5) / 2)
    index = 0
    matrix = [ [0.0 for ii in range(N)] for jj in range(N) ]
    for i in range(1,N+1):
        for j in range(i,N+1):
            matrix[i-1][j-1] = upptriangle[index]
            matrix[j-1][i-1] = upptriangle[index]
            index += 1
    return matrix
#-----------------------------------------------#
def lowt2matrix(lowtriangle):
    if lowtriangle is None: return None
    l = len(lowtriangle)
    N = int( (-1 + (1+8*l)**0.5) / 2)
    index = 0
    matrix = [ [0.0 for ii in range(N)] for jj in range(N) ]
    for i in range(1,N+1):
        for j in range(1,i+1):
            matrix[i-1][j-1] = lowtriangle[index]
            matrix[j-1][i-1] = lowtriangle[index]
            index += 1
    return matrix
#-----------------------------------------------#
def matrix2lowt(matrix):
    if matrix is None: return None
    nrows = len(matrix)
    ncols = len(matrix[0])
    low   = []
    for row in range(nrows):
        for col in range(0,row+1):
            Fij = matrix[row][col]
            low.append(Fij)
    return low
#-----------------------------------------------#
def time2human(t, units="secs"):
    ''' units --> "msecs", "secs" ,"mins" ,"hours" or "days" '''
    for ii in range(10):
        # deal with miliseconds
        if   units == "msecs":
           if   t > 1000: t,units = t/1000 , "secs"
        # deal with seconds
        elif units == "secs":
           if   t <  1  : t,units = t*1000 , "msecs"
           elif t > 60  : t,units = t/60   , "mins"
        # deal with minutes
        elif units == "mins":
           if   t <  1  : t,units = t*60   , "secs"
           elif t > 60  : t,units = t/60   , "hours"
        # deal with hours
        elif units == "hours":
           if   t <  1  : t,units = t*60   , "mins"
           elif t > 24  : t,units = t/24   , "days"
        # deal with days
        elif units == "days":
           if   t <  1  : t,units = t*24   , "hours"
        # other case, break
        else: break
    return (t, units)
#===============================================#




#===============================================#
# Functions related to molecules                #
#===============================================#
def vol2pressure(V,T): return pc.KB*T/V
#---------------------------------------------#
def pressure2vol(P,T): return pc.KB*T/P
#---------------------------------------------#
def symbol_and_atonum(symbol_or_atonum):
    try:
        atonum = int(symbol_or_atonum)
        symbol = dpt_z2s[atonum].strip()
    except:
        symbol = correct_symbol(symbol_or_atonum)
        atonum = dpt_s2z[symbol]
    return symbol,atonum
#---------------------------------------------#
def symbols_and_atonums(symbols_or_atonums):
    symbols = []
    atonums = []
    for atsym in symbols_or_atonums:
        symbol, atonum = symbol_and_atonum(atsym)
        symbols.append(symbol)
        atonums.append(atonum)
    return symbols, atonums
#---------------------------------------------#
def correct_symbol(symbol):
    symbol = symbol.strip()
    if symbol.upper() in "XX,X,DA": return "XX"
    return symbol[0].upper()+symbol[1:].lower()
#---------------------------------------------#
def clean_dummies(symbols,xcc=None,masses=None,gcc=None,Fcc=None):
    # initialize lists
    symbols_wo = []
    masses_wo  = []
    xcc_wo     = []
    gcc_wo     = []
    # first: clean only symbols
    indices_dummies = []
    for idx,symbol in enumerate(symbols):
        if str(symbol) == "0" or str(symbol).upper() == "XX":
           indices_dummies.append(idx)
           continue
        symbols_wo.append(symbol)
    # clean the rest of lists, if needed
    num0 =   len(symbols_wo)
    num1 = 3*len(symbols_wo)
    num2 = 3*len(symbols_wo)*(3*len(symbols_wo)+1)//2
    if masses is not None and len(masses) == num0: return symbols_wo,masses
    if xcc    is not None and len(xcc)    == num1: return symbols_wo,xcc
    if gcc    is not None and len(gcc)    == num1: return symbols_wo,gcc
    if Fcc    is not None and len(Fcc)    == num2: return symbols_wo,Fcc

    for idx,symbol in enumerate(symbols):
        if idx in indices_dummies: continue
        if   xcc    is not None: xcc_wo += xcc[3*idx:3*idx+3]
        elif gcc    is not None: gcc_wo += gcc[3*idx:3*idx+3]
        elif masses is not None: masses_wo.append(masses[idx])
    if   xcc    is not None: return symbols_wo,xcc_wo
    elif gcc    is not None: return symbols_wo,gcc_wo
    elif masses is not None: return symbols_wo,masses_wo
    # cleaning Fcc
    if Fcc is not None:
        # (a) in matrix format
        Fcc_matrix = lowt2matrix(Fcc)
        # (b) clean Fcc
        Fcc_wo = []
        for row,symbol1 in enumerate(symbols):
            the_row = []
            dummy1 = str(symbol1) == "0" or str(symbol1).upper() == "XX"
            if dummy1: continue
            for col,symbol2 in enumerate(symbols):
                dummy2 = str(symbol2) == "0" or str(symbol2).upper() == "XX"
                if dummy2: continue
                the_row.append( Fcc_matrix[row][col] )
            Fcc_wo.append(the_row)
        # (c) return to triangular list
        Fcc_wo = matrix2lowt(Fcc_wo)
        return symbols_wo,Fcc_wo
    # if nothing more, it just returns symbols_wo
    return symbols_wo
#---------------------------------------------#
def correct_symbols(symbols):
   '''
   for each symbol, first character in upper, the rest in lower
   '''
   return [correct_symbol(symbol) for symbol in symbols]
#---------------------------------------------#
def get_molformula(symbols):
    '''
    Returns the molecular formula of a given molecule
    '''
    formula_dict = {symbol:symbols.count(symbol) for symbol in symbols}
    molform = ""
    for key,value in sorted(formula_dict.items()):
        if value != 1: molform += "%s(%i)"%(key,value)
        if value == 1: molform += "%s"%(key)
    return molform
#---------------------------------------------#
def molformula2mass(mformu):
    digits  = "0123456789"
    symbols = ""
    nc      = len(mformu)
    for idx,character in enumerate(mformu):
        if character == ")":
           number   = int(number)
           current  = symbols.split()[-1]
           symbols += " ".join([current for ii in range(number-1)])
           symbols += " "
           continue
        if   character in digits:
           number += character
           continue
        if character in "(":
           number = ""
           continue
        symbols += character
        if idx+1 < nc:
           next_upper = (mformu[idx+1].upper() == mformu[idx+1])
           if next_upper: symbols += " "
    symbols = symbols.split()
    masses  = symbols2masses(symbols)
    totmass = sum(masses)
    return totmass, symbols
#---------------------------------------------#
def symbols2atonums(symbols):
    return [dpt_s2z[s] for s in symbols]
#---------------------------------------------#
def atonums2symbols(atonums):
    return [dpt_z2s[z].strip() for z in atonums]
#---------------------------------------------#
def get_symbols(atonums):
    return atonums2symbols(atonums)
#---------------------------------------------#
def atonums2masses(atonums):
    return [dpt_z2m[z] for z in atonums]
#---------------------------------------------#
def symbols2masses(symbols):
    return [dpt_s2m[s] for s in symbols]
#---------------------------------------------#
def get_atonums(symbols):
    return [dpt_s2z[s] for s in symbols]
#---------------------------------------------#
def howmanyatoms(xcc): return len(xcc)//3
#---------------------------------------------#
def xyz(xcc,at): return xcc[3*at:3*at+3]
#---------------------------------------------#
def x(xcc,at)  : return xcc[3*at+0]
#---------------------------------------------#
def y(xcc,at)  : return xcc[3*at+1]
#---------------------------------------------#
def z(xcc,at)  : return xcc[3*at+2]
#---------------------------------------------#
def get_centroid(xcc,indices=None):
    if indices is None: indices = range(len(xcc)//3)
    centroid_x = sum([x(xcc,idx) for idx in indices])/len(indices)
    centroid_y = sum([y(xcc,idx) for idx in indices])/len(indices)
    centroid_z = sum([z(xcc,idx) for idx in indices])/len(indices)
    return [centroid_x,centroid_y,centroid_z]
#---------------------------------------------#
def get_com(xcc,masses,indices=None):
    '''
    Returns the centre of mass of the selected atoms (indices)
    If no indices given, all atoms are considered
    '''
    if indices is None: indices = range(len(masses))
    tmass = sum([           masses[idx] for idx in indices])
    com_x = sum([x(xcc,idx)*masses[idx] for idx in indices])/tmass
    com_y = sum([y(xcc,idx)*masses[idx] for idx in indices])/tmass
    com_z = sum([z(xcc,idx)*masses[idx] for idx in indices])/tmass
    return [com_x,com_y,com_z]
#---------------------------------------------#
def get_distmatrix(xcc):
    nat = howmanyatoms(xcc)
    dmatrix = np.zeros( (nat,nat) )
    for ii in range(nat):
        xii = xyz(xcc,ii)
        for jj in range(ii+1,nat):
            xjj = xyz(xcc,jj)
            dmatrix[ii][jj] = distance(xii,xjj)
            dmatrix[jj][ii] = dmatrix[ii][jj]
    return dmatrix
#---------------------------------------------#
def set_origin(xcc,x0):
    '''
    returns xcc with x0 as origin
    '''
    nat = howmanyatoms(xcc)
    return [xi-xj for xi,xj in zip(xcc,nat*x0)]
#---------------------------------------------#
def shift2com(xcc,masses):
    '''
    Function to shift to center of mass
    '''
    com = get_com(xcc,masses)
    xcc = set_origin(xcc,com)
    return xcc
#---------------------------------------------#
def islinear(xcc):
    nats    = howmanyatoms(xcc)
    masses  = [1.0 for at in range(nats)]
    itensor = get_itensor_matrix(xcc,masses)
    linear  = get_itensor_evals(itensor)[3]
    return linear
#---------------------------------------------#
def same_geom(x1,x2,eps=EPS_GEOM):
    diff = [xi-xj for xi,xj in zip(x1,x2)]
    #value = norm(diff) / len(diff)
    value = max(diff)
    if value < eps: return True
    else          : return False
#---------------------------------------------#
def center_and_orient(xcc,gcc,Fcc,masses):
    # number of atoms?
    nats = howmanyatoms(xcc)
    # in center of mass
    xcc = shift2com(xcc,masses)
    if nats == 1: return xcc, gcc, Fcc
    # Is it linear?
    itensor = get_itensor_matrix(xcc,masses)
    evalsI, rotTs, rtype, linear = get_itensor_evals(itensor)
    if not linear: return xcc, gcc, Fcc
    # Get moments and axis of inertia
    itensor = np.matrix(itensor)
    evalsI, evecsI = np.linalg.eigh(itensor)
    evecsI = np.matrix(evecsI)
    # Rotation matrix 3N x 3N
    zeros = np.zeros((3, 3))
    R= []
    for at in range(nats):
        row = []
        count = 0
        while count < at:
              row.append(zeros)
              count += 1
        row.append(evecsI)
        while count < nats-1:
              row.append(zeros)
              count += 1
        R.append(row)
    R = np.matrix(np.block(R))
    # Rotate coordinates
    xcc = np.matrix(xcc) * R
    xcc = xcc.tolist()[0]
    # Rotate gradient
    if gcc is not None and len(gcc) != 0:
       gcc = np.matrix(gcc) * R
       gcc = gcc.tolist()[0]
    # rotate Fcc
    if Fcc is not None and len(Fcc) != 0:
       nr,nc = np.matrix(Fcc).shape
       if nr != nc: Fcc = lowt2matrix(Fcc)
       Fcc = R.transpose() * np.matrix(Fcc) * R
       Fcc = Fcc.tolist()
    return xcc, gcc, Fcc
#---------------------------------------------#
def cc2ms_x(xcc,masses,mu=1.0/pc.AMU):
    ''' cartesian --> mass-scaled; x'''
    if xcc is None or len(xcc) == 0: return xcc
    nat  = len(masses)
    cfs  = [ (masses[idx]/mu)**0.5 for idx in range(nat)]
    xms  = [ [x(xcc,at)*cfs[at],y(xcc,at)*cfs[at],z(xcc,at)*cfs[at]] for at in range(nat)]
    xms  = flatten_llist(xms)
    return xms
#---------------------------------------------#
def cc2ms_g(gcc,masses,mu=1.0/pc.AMU):
    ''' cartesian --> mass-scaled; gradient'''
    if gcc is None or len(gcc) == 0: return gcc
    nat  = len(masses)
    cfs  = [ (masses[idx]/mu)**0.5 for idx in range(nat)]
    gms  = [ [x(gcc,at)/cfs[at],y(gcc,at)/cfs[at],z(gcc,at)/cfs[at]] for at in range(nat)]
    gms  = flatten_llist(gms)
    return gms
#---------------------------------------------#
def cc2ms_F(Fcc,masses,mu=1.0/pc.AMU):
    ''' cartesian --> mass-scaled; force constant matrix'''
    if Fcc is None or len(Fcc) == 0: return Fcc
    nat    = len(masses)
    Fms    = [ [ 0.0 for at1 in range(3*nat) ] for at2 in range(3*nat)]
    for i in range(3*nat):
        mi = masses[int(i/3)]
        for j in range(3*nat):
            mj = masses[int(j/3)]
            f = mu / ((mi*mj)**0.5)
            Fms[i][j] = Fcc[i][j] * f
    return Fms
#---------------------------------------------#
def ms2cc_x(xms,masses,mu=1.0/pc.AMU):
    ''' mass-scaled --> cartesian; x'''
    if xms is None or len(xms) == 0: return xms
    nat  = len(masses)
    cfs  = [ (masses[idx]/mu)**0.5 for idx in range(nat)]
    xcc  = [ [x(xms,at)/cfs[at],y(xms,at)/cfs[at],z(xms,at)/cfs[at]] for at in range(nat)]
    xcc  = flatten_llist(xcc)
    return xcc
#---------------------------------------------#
def ms2cc_g(gms,masses,mu=1.0/pc.AMU):
    ''' mass-scaled --> cartesian; gradient'''
    if gms is None or len(gms) == 0: return gms
    nat  = len(masses)
    cfs  = [ (masses[idx]/mu)**0.5 for idx in range(nat)]
    gcc  = [ [x(gms,at)*cfs[at],y(gms,at)*cfs[at],z(gms,at)*cfs[at]] for at in range(nat)]
    gcc  = flatten_llist(gcc)
    return gcc
#---------------------------------------------#
def ms2cc_F(Fms,masses,mu=1.0/pc.AMU):
    ''' mass-scaled --> cartesian; force constant matrix'''
    if Fms is None or len(Fms) == 0: return Fms
    nat    = len(masses)
    Fcc    = [ [ 0.0 for at1 in range(3*nat) ] for at2 in range(3*nat)]
    for i in range(3*nat):
        mi = masses[int(i/3)]
        for j in range(3*nat):
            mj = masses[int(j/3)]
            f = mu / ((mi*mj)**0.5)
            Fcc[i][j] = Fms[i][j] / f
    return Fcc
#===============================================#



#===============================================#
# Functions related to rotations                #
#===============================================#
def gen_rotmatrix(axis,theta):
    '''
    Generates the rotation matrix around axis.
    Input:
      * axis: the (x,y,z) coordinates of the
              direction vector of the axis
      * theta: the rotation angle (radians)
    Returns:
      * rot_matrix: a 3x3 numpy matrix

    Note: The rotation considers that the axis is
    situated at the origin

    Direction of rotation given by right-hand rule
    '''

    # Theta in [0,2*pi]
    while theta < 0.0:
          theta = theta + 2.0*np.pi
    while theta > 2*np.pi:
          theta = theta - 2.0*np.pi
    # Theta in [-pi,pi]
    if theta > np.pi:
       theta = -(2.0*np.pi - theta)

    # Normalize axis length
    ux, uy, uz = normalize_vec(axis)
    st = np.sin(theta)
    ct = np.cos(theta)

    R11 = ux*ux*(1.- ct) +    ct
    R12 = ux*uy*(1.- ct) - uz*st
    R13 = ux*uz*(1.- ct) + uy*st

    R21 = uy*ux*(1.- ct) + uz*st
    R22 = uy*uy*(1.- ct) +    ct
    R23 = uy*uz*(1.- ct) - ux*st

    R31 = uz*ux*(1.- ct) - uy*st
    R32 = uz*uy*(1.- ct) + ux*st
    R33 = uz*uz*(1.- ct) +    ct

    rot_matrix = np.matrix( [ [R11,R12,R13],[R21,R22,R23],[R31,R32,R33] ] )
    return rot_matrix
#-----------------------------------------------#
def rotate_point(xyz,RotMat):
    RotMat = np.matrix(RotMat)
    # assert right orientation of new coordinate system
    if np.linalg.det(RotMat) < 0.0: RotMat = - RotMat
    # Rotate atoms in fragment
    rotated_xyz = []
    xyz = RotMat * np.matrix(xyz).transpose()
    xyz = np.array((xyz.transpose()).tolist()[0])
    return xyz
#-----------------------------------------------#
def rotate_molecule(xcc,RotMat):
    # assert right orientation of new coordinate system
    RotMat = np.matrix(RotMat)
    if np.linalg.det(RotMat) < 0.0: RotMat = - RotMat
    # rotate atom by atom
    natoms = len(xcc)//3
    final_xcc = []
    for atom in range(natoms):
        xyz = xcc[3*atom:3*atom+3]
        xyz = list(rotate_point(xyz,RotMat))
        final_xcc += xyz
    return final_xcc
#-----------------------------------------------#
def rotate_internal(xcc,ugraph,tbond,theta):
    if xcc is None: return None
    natoms = int(round(len(xcc)/3))
    # convert to numpy array
    xvector = np.array(xcc,copy=True)
    # Get two fragments
    idxA, idxB = tbond
    A_frag = set(ugraph.bfsearch1d(idxB,idxA))
    B_frag = set(ugraph.bfsearch1d(idxA,idxB))

    # Compare fragments. They may be equal in case of cyclic systems
    if (B_frag is None) or (A_frag is None) or (B_frag == A_frag):
       return None

    # Choose smaller fragment
    if len(A_frag) > len(B_frag):
       # if B_frag, rotation around A-->B
       x0   = xvector[3*idxA:3*idxA+3]
       axis = xvector[3*idxB:3*idxB+3] - x0
       target_fragment = B_frag.copy()
    else:
       # if A_frag, rotation around B-->A
       x0   = xvector[3*idxB:3*idxB+3]
       axis = xvector[3*idxA:3*idxA+3] - x0
       target_fragment = A_frag.copy()
    axis = axis / np.linalg.norm(axis)

    # Remove indices of the bond
    target_fragment.discard(idxA)
    target_fragment.discard(idxB)

    # Get rotation matrix
    R = gen_rotmatrix(axis,theta)

    # Rotate atoms in fragment
    rotated_xyz = []
    for idx in range(natoms):
        xyz    = xvector[3*idx:3*idx+3]
        if idx in target_fragment:
            xyz = xyz - x0
            xyz = R * np.matrix(xyz).transpose()
            xyz = np.array((xyz.transpose()).tolist()[0])
            xyz = xyz + x0
        rotated_xyz += xyz.tolist()
    return rotated_xyz
#===============================================#



#===============================================#
# Functions related to frequencies              #
#===============================================#
def afreq2wnum(angfreq):
    return angfreq / pc.TWOPI / pc.C0
#-----------------------------------------------#
def wnum2afreq(wavenum):
    return wavenum * pc.TWOPI * pc.C0
#-----------------------------------------------#
def eval2afreq(evalue,mu=1.0/pc.AMU):
    return sign(evalue) * (abs(evalue)/mu)**0.5
#-----------------------------------------------#
def eval2wnum(evalue,mu=1.0/pc.AMU):
    return wnum2afreq(eval2afreq(evalue,mu))
#-----------------------------------------------#
def afreq2zpe(angfreq):
    if angfreq < 0.0: return 0.0
    return pc.HBAR * angfreq / 2.0
#-----------------------------------------------#
def afreq2turnpoint(angfreq,mu):
    if angfreq < 0.0: return 1e10
    return np.sqrt( pc.HBAR / angfreq / mu )
#-----------------------------------------------#
def wnum2zpe(wavenum):
    angfreq = wnum2afreq(wavenum)
    if angfreq < 0.0: return 0.0
    return pc.HBAR * angfreq / 2.0
#-----------------------------------------------#
def eval2cm(evalue,mu=1.0/pc.AMU):
    return eval2wnum(evalue,mu)/pc.CM
#-----------------------------------------------#
def afreq2cm(angfreq):
    return afreq2wnum(angfreq)/pc.CM
#-----------------------------------------------#
def cm2afreq(cm):
    return wnum2afreq(cm * pc.CM)
#-----------------------------------------------#
def afreq2eV(angfreq):
    return afreq2wnum(angfreq)/pc.CM
#-----------------------------------------------#
def same_freqs(ccfreqs,icfreqs,epsilon=EPS_CCIC):
    # compare lengths
    if len(ccfreqs) != len(icfreqs):
        return False
    # compare freqs
    for ccf, icf in zip(ccfreqs,icfreqs):
        ccf = afreq2cm(ccf)
        icf = afreq2cm(icf)
        if abs(ccf-icf) > epsilon: return False
    return True
#-----------------------------------------------#
def numimag(freqs):
    return [freq<0.0 for freq in freqs].count(True)
#===============================================#




#===============================================#
# Rotations/Vibrations                          #
#===============================================#
def get_itensor_matrix(xcc,masses):
    ''' returns inertia tensor (au)'''
    nat = howmanyatoms(xcc)
    inertia = [[0.0 for i in range(3)] for j in range(3)]
    for i in range(nat):
        # Diagonal elements
        inertia[0][0] += masses[i] * (y(xcc,i)**2 + z(xcc,i)**2)
        inertia[1][1] += masses[i] * (z(xcc,i)**2 + x(xcc,i)**2)
        inertia[2][2] += masses[i] * (x(xcc,i)**2 + y(xcc,i)**2)
        # Offdiagonal elements
        inertia[0][1] -= masses[i] * x(xcc,i) * y(xcc,i)
        inertia[0][2] -= masses[i] * z(xcc,i) * x(xcc,i)
        inertia[1][2] -= masses[i] * y(xcc,i) * z(xcc,i)
    inertia[1][0] = inertia[0][1]
    inertia[2][0] = inertia[0][2]
    inertia[2][1] = inertia[1][2]
    return inertia
#---------------------------------------------#
def get_itensor_evals(itensor):
    evalsI, evecsI = np.linalg.eigh(itensor)

    Ia, Ib, Ic = evalsI

    bool_a  = abs(Ia)        < EPS_INERTIA # Ia = 0
    bool_ab = abs(Ia/Ib-1.0) < EPS_FLOAT   # Ia = Ib
    bool_bc = abs(Ib/Ic-1.0) < EPS_FLOAT   # Ib = Ic
    bool_abc = bool_ab and bool_bc

    if   bool_abc  : rtype = "spherical top"
    elif bool_ab   : rtype = "oblate symmetric top"
    elif bool_bc   :
         if  bool_a: rtype = "linear rotor"
         else      : rtype = "prolate symmetric top"
    else           : rtype = "asymmetric top"

    if rtype == "linear rotor":
       linear = True
       evalsI = [evalsI[1]]
    else:
       linear = False

    # rotational temperature
    rotTs = [pc.HBAR**2 / (2*Ii*pc.KB) for Ii in evalsI]
    return evalsI, rotTs, rtype, linear
#-----------------------------------------------#
def get_projectionmatrix(xcc,masses,v0=None):
    '''
    Generates matrix to project translation
    and rotation coordinates (mass scaled/weighted)
    Other coordinate can be projected by introducing it
    using v0 (in mass-scaled)
    '''
    vecs = []
    nat  = len(masses)
    # PROJECT TRA IN HESS FOR FREQS
    if PROJECT_TRA:
       # translation
       sqrtmasses = [np.sqrt(mass) for mass in masses]
       b1 = [term if ii==0 else 0.0 for term in sqrtmasses for ii in range(3)]
       b2 = [term if ii==1 else 0.0 for term in sqrtmasses for ii in range(3)]
       b3 = [term if ii==2 else 0.0 for term in sqrtmasses for ii in range(3)]
       norm1 = np.linalg.norm(b1)
       norm2 = np.linalg.norm(b2)
       norm3 = np.linalg.norm(b3)
       b1 /= norm1
       b2 /= norm2
       b3 /= norm3
       vecs += [b1,b2,b3]
    # PROJECT ROT IN HESS FOR FREQS
    if PROJECT_ROT:
       # rotation
       b4 = np.zeros(len(xcc))
       b5 = np.zeros(len(xcc))
       b6 = np.zeros(len(xcc))
       for i in range(nat):
           b4[3*i + 1] =   np.sqrt(masses[i]) * z(xcc,i)
           b4[3*i + 2] = - np.sqrt(masses[i]) * y(xcc,i)
           b5[3*i + 0] = - np.sqrt(masses[i]) * z(xcc,i)
           b5[3*i + 2] =   np.sqrt(masses[i]) * x(xcc,i)
           b6[3*i + 0] =   np.sqrt(masses[i]) * y(xcc,i)
           b6[3*i + 1] = - np.sqrt(masses[i]) * x(xcc,i)
       norm4 = np.linalg.norm(b4)
       norm5 = np.linalg.norm(b5)
       norm6 = np.linalg.norm(b6)
       if norm4 > EPS_NORM: b4 /= norm4; vecs.append(b4)
       if norm5 > EPS_NORM: b5 /= norm5; vecs.append(b5)
       if norm6 > EPS_NORM: b6 /= norm6; vecs.append(b6)
    # Gram Schmidt
    if len(vecs) != 0:
       X = np.matrix(vecs).transpose()
       X_gs, R = np.linalg.qr(X)
       projmatrix = X_gs * X_gs.H
    else:
       projmatrix = np.zeros( (3*nat,3*nat) )
    # PROJECT GRADIENT
    if v0 is not None:
       normv0 = np.linalg.norm(v0)
       if normv0 > EPS_NORM:
          v0 = np.matrix( v0 ) / normv0
          projmatrix += v0.transpose() * v0
    return projmatrix
#-----------------------------------------------#
def project_hessian(Fms,natoms,proj_matrix):
    ''' Fms: hessian in mass-scaled '''
    # identity matrix
    I = np.identity(3*natoms)
    # Calculate projected hessian matrix
    Fms_proj = (I - proj_matrix) * Fms * (I - proj_matrix)
    return Fms_proj
#-----------------------------------------------#
def diagonalize_hessian(Fms,mu=1.0/pc.AMU):
    # as Fms is symmetric --> diagonalize with eigh
    evalsF, evecsF = np.linalg.eigh(Fms)
    # Convert evals to angular frequencies
    ccfreqs = [eval2afreq(v,mu) for v  in evalsF]
    # evecsF to list of eigenvectors
    evecsF  = evecsF.transpose().tolist()
    # return data
    return ccfreqs, evalsF, evecsF
#-----------------------------------------------#
def detect_frozen(Fcc,nat):
    '''
    Fcc --> 3Nx3N
    '''
    frozen =  []
    if Fcc is None or len(Fcc) == 0: return frozen
    Fcc = np.matrix(Fcc)
    for at in range(nat):
        # get columns (rows are equivalent; symmetric matrix)
        colx = Fcc[:][3*at+0].tolist()[0]
        coly = Fcc[:][3*at+1].tolist()[0]
        colz = Fcc[:][3*at+2].tolist()[0]
        # Get norm of column
        normx = np.linalg.norm(colx)
        normy = np.linalg.norm(coly)
        normz = np.linalg.norm(colz)
        # are columns made of zeros?
        if normx != 0.0: continue
        if normy != 0.0: continue
        if normz != 0.0: continue
        frozen.append(at)
    return frozen
#-----------------------------------------------#
def calc_ccfreqs(Fcc,masses,xcc,mu=1.0/pc.AMU,v0=None):
    '''
    xcc has to be centered at com
    v0 in case gradient has to be removed
    '''
    # num of atoms and Fcc in matrix format
    nat = len(masses)
    if len(Fcc) != 3*nat: Fcc = lowt2matrix(Fcc)
    # Frozen atoms?
    frozen  = detect_frozen(Fcc,nat)
    boolN   = np.array([at not in frozen for at in range(nat)])
    bool3N  = np.array([at not in frozen for at in range(nat) for ii in range(3)])
    # if FROZEN atoms, then reduce system!!
    if len(frozen) != 0:
       masses = np.array(masses)[boolN]
       xcc    = np.array(xcc)[bool3N]
       if v0 is not None: v0 = np.array(v0)[bool3N]
       # force constant matrix
       Fcc = [[Fcc[idx1][idx2] for idx1 in range(3*nat) if bool3N[idx1]]\
                               for idx2 in range(3*nat) if bool3N[idx2]]
       # num atoms
       nat    = len(masses)
       # re-center system
       xcc = shift2com(xcc,masses)
    # Analyze system
    linear = islinear(xcc)
    if linear: nvdof = 3*nat-5
    else     : nvdof = 3*nat-6
    if v0 is not None: nvdof -= 1
    # mass-scaled hessian
    Fms = cc2ms_F(Fcc,masses,mu=mu)
    Fms = np.matrix(Fms)
    # projection matrix
    pmatrix = get_projectionmatrix(xcc,masses,v0)
    # projected hessian
    Fms = project_hessian(Fms,nat,pmatrix)
    # Diagonalization
    ccfreqs, evalsF, evecsF = diagonalize_hessian(Fms,mu)
    # remove extra frequencies
    idxs = sorted([(abs(fq),fq,idx) for idx,fq in enumerate(ccfreqs)])
    idxs.reverse()
    idxs = idxs[:nvdof]
    idxs = sorted([(fq,idx) for absfq,fq,idx in idxs])
    idxs = [idx for fq,idx in idxs]
    ccfreqs = [ccfreqs[idx] for idx in idxs]
    evalsF  = [ evalsF[idx] for idx in idxs]
    evecsF  = [ evecsF[idx] for idx in idxs]
    # Consider the removed atoms in the eigenvectors
    if len(frozen) != 0:
       for idx,evecF in enumerate(evecsF):
           evecsF[idx] = [evecF.pop(0) if booli else 0.0 for booli in bool3N]
    return ccfreqs, evalsF, evecsF
#-----------------------------------------------#
def scale_freqs(freqs,fscal):
    return [freq*fscal for freq in freqs]
#===============================================#


#=============================================#
#        Functions related to extrema         #
#=============================================#
def minima_in_list(lx,ly):
    '''
    in the list, it find points
    that may be local minima
    Returns the list of x-guesses
    '''
    np = len(lx)
    # initialize guesses
    guesses = []
    # initial point
    if ly[0] <= ly[1]:
       guesses.append(lx[0])
    # mid points
    for idx in range(1,np-1):
        xi, yi = lx[idx-1],ly[idx-1]
        xj, yj = lx[idx  ],ly[idx  ]
        xk, yk = lx[idx+1],ly[idx+1]
        if yj <= yi and yj <= yk: guesses.append(xj)
    # final points
    if ly[-1] <= ly[-2]:
       guesses.append(lx[-1])
    return guesses
#=============================================#

#=============================================#
# Functions to print certains variables       #
#=============================================#
def print_matrix(hessian,f="%+6.3f"):
    nrows, ncols = hessian.shape
    l       = len(f%1.0)
    sint    = "%%%ii"%l
    STRING  = "     * shape  = %ix%i\n"%(nrows,ncols)
    STRING += " "*8 +"  ".join([sint%(ii+1) for ii in range(ncols)])+"\n"
    for row in range(nrows):
        STRING += "%6i  "%(row+1) + "  ".join([f%value for value in hessian[row,:].tolist()[0]])+"\n" 
        STRING += "\n"
    print(STRING)
#---------------------------------------------#
def print_string(string,nbs=0):
    if string ==  "" : return
    if string == "\n": return
    for line in string.split("\n"): print(" "*nbs+line)
#=============================================#

