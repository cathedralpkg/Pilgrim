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
| Sub-module :  interpolate        |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#=============================================#
from   common.Spline  import Spline
#=============================================#


#=============================================#
#     Functions related to interpolations     #
#=============================================#
def interpolate(xvalues,yvalues,x,d=0):
    spl = Spline(xvalues,yvalues,tension=0.0)
    if d == 0: return spl(x)
    # get derivative
    return spl(x), spl.derivative(x)
#---------------------------------------------#
def interpolate_nones(xx,yy,mode="cubic"):
    '''
    element to correct in yy is given as None
    '''
    npts = len(yy)
    # correct Nones with spline
    if mode == "cubic":
       # create spline
       copy_xx = [x for idx,x in enumerate(xx) if yy[idx] is not None]
       copy_yy = [y for idx,y in enumerate(yy) if yy[idx] is not None]
       spl     = Spline(copy_xx,copy_yy,tension=0.0)
       # interpolate
       for idx in range(npts):
           x,y = xx[idx],yy[idx]
           if y is not None: continue
           yy[idx] = spl(x)
    # correct Nones with linear interpolation
    if mode == "linear":
       for ii in range(npts):
           x,y = xx[ii], yy[ii]
           if y is not None: continue
           idx0 = ii -1
           idx2 = None
           for jj in range(ii,npts):
               if yy[jj] is not None: idx2 = jj; break
           if idx2 is None: continue
           # interpolate
           x0 , y0 = xx[idx0], yy[idx0]
           x2 , y2 = xx[idx2], yy[idx2]
           tangent = (y2-y0)/(x2-x0)
           yy[ii]  = y0 + tangent * (x-x0)
    return yy
#=============================================#

