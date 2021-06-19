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
| Module     :  modpilgrim         |
| Sub-module :  names              |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the names of
files and folders used in Pilgrim.

There are also functions to get
such a names.
'''


#===============================================================#
UFOLDER = "UDATA/"                                              #
ANHDIR  = "ANHAR/"                                              #
TMP     = "TMP/"                                                #
#===============================================================#
TMPi   = TMP + "MEPcalcs_%s/"                                   #
TMPHLi = TMP + "HLcalcs_%s/"                                    #
DIR1   = "1-GTS/"                                               #
DIR2   = "2-PLG_DATA/"                                          #
DIR3   = "3-PLG_OUTPUT/"                                        #
DIR4   = "4-PLG_RST/"                                           #
DIR5   = "5-MOLDEN/"                                            #
DIR6   = "6-PLOTFILES/"                                         #
DIR7   = "7-SUMMARIES/"                                         #
#===============================================================#
IFILE0 = "tracking"                                             #
IFILE1 = "pif.struc"                                            #
IFILE2 = "pif.temp"                                             #
IFILE3 = "pif.path"                                             #
IFILE4 = "pif.calcs"                                            #
IFILE5 = "pif.chem"                                             #
IFILE6 = "pif.kmc"                                              #
IFILE7 = "pif.dlevel"                                           #
CFILE  = DIR4+"convergence.txt"                                 #
#===============================================================#


#===============================================================#
def name2data(name):                                  
    if name is None: return None, None
    if "." in name: ctc, itc = name.split(".")        
    else          : ctc, itc = name, None             
    return ctc, itc                                   
#---------------------------------------------------------------#
def struckey(ctc,itc=None):
    name = ctc                                        
    if itc is not None: name += "."+itc               
    return name                                       
#===============================================================#


#===============================================================#
def get_gts(ctc,itc):
    return DIR1 + "%s.gts"%struckey(ctc,itc)
#---------------------------------------------------------------#
def get_rst(ctc,itc):
    name = struckey(ctc,itc)
    return DIR4 + "%s.rst"%name
#---------------------------------------------------------------#
def get_dof(dlevel):
    if dlevel: fname = "data.dlevel"
    else     : fname = "data.slevel"
    return DIR2 + fname
#---------------------------------------------------------------#
def get_hlf():
    fname = "highlevel.txt"
    return DIR2 + fname
#---------------------------------------------------------------#
def get_pof(dlevel,case,target=None):
    fname = case
    if target is not None: fname += ".%s"%target
    if dlevel            : fname += ".dlevel.txt"
    else                 : fname += ".slevel.txt"
    return DIR3 + fname
#---------------------------------------------------------------#
def get_gtsmolden(ctc,itc):
    return DIR5 + "sp.%s.molden"%struckey(ctc,itc)
#---------------------------------------------------------------#
def get_rstxyz(ctc,itc):
    name = struckey(ctc,itc)
    return DIR5 + "path.%s.xyz"%name
#---------------------------------------------------------------#
def get_plf(dlevel=False):
    if dlevel: fname = "plots.dlevel.txt"
    else     : fname = "plots.slevel.txt"
    return DIR6 + fname
#===============================================================#


