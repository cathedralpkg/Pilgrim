'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.1
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
| Sub-module :  pilesso            |
| Last Update:  2020/03/01 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module is in charge of returning
the Python function used for the
calculation with selected software.
'''


#---------------------------------------------------------------#
import os
from   common.files      import read_file
from   common.fncs       import clean_lines
from   common.Exceptions import UnknownSoft
#---------------------------------------------------------------#


#===============================================================#
def get_dsoft():
    '''
    This functions is in charge of reading the file
    mesc.txt where the available softwares are listed
    '''
    # defined here:
    dsoft   = {}
    if True:
       dsoft['gaussian'] = 'modpilgrim.pilgau'
       dsoft['orca'    ] = 'modpilgrim.pilorca'
      #dsoft['clhbrpot'] = 'modpilgrim.pilasurf'

    # read from file
    if False:
       TXTFILE = os.path.dirname(os.path.realpath(__file__))+"/mesc.txt"
       # read file
       lines = read_file(TXTFILE)
       lines = clean_lines(lines,"#",True)
       # get info from lines
       for line in lines:
           if line == "\n": continue
           software, module = line.split()
           dsoft[software] = module
    # return dictionary
    return dsoft
#---------------------------------------------------------------#
def get_templates(ch,mtp,case="LL"):
    templates = {}
    # get info from file
    dsoft = get_dsoft()
    for software in dsoft.keys():
        module = dsoft[software]
        lib = __import__(module,fromlist=[''])
        globals()["mes"] = lib
        if   case == "LL": string = mes.pilgrim_template(ch,mtp)
        elif case == "HL": string = mes.pilgrim_templateHL(ch,mtp)
        else             : string = ""
        if string is None: string = software+"\n" # for asurf
        templates[software] = string
    return templates
#---------------------------------------------------------------#
def get_spc_fnc(software):
    '''
    Returns function to performs a single-point calculation
    using the corresponding interface
    '''
    # get info from file
    dsoft = get_dsoft()
    # valid software?
    if software not in dsoft.keys():
        raise UnknownSoft(Exception)
    # import module
    module = dsoft[software]
    lib = __import__(module,fromlist=[''])
    globals()["mesc"] = lib
    # call function from module and return data
    return mesc.pilgrim_spc
#===============================================================#


