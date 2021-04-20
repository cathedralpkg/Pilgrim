'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.2
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
| Sub-module :  optHARVEST         |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --harvest option
of Pilgrim
'''

#--------------------------------------------------#
import datetime
import time
import os
import sys
#--------------------------------------------------#
import modpilgrim.pilrw         as     RW
import modpilgrim.names             as     PN
import modpilgrim.strings           as     PS
#--------------------------------------------------#
from   common.fncs       import print_string
from   common.fncs       import time2human
#--------------------------------------------------#
from   common.Logger     import Logger
from   modpilgrim.diverse           import status_check
from   modpilgrim.diverse           import ffchecking
#--------------------------------------------------#


def main(idata,status,case,targets):

    stat2check = [5]
    mustexist  = []
    tocreate   = [PN.DIR7]
    #-------------------------------------------------------#
    # Read Pilgrim input files, check file/folder status    #
    # and expand tuple 'case'                               #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, dkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,stat2check)
    if fstatus == -1: exit()
    # existency of folders
    fstatus = ffchecking(mustexist,tocreate)
    if fstatus == -1: exit()
    # expand case
    (dof,hlf,plotfile),dlevel,software = case
    #-------------------------------------------------------#


    #---------------------------------#
    # files with data and output file #
    #---------------------------------#
    if len(targets) == 0 and dchem is not None: targets = dchem.keys()

    for rcname in sorted(targets):
        # pof file
        t1 = time.time()
        if dlevel: pof = PN.DIR7 + "summary.%s.dlevel.txt"%(rcname)
        else     : pof = PN.DIR7 + "summary.%s.slevel.txt"%(rcname)
        print("   --> %s"%rcname)
        if rcname not in dchem.keys():
           print("       reaction NOT FOUND!")
           print("")
           continue

        print("       * Summary file: %s"%pof)
        sys.stdout = Logger(pof,"w",True)
        #sys.stdout.writeinfile(PS.init_txt())
        sinit = "# DATA FILES FOR REACTION: %s #"%rcname
        sinit = "#"*len(sinit)+"\n"+sinit+"\n"+"#"*len(sinit)+"\n\n"
        sys.stdout.writeinfile(sinit)

        Rs, TS, Ps = dchem[rcname]
        string = ""
        # pfn for reactants, ts and PS
        for ctc in Rs+Ps+[TS]:
            if ctc is None: continue
            pof_i = PN.get_pof(dlevel,"pfn",ctc)
            if not os.path.exists(pof_i):
               print("       * WARNING! NOT FOUND: %s"%pof_i)
               continue
            try:
               with open(pof_i,'r') as asdf: lines = "".join(asdf.readlines())
               string += lines
            except:
               print("       * WARNING! Problems with file: %s"%pof_i)
        # path files
        if TS is not None:
           ctc, itc = PN.name2data(TS)
           if itc is not None: itcs = [itc]
           else              : itcs = sorted([itc for itc,w in dctc[ctc]._itcs])
           for itc in itcs:
               target = PN.struckey(ctc,itc)
               pof_i = PN.get_pof(dlevel,"path",target)
               if not os.path.exists(pof_i):
                  print("       * WARNING! NOT FOUND: %s"%pof_i)
                  continue
               try:
                  with open(pof_i,'r') as asdf: lines = "".join(asdf.readlines())
                  string += lines
               except:
                  print("       * WARNING! Problems with file: %s"%pof_i)
        # rcons files
        pof_i = PN.get_pof(dlevel,"rcons",rcname)
        if not os.path.exists(pof_i):
            print("       * WARNING! NOT FOUND: %s"%pof_i)
        else:
           try:
               with open(pof_i,'r') as asdf: lines = "".join(asdf.readlines())
               string += lines
           except:
               print("       * WARNING! Problems with file: %s"%pof_i)
        # write file
        sys.stdout.writeinfile(string)
        # print end of file
        t2    = time.time()
        etime = time2human(t2-t1,"secs")
        sys.stdout.writeinfile(PS.end_txt(*etime))
        # End logger
        sys.stdout = Logger(None)
        print("")


#---------------------------------------------------------------#


