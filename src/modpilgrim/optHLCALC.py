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
| Sub-module :  optHLCALC          |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#--------------------------------------------------#
import os
#--------------------------------------------------#
import modpilgrim.names        as PN
import modpilgrim.pilrw         as RW
import modpilgrim.strings      as PS
#--------------------------------------------------#
from   modpilgrim.pilesso          import get_spc_fnc
from   modpilgrim.diverse      import ffchecking
from   modpilgrim.diverse      import status_check
from   modpilgrim.diverse      import find_label_in_rst
from   modpilgrim.steepdesc    import sorted_points
#--------------------------------------------------#
import common.files      as ff
import common.fncs       as fncs
import common.Exceptions as Exc
#--------------------------------------------------#
from   common.files     import read_gtsfile
from   common.Molecule  import Molecule
#--------------------------------------------------#


#===============================================================#
def get_targets(targets,dctc):
    # No targets selected
    if targets == [] or targets == "*":
       targets = [PN.struckey(ctc,itc) for ctc in dctc.keys() for itc,weight in dctc[ctc]._itcs ]
       targets.sort()
       return targets
    # Check targets
    itargets = []
    for target in targets:
        np = target.count(".")
        # a whole ctc is selected
        if np == 0:
           ctc = target
           if ctc not in dctc.keys(): continue
           itargets += [PN.struckey(ctc,itc) for itc,weight in dctc[ctc]._itcs]
        # a specific itc
        elif np == 1:
            ctc,itc = target.split(".")
            if ctc not in dctc.keys(): continue
            itcs = [itc_i for itc_i,weight_i in dctc[ctc]._itcs]
            if itc not in itcs: continue
            itargets += [target]
    return itargets
#---------------------------------------------------------------#
def highlevel_sp(ctc,itc,keyHL,gtsfile,software,dtesHL,dhighlvl):
    # is calculation needed?
    if keyHL in dhighlvl.keys(): return None
    # folder for calculation and file name
    TMP   = PN.TMPHLi%keyHL
    mname = "%s.%s"%(keyHL,"sp")
    # Function for calculation and template
    spc_fnc = get_spc_fnc(software)
    tes     = dtesHL.get(software,{}).get(ctc,None)
    if tes is None: raise Exc.NoTemplateGiven(Exception)
    # extra-arguments for spc_fnc
    clean = False # do not delete files
    bhess = False # do not calculate hessian (no needed actually)
    eargs = (tes, TMP, clean)
    # Read gts file (Molecule instance)
    molecule = Molecule()
    molecule.set_from_gts(gtsfile)
    # HL-calculation
    out_spc = spc_fnc(molecule._xcc,molecule._symbols,bhess,mname,eargs)
    return out_spc[4]
#---------------------------------------------------------------#
def auto_points(drst,allpoints,nbw=3,nfw=3):
    '''
    allpoints must be sorted!!
    '''
    # extension of path
    sbw = drst[allpoints[ 0]][0]
    sfw = drst[allpoints[-1]][0]
    # select points
    deltabw = sbw/float(nbw)
    deltafw = sfw/float(nfw)
    svals = [0.0]
    if deltabw != 0.0: svals += [(ii+1)*deltabw for ii in range(nbw)]
    if deltafw != 0.0: svals += [(ii+1)*deltafw for ii in range(nfw)]
    svals.sort()
    # get closer value in drst
    selected = []
    idx1     = 0
    for s_i in svals:
        mdiff = float("inf")
        l_i   = None
        for idx2 in range(idx1,len(allpoints)):
            l_j  = allpoints[idx2]
            s_j  = drst[l_j][0]
            diff = abs(s_j-s_i)
            if diff < mdiff:
               l_i   = l_j
               idx1  = idx2-1
               mdiff = diff
            if diff > mdiff: break
        selected.append(l_i)
    return selected
#---------------------------------------------------------------#
def highlevel_mep(ctc,itc,keyHL_sp,software,dtesHL,dhighlvl,points):
    # key for dhighlvl
    keyHL_path = "%s.%s.path"%(ctc,itc)
    # initialize dictE
    dictE = {}
    # folder for calculation and file name
    TMP   = PN.TMPHLi%PN.struckey(ctc,itc)
    # template and function for calculation
    tes     = dtesHL.get(software,{}).get(ctc,None)
    spc_fnc = get_spc_fnc(software)
    # extra-arguments for spc_fnc
    clean = False # do not delete files
    bhess = False # do not calculate hessian (no needed actually)
    eargs = (tes, TMP, clean)
    # rst file
    rstfile = PN.get_rst(ctc,itc)
    tpath2, tcommon2, drst = ff.read_rst(rstfile)
    if drst == {}:
       yield None,keyHL_path,(None,None),None,"emptyrst"
       return
    ch,mtp,atonums,masses,mu = tcommon2
    symbols = fncs.atonums2symbols(atonums)
    # define points
    allpoints = sorted_points(drst,hess=False)
    for point in points:
        if point.startswith("auto"):
           auto,nbw,nfw = point.split("_")
           points = auto_points(drst,allpoints,int(nbw),int(nfw))
           break
    # loop in points
    for point_i in points:
        # name for files
        mname = "%s.%s.%s"%(ctc,itc,point_i)
        # get label
        label_i, s_i = find_label_in_rst(point_i,drst)
        # tuple for identify
        tuple_ID = (s_i,label_i)
        # not in drst?
        if label_i is None:
           yield point_i,keyHL_path,tuple_ID,None,"noinrst"
        # already calculated
        elif label_i in dhighlvl.get(keyHL_path,{}).keys(): 
           yield point_i,keyHL_path,tuple_ID,dhighlvl[keyHL_path][label_i],"already"
        # HL-calculation
        else:
           if s_i == 0.0:
              Ehl_sp   = dhighlvl.get(keyHL_sp,None)
              if Ehl_sp is not None:
                 yield point_i,keyHL_path,tuple_ID,Ehl_sp,"saddlefromsp"
                 continue
           # get xcc from rst
           xms = drst[label_i][2]
           xcc = fncs.ms2cc_x(xms,masses,mu)
           # calculation
           try:
              out_spc = spc_fnc(xcc,symbols,bhess,mname,eargs)
              yield point_i,keyHL_path,tuple_ID,out_spc[4],"calculated"
           except:
              yield point_i,keyHL_path,tuple_ID,None,"fail"
#===============================================================#



#===============================================================#
def main(idata,status,case,targets="*"):

    stat2check = [1,4,7]
    mustexist  = [PN.DIR1]
    tocreate   = [PN.DIR2,PN.TMP]
    #-------------------------------------------------------#
    # Read Pilgrim input files, check file/folder status    #
    # and expand tuple 'case'                               #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, tkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,stat2check)
    if fstatus == -1: exit()
    # existency of folders
    fstatus = ffchecking(mustexist,tocreate)
    if fstatus == -1: exit()
    # expand case
    (dof,hlf,plotfile),dlevel,software = case
    #-------------------------------------------------------#


    print("   Selected software: %s"%software)
    print("")

    targets = get_targets(targets,dctc)
    # Print targets
    if len(targets) != 0:
       print("   High-level (HL) calculations will be carried out for:")
       ml = max([len(target) for target in targets]+[1])
       for target in targets:
           ctc,itc = PN.name2data(target)
           TMP     = PN.TMPHLi%PN.struckey(ctc,itc)
           print("       %s (TMP folder: %s)"%("%%-%is"%ml%target,TMP))
       print("")
    else:
       print("   No valid target(s) for High-level (HL) calculations")
       print("")
       return

    # read high-level output file
    dhighlvl = RW.read_highlevelfile(hlf)

    # loop over each target
    print("   Carrying out HL-calculations:")
    print("")
    ml   = max([len(target) for target in targets])
    cols = ("target","point identifier","s value (bohr)","E_HL (hartree)")
    lineformat = " %%%is | %%%is | %%%is | %%%is "%(ml,16,14,16)
    thead = lineformat%cols
    tdivi = "-"*len(thead)
    print("       " + tdivi)
    print("       " + thead)
    print("       " + tdivi)
    for target in targets:
        ctc,itc = PN.name2data(target)
        sptype  = dctc[ctc]._type
        perform = False
        # the important files
        gtsfile = dctc[ctc].gtsfile(itc)
        rstfile = PN.get_rst(ctc,itc)
        # key for dlevel
        key = PN.struckey(ctc,itc)
        if key in ddlevel.keys():
           perform,points = True, ddlevel[key]
        else:
           perform,points = False, None
        if not perform:
           dataline = (target,"","","")
           print("       " + lineformat%dataline)
           continue
        #----------------------------#
        # carry on calculation of SP #
        #----------------------------#
        if os.path.exists(gtsfile):
           # key
           keyHL_sp = PN.struckey(dctc[ctc]._root,itc)
           try:
              Ehl_sp   = highlevel_sp(ctc,itc,keyHL_sp,gtsfile,software,dtesHL,dhighlvl)
              # calculation was performed
              if Ehl_sp is not None:
                 dhighlvl[keyHL_sp] = Ehl_sp
                 dataline = (target,"sp","","%+14.8f"%Ehl_sp)
                 extra    = ""
              # calculation was previously carried out
              else:
                 dataline = (target,"sp","","%+14.8f"%dhighlvl[keyHL_sp])
                 extra    = " *"
           # calculation failed
           except:
              dataline = (target,"sp","","calc. failed")
              extra    = ""
        else: dataline = (target,"sp","","NO GTS")
        print("       " + lineformat%dataline + extra)
        #------------------------------------#
        # carry on calculation of MEP points #
        #------------------------------------#
        if sptype == 0 or len(points) == 0: continue
        if os.path.exists(rstfile):
           idata = (ctc,itc,keyHL_sp,software,dtesHL,dhighlvl,points)
           dictE = {}
           for data_mep_point in highlevel_mep(*idata):
               point_i,keyHL,(s_i,label_i),Ehl,case = data_mep_point
               if case == "calculated":
                  dataline = (target,label_i,"%8.3f"%s_i,"%+14.8f"%Ehl)
                  extra    = ""
                  dictE[label_i] = Ehl
               elif case == "saddlefromsp":
                  dataline = (target,label_i,"%8.3f"%s_i,"%+14.8f"%Ehl)
                  extra    = " (from sp)"
                  dictE[label_i] = Ehl
               elif case == "already":
                  dataline = (target,label_i,"%8.3f"%s_i,"%+14.8f"%Ehl)
                  extra    = " *"
               elif case == "noinrst":
                  # point_i is a value
                  try   : dataline = (target,"","%8.3f"%float(point_i),"NOT IN RST")
                  # point_i is a label
                  except: dataline = (target,"%s"%point_i,"","NOT IN RST")
                  extra    = ""
               elif case == "emptyrst":
                  dataline = (target,   "","", "EMPTY RST")
                  extra    = ""
               elif case == "fail":
                  dataline = (target,label_i,"%8.3f"%s_i,"calc. failed")
                  extra    = ""
               print("       " + lineformat%dataline + extra)
           # update dhighlvl
           dhighlvl[keyHL] = dhighlvl.get(keyHL,{})
           dhighlvl[keyHL].update(dictE)
        else:
           dataline = (target,"path","","NO RST")
           print("       " + lineformat%dataline)
    print("       " + tdivi)
    print("       " + "lines with * ==> not calculated: in file")
    print("")
    print("   (Re)writing '%s' file..."%hlf)
    RW.write_highlevelfile(hlf,dhighlvl)
    print("")
#===============================================================#


