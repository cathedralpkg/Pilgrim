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
| Sub-module :  optPFN             |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --pfn option
of Pilgrim
'''

#--------------------------------------------------#
import time
import os
import sys
import numpy   as np
#--------------------------------------------------#
import modpilgrim.names   as PN
import modpilgrim.pilrw    as RW
import modpilgrim.strings as PS
#--------------------------------------------------#
from   common.fncs       import exp128
from   common.fncs       import print_string
from   common.fncs       import time2human
from   common.fncs       import add_iblank
from   common.fncs       import pressure2vol
from   common.physcons   import KB, AMU, KCALMOL, ML, KB_SI
from   common.physcons   import PRE0, VOL0
from   common.Logger     import Logger
from   common.Molecule   import Molecule
from   common.criteria   import EPS_AMU
from   common.criteria   import EPS_TEMP

import common.Exceptions as Exc
from   common.files      import read_file
from   common.files      import read_q2dtorout
from   common.files      import read_mstorout
#--------------------------------------------------#
from   modpilgrim.diverse         import ffchecking
from   modpilgrim.diverse         import status_check
from   modpilgrim.diverse         import get_input_data
#--------------------------------------------------#
from   modpilgrim.plotting   import manage_data_for_plot_weights
from   modpilgrim.plotting   import write_plotfile
#--------------------------------------------------#


WARNINGS = []

#===============================================================#
def check_ctc(targets,dctc):
    updated = []
    for ctc in targets:
        if ctc not in dctc.keys():
           print_string("Selected target (%s) not found as a STRUC. Skipped!"%ctc,2)
           print("")
           continue
        updated.append(ctc)
    return updated
#---------------------------------------------------------------#
def molecule_prep(gts,eslist=[],tiso=None,fscal=1.0):
    # generate instance of Molecule
    molecule = Molecule()
    # read gts
    molecule.set_from_gts(gts)
    # masses
    if tiso is not None:
       imods,imasses = tiso
       molecule.apply_imods(imods,imasses)
    # electronic states
    molecule.setvar(les=eslist)
    # set frequency scaling
    molecule.setvar(fscal=fscal)
    # setup
    molecule.setup()
    molecule.ana_freqs("cc")
    return molecule
#---------------------------------------------------------------#
def molecule_pfns(molecule,ltemp,sptype=0):
    # reference energies
    V0  = molecule._V0
    zpe = molecule._cczpe
    # calculate partition functions
    pf_tot, V1, (pf_PIB,pf_RR,pf_HO,pf_ele) = molecule.calc_pfns(ltemp,"cc",fmode=-sptype)
    # return data
    return V0, V1, pf_tot, (pf_PIB,pf_RR,pf_HO,pf_ele)
#---------------------------------------------------------------#
def conformer_contributions(ltemp,dpfn,ctc,itcs):
    # get contributions
    dchi = {}
    for itc,weight in itcs:
        V0i,V1i,Qi = dpfn[PN.struckey(ctc,itc)]
        V0 ,V1 ,Q  = dpfn[PN.struckey(ctc,"msho")]
        # calculate contribution
        expvec = [exp128(-(V1i-V1)/KB/T) for T in ltemp]
        chi    = weight * np.array(Qi) / np.array(Q) * np.array(expvec)
        # save data
        dchi[itc] = chi
    return dchi
#---------------------------------------------------------------#
def get_V0HL_list(ctc,itcs,dhighlvl):
    list_V0HL = []
    for itc,weight in itcs:
        key  = PN.struckey(ctc,itc)
        V0HL = dhighlvl.get(key,None)
        list_V0HL.append(V0HL)
    return list_V0HL
#---------------------------------------------------------------#
def calculate_QMSHO_for_ctc(Cluster,ltemp,dimasses={},tupleHL=None):
    global WARNINGS
    # expand data
    ctc          = Cluster._ctc
    itcs         = Cluster._itcs
    eslist       = Cluster._es
    sptype       = Cluster._type
    fscal        = Cluster._fscal
    diso         = Cluster._diso
    if tupleHL is None: dlevel, dhighlvl = False, {}
    else              : dlevel, dhighlvl = tupleHL
    # get list of gts files and check existency
    gtsfiles = Cluster.gtsfiles()
    exception = Exc.LackOfGts(Exception)
    for gts in gtsfiles:
        if not os.path.exists(gts):
           exception._var = gts
           raise exception
    # DLEVEL
    if dlevel: V0HLs = get_V0HL_list(Cluster._root,itcs,dhighlvl)
    else     : V0HLs = [None for itc,weight in itcs]
    # Generate Molecule instances
    molecules = []
    mass = None
    V0LLs = []

    for idx,gts in enumerate(gtsfiles):
        itc,weight = itcs[idx]
        # tuple with isotopic
        if   itc in diso.keys(): iso = diso[itc]
        elif "*" in diso.keys(): iso = diso["*"]
        else                   : iso = None
        tiso = (iso,dimasses)
        # prepare and add molecule
        molecule = molecule_prep(gts,eslist,tiso,fscal)
        # save low-level
        V0LLs.append( float(molecule._V0) )
        # apply dlevel
        if V0HLs[idx] is not None:
           molecule._V0 = V0HLs[idx]
           molecule._ccV1 = molecule._V0 + molecule._cczpe
        # save molecule
        molecules.append(molecule)
    # Get min of V0 and V1
    V0 = min([molecule._V0   for molecule in molecules])
    V1 = min([molecule._ccV1 for molecule in molecules])
    # String
    string = ""
    if dlevel:
       string += "    keyword --dlevel activated: applying High-Level energies to this STRUC\n"
       string += "\n"
       string += "        * Low-level (LL) and high-level (HL) energies (in hartree):\n"
       string += "\n"
       string += "          ---------------------------------------------------\n"
       string += "           conformer |     LL energy     |     HL energy     \n"
       string += "          ---------------------------------------------------\n"
       for idx,(itc,weight) in enumerate(itcs):
           if V0HLs[idx] is None: V0HL = "        -        "
           else                 : V0HL = "%+17.8f"%V0HLs[idx]
           string += "           %-9s | %+17.8f | %s \n"%(itc,V0LLs[idx],V0HL)
       string += "          ---------------------------------------------------\n"
       string += "\n"
       if None in V0HLs:
          string += "        * FAIL!! Not enough data... DUAL-LEVEL will be omitted!\n"
          string += "\n"
          WARNINGS.append("Dual-level was not applied to %s..."%ctc)
    string += add_iblank(PS.getstring_ctc(sptype,molecules,itcs,V0,V1),4)
    # Compare masses
    for molecule in molecules:
        if mass is None: mass = molecule._mass
        # compare masses
        if abs(mass-molecule._mass)*AMU > EPS_AMU:
           WARNINGS.append("Mass inconsistences in %s"%ctc)
           exception =  Exc.DiffMassConf(Exception)
           exception._var = string
           raise exception
    # get partition functions for each molecule (and for CTC)
    dctr       = {}
    dpfn       = {}
    dgibbs1cm3 = {}
    dgibbs1bar = {}
    pfn_tot    = [0.0 for T in ltemp]
    #pfn_gibbs  = [0.0 for T in ltemp] # goes without sigma_rot
    for idx,molecule in enumerate(molecules):
        itc,weight = itcs[idx]
        V0_i, V1_i, pf_i, tuple_i = molecule_pfns(molecule,ltemp,sptype)
        key = PN.struckey(ctc,itc)
        # save individual pfn and gibbs free energy
        rotsigma    = molecule._rotsigma
        # Considering volume --> 1 particle per cc
        gibbs1cm3_i  = [V1_i-KB*T*np.log(q_j*VOL0                ) for T,q_j in zip(ltemp,pf_i)]
        # Considering pression --> 1bar
        gibbs1bar_i  = [V1_i-KB*T*np.log(q_j*pressure2vol(PRE0,T)) for T,q_j in zip(ltemp,pf_i)]
        # save data
        dpfn[key]   = (V0_i,V1_i,pf_i)
        dgibbs1cm3[key] = gibbs1cm3_i
        dgibbs1bar[key] = gibbs1bar_i
        # update MSHO values
        rel_V1 = V1_i-V1
        for idx,T in enumerate(ltemp):
            beta = 1.0/(T*KB)
            pfn_tot[idx]   += weight*pf_i[idx]*exp128(-rel_V1*beta)
            #pfn_gibbs[idx] += weight*pf_i[idx]*exp128(-rel_V1*beta)*rotsigma
        # add to string
        pf_PIB,pf_RR,pf_HO,pf_ele = tuple_i
        pf_RRHO = [aa*bb for aa,bb in zip(pf_RR,pf_HO)]
        dctr[key] = (pf_PIB,pf_RRHO,pf_ele)
        string += "    "+"-"*len("Conformation: %s "%itc)+"\n"
        string += "    Conformation: %s\n"%itc
        string += "    "+"-"*len("Conformation: %s "%itc)+"\n"
        for line in molecule.info_string().split("\n"):
            string += "     | "+line+"\n"
        #string += add_iblank(molecule.info_string(),5)
        string += "     | Partition functions (pfns):\n"
        for line in PS.getstring_pfn1(ltemp,pf_i,pf_PIB,pf_RR,pf_HO,pf_ele).split("\n"):
            string += "     | "+line+"\n"
       #string += add_iblank(PS.getstring_pfn1(ltemp,pf_i,pf_PIB,pf_RR,pf_HO,pf_ele),7)
        string += "     | Gibbs free energy (hartree):\n"
        for line in PS.spfn_igibbs(ltemp,gibbs1cm3_i,gibbs1bar_i).split("\n"):
            string += "     | "+line+"\n"
       #string += add_iblank(PS.spfn_igibbs(ltemp,gibbs1cm3_i,gibbs1bar_i),7)
        # generate molden file
        molden  = PN.get_gtsmolden(ctc,itc)
        if not os.path.exists(molden):
           if not os.path.exists(PN.DIR5): os.mkdir(PN.DIR5)
           #molecule.genfile_xyz(molden)
           molecule.genfile_molden(molden)
           string += "     | file %s was generated\n"%molden
           string += "\n"
        string += "\n"
    # save Q_MSHO
    dpfn[  PN.struckey(ctc,"msho")] = (V0,V1,pfn_tot)
    pf_PIBele = [aa*bb for aa,bb in zip(pf_PIB,pf_ele)]
    pf_RRHO   = [aa/bb for aa,bb in zip(pfn_tot,pf_PIBele)]
    dctr[ctc] = (pf_PIB,pf_RRHO,pf_ele)
    # save Gibbs_MSHO
    gibbs1cm3_tot = [V1-KB*T*np.log(q_j*VOL0                ) for T,q_j in zip(ltemp,pfn_tot)]
    gibbs1bar_tot = [V1-KB*T*np.log(q_j*pressure2vol(PRE0,T)) for T,q_j in zip(ltemp,pfn_tot)]
    dgibbs1cm3[PN.struckey(ctc,"msho")] = gibbs1cm3_tot
    dgibbs1bar[PN.struckey(ctc,"msho")] = gibbs1bar_tot
    # Add to string
    nconfs = sum([weight for name,weight in itcs])
    if nconfs > 1:
       dchi = conformer_contributions(ltemp,dpfn,ctc,itcs)
       string += "\n"
       string += "    Total multi-structural HO pfn (QMS_HO) and Gibbs free energies (GFE):\n"
      #string += "    Total multi-structural HO pfn (QMS_HO) relative to min(V1):\n"
       string += add_iblank(PS.getstring_pfn2(ltemp,pfn_tot,gibbs1cm3_tot,gibbs1bar_tot),10)
       string += "    Individual contributions to the partition function:\n"
       string += add_iblank(PS.string_contributions(itcs,ltemp,dchi),10)
    else:
       dchi = {}
    return dpfn, dgibbs1cm3, dgibbs1bar, string, dchi, dctr
#---------------------------------------------------------------#
def deal_with_anh(anhfile,ltemp,ZPEMSHO,QMSHO):
    RATIO  = None
    DZPE   = None
    string = ""
    while True:
       if anhfile is None: break
       anhfile = PN.ANHDIR+anhfile
       if not os.path.exists(anhfile):
          string += "Anharmonicity file '%s' NOT FOUND!\n"%anhfile
          string += "\n"
          break
       string += "Anharmonicity file '%s' FOUND!\n"%anhfile
       string += "\n"
       # read lines in file
       lines = read_file(anhfile)
       # type of file??
       itisq2dtor = False
       itismstor  = False
       for line in lines[::-1]:
           if "End of Q2DTor output" in line: itisq2dtor = True; break
           if "End of MSTor output"  in line: itismstor  = True; break
       # get data
       if   itisq2dtor: Tlist, ZPE_MSHO,Q_MSHO,ZPE_ANH, Q_ANH = read_q2dtorout(anhfile)
       elif itismstor : Tlist, ZPE_MSHO,Q_MSHO,ZPE_ANH, Q_ANH = read_mstorout(anhfile)
       else:
          string += "  * unable to identify the program where this file come from...\n"
          string += "\n"
          break
       # print info
       if None in (ZPE_MSHO,Q_MSHO,ZPE_ANH, Q_ANH):
          string += "  * unable to read the file properly...\n"
          string += "\n"
          break
       # print zpe
       string += "  * Data in file:\n"
       string += "\n"
       string += "      ZPE_MSHO: %.8f hartree\n"%ZPE_MSHO
       string += "      ZPE_ANH : %.8f hartree\n"%ZPE_ANH
       string += "\n"
       # compare list of temperatures
       string += "  * checking list of temperatures\n"
       common = []
       for idxT1,T1 in enumerate(Tlist):
           for idxT2,T2 in enumerate(ltemp):
               if abs(T1-T2) < EPS_TEMP: common.append(idxT1)
       string += "    num_temps(%s) = %i\n"%(PN.IFILE2,len(ltemp))
       string += "    num_temps(%s) = %i\n"%(anhfile,len(Tlist))
       string += "    in common --> %i\n"%(len(common))
       if len(common) != len(ltemp):
           string += "    differences in temperatures detected [FAIL]\n"
           break
       elif len(common) != len(Tlist):
           Tlist  = [Tlist[idxT]  for idxT in common]
           Q_MSHO = [Q_MSHO[idxT] for idxT in common]
           Q_ANH  = [Q_ANH[idxT]  for idxT in common]
       string += "    temperatures match! [OK]\n"
       string += "\n"
       # print partition functions
       thead = "  T  (K)  |   Q_MSHO   |   Q_ANH    "
       tdivi = "-"*len(thead)
       string += "      "+tdivi+"\n"
       string += "      "+thead+"\n"
       string += "      "+tdivi+"\n"
       for idx,T in enumerate(Tlist):
           tline = " %8.2f | %10.3E | %10.3E "%(T, Q_MSHO[idx], Q_ANH[idx])
           string += "      "+tline+"\n"
       string += "      "+tdivi+"\n"
       string += "\n"
       # compare ZPE of MSHO
       string += "  * checking MSHO data...\n"
       maxdiff = abs(ZPE_MSHO - ZPEMSHO)*KCALMOL
       if maxdiff < 0.1:
          string += "    abs. diff. in ZPE is %.1E kcal/mol [OK]\n"%maxdiff
       else:
          string += "    abs. diff. in ZPE is %.1E kcal/mol [FAIL]\n"%maxdiff
          break
       # compare pfn of MSHO
       maxdiff = 100*max( [abs(pi-pj)/pi for pi,pj in zip(QMSHO,Q_MSHO)] )
       if maxdiff < 2.0:
          string += "    max. rel. diff in MSHO partition function is %.1E%% [OK]\n"%maxdiff
       else:
          string += "    max. rel. diff in MSHO partition function is %.1E%% [FAIL]\n"%maxdiff
          break
       string += "\n"
       # calculate ratio
       string += "  * Calculating anh. ratio = Q_ANH/Q_MSHO * exp(-(ZPE_ANH-ZPE_MSHO)*beta)\n"
       string += "\n"
       DZPE  = ZPE_ANH - ZPE_MSHO
       RATIO = [(Q_ANH[idx]/Q_MSHO[idx])*exp128(-DZPE/KB/T) for idx,T in enumerate(ltemp)]
       # print partition functions
       thead = "  T  (K)  | ANH. RATIO "
       tdivi = "-"*len(thead)
       string += "      "+tdivi+"\n"
       string += "      "+thead+"\n"
       string += "      "+tdivi+"\n"
       for idx,T in enumerate(ltemp):
           tline = " %8.2f | %10.3E "%(T, RATIO[idx])
           string += "      "+tline+"\n"
       string += "      "+tdivi+"\n"
       string += "\n"
       break
    return RATIO, string

#===============================================================#



#===============================================================#
def main(idata,status,case,targets="*"):
    global WARNINGS

    stat2check = [1,2]
    mustexist  = [PN.DIR1]
    tocreate   = [PN.DIR2,PN.DIR3,PN.DIR6]

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

    # read high level file
    if dlevel: dhighlvl = RW.read_highlevelfile(hlf)
    else     : dhighlvl = {}
    tupleHL = (dlevel,dhighlvl)

    # targets
    if targets == [] or targets == "*": targets = dctc.keys()

    # valid ctc?
    targets = check_ctc(targets,dctc)

    if targets == []: return

    #-------------------------------#
    # Calculation for each compound #
    #-------------------------------#
    dpfn      = {}
    dgibbs1cm3 = {}
    dgibbs1bar = {}
    danh      = {}
    dctr      = {}
    # read dall (raises exception if problem with temperatures)
    dall = RW.read_alldata(dof,ltemp)[0]

    for ctc in sorted(targets):
        NBS = 3
        # write data also in file
        pof = PN.get_pof(dlevel,"pfn",ctc)
        sys.stdout = Logger(pof,"w",True)
        t1 = time.time()
        sys.stdout.writeinfile(PS.init_txt())
        # write title
        print_string(PS.title_pfn(ctc,pof),NBS)
        # Get ClusterConf instance
        Cluster = dctc[ctc]
        # calculation
        dchi = {}
        try:
          idata = (Cluster,ltemp,dimasses,tupleHL)
          dpfn_i, dgibbs1cm3_i, dgibbs1bar_i, string, dchi, dctr_i = calculate_QMSHO_for_ctc(*idata)
        # any error?
        except Exception as exception:
            WARNINGS.append("Problem with %s"%ctc)
            if type(exception) == Exc.DiffMassConf:
               print_string(exception._var,NBS)
               continue
            if type(exception) == Exc.LackOfGts:
               error_line =  "ERROR: Unable to find '%s'"%exception._var
               print_string(error_line,7)
               print("")
               continue
            raise exception
        # print info
        print_string(string,NBS)
        # anharmonicity??
        anhfile = dctc[ctc]._anh
        V0, V1, QMSHO = dpfn_i[PN.struckey(ctc,"msho")]
        ZPEMSHO = V1-V0
        RATIO, string = deal_with_anh(anhfile,ltemp,ZPEMSHO,QMSHO)
        if string != "":
          print_string(string,7)
          if RATIO is None:
             WARNINGS.append("Anharmonicity IGNORED for %s"%ctc)
        # add data of weights to plot
        if plotfile is not None and dchi != {}:
           plotdata = {}
           plotdata.update(manage_data_for_plot_weights(ctc,ltemp,dchi))
           write_plotfile(plotfile,plotdata)
        # update dictionary
        dpfn.update(dpfn_i)
        dgibbs1cm3.update(dgibbs1cm3_i)
        dgibbs1bar.update(dgibbs1bar_i)
        dctr.update(dctr_i)
        if RATIO is not None: danh[ctc] = RATIO
        # print end of file
        t2    = time.time()
        etime = time2human(t2-t1,"secs")
        sys.stdout.writeinfile(PS.end_txt(*etime))
    # write file (read, check ltemp, update, write)
    sys.stdout = Logger(None)
    print_string("Updating data file: %s"%dof,3)
    dall = RW.read_alldata(dof,ltemp)[0]
    dall["pfn"].update(dpfn)
    dall["gibbs1cm3"].update(dgibbs1cm3)
    dall["gibbs1bar"].update(dgibbs1bar)
    dall["anh"].update(danh)
    dall["ctr"].update(dctr)
    RW.write_alldata(dof,ltemp,dall)
    print("")

    # print WARNINGS
    if len(WARNINGS) != 0:
       print_string("WARNINGS:",3)
       for warning in WARNINGS:
          print("      *",warning)
       print("")
#===============================================================#



