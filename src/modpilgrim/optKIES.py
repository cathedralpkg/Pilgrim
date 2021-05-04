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
| Module     :  modpilgrim         |
| Sub-module :  optKIES            |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --kies option
of Pilgrim
'''

#--------------------------------------------------#
#import datetime
#import time
#import os
#import sys
#--------------------------------------------------#
import numpy             as     np
import modpilgrim.pilrw   as     RW
import modpilgrim.names  as     PN
#--------------------------------------------------#
#import common.Exceptions as Exc
#--------------------------------------------------#
from   common.physcons   import ML
from   common.physcons   import SECOND
from   common.physcons   import NA
from   common.physcons   import KB
from   common.fncs       import print_string
from   common.fncs       import exp128
#--------------------------------------------------#
from   modpilgrim.diverse           import ffchecking
from   modpilgrim.diverse           import dlevel_to_files
from   modpilgrim.diverse           import status_check
from   modpilgrim.diverse           import get_contributions
from   modpilgrim.diverse           import get_transmissioncoeffs
from   modpilgrim.ChemReaction      import ChemReaction
import modpilgrim.strings as PS
#--------------------------------------------------#


#==================================================#
def ask_for_reactions(dchem):
    bool_end   = False
    bool_cont  = False
    reactionH = None
    reactionD = None

    sinit = " "*5

    # ask for reactions
    print(sinit+"Introduce reactions without/with isotopic substitution(s)")
    print(sinit+"Type 'end()' or 'exit()' to finish")

    for question in range(2):

        if question == 0: answer = input(sinit+" >> without: ").strip()
        if question == 1: answer = input(sinit+" >> with   : ").strip()
       #if question == 0: answer = "R1h.fw"
       #if question == 1: answer = "R1d.fw"
        nr     = len(answer.split())

        # User end program?
        if "end()" in answer or "exit()" in answer:
            bool_end = True
            return bool_end,bool_cont,reactionH,reactionD
        # assert only one reaction
        elif nr == 0:
            print("    A reaction must be introduced!")
            bool_cont = True
            return bool_end,bool_cont,reactionH,reactionD
        elif nr != 1:
            print("    Only one reaction must be introduced!")
            bool_cont = True
            return bool_end,bool_cont,reactionH,reactionD
        # Reactions present the dot?
        elif "." not in answer:
            print("    Do not forget to include the dot (.)!")
            bool_cont = True
            return bool_end,bool_cont,reactionH,reactionD
        # Reaction defined?
        else:
           rcname, direction = answer.split(".")[0:2]
           if rcname not in dchem.keys():
              print("    Reaction %s is not defined!"%answer)
              bool_cont = True
              return bool_end,bool_cont,reactionH,reactionD

        # Save reaction
        if question == 0: reactionH = answer
        if question == 1: reactionD = answer

    # Return data
    return bool_end,bool_cont,reactionH,reactionD
#==================================================#



#==================================================#
def get_data_reactants(chemreac,direction):
    if direction == "fw": return chemreac.obtain_pfnparts("reactants")
    if direction == "bw": return chemreac.obtain_pfnparts("products")
#--------------------------------------------------#
def get_anh_and_k(chemreac,direction):
    if direction == "fw":
       anh_k = chemreac._ANHkfw
       dk    = chemreac._kfw
    if direction == "bw":
       anh_k = chemreac._ANHkbw
       dk    = chemreac._kbw
    return anh_k, dk
#--------------------------------------------------#
def main_data(chemreacH,chemreacD,dirH,dirD):
    QtrH_R, QrvH_R, QelH_R, V1H_R = get_data_reactants(chemreacH,dirH)
    QtrD_R, QrvD_R, QelD_R, V1D_R = get_data_reactants(chemreacD,dirD)

    anh_kH, dkH = get_anh_and_k(chemreacH,dirH)
    anh_kD, dkD = get_anh_and_k(chemreacD,dirD)

    QtrH_TS, QrvH_TS, QelH_TS, V1H_TS = chemreacH.obtain_pfnparts("ts")
    QtrD_TS, QrvD_TS, QelD_TS, V1D_TS = chemreacD.obtain_pfnparts("ts")

    return QtrH_R, QrvH_R, QelH_R, V1H_R,\
           QtrD_R, QrvD_R, QelD_R, V1D_R,\
           anh_kH, dkH,\
           anh_kD, dkD,\
           QtrH_TS, QrvH_TS, QelH_TS, V1H_TS,\
           QtrD_TS, QrvD_TS, QelD_TS, V1D_TS
#--------------------------------------------------#
def contrib_tr_rv_el_tor(chemreacH,chemreacD,dirH,dirD):
    # get data
    QtrH_R, QrvH_R, QelH_R, V1H_R,\
    QtrD_R, QrvD_R, QelD_R, V1D_R,\
    anh_kH, dkH,\
    anh_kD, dkD,\
    QtrH_TS, QrvH_TS, QelH_TS, V1H_TS,\
    QtrD_TS, QrvD_TS, QelD_TS, V1D_TS = main_data(chemreacH,chemreacD,dirH,dirD)

    # calculate contributions
    ltemp = chemreacH._ltemp
    dE = (V1H_TS-V1H_R)-(V1D_TS-V1D_R)
    vE = np.array([exp128(-dE/KB/T) for T in ltemp])
    kies_tr = (QtrH_TS/QtrH_R) / (QtrD_TS/QtrD_R)
    kies_rv = (QrvH_TS/QrvH_R) / (QrvD_TS/QrvD_R) * vE
    kies_el = (QelH_TS/QelH_R) / (QelD_TS/QelD_R)
    kies_tor = anh_kH/anh_kD
    return kies_tr,kies_rv,kies_el,kies_tor
#--------------------------------------------------#
def contrib_vtun_tot(chemreacH,chemreacD,dirH,dirD,kies_tr,kies_rv,kies_tor):
    # get data
    QtrH_R, QrvH_R, QelH_R, V1H_R,\
    QtrD_R, QrvD_R, QelD_R, V1D_R,\
    anh_kH, dkH,\
    anh_kD, dkD,\
    QtrH_TS, QrvH_TS, QelH_TS, V1H_TS,\
    QtrD_TS, QrvD_TS, QelD_TS, V1D_TS = main_data(chemreacH,chemreacD,dirH,dirD)

    # calculate contributions
    kies_vtun    = {}
    kies_tot     = {}
    ltemp = chemreacH._ltemp
    for X in dkH.keys():
        if dkH[X] is None: continue
        elif X == "tst": kies_vtun[X] = np.array([1.0 for T in ltemp])
        else           : kies_vtun[X] = (dkH[X]/dkH["tst"]) / (dkD[X]/dkD["tst"])
        kies_tot[X] = kies_tr*kies_rv*kies_vtun[X]*kies_tor
    return kies_vtun,kies_tot
#--------------------------------------------------#
def contrib_ind_rv_vtun(chemreacH,chemreacD,dirH,dirD,kies_rv):
    # get data
    QtrH_R, QrvH_R, QelH_R, V1H_R,\
    QtrD_R, QrvD_R, QelD_R, V1D_R,\
    anh_kH, dkH,\
    anh_kD, dkD,\
    QtrH_TS, QrvH_TS, QelH_TS, V1H_TS,\
    QtrD_TS, QrvD_TS, QelD_TS, V1D_TS = main_data(chemreacH,chemreacD,dirH,dirD)

    # get more data
    dchiH = chemreacH._tschi["tst"]
    dchiD = chemreacD._tschi["tst"]
    dtcH  = chemreacH._dtcoef
    dtcD  = chemreacD._dtcoef
    ltemp = chemreacH._ltemp
    itc0  = chemreacH._tsitc0

    # calculate individual rovib
    ikies_rv = {}
    for itc in sorted(dchiH.keys()):
        chiH_j = np.array(dchiH[itc])
        chiD_j = np.array(dchiD[itc])
        ikies_rv[itc] = chiH_j/chiD_j * np.array(kies_rv)

    # calculate individual vtun
    ikies_vtun = {}
    for X in dkH.keys():
        if X == "tst": continue
        ikies_vtun[X] = {}
        for itc in sorted(dchiH.keys()):
            if   X.startswith("ms"): XG,itcG = X[2:], itc0
            elif X.startswith("mp"): XG,itcG = X[2:], itc
            else                   : continue
            gammaH_j = np.array(dtcH[XG][itcG])
            gammaD_j = np.array(dtcD[XG][itcG])
            ikies_vtun[X][itc] = np.array(gammaH_j) / np.array(gammaD_j)
    ikies_vtun["tst"] = {itc : np.array([1.0 for T in ltemp]) for itc in dchiH.keys()}

    return ikies_rv, ikies_vtun
#--------------------------------------------------#
def contrib_pjd_pjh(chemreacH,chemreacD,dirH,dirD,kies_rv,kies_vtun,ikies_rv,ikies_vtun):

    # get data
    QtrH_R, QrvH_R, QelH_R, V1H_R,\
    QtrD_R, QrvD_R, QelD_R, V1D_R,\
    anh_kH, dkH,\
    anh_kD, dkD,\
    QtrH_TS, QrvH_TS, QelH_TS, V1H_TS,\
    QtrD_TS, QrvD_TS, QelD_TS, V1D_TS = main_data(chemreacH,chemreacD,dirH,dirD)

    # get more data
    dchiH = chemreacH._tschi["tst"]
    dchiD = chemreacD._tschi["tst"]
    dtcH  = chemreacH._dtcoef
    dtcD  = chemreacD._dtcoef
    ltemp = chemreacH._ltemp
    itc0  = chemreacH._tsitc0

    pjh   = {X:{} for X in dkH.keys()}
    pjd   = {X:{} for X in dkD.keys()}
    for X in dkD.keys():
        # average gamma_D
        if   X.startswith("ms"): gammaD_ave = np.array(dtcD[X[2:]][itc0])
        elif X.startswith("mp"): gammaD_ave = np.array(dtcD[X[2:]]["averaged"])
        else                   : gammaD_ave = np.array([1.0 for T in ltemp])
        # calculate P_j,D & P_j,H
        for itc in sorted(dchiD.keys()):
            chiD_j   = np.array(dchiD[itc])
            # gammaD_j
            if   X.startswith("ms"): gammaD_j = np.array(dtcD[X[2:]][itc0])
            elif X.startswith("mp"): gammaD_j = np.array(dtcD[X[2:]][itc])
            else                   : gammaD_j = np.array([1.0 for T in ltemp])
            # pjD
            if   X.startswith("ms"): pjd[X][itc] = chiD_j
            elif X.startswith("mp"): pjd[X][itc] = chiD_j * gammaD_j / gammaD_ave
            else                   : pjd[X][itc] = chiD_j
            # pjH
            if   X.startswith("ms"): vtun = np.array([1.0 for T in ltemp])
            elif X.startswith("mp"): vtun = ikies_vtun[X][itc] / np.array(kies_vtun[X])
            else                   : vtun = np.array([1.0 for T in ltemp])
            rv = ikies_rv[itc] / np.array(kies_rv)
            pjh[X][itc] = pjd[X][itc] * rv * vtun
    return pjh,pjd
#--------------------------------------------------#
def contrib_totikies(pjh,pjd,kies_tot):
    lX = pjh.keys()
    totkies_j   = {X:{} for X in lX}
    totkiest_j  = {X:{} for X in lX}
    for X in lX:
        for itc in sorted(pjh[X].keys()):
            totkiest_j[X][itc] = pjh[X][itc]*kies_tot[X]
            totkies_j[ X][itc] = totkiest_j[X][itc] / pjd[X][itc]
    return totkies_j,totkiest_j
#==================================================#


#==================================================#
def kies_from_pair_of_reactions(reactionH,reactionD,ltemp,dchem,dctc,dall):

    # a) Deal with reaction without isotopes
    rcnameH, dirH = reactionH.split(".")[0:2]
    RsH,TSH,PsH = dchem[rcnameH]
    chemreacH = ChemReaction(rcnameH,ltemp,dctc)
    chemreacH.external_data(dall)
    chemreacH.add_reactant(RsH)
    chemreacH.add_products(PsH)
    chemreacH.add_ts(TSH)


    # b) Deal with reaction with isotopes
    rcnameD, dirD = reactionD.split(".")[0:2]
    RsD,TSD,PsD = dchem[rcnameD]
    chemreacD = ChemReaction(rcnameD,ltemp,dctc)
    chemreacD.external_data(dall)
    chemreacD.external_data(dall)
    chemreacD.add_reactant(RsD)
    chemreacD.add_products(PsD)
    chemreacD.add_ts(TSD)

    # c) chi_i's, transmission coefficients, anharmonicity, rate constants
    chemreacH.obtain_pfns()
    chemreacD.obtain_pfns()
    chemreacH.calculate_transcoeffs()
    chemreacD.calculate_transcoeffs()
    chemreacH.calculate_anharmonicity()
    chemreacD.calculate_anharmonicity()
    chemreacH.calculate_rateconstants()
    chemreacD.calculate_rateconstants()

    #lXH = chemreacH._tschi.keys()
    #lXD = chemreacD._tschi.keys()

    # d) translational, rovibrational and torsional KIES
    print_string(PS.skies_definitions(),5)
    kies_tr,kies_rv,kies_el,kies_tor = contrib_tr_rv_el_tor(chemreacH,chemreacD,dirH,dirD)
    print_string(PS.skies_basiccontris(ltemp,kies_tr,kies_rv,kies_tor),5)

    # e) vtun, total and total with anharmonicity
    kies_vtun, kies_tot = contrib_vtun_tot(chemreacH,chemreacD,dirH,dirD,kies_tr,kies_rv,kies_tor)
    print_string(PS.skies_vtun_tot(ltemp,kies_vtun,kies_tot),5)

    # f) individual contributions
    ikies_rv, ikies_vtun = contrib_ind_rv_vtun(chemreacH,chemreacD,dirH,dirD,kies_rv)
    pjh, pjd             = contrib_pjd_pjh(chemreacH,chemreacD,dirH,dirD,\
                           kies_rv,kies_vtun,ikies_rv,ikies_vtun)
    totkies_j,totkiest_j = contrib_totikies(pjh,pjd,kies_tot)

    print_string(PS.skies_indkies(ltemp,ikies_rv,ikies_vtun,pjh,pjd,totkies_j,totkiest_j),5)
            
#==================================================#
      

def main(idata,status,case):

    stat2check = [2,5]
    mustexist  = []
    tocreate   = []
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

    # read dofs
    dall, ltemp = RW.read_alldata(dof)
    # get saved reactions
    print_string(PS.skies_summary(dchem,dall),5)

    # Ask user which reactions
    count = 0
    while True:
         # Ask user for reactions
         bool_end,bool_cont,reaction1,reaction2 = ask_for_reactions(dchem)
         print("")
         print("")
         if bool_end : break
         if bool_cont: continue
         kies_from_pair_of_reactions(reaction1,reaction2,ltemp,dchem,dctc,dall)
         #exit()


