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
| Sub-module :  optRCONS           |
| Last Update:  2020/03/01 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#--------------------------------------------------#
import time
import datetime
import os
import sys
import numpy   as np
#--------------------------------------------------#
import common.Exceptions as Exc
import common.partfns    as partfns
#--------------------------------------------------#
from   common.fncs       import time2human
from   common.fncs       import prod_list
from   common.fncs       import print_string
import common.physcons as pc
from   common.physcons   import KB, ML, H, AMU
from   common.Logger     import Logger
from   common.Molecule   import Molecule
#--------------------------------------------------#
from modpilgrim.ChemReaction import ChemReaction
import modpilgrim.names   as PN
import modpilgrim.pilrw    as RW
import modpilgrim.strings as PS
#--------------------------------------------------#
from   modpilgrim.fit2anarc         import fit2anarc
from   modpilgrim.diverse           import ffchecking
from   modpilgrim.diverse           import get_input_data
from   modpilgrim.diverse           import status_check
from   modpilgrim.diverse           import get_transmissioncoeffs
from   modpilgrim.plotting          import manage_data_for_plot_rcons
from   modpilgrim.plotting          import write_plotfile
#--------------------------------------------------#



#==========================================================#
def prepare_plot_data(chemreac):
    plotdata = {}
    for direc in "fw,bw".split(","):
        if direc == "fw": nR, ks = chemreac._nR, chemreac._kfw
        if direc == "bw": nR, ks = chemreac._nP, chemreac._kbw
        # human units: string and value
        hu = pc.ML**(nR-1.0) / pc.SECOND
        if nR != 1: units = "s^{-1} (cm^3)^{%i} molecule^{-%i}"%(nR-1,nR-1)
        else      : units = "s^{-1}"
        # generate plotdata
        d4plot = {X:([ki*hu for ki in k],{}) for X,k in ks.items() if k is not None}
        args   = (chemreac._rcname,direc,chemreac._ltemp,d4plot,units)
        # update plotdata dictionary
        plotdata.update( manage_data_for_plot_rcons(*args) )
    return plotdata
#==========================================================#



#===============================================================#
def main(idata,status,case,targets="*"):

    stat2check = [1,2,5]
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



    # no specific target selected
    if "*" in targets or len(targets) == 0: targets = dchem.keys()

    # clean targets
    unknown = sorted([target for target in targets if target not in dchem.keys()])
    targets = sorted([target for target in targets if target in dchem.keys()])

    if unknown != []:
       print_string("The following reactions are not defined:",3)
       for unkn in unknown: print_string("* %s"%unkn,6)
       print("")

    # read dof
    dall = RW.read_alldata(dof,ltemp)[0]

    #--------------------------------#
    # Calculations for each reaction #
    #--------------------------------#
    for rcname in targets:
        # pof file
        t1 = time.time()
        pof = PN.get_pof(dlevel,"rcons",rcname)
        sys.stdout = Logger(pof,"w",True)
        sys.stdout.writeinfile(PS.init_txt())

        # preparation of reaction
        Rs, TS, Ps = dchem[rcname]
        print_string(PS.srcons_init(rcname,pof,Rs,TS,Ps),3)

        # Create instante of ChemReaction
        chemreac = ChemReaction(rcname,ltemp,dctc,dimasses)
        chemreac.external_data(dall)
        chemreac.add_reactant(Rs)
        chemreac.add_products(Ps)
        chemreac.add_ts(TS)
        chemreac.read_gtsfiles()
        # conservation of mass and charge
        chemreac.check_conservation()
        print_string(PS.srcons_conservation(chemreac),5)
        if chemreac._problem:
           print_string("ERROR: problem with conservation of properties!",5)
           print("")
           continue

        # print energies
        chemreac.obtain_pfns()
        print_string(PS.srcons_relenergies(chemreac),5)

        # get partition functions and anharmonicity
        chemreac.calculate_anharmonicity()
        print_string(PS.srcons_anhtable(chemreac),5)

        # calculate equilibrium constant
        chemreac.calculate_eqconstant()
        print_string(PS.srcons_keq(chemreac),5)

        # calculate transmission coefficients
        chemreac.calculate_transcoeffs()
        print_string(PS.srcons_indivtcoefs(chemreac),5)
        print_string(PS.srcons_avertcoefs(chemreac),5)
        print_string(PS.srcons_tscontributions(chemreac),5)

        # (a) forward  rate constants
        if chemreac._ts is None: continue
        chemreac.calculate_rateconstants()
        print_string(PS.srcons_rateconstants(chemreac,"fw"),5)
        print_string(PS.srcons_actgibbs(chemreac,"fw"),5)

        # (b) backward rate constants
        print_string(PS.srcons_rateconstants(chemreac,"bw"),5)
        print_string(PS.srcons_actgibbs(chemreac,"bw"),5)

        # Update dall
        dall = RW.read_alldata(dof,ltemp)[0]

        for X,ks in chemreac._kfw.items():
            if ks is None: continue
            dall["rcons"]["%s.%s.fw"%(X,rcname)] = ks
        for X,ks in chemreac._kbw.items():
            if ks is None: continue
            dall["rcons"]["%s.%s.bw"%(X,rcname)] = ks

        for X in chemreac._dtcoef.keys():
            av_tcoef = chemreac._dtcoef[X]["averaged"]
            if av_tcoef is None: continue
            dall["cfs"]["%s.%s"%(X,rcname)] = av_tcoef
        RW.write_alldata(dof,ltemp,dall)

        # Prepare plotdata
        plotdata = prepare_plot_data(chemreac)
        if plotfile is not None and plotdata != {}:
            print_string("Updating plot file: %s"%plotfile,3)
            print("")
            write_plotfile(plotfile,plotdata)


        # Print special case, when Rs == Ps
        if (sorted(chemreac._reacts) == sorted(chemreac._prods)) and (chemreac._ts is not None):
           print_string(PS.srcons_rconsX2(chemreac),5)

        # print end of file and finish logger
        print("")
        print("")
        t2    = time.time()
        etime = time2human(t2-t1,"secs")
        sys.stdout.writeinfile(PS.end_txt(*etime))
        sys.stdout = Logger(None)
#===============================================================#

