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
| Sub-module :  optFIT             |
| Last Update:  2020/03/01 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --fit option
of Pilgrim
'''

#--------------------------------------------------#
import common.Exceptions as Exc
import common.partfns    as partfns
#--------------------------------------------------#
from   common.fncs       import print_string
import common.physcons as pc
#--------------------------------------------------#
import modpilgrim.names   as PN
import modpilgrim.pilrw    as RW
import modpilgrim.strings as PS
#--------------------------------------------------#
from   modpilgrim.fit2anarc         import fit2anarc
from   modpilgrim.ChemReaction      import ChemReaction
from   modpilgrim.diverse           import ffchecking
from   modpilgrim.diverse           import status_check
from   modpilgrim.plotting          import manage_data_for_plot_rcons
from   modpilgrim.plotting          import write_plotfile
#--------------------------------------------------#

METHODS = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")

NINIT = 3

#===============================================================#
def fit2analytic(ltemp,rcons,nR):
    # fit to analytic
    if len(ltemp) > 1:
       human_units = pc.ML**(nR-1.0) / pc.SECOND
       k    = [val*human_units for val in rcons]
       dfit = fit2anarc(ltemp,k,log=True)
       # save data
       return (k,dfit)
    # return data
    return (None,None)
#---------------------------------------------------------------#
def fit_reaction(rcname,nR,nP,dall,ltemp):
    fitting       = {}
    fitting["fw"] = {}
    fitting["bw"] = {}
    # perform fittings
    for direc in "fw,bw".split(","):
        if direc == "fw": num = nR
        if direc == "bw": num = nP
        for X in METHODS:
           key = "%s.%s.%s"%(X,rcname,direc)
           ks  = dall["rcons"].get(key,None)
           if ks is None: continue
           ks_human, dfit = fit2analytic(ltemp,ks,num)
           fitting[direc][X] = (key,dfit,ks_human)
    return fitting
#---------------------------------------------------------------#
def prepare_plotdata(rcname,nR,nP,ltemp,fitting):
    plotdata = {}
    nones    = (None,None,None)
    # Deal with plots
    for direc in "fw,bw".split(","):
        if direc == "fw": num = nR
        if direc == "bw": num = nP
        # string for units (for plotting)
        if   num == 1: units = "s^{-1}"
        elif num == 2: units = "cm^3 / molecule / s"
        else         : units = "(cm^3)^%i / s /molecule^%i"%(num-1,num-1)
        d4plot = {}
        for X in METHODS:
           key1, dfit, ks_human = fitting[direc].get(X,nones)
           if key1 is None: continue
           # deal with plot data
           d4plot[X] = (ks_human,dfit)
        plotdata.update( manage_data_for_plot_rcons(rcname,direc,ltemp,d4plot,units) )
    return plotdata
#---------------------------------------------------------------#
def get_summary_lines(rcname,fitting,ltemp,ml=4):
    nones = (None,None,None)
    LINES = {}

    # forward and backward
    for direc in ["fw","bw"]:
        if len(fitting[direc].keys()) == 0: continue
        key0 = ("%%-%is"%ml)%("%s.%s"%(rcname,direc))
        # methods
        for X in METHODS:
            key1, dfit, ks_human = fitting[direc].get(X,nones)
            if dfit is None: continue
            # type of fitting
            for anatype in dfit.keys():
                coefs, r2 = dfit[anatype]
                if anatype in [4,5]: coefs = coefs[0:3]+[coefs[4],coefs[3]]
                coefs = " ".join(["%11.4E"%coef for coef in coefs])
                key2 = (anatype,X,rcname,direc)
                line  = "k(%s) analytic%i %-55s # r^2 = %.8f\n"%(key0,anatype,coefs,r2)
                LINES[key2] = line
    return LINES
#---------------------------------------------------------------#
def print_summary(targets,LINES):
    # Print lines for KMC
    string  = ""
    string += "=================\n"
    string += " FITTING SUMMARY \n"
    string += "=================\n"
    string += "\n"
    for X in METHODS:
        bloque_X = ""
        for anatype in [1,2,3,4,5]:
            bloque_at = ""
            for rcname in targets:
                for direc in "fw,bw".split(","):
                    key = (anatype,X,rcname,direc)
                    line = LINES.get(key,None)
                    if line is None: continue
                    bloque_at += line
                    #bloque_at += line.split("#")[0]+"\n"
            if bloque_at == "": continue
            bloque_X += bloque_at+"\n"
        if bloque_X == "": continue
        # add to string
        string += "%s\n"%PS.KEYNICE[X]
        string += "\n"+bloque_X+"\n"
    print_string(string,NINIT)
#===============================================================#



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
    targets = sorted([target for target in targets if target in dchem.keys()])

    # read dof
    dall = RW.read_alldata(dof,ltemp)[0]
    
    #--------------------------------#
    # Calculations for each reaction #
    #--------------------------------#
    print(" "*NINIT+"Analytic expressions:")
    print("")
    print(" "*NINIT+"(1) k = A exp(-B/T)              ")
    print(" "*NINIT+"(2) k = A*T^n*exp(-B/T)          ")
    print(" "*NINIT+"(3) k = A*(T/Tr)^n*exp(-B/T)     ")
    print(" "*NINIT+"(4) k = A*(T/Tr)^n*exp(-B*(T+T0)/(T^2+T0^2))")
    print(" "*NINIT+"(5) k = A*((T+T0)/Tr)^n*exp(-B*(T+T0)/(T^2+T0^2))")
    print("")

    LINES, plotdata = {}, {}
    ml = max(len(rcname) for rcname in targets)+3
    for rcname in targets:
        print_string(PS.sfit_reaction(rcname),3)
        Rs, TS, Ps = dchem[rcname]
        nR, nP     = len(Rs), len(Ps)
        # perform fitting
        fitting_i = fit_reaction(rcname,nR,nP,dall,ltemp)
        # generate plotdata
        plotdata_i = prepare_plotdata(rcname,nR,nP,ltemp,fitting_i)
        # generate summary lines
        LINES_i = get_summary_lines(rcname,fitting_i,ltemp,ml)
        # update dictionaries
        plotdata.update(plotdata_i)
        LINES.update(LINES_i)
        print_string(PS.sfit_fitting(rcname,fitting_i,ltemp),5)

    # update file with plots
    if plotfile is not None and plotdata != {}:
       write_plotfile(plotfile,plotdata)

    # print KMC lines
    if LINES != {}: print_summary(targets,LINES)

#===============================================================#

