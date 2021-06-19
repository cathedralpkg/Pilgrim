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
| Sub-module :  optKMC             |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --kmc option
of Pilgrim
'''

#--------------------------------------------------#
import datetime
import time
import os
import sys
#--------------------------------------------------#
import numpy             as     np
import modpilgrim.pilrw         as     RW
import modpilgrim.names             as     PN
import modpilgrim.strings           as     PS
#--------------------------------------------------#
from   common.fncs       import print_string
from   common.fncs       import time2human
from   common.fncs       import exp128
#--------------------------------------------------#
#--------------------------------------------------#
import common.Exceptions as Exc
#--------------------------------------------------#
from   common.Logger     import Logger
#--------------------------------------------------#
from   common.physcons   import ML
from   common.physcons   import SECOND
from   common.physcons   import NA
#--------------------------------------------------#
from   modpilgrim.fit2anarc  import log_anarc1
from   modpilgrim.fit2anarc  import log_anarc2
from   modpilgrim.fit2anarc  import log_anarc3
from   modpilgrim.fit2anarc  import log_anarc4
from   modpilgrim.fit2anarc  import log_anarc5
from   modpilgrim.plotting   import manage_data_for_plot_kmc
from   modpilgrim.plotting   import write_plotfile
from   modpilgrim.diverse    import get_input_data
from   modpilgrim.diverse    import ffchecking
from   modpilgrim.diverse    import status_check
from   modpilgrim.kmc        import kmc
#--------------------------------------------------#



#---------------------------------------------------------------#
def get_ratecons(rcs,dchem,dall,idx,temp):
    drcons    = dall.get("rcons",{})
    processes = []
    for key,(rctype,weight,coefs) in sorted(rcs.items()):
        rctype = rctype.lower()
        # reaction name and direction
        if "." in key: rcname, direction = key.split(".")
        else         : rcname, direction = key, "both"
        # elements in reaction
        Rs,TS,Ps = dchem[rcname]
        nR,nP = len(Rs),len(Ps)
        # read/calculate rate constant
        if "analytic" not in rctype:
           key_fw = "%s.%s.%s"%(rctype,rcname,"fw")
           key_bw = "%s.%s.%s"%(rctype,rcname,"bw")
           # get rate constants
           kfw = drcons.get(key_fw,None)
           kbw = drcons.get(key_bw,None)
           if kfw is not None: kfw = kfw[idx]
           if kbw is not None: kbw = kbw[idx]
           # any non-desired rate constants
           if direction == "fw": kbw = None
           if direction == "bw": kfw = None
           # none
           if kbw is None and kfw is None:
              exception =Exc.NoRateCons(Exception)
              exception._var = (rctype,key)
              raise exception
        else:
           if   rctype.lower() == "analytic1": k = log_anarc1(temp,*coefs)
           elif rctype.lower() == "analytic2": k = log_anarc2(temp,*coefs)
           elif rctype.lower() == "analytic3": k = log_anarc3(temp,*coefs)
           elif rctype.lower() == "analytic4": k = log_anarc4(temp,*coefs)
           elif rctype.lower() == "analytic5": k = log_anarc5(temp,*coefs)
           else                              : k = None
           # log --> exp
           if k is not None: k = exp128(k)
           # save data properly
           if   direction in ["fw","both"]: kfw, kbw = k   , None
           elif direction in ["bw"       ]: kfw, kbw = None, k
           else                           : kfw, kbw = None, None
           # in atomic units
           hunitsFW = ML**(nR-1.0) / SECOND
           hunitsBW = ML**(nP-1.0) / SECOND
           if kfw is not None: kfw /= hunitsFW
           if kbw is not None: kbw /= hunitsBW
        # ignore reactions giving rise to bimolecular products
        if len(Ps) > 1 and direction != "bw": kbw = None
        # save in processes
        if kfw is not None: processes.append( (Rs,Ps,weight*kfw) )
        if kbw is not None: processes.append( (Ps,Rs,weight*kbw) )
    return processes
#---------------------------------------------------------------#
#def kmc_savedata(kmcof,data):
#    # molecules and same length
#    molecules = sorted(data[0][2].keys())
#    len1 = max([len(string) for string in molecules+["time","-"*10]])
#    ss   = "%%%is"%len1
#    # string for gnuplot data file
#    string = ""
#    string += "# Time in ps; populations in molecules\n"
#    string += "\n\n"
#    # string for gnuplot string
#    gpstr  = "datafile='%s'\n"%kmcof
#    gpstr += "plot\\\n"
#    # data for each temperature
#    for ii,(stemp,xvalues,yvalues) in enumerate(data):
#         # add line in gnuplot script
#         gpstr += "    datafile  index %i  using 1:2 with points,\\\n"%ii
#         # add data in gnuplot file
#         string += "# T = %s K\n"%stemp
#         string += "#"+ss%"time" + " ".join([ss%mol for mol in molecules]) + "\n"
#         for jj,time in enumerate(xvalues):
#             time = 1e12*SECOND*time
#             string += "  %10.3E "%time
#             for molecule in molecules:
#                 pop = yvalues[molecule][jj]/NA
#                 string += " %10.3E "%pop
#             string += "\n"
#         string += "\n\n"
#    # write files
#    with open(kmcof,'w') as asdf: asdf.write(string)
#    with open("plotkmc.gp",'w') as asdf: asdf.write(gpstr)
#---------------------------------------------------------------#



def main(idata,status,case,targets):

    stat2check = [2,5,6]
    mustexist  = []
    tocreate   = [PN.DIR3,PN.DIR6]
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
    if len(targets) == 0 and dkmc is not None: targets = dkmc.keys()

    for target in sorted(targets):
        print("   --> %s"%target)
        print("")
        if target not in dkmc.keys():
           print("       KMC not defined!")
           print("")
           continue
        t1 = time.time()
        pof = PN.get_pof(dlevel,"kmc.%s"%target)
        print("   Pilgrim output file: %s"%pof)
        print("")

        sys.stdout = Logger(pof,"w",True)
        sys.stdout.writeinfile(PS.init_txt())

        # expand KMC tuple
        ipops,rcs,psteps,volume,timeunits,excess = dkmc[target]
        valid_tu = "fs,ps,mcs,ms,s,min,hr".split(",")
        if timeunits not in valid_tu: timeunits = "ps"

        # reactivo limitante
        try   : POP0 = min([ipop for species,ipop in sorted(ipops.items()) if ipop != 0.0])
        except: POP0 = 0

        print_string(PS.skmc_init(ipops,POP0,excess,rcs,psteps,volume,timeunits,dof),5)

        # continue?
        if POP0 == 0:
            print("   All initial populations are zero...")
            print("")
            return
        if rcs == {}:
            print("   No reactions considered in %s"%PN.IFILE6)
            print("")
            return

        # read dof
        dall = RW.read_alldata(dof,ltemp)[0]

        # perform KMC for each temperature
        data = []
        fratios = {}
        ftimes  = []
        for idx,temp in enumerate(ltemp):
            # title
            title    = " T = %.3f Kelvin "%temp
            division ="-"*len(title)
            string  = "     "+division+"\n"
            string += "     "+title   +"\n"
            string += "     "+division+"\n"
            print(string)
            # get rate constants
            processes = get_ratecons(rcs,dchem,dall,idx,temp)
            # print initial information before simulation
            print_string(PS.skmc_processes(temp,processes),9)
            # perform kmc simulation
            xvalues, yvalues = kmc(ipops,processes,excess,volume,psteps)
            fratios["%7.2f"%temp] = {}
            for species in yvalues.keys():
                fratios["%7.2f"%temp][species] = yvalues[species][-1] / POP0
            # print data from simulation
            print_string(PS.skmc_results(xvalues,yvalues,timeunits),9)
            # save data needed for txt and pdf files
            data += [ ("%.2f"%temp,xvalues,yvalues) ]
            ftimes.append(xvalues[-1])

        # print final ratios
        species = sorted(yvalues.keys())
        print_string(PS.skmc_finaltimes(ltemp,ftimes,timeunits),5)
        print_string(PS.skmc_finalratio(ltemp,species,fratios),5)

        # save data for plotting
        if plotfile is not None:
           plotdata = {}
           plotdata.update(manage_data_for_plot_kmc(target,data,fratios,volume,timeunits))
           write_plotfile(plotfile,plotdata)
        #kmc_savedata(kmcof,data)

        # print end of file
        t2    = time.time()
        etime = time2human(t2-t1,"secs")
        sys.stdout.writeinfile(PS.end_txt(*etime))
        # End logger
        sys.stdout = Logger(None)

#---------------------------------------------------------------#


