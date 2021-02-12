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
| Sub-module :  diverse            |
| Last Update:  2020/03/01 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
functions used in Pilgrim
'''

#---------------------------------------------------------------#
import os
import numpy                  as     np
import modpilgrim.pilrw        as     RW
import modpilgrim.names       as     PN
import modpilgrim.strings     as     PS
#---------------------------------------------------------------#
from   modpilgrim.PathVars    import PathVars
from   modpilgrim.steepdesc   import sorted_points
#---------------------------------------------------------------#
import common.fncs            as     fncs
import common.Exceptions      as     Exc
import common.physcons        as     pc
#---------------------------------------------------------------#
from   common.files      import read_gtsfile
from   common.internal   import ic2string
from   common.Molecule   import Molecule
from   common.criteria   import EPS_MEPS
#---------------------------------------------------------------#


#===============================================================#
# Reading of all Pilgrim input files at once                    #
#===============================================================#
def get_input_data():
    '''
    * reads all Pilgrim input files
    * it also returns a status list where:
        * -1 : file does not exists
        *  0 : file exists but it is empty
        *  1 : file exists and contains data
    '''
    # read data
    idata  = []
    idata += [ RW.read_ctc()   ] # idx+1 = 1
    idata += [ RW.read_temp()  ] # idx+1 = 2
    idata += [ RW.read_path()  ] # idx+1 = 3
    idata += [ RW.read_tes()   ] # idx+1 = 4
    idata += [ RW.read_chem()  ] # idx+1 = 5
    idata += [ RW.read_kmc()   ] # idx+1 = 6
    idata += [ RW.read_dlevel()] # idx+1 = 7
    # Prepare dpath
    dpath, (filename,status) = idata[2]
    dpath = prepare_dpath(dpath)
    idata[2] = (dpath, (filename,status))
    # string with status
    string  = ""
    string += "Status of input files\n"
    string += "\n"
    string += "   ---------------------------------\n"
    string += "          input file       | status \n" 
    string += "   ---------------------------------\n"
    idx = 0
    for data,(fname,status) in idata: 
        idx += 1
        string += "    #%i : %-17s |   %2i   \n"%(idx,fname,status)
    string += "   ---------------------------------\n"
    string += "   status = -1 ==> file does not exist\n"
    string += "   status =  0 ==> file exists but it is empty\n"
    string += "   status =  1 ==> file exists and contains data\n"
    # split data
    the_data   = [data   for data,(fname,status) in idata]
    the_status = [status for data,(fname,status) in idata]
    # return data and string
    return the_data, the_status, string
#---------------------------------------------------------------#
def prepare_dpath(dpath):
    # Read data for each ctc
    for ctc, (pathtype,lines) in dpath.items():
        pathvars = PathVars(pathtype)
        # Read info in lines
        for line in lines:
            line  = line.lower()
            key   = line.split()[0]
            value = " ".join(line.split()[1:]).lower()
            try   :
               pathvars.setvar(key,value)
            except:
               exception = Exc.ReadProblem(Exception)
               exception._file = filename
               exception._var  = line.split("\n")[0]
               raise exception
        dpath[ctc] = pathvars
    return dpath
#---------------------------------------------------------------#
def dlevel_to_files(dlevel):
    # files according to case
    dof      = PN.get_dof(dlevel)
    plotfile = PN.get_plf(dlevel)
    hlf      = PN.get_hlf() # high-level file
    # string with information
    string  = "Files of importance:\n"
    string += "\n"
    if dlevel     : string += "  - dlevel --> yes\n"
    else          : string += "  - dlevel --> no\n"
    string += "\n"
    string += "  * file for data       storing: '%s'\n"%dof
    string += "  * file for high-level storing: '%s'\n"%hlf
    string += "  * file for plot       storing: '%s'\n"%plotfile
    string += "\n"
    return (dof,hlf,plotfile),string
#===============================================================#


#===============================================================#
# Functions where some things are checked                       #
#===============================================================#
def ffchecking(mustexist=[],tocreate=[]):
    '''
    * checks the existence of folders in 'mustexist'
    * creates folders in 'tocreate'
    '''
    # folders that must exist
    for folder in mustexist:
        if os.path.exists(folder): continue
        return -1
    # folders to create
    for folder in tocreate:
        if not os.path.exists(folder):
           os.mkdir(folder)
    return 0
#---------------------------------------------------------------#
def status_check(lstatus,indices):
    '''
    * checks the status in 'lstatus' associated
      with the indices in filesIDX
    * idx in indices, from 1 to 8
    '''
    final_status = 0
    for idx in indices:
        if lstatus[idx-1] != 1:
           print("     - status of input file #%i is not 1\n"%(idx))
           final_status = -1
    return final_status
#===============================================================#


#===============================================================#
# Functions to do some useful things                            #
#===============================================================#
def calc_fwdir(gtsfile):
    '''
    * calculates which internal coordinate changes the most
      in the transition state
    * it also gives the sign of the variation in the forward
      direction of the MEP
    '''
    if (gtsfile is None) or (not os.path.exists(gtsfile)): return (None,None)
    # Generate Molecule from gts file
    molecule = Molecule()
    molecule.set_from_gts(gtsfile)
    # setup (with frequency calculation)
    molecule.setup()
    ic, fwsign = molecule.get_imag_main_dir()
    # return data
    return (ic2string(ic),fwsign)
#---------------------------------------------------------------#
def find_label_in_rst(s_i,drst):
    # initialize output label
    label  = None
    svalue = None
    # get all points
    allpoints = sorted_points(drst,hess=False)
    # as_i is indeed a label!
    if s_i in allpoints: return s_i, drst[s_i][0]
    # s_i is not a float number
    try   : s_i = float(s_i)
    except: return label, s_i
    # go one by one
    for l_j in allpoints:
        s_j  = drst[l_j][0]
        if abs(s_j-s_i) < EPS_MEPS:
           label  = l_j
           svalue = s_j
           break
    return label, svalue
#===============================================================#


#===============================================================#
def get_contributions(ctc,dall,dctc,ltemp):
    if ctc not in dctc.keys()   : return None
    # Only one conformer?
    if len(dctc[ctc]._itcs) == 1:
       itc,weight = dctc[ctc]._itcs[0]
       dchi = {itc : [1.0 for T in ltemp]}
       return dchi
    # Get total partition function
    V0, V1, PFN  = dall["pfn"]["%s.msho"%ctc]
    # Get ratios
    dchi = {}
    for itc,weight in dctc[ctc]._itcs:
        # Get individual partition functions
        try:
           V0i,V1i,PFNi = dall["pfn"][PN.struckey(ctc,itc)]
           # Calculate contribution
           dE        = (V1i-V1)
           exp_arg   = [-dE/pc.KB/T for T in ltemp]
           ratio_pfn = [weight*pfi/pftot for pfi,pftot in zip(PFNi,PFN)]
           chi_i     = [aa*fncs.exp128(bb) for aa,bb in zip(ratio_pfn,exp_arg)]
           dchi[itc] = np.array(chi_i)
        except:
           exception = Exc.LostConformer(Exception)
           exception._var = PN.struckey(ctc,itc)
           raise exception
    # Return data
    return dchi
#---------------------------------------------------------------#
def get_transmissioncoeffs(ctc,dall,dctc,ltemp):
    # Get contributions
    dchi = get_contributions(ctc,dall,dctc,ltemp)
    if ctc not in dctc.keys(): return None, dchi
    # Get individual transmission coefficients
    dtc  = {X:{} for X in "tst,tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")}
    # Get transmission coefficients
    for itc,weight in dctc[ctc]._itcs:
        tsname = PN.struckey(ctc,itc)
        cvt    = dall.get("cvt"   ,{}).get(tsname,None)
        zct    = dall.get("zct"   ,{}).get(tsname,None)
        sct    = dall.get("sct"   ,{}).get(tsname,None)
        cagtst = dall.get("cagtst",{}).get(tsname,None)
        cagcvt = dall.get("cagcvt",{}).get(tsname,None)
        # now total transmission coefficients
        try   : tc_tstzct = np.array(fncs.prod_list((zct,    cagtst)))
        except: tc_tstzct = None
        try   : tc_tstsct = np.array(fncs.prod_list((sct,    cagtst)))
        except: tc_tstsct = None
        try   : tc_cvt    = np.array(fncs.prod_list((    cvt,      )))
        except: tc_cvt    = None
        try   : tc_cvtzct = np.array(fncs.prod_list((zct,cvt,cagcvt)))
        except: tc_cvtzct = None
        try   : tc_cvtsct = np.array(fncs.prod_list((sct,cvt,cagcvt)))
        except: tc_cvtsct = None
        # save
        dtc["tst"   ][itc] = np.array([1.0 for T in ltemp])
        dtc["tstzct"][itc] = tc_tstzct
        dtc["tstsct"][itc] = tc_tstsct
        dtc["cvt"   ][itc] = tc_cvt
        dtc["cvtzct"][itc] = tc_cvtzct
        dtc["cvtsct"][itc] = tc_cvtsct
    # Get averaged transmission coefficient
    for X in dtc.keys():
        # Initializing averaged value
        dtc[X]["averaged"] = np.array([0.0 for T in ltemp])
        # Calculating averaged coef
        for itc,weight in dctc[ctc]._itcs:
            # No data for this method?
            if dtc[X][itc] is None:
               dtc[X]["averaged"] = None
               break
            # update averaged transmission coeff
            else: dtc[X]["averaged"] += dchi[itc] * dtc[X][itc]
    return dtc, dchi
#===============================================================#

