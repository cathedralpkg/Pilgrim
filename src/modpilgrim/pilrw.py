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
| Sub-module :  pilrw              |
| Last Update:  2020/03/01 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#--------------------------------------------------#
import os
import sys
import numpy as np
import modpilgrim.names              as     PN
#--------------------------------------------------#
from   modpilgrim.ClusterConf  import ClusterConf
import modpilgrim.helps   as PH
#--------------------------------------------------#
import common.Exceptions  as     Exc
from   common.files       import read_file
from   common.fncs        import clean_line
from   common.fncs        import clean_lines
from   common.fncs        import afreq2cm
from   common.fncs        import xyz
from   common.fncs        import same_lfloats
from   common.files       import read_gtsfile
from   common.dicts       import dpt_im
from   common.internal    import unmerge_ics, string2ic
from   common.physcons    import AMU, KCALMOL, ML
import common.fncs        as     fncs
#--------------------------------------------------#
from   common.criteria   import EPS_TEMP
#--------------------------------------------------#


#====================================================#
def file2lines(filename):
    if not os.path.exists(filename):
        lines  = []
        status = -1
    else:
        lines = read_file(filename)
        lines = clean_lines(lines,"#",True)
        lines = [line for line in lines if line != "\n"]
        if lines == []: status = 0
        else          : status = 1
    return lines, status
#====================================================#


#====================================================#
#         Functions for: READING INPUT FILES         #
#====================================================#
def read_track(filename=PN.IFILE0):
    lines, status = file2lines(filename)
    # Get info from lines
    dict_track = {}
    for line in lines:
        ifile, arrow, ffile = line.split()
        dict_track[ifile] = ffile
    return dict_track, (filename,status)
#----------------------------------------------------#
def read_ctc(filename=PN.IFILE1):
    '''
    - Returns dictionaries:
          dctc[ctc]       = ClusterConf instance
          dimasses[label] = mass
    '''
    lines, status = file2lines(filename)
    # Get info for lines
    record1  = False
    record2  = False
    dctc     = {}
    dimasses = {}
    for line in lines:
        #---------------#
        # isomass block #
        #---------------#
        if "end_isomass"   in line: record1 = False
        if record1:
           if "=" in line: label,mass = line.split("=")
           else          : label,mass = line.split()
           label = label.strip()
           mass  = mass.strip()
           dimasses[label] = float(mass)/AMU
        if "start_isomass" in line: record1 = True
        #-----------#
        # ctc block #
        #-----------#
        if "end_struc"   in line or "end_ctc" in line:
           record2  = False
           # set from lines
           CTC = ClusterConf(name)
           CTC.set_from_piflines("".join(ctclines))
           # save ClusterConf instance
           dctc[name] = CTC
        if record2: ctclines.append(line)
        if "start_struc" in line or "start_ctc" in line:
           record2  = True
           name     = line.split()[1]
           ctclines = []
    return (dctc,dimasses),(filename,status)
#----------------------------------------------------#
def read_temp(filename=PN.IFILE2):
    lines, status = file2lines(filename)
    # Get info for lines
    Tlist = []
    for line in lines:
        # special case - not documented
        if "range" in line:
           dummy, T0, Tf, dT = line.split()
           T0 = float(T0)
           Tf = float(Tf)
           dT = float(dT)
           nsteps = int(round((Tf-T0)/dT))+1
           Tlist += [T0+k*dT for k in range(nsteps)]
           continue
        Tlist += [float(T) for T in line.split()]
    return sorted(Tlist), (filename,status)
#----------------------------------------------------#
def read_path(filename=PN.IFILE3):
    lines, status = file2lines(filename)
    # Get info for lines
    dpath = {}
    ctc_lines    = None
    bool_mep     = False
    bool_drp     = False
    saveline     = False
    for line in lines:
        #-- MEP --#
        if line.startswith("end_mep"):
             bool_mep  = False
             saveline  = False
             dpath[ctc] = ("mep",pathlines)
             continue
        #-- DRP --#
        if line.startswith("end_drp"):
             bool_drp  = False
             saveline  = False
             dpath[ctc] = ("drp",pathlines)
             continue
        # save line
        if saveline: pathlines.append(line)
        #-- MEP --#
        if line.startswith("start_mep"):
             bool_mep  = True
             saveline  = True
             ctc       = line.split()[1]
             pathlines = []
             continue
        #-- DRP --#
        if line.startswith("start_drp"):
             bool_drp  = True
             saveline  = True
             ctc       = line.split()[1]
             pathlines = []
             continue
    return dpath, (filename,status)
#----------------------------------------------------#
def read_tes(filename=PN.IFILE4):
    lines  = read_file(filename)
    dtesLL = {}
    dtesHL = {}
    bkeep = False
    sk1 = "start_meppoint"
    sk2 = "start_highlevel"
    ek1 = "end_meppoint"
    ek2 = "end_highlevel"
    for line in lines:
        # save LL data
        if   line.startswith(ek1):
             dtesLL[software] = dtesLL.get(software,{})
             dtesLL[software].update({ctc:string})
             bkeep = False
        # save HL data
        elif line.startswith(ek2):
             dtesHL[software] = dtesHL.get(software,{})
             dtesHL[software].update({ctc:string})
             bkeep = False
        elif bkeep:
             string += line
        elif line.startswith(sk1) or line.startswith(sk2):
             ctc, software = line.split()[1:]
             string = ""
             bkeep  = True
        else: continue
    if not os.path.exists(filename)   : status = -1
    elif dtesLL == {} and dtesHL == {}: status =  0
    else                              : status =  1
    return (dtesLL,dtesHL), (filename,status)
#----------------------------------------------------#
def read_chem(filename=PN.IFILE5):
    lines, status = file2lines(filename)
    # Get info for lines
    dchem = {}
    for line in lines:
        # read line and split
        if line.count(":") == 1:
           rcname, reaction = line.split(":")
        else: continue
        # reaction formula
        Rs, TS, Ps = reaction.split("-->")
        # some strips
        rcname = rcname.strip()
        Rs     = [r.strip() for r in Rs.split("+")]
        TS = TS.strip()
        Ps     = [p.strip() for p in Ps.split("+")]
        if Rs == [""]: Rs = []
        if TS ==  "" : TS = None
        if Ps == [""]: Ps = []
        dchem[rcname] = (Rs,TS,Ps)
    return dchem, (filename,status)
#----------------------------------------------------#
def data_in_kmc_block(kmc_block):
    # first line --> name for this kmc
    fline = kmc_block.split("\n")[0]
    if not fline.startswith("kmc"): return None, None
    # initialize variables
    kmc_name        = fline.split()[1]
    ipopulations    = {}
    numbersteps     = 10000
    volume          = 1.0 / ML
    timeunits       = "ps"
    excess          = []
    rateconstants   = {}
    # Get info for lines
    for line in kmc_block.split("\n"):
        try:
           ldata = line.split()
           variable = ldata[0]
           if   variable == "psteps"   : numbersteps =   int(ldata[1])
           elif variable == "volume"   : volume      = float(ldata[1]) / ML
           elif variable == "timeunits": timeunits   = ldata[1]
           elif "pop0(" in variable:
               species = variable.split("pop0(")[1].split(")")[0].strip()
               ipop   = ldata[1]
               if len(ldata) == 2 and ldata[1].lower() == "excess":
                   excess.append(species)
               ipopulations[species] = float(ipop)
           elif "k(" in variable:
               reaction = line.split("(")[1].split(")")[0].strip()
               if "*" in line:
                  weight = int(line.split("*")[1].split()[0])
                  rcline = " ".join(line.split("*")[1].strip().split()[1:])
               else:
                  weight = 1
                  rcline = line.split(")")[1].strip()
               rctype = rcline.split()[0].lower()
               coefs  = None
               if   rctype == "analytic1":
                   A,B = rcline.split()[1:3]
                   coefs = [float(A),float(B)]
               elif rctype == "analytic2":
                   A,B,n = rcline.split()[1:4]
                   coefs = [float(A),float(B),float(n)]
               elif rctype == "analytic3":
                   A,B,n,Tr = rcline.split()[1:5]
                   coefs = [float(A),float(B),float(n),float(Tr)]
               elif rctype == "analytic4":
                   A,B,n,Tr,T0 = rcline.split()[1:6]
                   coefs = [float(A),float(B),float(n),float(T0),float(Tr)]
               elif rctype == "analytic5":
                   A,B,n,Tr,T0 = rcline.split()[1:6]
                   coefs = [float(A),float(B),float(n),float(T0),float(Tr)]
               # in case user is not using ms or mp prefix
               elif rctype != "tst" and not rctype.startswith("m") and not rctype.startswith("ana"):
                    rctype = "ms"+rctype
               rateconstants[reaction] = (rctype, weight, coefs)
        except: pass
    kmc_tuple = (ipopulations, rateconstants, numbersteps, volume, timeunits, excess)
    return kmc_name,kmc_tuple
#----------------------------------------------------#
def read_kmc(filename=PN.IFILE6):
    # read lines
    lines, status = file2lines(filename)
    # split into blocks
    blocks = [block.split("end_kmc")[0] for block in "".join(lines).split("start_")]
    # Analyze each block
    dkmc = {}
    for block in blocks:
        kmc_name, kmc_tuple = data_in_kmc_block(block)
        if kmc_name is None: continue
        dkmc[kmc_name] = kmc_tuple
    # Merge data
    if status in [0,-1]: dkmc = None
    return dkmc, (filename,status)
#----------------------------------------------------#
def read_dlevel(filename=PN.IFILE7):
    lines, status = file2lines(filename)
    ddlevel = {}
    for line in lines:
        line = line.strip()
        if line == "": continue
        key     = line.split()[0]
        other   = " ".join(line.split()[1:])
        spoints = []
        if "{" in other and "}" in other:
          spoints += [si for si in other.split("{")[1].split("}")[0].split(",")]
        ddlevel[key] = spoints
    return ddlevel, (filename,status)
#====================================================#


#====================================================#
#         Functions for: WRITING INPUT FILES         #
#====================================================#
def write_track(dict_track,filename=PN.IFILE0):
    if dict_track == {}:
       string = "# NO DATA INSIDE %s\n"%PN.UFOLDER
    else:
       # lengths
       mlkey = max([len(key) for key in dict_track.keys()  ])
       mlval = max([len(val) for val in dict_track.values()])
       # format of each line
       lformat = "%%%is   -->  %%%is\n"%(mlkey,mlval)
       string  = "".join([lformat%(ifile,ffile) for ifile,ffile in sorted(dict_track.items())])
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_ctc(dctc,dimasses={},filename=PN.IFILE1):
    string = ""
    #------------------------#
    # write info in dimasses #
    #------------------------#
    if dimasses != {}:
       string += "#------------------------------------------#\n"
       string += "#    Definition of the isotopic masses     #\n"
       string += "#------------------------------------------#\n"
       string += "start_isomass\n"
       toloop = sorted([(imass,label) for label,imass in dimasses.items()])
       for imass,label in toloop:
           imass *= AMU
           string += "   %-5s = %12.7f\n"%(label,imass)
       string += "end_isomass\n"
       string += "\n"
    #------------------------#
    # write info in dctc     #
    #------------------------#
    string0 = ""
    string1 = ""
    stringn = ""
    for name,ctc in sorted(dctc.items()):
        string_ctc  = "start_ctc %s\n"%name
        for line in ctc.get_piflines().split("\n"):
            if line == "": continue
            string_ctc+= "   "+line+"\n"
        string_ctc += "end_ctc\n\n"
        if   int(ctc._type) == 0: string0 += string_ctc
        elif int(ctc._type) == 1: string1 += string_ctc
        else                    : stringn += string_ctc
    # Generate string with all data
    if string0 != "":
       string += "#------------------------------------------#\n"
       string += "# Cluster of torsional conformers:  MINIMA #\n"
       string += "#------------------------------------------#\n"
       string += string0
    if string1 != "":
       string += "#------------------------------------------#\n"
       string += "# Cluster of torsional conformers: SADDLES #\n"
       string += "#------------------------------------------#\n"
       string += string1
    if stringn != "":
       string += "#------------------------------------------#\n"
       string += "# Cluster of torsional conformers: OTHERS  #\n"
       string += "#------------------------------------------#\n"
       string += stringn
    string += PH.INFO_pifstruc
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_temp(Tlist,filename=PN.IFILE2):
    string = ""
    for idx in range(0,len(Tlist),7):
        line = "  ".join(["%7.2f"%T for T in Tlist[idx:idx+7]])
        string += "   " + line + "\n"
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_path(dpath,filename=PN.IFILE3):
    string = "\n"
    for target,pathvars in sorted(dpath.items()):
        string += pathvars.string4pif(target)
        string += "\n"
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_tes(dtesLL,dtesHL,filename=PN.IFILE4):
    dstring = {}
    for case,dtes in [("LL",dtesLL),("HL",dtesHL)]:
        # loop for each software
        for software in sorted(dtes.keys()):
            # initialize dictionary
            dstring[software]       = dstring.get(software,{})
            dstring[software][case] = dstring[software].get(case,"")
            # block case?
            if case == "LL": blockkey = "meppoint"
            if case == "HL": blockkey = "highlevel"
            # generate string
            string  = ""
            for ctc in sorted(dtes[software].keys()):
                string += "#- - - - - - - - - - - - - -#\n"
                string += "start_%s %s %s\n"%(blockkey,ctc,software)
                string += dtes[software][ctc]
                string += "end_%s\n"%blockkey
                string += "#- - - - - - - - - - - - - -#\n"
                string += "\n"
            dstring[software][case] += string

    # Final string
    final_string = ""
    for software in dstring.keys():
        final_string += "#===========================%-10s=#\n"%("="*10)
        final_string += "# INPUT FILE TEMPLATES FOR: %-10s #\n"%(software.upper())
        final_string += "#===========================%-10s=#\n"%("="*10)
        final_string += "\n"
        for case in dstring[software].keys():
            if case == "LL":
               final_string += "#-------------------------#\n"
               final_string += "# MEP POINT CALCULATIONS  #\n"
               final_string += "#-------------------------#\n"
            if case == "HL":
               final_string += "#-------------------------#\n"
               final_string += "# HIGH-LEVEL CALCULATIONS #\n"
               final_string += "#-------------------------#\n"
            final_string += "\n"
            final_string += dstring[software][case]
        final_string += "\n"

    with open(filename,'w') as asdf: asdf.write(final_string)
#----------------------------------------------------#
def write_chem(dchem,filename=PN.IFILE5):
    string  = PH.INFO_pifchem+"\n"
    # nice formatting
    ml1,ml2,ml3,ml4 = 3,3,3,3
    for rcname,(Rs,TS,Ps) in dchem.items():
        if TS is None: TS = ""
        l1 = len(rcname)
        l2 = len(" + ".join(Rs))
        l3 = len(TS)
        l4 = len(" + ".join(Ps))
        if l1 > ml1: ml1 = l1
        if l2 > ml2: ml2 = l2
        if l3 > ml3: ml3 = l3
        if l4 > ml4: ml4 = l4
    # remove blanks in string
    string = "\n".join([line.strip() for line in string.split("\n")])
    # add to string
    for label,(Rs,TS,Ps) in sorted(dchem.items()):
        Rs = " + ".join(Rs)
        Ps = " + ".join(Ps)
        if TS is None: TS = " "
        # nice formating
        label = ("%%-%is"%ml1)%label
        Rs    = ("%%-%is"%ml2)%Rs
        TS    = ("%%-%is"%ml3)%TS
        Ps    = ("%%-%is"%ml4)%Ps
        # add to string
        string += "%s : %s --> %s --> %s\n"%(label,Rs,TS,Ps)
    string += "\n"
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_kmc(dkmc,filename=PN.IFILE6):
    string = ""
    for kmc_name,kmc_tuple in dkmc.items():
        ipopulations, rateconstants, numbersteps, volume, timeunits, excess = kmc_tuple
        if timeunits not in "fs,ps,mcs,ms,s,min,hr".split(","): timeunits = "ps"
        # get l1 and l2 for nice formatting
        try   : l1 = max(len(ctc)  for ctc  in ipopulations.keys())+7
        except: l1 = 7
        try   : l2 = max(len(reac) for reac in rateconstants.keys())+3
        except: l2 = 3
        # Write step for printing
        if kmc_name is not None: string += "start_kmc %s\n"%kmc_name
        string += "   # KMC Parameters\n"
        string += "   psteps     %-5i      # print data each nstp steps\n"%numbersteps
        string += "   volume     %8.2E   # simulation volume (mL)\n"%(volume*ML)
        string += "   timeunits  %-5s      # units for time variable\n"%timeunits
        string += "\n"
        # Concentration of each species
        string += "   # Initial (non-zero) populations (number of molecules)\n"
        for specie,ipop in sorted(ipopulations.items()):
            specie = "pop0(%s)"%specie
            string += ("   %%-%is   %%.2e"%l1)%(specie,ipop)
            if specie in excess: string += "  excess"
            string += "\n"
        string += "\n"
        # Type of rate constant to use
        string += "   # Selection of the rate constant to use\n"
        for reaction,(rctype,weight,rcdata) in sorted(rateconstants.items()):
            reaction = "k(%s)"%reaction
            if weight == 1: string += ("   %%-%is  %%s  "%l2)%(reaction,rctype.lower())
            else          : string += ("   %%-%is*%%i  %%s  "%l2)%(reaction,weight,rctype.lower())
            if rcdata is not None:
               string += " ".join(["%.3e"%coef for coef in rcdata])
            string += "\n"
        if kmc_name is not None: string += "end_kmc\n"
        string += "\n"
    # Write string
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def write_dlevel(ddlevel,filename=PN.IFILE7):
    string = ""
    ml = max([len(ctc) for ctc in ddlevel.keys()]+[1])
    for key,points in sorted(ddlevel.items()):
        string += "%%-%is "%ml%key
        if len(points) != 0:
           string += "{%s}"%",".join(points)
        string += "\n"
    with open(filename,'w') as asdf: asdf.write(string)
#====================================================#


#====================================================#
#        Functions for:  READING/WRITING DATA        #
#====================================================#
def read_highlevelfile(filename):
    '''
    dhighlvl: key ==> value

    Formats :
         ctc.itc ==> Energy (for sp)
         ctc.itc ==> {s:E}  (for path points)
    '''
    dhighlvl = {}
    # read file
    lines, status = file2lines(filename)
    for line in lines:
        key, slabel, Ehl = line.split()
        # save to dictionary
        if slabel == "sp":
           dhighlvl[key] = float(Ehl)
        else:
           dhighlvl[key] = dhighlvl.get(key,{})
           dhighlvl[key][slabel] = float(Ehl)
    return dhighlvl
#----------------------------------------------------#
def write_highlevelfile(filename,dhighlvl):
    # if exists, read it before
    dhighlvl_0 = read_highlevelfile(filename)
    # update with new
    for key,value in dhighlvl.items():
        # new key
        if key not in dhighlvl_0.keys():
           dhighlvl_0[key] = value
           continue
        # old key and "sp"
        if type(value) == float:
           dhighlvl_0[key] = value
           continue
        # old key and path points
        for slabel,E in value.items():
           dhighlvl_0[key][slabel] = E
    # save with correct name and del
    dhighlvl = dhighlvl_0
    del dhighlvl_0
    # for nice formatting
    ml = max([len(key) for key in dhighlvl.keys()])
    # write string (two parts: one for sp, one for MEP points)
    string_sp   = ""
    string_path = ""
    for key,value in sorted(dhighlvl.items()):
        key = "%%-%is"%ml%key
        # stationary points
        if type(value) == float:
           string_sp   += "%s  %-6s  %-+14.8f\n"%(key,"sp",value)
        # path points
        if type(value) == dict:
           for slabel,E in sorted(value.items()):
               string_path += "%s  %-6s  %-+14.8f\n"%(key,slabel,E)
           string_path += "\n"
    string = string_sp + "\n" + string_path
    with open(filename,'w') as asdf: asdf.write(string)
#----------------------------------------------------#
def read_alldata(filename,ltemp=None):
    #initialize dictionary
    variables = ["gibbs1cm3","gibbs1bar","pfn","ctr","anh","zct","sct","cvt","cagtst","cagcvt","cfs","rcons"]
    dall  = { variable:{} for variable in variables}
    # read file
    lines, status = file2lines(filename)
    # divide into sections
    lines  = "".join(lines)
    blocks = lines.split("start_")
    # temperatures in pif.temp
    if ltemp is None:
       ltemp,(dummy,status) = read_temp(filename=PN.IFILE2)
       if status != 1: ltemp = None
    # get data
    for block in blocks:
        # get clean block
        block = block.split("end_")[0]
        # divide into lines
        lines = block.split("\n")
        lines = [line for line in lines if line != ""]
        if lines == []: continue
        # variable?
        variable = lines[0].split()[0]
        key      = lines[0].split()[1]
        lines    = lines[1:]
        if variable not in variables: continue
        # V0 and V1
        if variable == "pfn":
           V0, V1 = float(lines[0]),float(lines[1])
           lines  = lines[2:]
        # lines to float
        data   = [line.split() for line in lines]
        if variable == "ctr":
           ltemp2 = [float(T)      for T,val1,val2,val3 in data]
           PIB    = [float(val1)   for T,val1,val2,val3 in data]
           RRHO   = [float(val2)   for T,val1,val2,val3 in data]
           ELE    = [float(val3)   for T,val1,val2,val3 in data]
        else:
           ltemp2 = [float(T)         for T,val in data]
           data2  = []
           for T,val in data:
               if "E" not in val                          : data2.append(float(val))
               elif -300 < float(val.split("E")[-1]) < 300: data2.append(float(val))
               else                                       : data2.append(np.float128(val))
           data = data2
        # save data
        if   variable == "pfn": dall["pfn"][key] = (V0,V1,data)
        elif variable == "ctr": dall["ctr"][key] = (PIB,RRHO,ELE)
        else                  : dall[variable][key] = data
        # check temperatures
        if ltemp is None:
           ltemp = ltemp2
        else:
           ok = fncs.same_lfloats(ltemp,ltemp2,EPS_TEMP)
           # raise exception if wrong
           exception = Exc.IncompData
           exception._var = filename
           if not ok: raise exception
    return dall, ltemp
#----------------------------------------------------#
def write_alldata(filename,ltemp,dall):
    variables = ["gibbs1cm3","gibbs1bar","pfn","ctr","anh","zct","sct","cvt","cagtst","cagcvt","cfs","rcons"]
    # write variables in order
    string = ""
    for variable in variables:
        if variable not in dall.keys(): continue
        string += "#-------------------------#\n"
        string += "# Adding data for: %6s #\n"%variable
        string += "#-------------------------#\n"
        # add info about variable
        if variable == "gibbs1cm3": string += "# Gibbs free energy in hartree (ref: 1 molec/cc)\n"
        if variable == "gibbs1bar": string += "# Gibbs free energy in hartree (ref: 1 bar)\n"
        if variable == "pfn"    : string += "# V0, V0+zpe and partition functions\n"
        if variable == "ctr"    : string += "# transl, rovib and elec partition functions\n"
        if variable == "anh"    : string += "# anharmonicity\n"
        if variable == "zct"    : string += "# ZCT transmission coefficient\n"
        if variable == "sct"    : string += "# SCT transmission coefficient\n"
        if variable == "cagtst" : string += "# CAG/TST coefficient\n"
        if variable == "cagcvt" : string += "# CAG/CVT coefficient\n"
        if variable == "cvt"    : string += "# CVT coefficient\n"
        if variable == "cfs"    : string += "# TST correction factors/coefficients\n"
        if variable == "rcons"  : string += "# Rate constants (in au)\n"
        for key,data in sorted(dall[variable].items()):
            if data is None: continue
            string += "start_%s %s\n"%(variable,key)
            if variable == "pfn":
               V0,V1   = data[0:2]
               data    = data[2]
               string += "   %.8f\n"%V0
               string += "   %.8f\n"%V1
               string += "".join("   %7.2f  %.5E\n"%(T,value) for T,value in zip(ltemp,data))
            elif variable == "ctr":
               PIB,RRHO,ELE = data
               for idx,T in enumerate(ltemp):
                   string += "   %7.2f  %.5E  %.5E  %.5E\n"%(T,PIB[idx],RRHO[idx],ELE[idx])
            elif variable in ["gibbs1cm3","gibbs1bar"]:
               string += "".join("   %7.2f  %.8f\n"%(T,value) for T,value in zip(ltemp,data))
            else:
               for T,val in zip(ltemp,data): string += "   %7.2f  %s\n"%(T,fncs.eformat(val,5))
              #string += "".join("   %7.2f  %.5E\n"%(T,val) for T,val in zip(ltemp,data))
            string += "end_%s\n\n"%variable
    with open(filename,'w') as asdf: asdf.write(string)
#====================================================#


#====================================================#
def read_frozen(filename):
    frozen_symbols = []
    frozen_xcc     = []
    if not os.path.exists(filename): return frozen_symbols, frozen_xcc
    with open(filename,'r') as asdf: lines = asdf.readlines()
    for line in lines:
        line = line.strip()
        if line == "": continue
        symbol,x,y,z = line.split()
        frozen_symbols.append(symbol)
        frozen_xcc += [float(x),float(y),float(z)]
    return frozen_symbols, frozen_xcc
#----------------------------------------------------#
def write_frozen(filename,frozen_xcc,frozen_symbols):
    if len(frozen_symbols) == 0: return
    string = ""
    for idx,symbol in enumerate(frozen_symbols):
        x,y,z = xyz(frozen_xcc,idx)
        string += " %-2s  %+9.5f  %+9.5f  %+9.5f\n"%(symbol,x,y,z)
    with open(filename,'w') as asdf: asdf.write(string)
#====================================================#



#====================================================#
def read_rstconv():
    drstconv = {}
    if not os.path.exists(PN.CFILE): return drstconv
    with open(PN.CFILE,'r') as asdf: lines = asdf.readlines()
    for line in lines:
        line = line.strip()
        if line == "": continue
        itc, lowtemp, useics, scterr, converged = line.split()
        drstconv[itc] = (float(lowtemp),useics,float(scterr),converged)
    return drstconv
#----------------------------------------------------#
def write_rstconv(drstconv):
    lformat = "%s  %.6f %s  %.6f  %s\n"
    string  = "".join([lformat%(itc,lowtemp,useics,scterr,converged)
                       for itc,(lowtemp,useics,scterr,converged)
                       in  sorted(drstconv.items())])
    with open(PN.CFILE,'w') as asdf: asdf.write(string)
#====================================================#

