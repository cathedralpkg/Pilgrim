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
| Sub-module :  optINPUT           |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#--------------------------------------------------#
import os
import copy
import readline
import sys
import time
#--------------------------------------------------#
from   common.dicts       import dpt_im
#--------------------------------------------------#
from   common.fncs        import print_string
from   common.fncs        import remove_float
from   common.fncs        import uniquify_flist
from   common.fncs        import fill_string
from   common.fncs        import is_string_valid
#--------------------------------------------------#
from   common.MyCompleter import MyCompleter
#--------------------------------------------------#
from   common.physcons    import AMU, ML
#--------------------------------------------------#
from   modpilgrim.optGATHER          import ls_struc
#--------------------------------------------------#
import modpilgrim.names    as     PN
import modpilgrim.pilrw     as     RW
#--------------------------------------------------#
from   modpilgrim.diverse            import get_input_data
from   modpilgrim.diverse            import status_check
from   modpilgrim.diverse            import ffchecking
from   modpilgrim.diverse            import calc_fwdir
#--------------------------------------------------#
from   modpilgrim.pilesso                import get_templates
#--------------------------------------------------#
from   modpilgrim.PathVars           import PathVars
#--------------------------------------------------#
from modpilgrim.helps import MENU
from modpilgrim.helps import HELP_input_temp    as HELP_temp
from modpilgrim.helps import HELP_input_chem    as HELP_chem
from modpilgrim.helps import HELP_input_path    as HELP_path
from modpilgrim.helps import HELP_input_kmc     as HELP_kmc
from modpilgrim.helps import HELP_input_dlevel  as HELP_dlevel
from modpilgrim.helps import HELP_input_isomass as HELP_isomass
from modpilgrim.helps import HELP_input_struc   as HELP_struc
#--------------------------------------------------#

END     = ["end","..","exit","exit()","end()"]
ILINE1  = "    > "
IBLANK1 = "      "
ILINE2  = "      >> "
IBLANK2 = "         "

#===============================================================#
# Activate autocomplete
#===============================================================#
def set_completer(case=0):
    if   case == 0: options = []
    elif case == 1: options = ["add","mod","rm","ls","help",]+END
    elif case == 2: options = END
    else          : exit("Unknown case for function 'set_completer'...")
    completer = MyCompleter(options)
    readline.set_completer(completer.complete)
    readline.parse_and_bind('tab: complete')
#===============================================================#



#===============================================================#
def ls_temp(ltemp):
    for idx in range(0,len(ltemp),8):
        print(IBLANK2+"  ".join("%7.2f"%T for T in ltemp[idx:idx+8]))
#---------------------------------------------------------------#
def add_temp(data,ltemp=[]):
    set_completer(2)
    copy_ltemp = list(ltemp)
    try:
      for string in data:
          if "range" in string:
             Ti,Tf,dT = string.split("range(")[1].split(")")[0].split(",")
             Ti,Tf,dT = float(Ti),float(Tf),float(dT)
             nsteps = int(round(1.0*(Tf-Ti)/dT))+1
             ltemp += [Ti+idx*dT for idx in range(nsteps)] 
          else:
             ltemp.append(float(string))
      ltemp = uniquify_flist(ltemp,eps=1e-3)
    except:
        print(IBLANK1+"invalid line...")
    if copy_ltemp != ltemp:
       print(IBLANK1+"list of temperatures was modified!")
       ls_temp(ltemp)
    return sorted(ltemp)
#---------------------------------------------------------------#
def rm_temp(data,ltemp=[]):
    copy_ltemp = list(ltemp)
    eps = 0.001 # kelvin
    if "all" in data:
        ltemp = []
    elif "from" in data:
        s_from, Ti, s_to, Tf = data
        Ti, Tf = float(Ti), float(Tf)
        ltemp = [T for T in ltemp if (T<Ti-eps) or (T>Tf+eps)]
    else:
        for T in data:
            ltemp = remove_float(float(T),ltemp,eps=1e-3)
    if copy_ltemp != ltemp:
       print(IBLANK1+"list of temperatures was modified!")
       ls_temp(ltemp)
    return sorted(ltemp)
#===============================================================#


#===============================================================#
def ls_path(dpath,targets=None):
    string = ""
    for target,pathvars in sorted(dpath.items()):
        if targets is not None:
           if target not in targets: continue
        string += pathvars.string4inp(target)
        string += "\n"
    if string.endswith("\n"): string = string[:-1]
    if string == "": return
    for line in (string).split("\n"): print(IBLANK2+line.strip())
#---------------------------------------------------------------#
def add_path(targets,dpath,dctc,dtesLL,lts,mod=False):
    set_completer(2)
    # defaults
    if "*" in targets or targets == []:
        if mod: targets = dpath.keys()
        else  : targets = lts
    targets = [target for target in targets if target in lts]
    # no targets? Finish
    if len(targets) == 0: return dpath, dtesLL
    # Generate PathVars object for all of them
    for target in targets:
        ctc, itc = PN.name2data(target)
        if ctc not in dctc.keys(): return dpath, dtesLL
        if mod and (ctc in dpath.keys()): continue
        if ctc in dpath.keys(): continue
        dpath[ctc] = PathVars("mep")
        # get fwdir (this may take some time...)
        gtsfile = dctc[ctc].gtsfiles()[0]
        try   : fwdir = calc_fwdir(gtsfile)
        except: continue
        dpath[ctc]._fwdir = fwdir
    # print path(s)
    ls_path(dpath,targets)
    # modify them
    while True:
          try:
            line = input(ILINE2).strip()
            if line.split()[0] in END: break
            if "=" not in line:
                line = line.split()
            else:
                line = [ string.strip() for string in line.split("=")]
            if len(line) != 2: raise Exception
            var, val = line
            var = var.lower()
            if var not in ["sbw","sfw","ds","hsteps","paral","scterr","sctmns"]: raise Exception
            for target in targets:
                dpath[target].setvar(var,val)
          except KeyboardInterrupt:
            print("")
            break
          except:
            print(IBLANK2+"invalid line...")
            print("")
            ls_path(dpath,targets)
            continue
          # print 
          ls_path(dpath,targets)
    # update dtesLL
    for ctc in dpath.keys():
        try   : ch,mtp = dctc[ctc]._ch, dctc[ctc]._mtp
        except: ch,mtp = 0, 1
        # templates by default
        default_templates = get_templates(ch,mtp,"LL")
        # add default?
        for software,string in default_templates.items():
            dtesLL[software] = dtesLL.get(software,{})
            if ctc in dtesLL[software].keys(): continue
            dtesLL[software][ctc] = string
    return dpath, dtesLL
#---------------------------------------------------------------#
def rm_path(data,dpath):
    ml = max([len(target) for target in data]+[1])
    for target in data:
        if target in dpath.keys():
           dpath.pop(target)
           print(IBLANK1+"%%-%is was removed"%ml%target)
    return dpath
#===============================================================#


#===============================================================#
def ls_chem(dchem):
    for rname in dchem.keys():
        Rs, TS, Ps = dchem[rname]
        if TS is None: TS = ""
        print(IBLANK2+"'%s' : '%s --> %s --> %s'"%(rname,"+".join(Rs),TS,"+".join(Ps)))
#---------------------------------------------------------------#
def add_chem(data,dchem={}):
    set_completer(2)
    try:
      line   = " ".join(data)
      if ":" not in line: raise Exception
      rname  = line.split(":")[0]
      rformu = line.split(":")[1]
      # strip strings
      rname  = rname.strip()
      rformu = rformu.strip()
      if rname  == "": raise Exception
      if rformu == "": raise Exception
      if not is_string_valid(rname,extra="_"):
         print(IBLANK1+"invalid name for the reaction (%s)!"%rname)
         return dchem
      # get reactant(s), TS and product(s)
      reac, ts, prod = rformu.split("-->")
      reac = [R.strip() for R in reac.split("+")]
      ts   = ts.strip()
      prod = [P.strip() for P in prod.split("+")]
      if "" in reac or "??" in reac: reac = []
      if "" in prod or "??" in prod: prod = []
      if "" == ts   or "??" == ts  : ts   = None
      dchem[rname] = (reac,ts,prod)
    except:
        print(IBLANK1+"invalid line...")
    return dchem
#---------------------------------------------------------------#
def rm_chem(data,dchem={}):
    ml = max([len(target) for target in data]+[1])
    for reaction in data:
       if reaction in dchem.keys():
           dchem.pop(reaction)
           print(IBLANK1+"%%-%is was removed"%ml%reaction)
    return dchem
#===============================================================#


#===============================================================#
def ls_kmc(dkmc,kmc_name=None):
    string = "\n"
    if dkmc is None:
        string += "no kmc defined"
    else:
        if kmc_name is None: tolist = sorted(dkmc.keys())
        else               : tolist = [kmc_name]
        for kmc_name in tolist:
            tkmc = dkmc[kmc_name]
            ipops, rcs, psteps, volume, timeunits, excess = tkmc
            # get l1 and l2 for nice formatting
            try   : l1 = max(len(ctc) for ctc in ipops.keys())
            except: l1 = 0
            try   : l2 = max(len(reac) for reac in rcs.keys())
            except: l2 = 0
            # generating string
            string += "--> kmc name: %s\n"%kmc_name
            string += "volume    = %.2e [in mL]\n"%(volume*ML)
            string += "psteps    = %i\n"%psteps
            string += "timeunits = %s\n"%timeunits
            string += "initial populations (number of particles):\n"
            for species,ipop in sorted(ipops.items()):
                string += ("   pop0(%%-%is)  = %%.3e"%l1)%(species,ipop)
                if species in excess: string += " (excess)"
                string += "\n"
            string += "Rate constants:\n"
            for reaction,(rctype,weight,rcdata) in sorted(rcs.items()):
                if weight == 1: string += ("   k(%%-%is) = %%-7s  "%l2)%(reaction,rctype.lower())
                else          : string += ("   k(%%-%is)*%%i = %%-7s  "%l2)%(reaction,weight,rctype.lower())
                if "analytic" in rctype.lower() and rcdata is not None:
                    string += " ".join(["%.3e"%coef for coef in rcdata])
                string += "\n"
            string += "\n"
    # print data
    for line in (string).split("\n"): print(IBLANK2+line)
#---------------------------------------------------------------#
def mod_kmc(val,dkmc,dchem):
    set_completer(2)
    # some initial checkings
    if dkmc is None :
       print(IBLANK2+"First initialize KMC with 'add kmc'...")
       return dkmc
    if len(val) == 0:
       print(IBLANK2+"A KMC must be selected...")
       return dkmc
    kmc_name = val[0]
    if kmc_name not in dkmc.keys():
       print(IBLANK2+"Selected KMC (%s) does not exist..."%kmc_name)
       return dkmc
    while True:
        ls_kmc(dkmc,kmc_name)
        (ipops,rcs,psteps,volume,timeunits,excess) = dkmc[kmc_name]
        # user introduces info
        line = input(ILINE2).strip()
        # exit?
        if line in END : return dkmc
        # proper format?
        if "=" not in line: continue
        try:
           # extract data from line
           command,data = [ss.strip() for ss in line.split("=")]
           data = data.split()
           # user asks for listing
           if   command == "ls"       : ls_kmc(dkmc,kmc_name); continue
           # modify data
           if   command == "psteps"   : psteps    = int(data[0])
           elif command == "volume"   : volume    = float(data[0])/ML
           elif command == "timeunits": timeunits = data[0]
           elif "pop0" in command:
                species = command.split("pop0(")[1].split(")")[0].strip()
                ipop = data[0]
                # remove species
                if ipop.lower() == "none":
                    if species in ipops.keys(): ipops.pop(species)
                    if species in excess      : excess.remove(species)
                # add species
                else:
                    ipops[species] = float(ipop)
                # no excess
                if len(data) == 1:
                    if species in excess: excess.remove(species)
                # in excess
                elif len(data) == 2 and data[1].lower() == "excess":
                    excess.append(species)
           elif "k(" in command:
               rcname = command.split("k(")[1].split(")")[0].strip()
               weight = 1
               if "*" in command:
                   weight = int(command.split("*")[1].split()[0])
               rctype = data[0].lower()
               rcdata = None
               # remove rctype?
               if (rctype == "none") and (rcname in rcs.keys()): rcs.pop(rcname)
               # analytic
               if   rctype == "analytic1":
                   rcdata = [float(coef) for coef in data[1:3]]
                   assert len(rcdata) == 2
               elif rctype == "analytic2":
                   rcdata = [float(coef) for coef in data[1:4]]
                   assert len(rcdata) == 3
               elif rctype == "analytic3":
                   rcdata = [float(coef) for coef in data[1:5]]
                   assert len(rcdata) == 4
               elif rctype == "analytic4":
                   rcdata = [float(coef) for coef in data[1:6]]
                   assert len(rcdata) == 5
               rcs[rcname] = (rctype,weight,rcdata)

           if timeunits not in "fs,ps,mcs,ms,s,min,hr".split(","): timeunits = "ps"
           dkmc[kmc_name] = (ipops,rcs,psteps,volume,timeunits,excess)
        except: continue
    return dkmc
#---------------------------------------------------------------#
def add_kmc(val,dkmc,dchem):
    set_completer(2)
    if dkmc is None: dkmc = {}
    if len(val) == 0: kmc_name = "mechanism"
    else            : kmc_name = val[0]
    # Already exists
    if kmc_name in dkmc.keys():
       print(IBLANK2+"KMC '%s' already exists... Try with 'mod kmc [kmc_name]'..."%kmc_name)
       return dkmc
    # Default values
    ipops     = {}
    rcs       = {}
    psteps    = 1000
    volume    = 1.0/ML
    timeunits = "ps"
    excess    = []
    # add compounds to ipops and reactions to rcs
    for reaction,(Rs,TS,Ps) in dchem.items():
        for compound in Rs+Ps:
            ipops[compound] = ipops.get(compound,0.0)
        # omit A<-->A reactions
        if sorted(Rs) == sorted(Ps): continue
        # add reactions to rcs
        rcnames = [key.split(".")[0] for key in rcs.keys()]
        if reaction in rcnames: continue
        rcs[reaction]  = rcs.get(reaction,("TST",1,None))
    # modify by user
    dkmc[kmc_name] = (ipops, rcs, psteps, volume, timeunits, excess)
    dkmc = mod_kmc([kmc_name],dkmc,dchem)
    return dkmc
#---------------------------------------------------------------#
def rm_kmc(data,dkmc):
    ml = max([len(target) for target in data]+[1])
    for target in data:
        if target in dkmc.keys():
           dkmc.pop(target)
           print(IBLANK1+"%%-%is was removed"%ml%target)
    print("")
    return dkmc
#===============================================================#


#===============================================================#
def ls_ind_struc(CTC,one=True):
    string  = "spname ==> %s\n"%CTC._ctc
    divi    = "-"*len(string)+"\n"
    string  = divi+string+divi
    if not one:
       string += "freqscal    = %.4f\n"%CTC._fscal
       for idx,(itc,weight) in enumerate(CTC._itcs):
           string += "weight(%3s) = %i # %s\n"%(itc,weight,str(CTC._lpg[idx]))
    else:
       string += "root        = %s\n"%CTC._root
       string += "freqscal    = %.4f\n"%CTC._fscal
       for idx,(itc,weight) in enumerate(CTC._itcs):
           string += "weight(%3s) = %i # %s\n"%(itc,weight,str(CTC._lpg[idx]))
       if CTC._anh is None:
          string += "anharfile   = -\n"
       else:
          string += "anharfile   = %s\n"%CTC._anh
       if CTC._diso == {}:
           string += "iso         = -\n"
       for key in CTC._diso.keys():
           if key == "*": xx = "    "
           else         : xx = ".%3s"%key
           string += "iso%s     = %s\n"%(xx," ".join(CTC._diso[key]))
    string += divi
    # print data
    for line in (string).split("\n"): print(IBLANK2+line)
#---------------------------------------------------------------#
def mod_struc(target,dctc,dimasses):
    if len(target) == 0:
       print(IBLANK1+"A structure has to be selected (mod struc SPNAME)!"%target)
       ls_struc(dctc)
       return dctc
    targets = []
    for t_i in target:
        if t_i == "*":
           targets = dctc.keys()
           break
        if t_i not in dctc.keys():
           print(IBLANK1+"'%s' is unknown!"%t_i)
           continue
        targets.append(t_i)
    if len(targets) == 0: return dctc
    # read gts files
    dC1 = {}
    for target in targets:
        dctc[target].read_gtsfiles(dimasses)
        if len(targets) == 1: ls_ind_struc(dctc[target])
        else                : ls_ind_struc(dctc[target],False)
        allC1 = [itc for (itc,w),pg in zip(dctc[target]._itcs,dctc[target]._lpg) if pg == "C1"]
        dC1[target] = allC1

    while True:
          try:
            line = input(ILINE2).strip()
            if line.split()[0] in END: break
            line = line.replace("="," ").split()
            if len(line) == 0: raise Exception
            var   = line[0]
            value = line[1:]
            for target in targets:
                # mod freqscal
                if var == "freqscal" : dctc[target]._fscal = float(value[0])
                # mod weights
                if var.startswith("weight"):
                   # get itcs
                   itcs = var.split("(")[1].split(")")[0]
                   if itcs == "all": itcs = dC1[target]
                   else:
                      if "," in itcs: itcs = itcs.split(",")
                      else          : itcs = [itcs]
                      itcs2 = []
                      for itc in itcs:
                          if "-" in itc:
                             itc1,itc2 = itc.split("-")
                             itcs2 += ["%003i"%itc for itc in range(int(itc1),int(itc2)+1,1)]
                          else:
                             itcs2 += [itc]
                      itcs = itcs2
                   # only those with C1 symmetry
                   itcs = [itc for itc in itcs if itc in dC1[target]]
                   # get weight
                   weight = int(value[0])
                   if not 1 <= weight <= 2: weight = 1
                   # apply
                   for itc_i in itcs:
                       for idx,(itc_j,w_j) in enumerate(dctc[target]._itcs):
                           if itc_i == itc_j: dctc[target]._itcs[idx] = (itc_i,weight)
            if len(targets) != 1:
                for target in targets: ls_ind_struc(dctc[target],False)
                continue
            # only one target
            if var == "root"     : dctc[target]._root  = value[0]
            if var == "anharfile":
                 if value == []: dctc[target]._anh   = None
                 else          : dctc[target]._anh   = value[0]
            target = targets[0]
            if var.startswith("iso"):
               if "." in var: var,xx = var.split(".")
               else         : var,xx = var,"*"
               if value == [] and xx in dctc[target]._diso.keys():
                    dctc[target]._diso.pop(xx)
               elif value != []:
                    dctc[target]._diso[xx] = value
            if var == "copywith":
               value = " ".join(value)
               isomods,new_ctc = value.split(" as ")
               new_ctc = new_ctc.strip()
               isomods  = isomods.split()
               dctc[new_ctc] = copy.deepcopy(dctc[target])
               dctc[new_ctc]._ctc = new_ctc
               dctc[new_ctc]._anh = None
               if len(isomods) == 0: dctc[new_ctc]._diso = {}
               else                : dctc[new_ctc]._diso = {"*":isomods}
               ls_struc(dctc)
               break
            ls_ind_struc(dctc[target])
          except KeyboardInterrupt:
            print("")
            break
          except:
            print(IBLANK2+"invalid line...")
            print("")
            continue
    return dctc
#---------------------------------------------------------------#
def rm_struc(data,dctc):
    ml = max([len(target) for target in data]+[1])
    for target in data:
        if target in dctc.keys():
           dctc.pop(target)
           print(IBLANK1+"%%-%is was removed"%ml%target)
    print("")
    return dctc
#===============================================================#

#===============================================================#
def ls_isomass(dimasses):
    string = ""
    if dimasses != {}:
        string += "Isotopic masses:\n"
        sorted_list = sorted([(imass,label) for label,imass in dimasses.items()])
        ml = max([len(key) for key in dimasses.keys()])
        for imass, label in sorted_list:
            string += "    %s = %10.4f amu\n"%("%%-%is"%ml%label,imass*AMU)
    for line in string.split("\n"): print(IBLANK2+line)
#---------------------------------------------------------------#
def add_isomass(value,dimasses):
    value = " ".join(value)
    # modify by user
    try:
       # isomass??
       if "=" in value: label, mass = value.split("=")
       else           : label, mass = value.split()
       dimasses[label.strip()] = float(mass)/AMU
    except:
       print(IBLANK2+"invalid line...")
       print("")
    ls_isomass(dimasses)
    return dimasses
#---------------------------------------------------------------#
def rm_isomass(data,dimasses):
    ml = max([len(target) for target in data]+[1])
    for target in data:
        if target in dimasses.keys():
           dimasses.pop(target)
           print(IBLANK1+"%%-%is was removed"%ml%target)
    print("")
    return dimasses
#===============================================================#

#===============================================================#
def ls_dlevel(ddlevel):
    string = ""
    ml = max([len(ctc) for ctc in ddlevel.keys()]+[1])
    for key,points in sorted(ddlevel.items()):
        string += "* %%-%is "%ml%key
        points = [str(pp) for pp in points if pp != "sp"]
        if len(points) != 0:
            string += "   s={%s}"%",".join(points)
        string += "\n"
    for line in string.split("\n"): print(IBLANK2+line)
#---------------------------------------------------------------#
def add_dlevel(data,ddlevel,dctc,dtesHL):
    set_completer(2)
   #print IBLANK2+"Example input file will be generated at the end!"
    ctc_sp0 = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 0])
    ctc_sp1 = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 1])
    for ctcs,sptype in [(ctc_sp0,0),(ctc_sp1,1)]:
        for ctc in ctcs:
            ch,mtp = dctc[ctc]._ch, dctc[ctc]._mtp
            itcs   = dctc[ctc]._itcs
            for itc,weight in sorted(itcs):
                if ctc+"."+itc in ddlevel.keys(): continue
                if sptype == 0: ddlevel[ctc+"."+itc] = []
                if sptype == 1: ddlevel[ctc+"."+itc] = ["auto_3_3"]
            # templates by default
            default_templates = get_templates(ch,mtp,"HL")
            # add default?
            for software,string in default_templates.items():
                dtesHL[software] = dtesHL.get(software,{})
                if ctc in dtesHL[software].keys(): continue
                dtesHL[software][ctc] = string
    ls_dlevel(ddlevel)
    return ddlevel, dtesHL
#---------------------------------------------------------------#
def rm_dlevel(data,ddlevel={}):
    ml = max([len(target) for target in data]+[1])
    for target in data:
        if target in ddlevel.keys():
           ddlevel.pop(target)
           print(IBLANK1+"%%-%is was removed"%ml%target)
    return ddlevel
#===============================================================#


#===============================================================#
def main(idata,status):

    stat2check = []
    mustexist  = []
    tocreate   = []

    #-------------------------------------------------------#
    # Read Pilgrim input files and check file/folder status #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, dkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,[])
    if fstatus == -1: exit()
    # special warning
    if status[0] == -1:
        print("     WARNING! '%s' file not found."%PN.IFILE1)
        print("     In order to generate it, run")
        print("        pilgrim --gather")
    # existency of folders
    fstatus = ffchecking(mustexist=[],tocreate=[])
    if fstatus == -1: exit()
    #-------------------------------------------------------#

    # print menu
    print_string(MENU,1)

    # setup dpath
    for ctc in dpath.keys():
        dpath[ctc].setup1()
        dpath[ctc].setup2()


    # interactive menu
    dbool = {}
    while True:
          set_completer(1)
          line = input(ILINE1).split()
          if len(line) == 0: continue
          command  = line[0].lower()
          try   : var = line[1].lower()
          except: var = None
          try   : val = line[2:]
          except: val = None
          # update list of transition states
          lts  = [ctc for ctc in dctc.keys() if dctc[ctc]._type == 1]
          # user ask for help
          if   command == "help":
             if   var is None     : print_string(MENU,3)
             elif var == "temp"   : print_string(HELP_temp,6)
             elif var == "path"   : print_string(HELP_path,6)
             elif var == "chem"   : print_string(HELP_chem,6)
             elif var == "kmc"    : print_string(HELP_kmc ,6)
             elif var == "dlevel" : print_string(HELP_dlevel,6)
             elif var == "isomass": print_string(HELP_isomass,6)
             elif var == "struc"  : print_string(HELP_struc,6)
          # user ask for list
          elif command == "ls":
             if   var == "struc"  : ls_struc(dctc)
             elif var == "temp"   : ls_temp(ltemp)
             elif var == "path"   : ls_path(dpath)
             elif var == "chem"   : ls_chem(dchem)
             elif var == "kmc"    : ls_kmc(dkmc)
             elif var == "dlevel" : ls_dlevel(ddlevel)
             elif var == "isomass": ls_isomass(dimasses)
             else                 : continue
          # add
          elif command == "add":
               if   var == "temp"   : ltemp          = add_temp(val,ltemp)
               elif var == "path"   : dpath,dtesLL   = add_path(val,dpath,dctc,dtesLL,lts)
               elif var == "chem"   : dchem          = add_chem(val,dchem)
               elif var == "kmc"    : dkmc           = add_kmc(val,dkmc,dchem)
               elif var == "dlevel" : ddlevel,dtesHL = add_dlevel(val,ddlevel,dctc,dtesHL)
               elif var == "isomass": dimasses       = add_isomass(val,dimasses)
               dbool[var] = True
          # mod
          elif command == "mod":
               if   var == "path"  : dpath,dtesLL   = add_path(val,dpath,dctc,dtesLL,lts,mod=True)
               elif var == "struc" : dctc           = mod_struc(val,dctc,dimasses)
               elif var == "kmc"   : dkmc           = mod_kmc(val,dkmc,dchem)     
               dbool[var] = True
          # rm
          elif command == "rm":
               if   var == "temp"   : ltemp    = rm_temp(val,ltemp)
               elif var == "path"   : dpath    = rm_path(val,dpath)
               elif var == "chem"   : dchem    = rm_chem(val,dchem)
               elif var == "dlevel" : ddlevel  = rm_dlevel(val,ddlevel)
               elif var == "struc"  : dctc     = rm_struc(val,dctc)
               elif var == "isomass": dimasses = rm_isomass(val,dimasses)
               elif var == "kmc"    : dkmc     = rm_kmc(val,dkmc)
               dbool[var] = True
          # user ask for finishing
          elif command in END:
             if dbool != {}:
                print("")
                print("  (Re)Writing files:")
                if dbool.get("isomass",False) or dbool.get("struc",False):
                   print("     %s"%PN.IFILE1)
                   RW.write_ctc(dctc,dimasses)
                if dbool.get("temp",False):
                   print("     %s"%PN.IFILE2)
                   RW.write_temp(ltemp)
                if dbool.get("path",False):
                   print("     %s"%PN.IFILE3)
                   RW.write_path(dpath)
                if dbool.get("path",False) or dbool.get("dlevel",False):
                   print("     %s"%PN.IFILE4)
                   RW.write_tes(dtesLL,dtesHL)
                if dbool.get("chem",False):
                   print("     %s"%PN.IFILE5)
                   RW.write_chem(dchem)
                if dbool.get("kmc",False):
                   print("     %s"%PN.IFILE6)
                   RW.write_kmc(dkmc)
                if dbool.get("dlevel",False):
                   print("     %s"%PN.IFILE7)
                   RW.write_dlevel(ddlevel)
                print("")
             break
          # user ask for nothing
          else: continue
#===============================================================#

