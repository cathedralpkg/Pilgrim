#!/usr/bin/python3.6
'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2020.2
License     : MIT/x11

Copyright (c) 2020, David Ferro Costas (david.ferro@usc.es) and
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
| Module     :  common             |
| Sub-module :  files              |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different functions
related to reading/writing diverse types
of files
''' 

#=============================================#
import random
import os
import sys
import fcntl
import common.Exceptions as Exc
import common.fncs     as     fncs
#---------------------------------------------#
from   common.pgs      import get_pgs
from   common.dicts    import dpt_s2m
from   common.physcons import ANGSTROM, AMU, KCALMOL, CM2H, KJMOL
from   common.physcons import METER, JOULE, KB
#=============================================#


#=============================================#
def random_filename(n=10):
    characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789."
    while True:
       fname = [ch for ch in characters]
       random.shuffle(fname)
       fname = "".join(fname[0:n])
       if os.path.exists(fname): continue
       break
    return fname
#=============================================#

#=============================================#
def mkdir_recursive(path):
    # already exists?
    if os.path.exists(path): return
    # try to create it
    try: os.mkdir(path)
    # anything failed?
    except:
      # try with sub-path
      sub_path = os.path.dirname(path)
      if not os.path.exists(sub_path): mkdir_recursive(sub_path)
      if not os.path.exists(    path): os.mkdir(path)
#=============================================#

#=============================================#
# Functions related to reading/writing files  #
#=============================================#
def read_file(filename):
    '''
    Just returns a list of lines
    '''
    lines = []
    if os.path.exists(filename):
       with open(filename,'r') as ff: lines = ff.readlines()
    return lines
#---------------------------------------------#
def write_file(filename,string,mode="w"):
    with open(filename,mode) as ff: lines = ff.write(string)
#=============================================#


#===============================================#
# Function(s) to deal with: XYZ file(s)         #
#===============================================#
def read_xyz(xyzfile):
    '''
    Read standard .xyz file
    Returns data in bohr
    '''
    xcc,symbs,masses = [], [], []
    lines   = read_file(xyzfile)
    natoms  = int(lines[0])
    comment =     lines[1]
    for line in lines[2:]:
        datainline = line.split()
        if len(datainline) == 4:
           symb, x, y, z = datainline
           mass = dpt_s2m[symb]
        elif len(datainline) == 5:
           symb, x, y, z, mass = datainline
           mass = float(mass) / AMU # mass in au
        else: continue
        masses.append(mass)
        symbs.append(symb)
        xcc += [float(x),float(y),float(z)]
    if "[BOHR]" not in comment: xcc = [xi/ANGSTROM for xi in xcc]
    return xcc, symbs, masses
#-----------------------------------------------#
def write_xyz(filename,xcc,symbs,comment="info line",mode="w",ext=".xyz"):
    '''
    Write standard .xyz file
    * xcc in bohr
    * mode = "a" or "w"
    '''
    if not filename.endswith(ext): filename += ext
    nat = len(symbs)
    if len(xcc) == nat: xcc = flatten_llist(xcc)
    # symbs is really a list of symbols or of atonums?
    if type(symbs[0]) == type(1):
        symbs = fncs.get_symbols(symbs)
    # Generate string
    string  = ""
    string += "%i\n"%nat
    string += "%s\n"%comment
    for idx in range(nat):
        symbol = symbs[idx]
        xx = fncs.x(xcc,idx)*ANGSTROM
        yy = fncs.y(xcc,idx)*ANGSTROM
        zz = fncs.z(xcc,idx)*ANGSTROM
        string += " %2s   %+10.5f   %+10.5f   %+10.5f\n"%(symbol,xx,yy,zz)
    # Write file
    write_file(filename,string,mode)
#===============================================#


#===============================================#
# Function(s) to deal with: MOLDEN file(s)      #
#===============================================#
def read_molden(filename):
    xcc, symbols, freqs, evecs = [], [], [], []
    lines = read_file(filename)
    bool_freq = False
    bool_geom = False
    bool_evec = False
    for line in lines:
        if bool_geom:
           symbol, xx, yy, zz = line.split()
           symbols.append(symbol)
           xcc += [float(xx),float(yy),float(zz)]
        if bool_freq:
           freqs.append(float(line))
        if bool_evec:
           if "vibration" in line:
              evec = []
           else:
              exx, eyy, ezz = line.split()
              evec += [float(exx),float(eyy),float(ezz)]
        if "[FR-COORD]"      in line:
           bool_geom = True
           bool_freq = False
           bool_evec = False
        if "[FREQ]"          in line:
           bool_geom = False
           bool_freq = True
           bool_evec = False
        if "[FR-NORM-COORD]" in line:
           bool_geom = False
           bool_freq = False
           bool_evec = True
    return xcc, symbols, freqs, evecs
#-----------------------------------------------#
def write_molden(filename,xcc,symbs,freqs,evecs):
    '''evecs NOT in mass-scaled!'''
    natoms = len(symbs)
    nfreqs = len(freqs)
    STRING  = "[Molden Format]\n"
    STRING += "[FR-COORD] # Coordinates in bohr\n"
    for at in range(natoms):
        symbol = symbs[at]
        x,y,z  = fncs.xyz(xcc,at)
        STRING += " %2s  %+11.6f  %+11.6f  %+11.6f \n"%(symbol,x,y,z)
    if freqs not in [None,[]]:
       STRING += "[FREQ] # Frequencies in cm^-1\n"
       for idx in range(nfreqs):
           freq = fncs.afreq2cm(freqs[idx])
           STRING += " %9.4f\n"%freq
    if evecs not in [None,[]]:
       STRING += "[FR-NORM-COORD] # Displacements in bohr\n"
       for idx in range(nfreqs):
           STRING += "vibration  %i\n"%(idx+1)
           evec = evecs[idx]
           for at in range(natoms):
               # It may fail with frozen atoms...
               try   : vx, vy, vz = fncs.xyz(evec,at)
               except: vx, vy, vz = 0.0, 0.0, 0.0
               STRING += "   %+9.3f  %+9.3f  %+9.3f\n"%(vx,vy,vz)
    write_file(filename,STRING)
#===============================================#


#===============================================#
# Function(s) to deal with: GTS file(s)         #
#===============================================#
def read_gtsfile(gtsfile):
    lines = read_file(gtsfile)
    # check extension
    end_gts = (gtsfile.lower()).endswith(".gts")
    if (not end_gts): raise Exc.FileType(Exception)
    # CHECK FILE IS GTS
    correct = False
    for line in lines:
        if line.startswith("start_cc"):
           correct = True
           break
    if not correct:
       exception = Exc.FileIsNotGTS(Exception)
       exception._var = gtsfile
       raise exception
    # Read cartesian coordinates
    info_cc = fncs.extract_lines(lines,"start_cc","end_cc")
    xcc,atonums = [],[]
    for line in info_cc:
        atonum,xx,yy,zz = line.split()
        xcc += [float(xx),float(yy),float(zz)]
        atonums.append(int(atonum))
    # Read basic information
    info_basic = fncs.extract_lines(lines,"start_basic","end_basic")
    ch, mtp, E, pgroup, rotsigma = None, None, None, None, None
    for line in info_basic:
        if line.startswith("pointgroup"):   pgroup   =       line.split()[1]
        if line.startswith("charge"):       ch       =   int(line.split()[1])
        if line.startswith("multiplicity"): mtp      =   int(line.split()[1])
        if line.startswith("rotsigma"):     rotsigma =   int(line.split()[1])
        if line.startswith("energy"):       E        = float(line.split()[1])
    # Read cartesian gradient
    info_grad = fncs.extract_lines(lines,"start_grad","end_grad")
    gcc = []
    for line in info_grad:
        gx , gy , gz = line.split()
        gcc += [float(gx),float(gy),float(gz)]
    # Read Hessian matrix
    info_hess = fncs.extract_lines(lines,"start_hess","end_hess")
    Fcc = []
    for line in info_hess:
        Fcc += [float(Fij) for Fij in line.split()]
    # Read frequencies (only special cases)
    info_freqs = fncs.extract_lines(lines,"start_freqs","end_freqs")
    freq_list = []
    for line in info_freqs:
        freq_list += [float(freqs)*CM2H for freqs in line.split()]
    # Return data
    masses = fncs.atonums2masses(atonums)
    return xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma,freq_list
#-----------------------------------------------#
def write_gtsfile(xcc,atonums,ch,mtp,E,pgroup,rotsigma,gcc,Fcc,gtsfile,freqs=None,level=""):
    if gcc   is None: gcc = []
    if Fcc   is None: Fcc = []
    if freqs is None: freqs = []
    nat      = len(atonums)
    str_gts  = ""
    # Write atomic numbers and cartesian coordinates
    if level != "": str_gts += "# level: %s\n"%level
    str_gts += "# Atomic number and non-scaled cartesian coordinates [bohr]\n"
    str_gts += "start_cc\n"
    for at in range(nat):
        atonum = atonums[at]
        xx,yy,zz = fncs.xyz(xcc,at)
        str_gts += "   %03i   %+15.8E  %+15.8E  %+15.8E\n"%(atonum,xx,yy,zz)
    str_gts += "end_cc\n\n"
    # Write basic data
    str_gts += "# Charge, multiplicity, energy [hartree],\n"
    str_gts += "# point group and rotational symmetry number\n"
    str_gts += "start_basic\n"
    str_gts += "   charge        %i\n"%ch
    str_gts += "   multiplicity  %i\n"%mtp
    str_gts += "   energy       %-+16.8f # Total energy in hartree\n"%E
    str_gts += "   pointgroup    %-5s           # Point group\n"%pgroup
    str_gts += "   rotsigma      %-5i           # Rotational sigma\n"%rotsigma
    str_gts += "end_basic\n\n"
    # Write cartesian gradiend
    if len(gcc) != 0:
       str_gts += "# Non-scaled cartesian gradient [hartree/bohr]\n"
       str_gts += "start_grad\n"
       for at in range(nat):
           gxx,gyy,gzz = fncs.xyz(gcc,at)
           str_gts += "   %+15.8E  %+15.8E  %+15.8E\n"%(gxx,gyy,gzz)
       str_gts += "end_grad\n\n"
    # Write force constant matrix (i.e. hessian matrix)
    if len(Fcc) != 0:
       # matrix in triangular form
       if len(Fcc) == 3*nat: Fcc = fncs.matrix2lowt(Fcc)
       str_gts += "# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]\n"
       str_gts += "# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...\n"
       str_gts += "start_hess\n"
       for idx in range(0,len(Fcc),5):
           line = "  ".join(["%+15.8E"%Fij for Fij in Fcc[idx:idx+5]])
           str_gts += "   %s\n"%line
       str_gts += "end_hess\n\n"
    # Write freqs
    elif len(freqs) != 0:
       str_gts += "start_freqs\n"
       for idx in range(0,len(freqs),5):
           line = "  ".join(["%7.2f"%(freq) for freq in freqs[idx:idx+5]])
           str_gts += "           %s\n"%line
       str_gts += "end_freqs\n\n"
    # write
    write_file(gtsfile,str_gts)
#===============================================#


#===============================================#
# Function(s) to deal with: RST file(s)         #
#===============================================#
def read_rst(rstfile):
    if not os.path.exists(rstfile):
        return None, None, {}
    lines = read_file(rstfile)
    lines = fncs.clean_lines(lines,"#",True)
    if len(lines) == 0:
       return None, None, {}
    # pathinfo lines
    lines_path   = fncs.extract_lines(lines,"start_pathinfo","end_pathinfo")
    for line in lines_path:
        if line.startswith("path ")  : path   = line.split()[1].lower()
        if line.startswith("mtype ") : path   = line.split()[1].lower()
        if line.startswith("mu ")    : mu     = float(line.split()[1]) / AMU
        if line.startswith("ds ")    : ds     = float(line.split()[1])
        if line.startswith("hsteps "): hsteps = int(line.split()[1])
        if line.startswith("cubic "):
            cubic = line.split()[1].lower()
            if cubic != "no": cubic = float(cubic)
            else            : cubic = None
    # structinfo lines
    lines_struct = fncs.extract_lines(lines,"start_structinfo","end_structinfo")
    bool1  = False
    bool2  = False
    atonums  = []
    masses   = []
    for line in lines_struct:
        if   line.startswith("atonums") :
             bool1,bool2 = True, False
             continue
        elif line.startswith("masslist") or line.startswith("masses"):
             bool1,bool2 = False, True
             continue
        elif line.startswith("ch "):
             bool1,bool2 = False, False
             ch  = int(line.split()[1])
        elif line.startswith("mtp "):
             bool1,bool2 = False, False
             mtp = int(line.split()[1])
        if bool1: atonums  += [int(atonum)     for atonum in line.split()]
        if bool2: masses   += [float(mass)/AMU for mass   in line.split()]
    # tuples of common data and of path
    tcommon = (ch,mtp,atonums,masses,mu)
    tpath   = (path,mu,ds,hsteps,cubic)
    # initialize structures
    drst   = {}
    mins   = +float("inf")
    maxs   = -float("inf")
    record = False
    lines  = iter(lines)
    # Read structures
    for line in lines:
        if "==>" in line:
           record = True
           t_value   = None
           data_line = line.split()[1:]
           meplabel  = data_line[0]
           if meplabel.endswith("bw"): meplabel = "bw"+meplabel.split("bw")[0]
           if meplabel.endswith("fw"): meplabel = "fw"+meplabel.split("fw")[0]
           variables = data_line[1].split("-")
           mepsvalue = float(data_line[2])
           Vmep      = float(data_line[3])
           if len(data_line) == 5: t_value  = float(data_line[4])
           xms,gms,Fms,v0,v1 = None,None,None,None,None
           # Read geometry
           if "x" in variables:
              xms = [] 
              while True:
                   line = next(lines)
                   if "---" in line: break
                   xms += line.split()
              xms = [float(xx) for xx in xms]
              #xcc = ms2cc_x(xms,masses,mu)
           # Read gradient
           if "g" in variables:
              gms = [] 
              while True:
                   line = next(lines)
                   if "---" in line: break
                   gms += line.split()
              gms = [float(gg) for gg in gms]
              #gcc = ms2cc_g(gms,masses,mu)
           # Read hessian
           if "F" in variables:
              Fms = [] 
              while True:
                   line = next(lines)
                   if "---" in line: break
                   Fms += line.split()
              Fms = [float(FF) for FF in Fms]
              Fms = fncs.lowt2matrix(Fms)
              #Fcc = ms2cc_F(Fms,masses,mu)
           # Read v0
           if "v0" in variables:
              v0  = [] 
              while True:
                   line = next(lines)
                   if "---" in line: break
                   v0 += line.split()
              v0 = [float(xx) for xx in v0]
           # Read v1
           if "v1" in variables:
              v1  = [] 
              while True:
                   line = next(lines)
                   if "---" in line: break
                   v1 += line.split()
              v1 = [float(xx) for xx in v1]
           # Generate Structure
           drst[meplabel] = (mepsvalue,Vmep,xms,gms,Fms,v0,v1,t_value)
           # limit values
           if mepsvalue < mins: mins = mepsvalue
           if mepsvalue > maxs: maxs = mepsvalue
    return tpath, tcommon, drst
#-----------------------------------------------#
def write_rst_head(rstfile,tpath,tcommon):
    NEPL = 5
    str_rst  = ""

    path,mu,ds,hsteps,cubic  = tpath
    ch,mtp,atonums,masses,mu = tcommon

    # Write structinfo-block
    str_rst += "start_structinfo\n"
    str_rst += "   ch         %i\n"%ch
    str_rst += "   mtp        %i\n"%mtp
    str_rst += "   atonums\n"
    for idx in range(0,len(atonums),2*NEPL):
        str_rst += " "*5 + "  ".join( "%2i"%atonum for atonum in atonums[idx:idx+2*NEPL] ) + "\n"
    str_rst += "   masses\n"
    for idx in range(0,len(atonums),NEPL):
        str_rst += " "*5 + "  ".join( "%14.8E"%(mass*AMU) for mass in masses[idx:idx+NEPL] ) + "\n"
    str_rst += "end_structinfo\n"
    str_rst += "\n"

    # Write pathinfo-block
    str_rst += "start_pathinfo\n"
    str_rst += "   mtype      %s\n"%path
    str_rst += "   mu         %f\n"%(mu*AMU)
    str_rst += "   ds         %7.6f\n"%ds
    if hsteps in [None,"no"]: str_rst += "   hsteps     no\n"
    else                    : str_rst += "   hsteps     %i\n"%hsteps
    if cubic in ["no",False,None]: str_rst += "   cubic      no\n"
    else                         : str_rst += "   cubic      %7.6f\n"%cubic
    str_rst += "end_pathinfo\n"
    str_rst += "\n"

    # Head for the structures
    str_rst += "STRUCTURES\n"
    str_rst += "\n"
    str_rst += "# For each structure: mass-scaled cartesian coordinates and atomic units.\n"
    str_rst += "# Data in ==>: label, written info, s value, energy, t value (only if path is pm).\n"
    str_rst += "\n"

    write_file(rstfile,str_rst)
#-----------------------------------------------#
def rst_string(label,data_structure):
    NEPL = 5
    # expand data
    s_i, E_i, xms_i,gms_i,Fms_i,v0_i,v1_i,t_i = data_structure

    # generate string for molecule
    info = "x"
    if (gms_i is not None) and (len(gms_i) != 0): info = info+"-g"
    if (Fms_i is not None) and (len(Fms_i) != 0): info = info+"-F"
    if (v0_i  is not None) and (len(v0_i ) != 0): info = info+"-v0"
    if (v1_i  is not None) and (len(v1_i ) != 0): info = info+"-v1"

    # first line
    if t_i is None: string = " ==> %7s   %11s   %+11.6f   %.10f      \n"%(label, info, s_i, E_i)
    else:           string = " ==> %7s   %11s   %+11.6f   %.10f  %.5E\n"%(label, info, s_i, E_i, t_i)

    # Data
    written_variables = info.split("-")
    for variable in written_variables:
        if variable == "x":  var2write = xms_i
        if variable == "g":  var2write = gms_i
        if variable == "F":  var2write = fncs.matrix2lowt(Fms_i)
        if variable == "v0": var2write = v0_i
        if variable == "v1": var2write = v1_i
        # Write list 
        for idx in range(0,len(var2write),NEPL):
            line = "  ".join(["%15.7E"%value for value in var2write[idx:idx+NEPL]])
            string += "  " + line + "\n"
        string += "    ---\n"
    return string
#-----------------------------------------------#
def write_rst_add(rstfile,label,data_structure):
    # get string for structure
    string = rst_string(label,data_structure)
    # Write file using fcntl (in case of multiple processes)
    with open(rstfile, "a") as g:
         fcntl.flock(g, fcntl.LOCK_EX)
         g.write(string)
         fcntl.flock(g, fcntl.LOCK_UN)
#-----------------------------------------------#
def write_rst(rstfile,tpath,tcommon,drst):
    # random name (no exists)
    for i in range(1000):
        provisional = "abcdefghijklmnopqrstuvwxyz"
        letters = [lett for lett in "abcdefghijklmnopqrstuvwxyz"]
        random.shuffle(letters)
        part1 = "".join(letters[:5])
        numbers = [num  for num  in "0123456789"]
        random.shuffle(numbers)
        part2 = "".join(numbers[:5])
        provisional = rstfile+"."+part1+"-"+part2
        if not os.path.exists(provisional): break
    if os.path.exists(provisional): provisional = rstfile
    # write head
    write_rst_head(provisional,tpath,tcommon)
    # string with geoms
    string = ""
    thelist = sorted([ (data[0],label) for label,data in drst.items()])
    for s_i,label in thelist:
        data_structure = drst[label]
        string += rst_string(label,data_structure)
    # write file
    write_file(provisional,string,mode="a")
    # rename to rstfile
    if provisional != rstfile: os.rename(provisional,rstfile)
#-----------------------------------------------#
def rst2xyz(rst,xyz,onlyhess=True,Eref=None):
    tpath, tcommon, drst = read_rst(rst)
    thelist = [(data[0],label) for (label,data) in drst.items()]
    thelist.sort()
    (ch,mtp,atonums,masses,mu) = tcommon
    symbols = fncs.atonums2symbols(atonums)
    natoms  = len(symbols)
    string = ""
    for s_i,label in thelist:
        s_i,E_i,xms_i,gms_i,Fms_i,v0_i,v1_i,t_i = drst[label]
        if onlyhess and Fms_i is None: continue
        if Eref is not None:
            E_i = (E_i - Eref) * KCALMOL
            energy = "%+.4f kcal/mol"%E_i
        else:
            energy = "%+.6f hartree"%E_i
        string += " %i\n"%natoms
        string += " * E = %s ; s = %7.4f (%s)\n"%(energy,s_i,label)
        xcc_i = fncs.ms2cc_x(xms_i,masses,mu)
        for idx,symbol in enumerate(symbols):
            x,y,z = fncs.xyz(xcc_i,idx)
            x *= ANGSTROM
            y *= ANGSTROM
            z *= ANGSTROM
            string += " %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,x,y,z)
    write_file(xyz,string)
#===============================================#

#===============================================#
# Function(s) to deal with: GAUSSIAN file(s)    #
#===============================================#
def read_gauout(gauout):
    '''
    Read Gaussian output (data in final message)
    and return important data
    '''
    # check extension
    end_out = (gauout.lower()).endswith(".out")
    end_log = (gauout.lower()).endswith(".log")
    if (not end_out) and (not end_log):
        raise Exc.FileType(Exception)
    # read lines
    lines = read_file(gauout)
    # CHECK FILE IS GAUSSIAN OUT
    correct = False
    for line in lines:
        if "Entering Gaussian System" in line:
           correct = True
           break
    if not correct: raise Exc.FileType(Exception)
    # ONIOM energy?
    try:
       E_ONIOM = None
       for line in lines[::-1]:
           if "ONIOM: extrapolated energy" in line:
               E_ONIOM = float(line.split()[-1])
               break
    except: E_ONIOM = None
    # Get Forces if exists
    key1 = "Forces (Hartrees/Bohr)"
    key2 = "Cartesian Forces:"
    try:
       gcc = []
       for line in fncs.extract_string(lines,key1,key2).split("\n")[3:-3]:
           dummy, dummy, gx, gy, gz = line.split()
           gcc += [float(gx),float(gy),float(gz)]
       # forces to gradient
       gcc = [-g_i for g_i in gcc]
    except: gcc = None
    if gcc == []: gcc = None
    # Read string from gaussian output
    DATA_LAST = []
    DATA_WFCC = []
    key1, key2 = "\GINC-","@"
    strings = fncs.extract_string(lines,key1,key2,accumulate=True)
    for string in strings:
        lines = "".join([line.strip() for line in string.split()])
        lines = lines.split("\\\\")
        if lines == [""]: return None, None, None, None, None, None

        # lines[3]: ch, mtp, symbols, xcc
        str_geom = lines[3].split("\\")
        ch , mtp = str_geom[0].split(",")
        xcc      = []
        symbols  = []
        ch  = int(ch)
        mtp = int(mtp)
        E   = None
        masses = []
        for line in str_geom[1:]:
            coords = line.split(",")
            if   len(coords) == 4:
               symbol, x, y, z = line.split(",")
            else:
               symbol, tmp, x, y, z = line.split(",")
            xcc += [float(x),float(y),float(z)]
            symbols.append(symbol)
        xcc = [xi/ANGSTROM for xi in xcc] # in bohr
        # lines[4]: energy
        E = float(lines[4].split("HF=")[1].split("\\")[0])
        # Gradient
        if gcc is None or len(gcc) == 0: gcc = [0.0 for x in xcc]
        # Hessian
        Fcc = []
        for idx in range(len(lines)):
            line = lines[idx]
            if "NImag" not in line: continue
            Fcc = [float(fij) for fij in lines[idx+1].split(",")]
        # other lists
        atonums = fncs.get_atonums(symbols)
        masses  = fncs.atonums2masses(atonums)
        calclevel = ""
        # Does it have hessian?
        if Fcc != []: DATA_WFCC = [xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel]
        # Save data
        DATA_LAST = [xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel]
    # Use ONIOM extrapolated energy if found!
    if E_ONIOM is not None:
       if len(DATA_WFCC) != 0: DATA_WFCC[4] = E_ONIOM
       if len(DATA_LAST) != 0: DATA_LAST[4] = E_ONIOM
    # return
    if DATA_WFCC != []: return DATA_WFCC
    else              : return DATA_LAST
#-----------------------------------------------#
def read_fchk(fchkfile):
    #------------------------#
    # The labels to look for #
    #------------------------#
    labels     = {}                              ; found_dict    = {}
    labels[0]  = "Number of atoms"               ; found_dict[0] = False
    labels[1]  = "Charge"                        ; found_dict[1] = False
    labels[2]  = "Multiplicity"                  ; found_dict[2] = False
    labels[3]  = "Total Energy"                  ; found_dict[3] = False
    labels[4]  = "Atomic numbers"                ; found_dict[4] = False
    labels[5]  = "Current cartesian coordinates" ; found_dict[5] = False
    labels[6]  = "Real atomic weights"           ; found_dict[6] = False
    labels[7]  = "Cartesian Gradient"            ; found_dict[7] = False
    labels[8]  = "Cartesian Force Constants"     ; found_dict[8] = False
    # mandatory labels in fchk
    idx_basic_labels = [0,1,2,3,4,5,6]
    # initialization
    natoms   = None ; ch     = None ; mtp  = None  ; E = None
    atonums  = []   ; masses = []
    xcc      = []   ; gcc    = []   ; Fcc  = []
    # read file
    lines = read_file(fchkfile)
    # check extension
    end_fchk = (fchkfile.lower()).endswith(".fchk")
    if (not end_fchk): raise Exc.FileType(Exception)
    # CHECK FILE IS FCHK
    correct = False
    for line in lines:
        if line.startswith("Current cartesian coordinates  "):
           correct = True
           break
    if not correct: raise Exc.FileType(Exception)
    # Get level of calculation
    calclevel = " ".join(lines[1].split()[1:])
    # get data from lines
    for idx in range(len(lines)):
        line = lines[idx]
        # Number of atoms
        if line.startswith(labels[0]):
            found_dict[0] = True
            natoms = int(line.split()[-1])
        # Charge
        elif line.startswith(labels[1]):
            found_dict[1] = True
            ch = int(line.split()[-1])
        # Spin multiplicity
        elif line.startswith(labels[2]):
            found_dict[2] = True
            mtp = int(line.split()[-1])
        # Total Energy
        elif line.startswith(labels[3]):
            found_dict[3] = True
            E = float(line.split()[-1])
        # Atomic Numbers
        elif line.startswith(labels[4]):
            found_dict[4] = True
            length = int(line.split()[-1])
            idx2   = idx+1
            while len(atonums) != length:
                  nextline = lines[idx2]
                  atonums += [int(i) for i in nextline.split()]
                  idx2    += 1
        # Cartesian Coordinates
        elif line.startswith(labels[5]):
            found_dict[5] = True
            length = int(line.split()[-1])
            idx2   = idx+1
            while len(xcc) != length:
                  nextline = lines[idx2]
                  xcc     += [float(i) for i in nextline.split()]
                  idx2    += 1
        # List of atomic masses
        elif line.startswith(labels[6]):
            found_dict[6] = True
            length = int(line.split()[-1])
            idx2   = idx+1
            while len(masses) != length:
                  nextline = lines[idx2]
                  masses  += [float(i) for i in nextline.split()]
                  idx2    += 1
        # Cartesian Gradient
        elif line.startswith(labels[7]) and natoms != 1:
            found_dict[7] = True
            length = int(line.split()[-1])
            idx2   = idx+1
            while len(gcc) != length:
                  nextline = lines[idx2]
                  gcc     += [float(i) for i in nextline.split()]
                  idx2    += 1
        # Cartesian Force Constant Matrix
        elif line.startswith(labels[8]) and natoms != 1:
            found_dict[8] = True
            length = int(line.split()[-1])
            idx2   = idx+1
            while len(Fcc) != length:
                  nextline = lines[idx2]
                  Fcc     += [float(fij) for fij in nextline.split()]
                  idx2    += 1
    # Return data
    #if gcc == []: gcc = None
    #if Fcc == []: Fcc = None
    return xcc,atonums,ch,mtp,E,gcc,Fcc,masses, calclevel.strip()
#===============================================#


#===============================================#
# Function(s) to deal with: ORCA file(s)        #
#===============================================#
def read_orca(orca_out):
    xcc, atonums, ch, mtp, E, calclevel = read_orcaout(orca_out)
    # Looking for files with gcc and Fcc
    orca_engrad = ".".join(orca_out.split(".")[:-1]) + ".engrad"
    orca_hess   = ".".join(orca_out.split(".")[:-1]) + ".hess"
    # Read them
    gcc, dummy = read_orcaengrad(orca_engrad)
    Fcc        = read_orcahess(orca_hess)
    # Get masses
    masses  = fncs.atonums2masses(atonums)
    # return data
    return xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel
#-----------------------------------------------#
def read_orcaout(orca_out):
    # Read lines
    lines = read_file(orca_out)
    # CHECK FILE IS ORCA OUT
    correct = False
    for line in lines:
        if "This ORCA versions uses" in line:
           correct = True
           break
    if not correct: raise Exc.FileType(Exception)
    # Initialize
    ch, mtp, E = None, None, None
    xcc, symbols, atonums = [], [], []
    # basis set
    basisset = ""
    try:
       for line in lines:
        if "utilizes the basis" in line:
            basisset = line.split(":")[1]
            break
    except: pass
    # hamiltonian
    hamiltonian = ""
    try:
        for idx,line in enumerate(lines):
            if "Hamiltonian:" in line:
                hamiltonian = lines[idx+1].split()[-1]
                break
    except: pass
    # level of calculation
    calclevel = hamiltonian+" "+basisset
    # Get ch, mtp, E
    pos = None
    for idx in range(len(lines)):
        line = lines[idx]
        if "Total Charge           Charge" in line: ch  = int(line.split()[-1])
        if "Multiplicity           Mult"   in line: mtp = int(line.split()[-1])
        if "FINAL SINGLE POINT ENERGY "    in line: E   = float(line.split()[-1])
        if "CARTESIAN COORDINATES (A.U.)"  in line: pos = idx
    # Get xcc, symbols, atonums
    if pos is None: sys.exit("Unable to find 'CARTESIAN COORDINATES (A.U.)' label in file!")
    pos = pos + 3
    for line in lines[pos:]:
        line = line.strip()
        if line == "": break
        idx, symbol, atonum, dummy, mass, x, y, z = line.split()
        xcc = xcc + [float(x),float(y),float(z)]
        symbols.append(symbol)
        atonums.append(int(atonum.split(".")[0]))
    # Return
    return xcc, atonums, ch, mtp, E, calclevel.strip()
#-----------------------------------------------#
def read_orcaengrad(orca_engrad):
    # Initialize
    natoms   = None
    atonums  = []   ; Etot = None
    x_cc     = []   ; g_cc = []
    # Read lines
    if not os.path.exists(orca_engrad): return g_cc, (x_cc, atonums, Etot)
    lines  = read_file(orca_engrad)
    # Read natoms, Etot, g_cc, atonums, xcc
    natoms = int(lines[3])
    Etot   = float(lines[7])
    g_cc   = [float(line)  for line in lines[11:11+3*natoms]]
    Zxyz   = [line.split() for line in lines[-natoms:]]
    for Z,x,y,z in Zxyz:
       x_cc = x_cc + [float(x),float(y),float(z)]
       atonums = atonums + [int(Z)]
    # Return data
    return g_cc, (x_cc, atonums, Etot)
#-----------------------------------------------#
def read_orcahess(orca_hess):
    if not os.path.exists(orca_hess): return []
    # Read lines
    lines  = read_file(orca_hess)
    # Read hessian
    for idx in range(len(lines)):
        line = lines[idx]
        if "$hessian" in line:
           nrows = int(lines[idx+1])
           idx_hessian = idx
    F_cc = [ [0.0]*nrows for i in range(nrows)]
    fcol = int(lines[idx_hessian+2].split()[ 0]) # first column
    lcol = int(lines[idx_hessian+2].split()[-1]) # last  column
    row  = 0
    for line in lines[idx_hessian+3:]:
        if row == nrows-1:
           fcol = int(line.split()[ 0]) # first column
           lcol = int(line.split()[-1]) # last  column
           row  = 0
           continue
        # Split data in line
        line_data = line.split()
        # Save data
        row = int(line_data[0])
        F_cc[row][fcol:lcol+1] = [float(value) for value in line_data[1:]]
        # Finished?
        if row == nrows -1 and lcol == nrows-1: break
    # Promediate hessian
    for i in range(1,nrows,1):
        for j in range(i):
           Fij = F_cc[i][j]
           Fji = F_cc[j][i]
           average = 0.5*(Fij + Fji)
           F_cc[i][j] = average
           F_cc[j][i] = average
    # Get lower triangle
    F_lt = []
    for i in range(nrows):
        for j in range(0,i+1):
            F_lt.append(F_cc[i][j])
    F_cc = F_lt
    # Return data
    return F_cc
#===============================================#

#===============================================#
# Function(s) to deal with: Q2DTor output file  #
#===============================================#
def read_q2dtorout(filename):
    # read file
    if type(filename) != type([]): lines = read_file(filename)
    else                         : lines = filename
    # some keywords
    key1 = "Energy of the lowest zero point level of the torsional PES:"
    key2 = "(a) Rovibrational partition functions"
    key3 = "Components of E2DT"
    key4 = "Translational and electronic partition functions"
    key5 = "Total partition functions, from rovibrational partition functions in (b)"
    bool1 = False
    bool2 = False
    bool3 = False
    # initialize variables
    ltemp, zpe_MSHO, Q_MSHO, zpe_E2DT, Q_E2DT = None, None, None, None, None
    Q_tr, Q_el = None, None
    # read file
    for idx,line in enumerate(lines):
        # get ZPE
        if bool1 and "* MSHO =>" in line:
            value = line.split()[3]
            units = line.split()[4]
            if "kcal" in units.lower(): zpe_MSHO = float(value) / KCALMOL
            else                      : zpe_MSHO = float(value) / KJMOL
        if bool1 and "* E2DT =>" in line:
            value = line.split()[3]
            units = line.split()[4]
            if "kcal" in units.lower(): zpe_E2DT = float(value) / KCALMOL
            else                      : zpe_E2DT = float(value) / KJMOL
            bool1 = False
        if key1 in line:
           bool1    = True
           zpe_MSHO = None
           zpe_E2DT = None
        # get rv partition function
        if key3 in line: bool2 = False
        if bool2:
            if "|" not in line: continue
            if "rv(" in line: continue
            T, q1WHO_i, qMSHO_i, qE2DT_i, ratio_i = line.split("|")
            ltemp.append(float(T))
            Q_MSHO.append(float(qMSHO_i))
            Q_E2DT.append(float(qE2DT_i))
        if key2 in line:
           bool2  = True
           ltemp  = []
           Q_MSHO = []
           Q_E2DT = []
        # get translational and electronic partition function
        if key5 in line: bool3 = False
        if bool3:
            if "|" not in line: continue
            if "Qtra" in line: continue
            T, qtr, qel = line.split("|")
            Q_tr.append(float(qtr))
            Q_el.append(float(qel))
        if key4 in line:
           bool3  = True
           Q_tr   = []
           Q_el   = []
    # Qtr to au
    if Q_tr is not None:
       p0     = 1E5 * (METER**3 / JOULE)
       Q_tr   = [Q_tr[idx]/(KB*T/p0) for idx,T in enumerate(ltemp)]
    # rovib to total
    if Q_MSHO is not None:
       Q_MSHO = [Q_tr[idx]*Q_MSHO[idx]*Q_el[idx] for idx in range(len(ltemp))]
    if Q_E2DT is not None:
       Q_E2DT = [Q_tr[idx]*Q_E2DT[idx]*Q_el[idx] for idx in range(len(ltemp))]
    return ltemp, zpe_MSHO, Q_MSHO, zpe_E2DT, Q_E2DT
#===============================================#


#===============================================#
# Function(s) to deal with: MSTor output file   #
#===============================================#
def read_mstorout(filename):
    # read file
    if type(filename) != type([]): lines = read_file(filename)
    else                         : lines = filename
    # some keywords
    key1 = "ZPE of the system = "
    key2 = "T(K)    Elec.        Trans."
    key3 = "Total Partition Function"
    key4 = "Partition Functions Including All Structures"
    key5 = "Zero of energy is set to ZPE"
    key6 = "============="
    # initialize variables
    ltemp, zpe_MSHO, Q_MSHO, zpe_MSTOR, Q_MSTOR = None, None, None, None, None
    Q_tr, Q_el = None, None
    bool1 = False
    bool2 = False
    bool3 = False
    # read file
    for idx,line in enumerate(lines):
        # get ZPE
        if key1 in line:
           zpe_MSHO  = float(line.split()[-1]) / KCALMOL
           zpe_MSTOR = float(zpe_MSHO)
        # get translational and electronic partition function
        if key3 in line: bool1 = False
        if bool1:
            data = line.split()
            if len(data) != 5: continue
            T, qel, qtr = data[0:3]
            ltemp.append(float(T))
            Q_tr.append(float(qtr))
            Q_el.append(float(qel))
        if key2 in line:
           bool1  = True
           ltemp  = []
           Q_tr   = []
           Q_el   = []
        # get rv partition function
        if key6 in line: bool3 = False
        if bool3:
            data = line.split()
            if len(data) != 5: continue
            T, qMSHO_i, qMSTOR_i = data[0:3]
            Q_MSHO.append(float(qMSHO_i))
            Q_MSTOR.append(float(qMSTOR_i))
        if bool2 and key5 in line: bool3 = True
        if key4 in line:
           bool2   = True
           Q_MSHO  = []
           Q_MSTOR = []
    # Qtr to au
    if Q_tr is not None:
       p0     = 1E5 * (METER**3 / JOULE)
       Q_tr   = [Q_tr[idx]/(KB*T/p0) for idx,T in enumerate(ltemp)]
    # rovib to total
    if Q_MSHO is not None:
       Q_MSHO = [Q_tr[idx]*Q_MSHO[idx]*Q_el[idx] for idx in range(len(ltemp))]
    if Q_MSTOR is not None:
       Q_MSTOR = [Q_tr[idx]*Q_MSTOR[idx]*Q_el[idx] for idx in range(len(ltemp))]
    return ltemp, zpe_MSHO, Q_MSHO, zpe_MSTOR, Q_MSTOR
#===============================================#


#===============================================#
# Function(s) to deal with: PDB file(s)         #
#===============================================#
def read_geoms_pdb(pdbfile):
    lines  = read_file(pdbfile)
    geoms  = {}
    record = False
    for line in lines:
        if "ENDMDL" in line:
           record = False
           geoms[number] = xvec
        if record and "ATOM " in line:
           x, y, z, dd, dd, symbol = line.split()[5:]
           xvec += [float(x),float(y),float(z)]
           symbols.append(symbol)
        if "MODEL " in line:
           record  = True
           number  = int(line.split()[1])
           xvec    = []
           symbols = []
    return geoms, symbols
#===============================================#

