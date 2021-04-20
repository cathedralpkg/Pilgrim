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
| Module     :  common             |
| Sub-module :  orca               |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Interface for the Electronic Structure calculation
using ORCA
'''

#=============================================#
#           TO BE MODIFIED BY USER            #
#=============================================#
#EXE="/home/david/Software/orca_4_0_1_2/orca" #
#=============================================#

#=============================================#
import os
import numpy as np
#---------------------------------------------#
import common.Exceptions as     Exc
from   common.fncs       import xyz
from   common.fncs       import clean_lines
from   common.fncs       import atonums2masses
from   common.files      import read_file
#=============================================#


#==========================================================#
def set_EXE():
    global EXE
    txt = os.path.dirname(os.path.realpath(__file__))+"/paths.txt"
    # Defined in this file
    if 'EXE' in globals():
        return
    # Try to export it from bashrc
    elif "OrcaExe" in os.environ:
        # in .bashrc: export OrcaExe="$MYHOME/Software/orca_4_0_1_2/orca"
        EXE = os.environ["OrcaExe"]
        return
    # Export it from file
    elif os.path.exists(txt):
        lines = read_file(txt)
        lines = clean_lines(lines,"#",True)
        for line in lines:
            if line == "\n": continue
            name, path = line.split()
            if name == "orca":
                EXE = path
                return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    if not os.path.exists(EXE): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def execute(ifile,ofile,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       ifile = folder + ifile
       ofile = folder + ofile
       err   = folder + err
    # Execution file?
    set_EXE()
    check_EXE()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = ofile
    # Execute command
    command = "%s %s 1>%s 2>%s"%(EXE,ifile,ofile,err)
    try   : status  = os.system(command)
    except: raise exception
    return status
#----------------------------------------------------------#
def normal_termination(ofile):
    with open(ofile,"r") as asdf: olines = asdf.readlines()
    if len(olines) == 0: return False
    if "ORCA TERMINATED NORMALLY" not in olines[-2]: return False
    return True
#----------------------------------------------------------#
def iofiles(name,folder=None):
    '''
    For a given name, it returns the name of all files
    for the Gaussian calculation
    '''

    if   folder is None          : folder = ""
    elif not folder.endswith("/"): folder = folder+"/"
    else                         : folder = folder
    wname  = folder + name # whole name
    ifile  = wname + ".inp"
    ofile  = wname + ".out"
    engrad = wname + ".engrad"
    hess   = wname + ".hess"
    gbw    = wname + ".gbw"
    err    = wname + ".err"
    return wname, ifile, ofile, engrad, hess, gbw, err
#==========================================================#


#==========================================================#
# Function(s) to deal with: ORCA file(s)                   #
#==========================================================#
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
    pos_cart = None
    pos_zmat = None
    for idx in range(len(lines)):
        line = lines[idx]
        if "Total Charge           Charge" in line: ch  = int(line.split()[-1])
        if "Multiplicity           Mult"   in line: mtp = int(line.split()[-1])
        if "FINAL SINGLE POINT ENERGY "    in line: E   = float(line.split()[-1])
        if "CARTESIAN COORDINATES (A.U.)"  in line: pos_cart = idx
        if "INTERNAL COORDINATES (ANGST"   in line: pos_zmat = idx
    # Get xcc, symbols, atonums
    if pos_cart is None: sys.exit("Unable to find 'CARTESIAN COORDINATES (A.U.)' label in file!")
    pos_cart += 3
    for line in lines[pos_cart:]:
        line = line.strip()
        if line == "": break
        idx, symbol, atonum, dummy, mass, x, y, z = line.split()
        if symbol == "-": symbol = "XX"
        xcc += [float(x),float(y),float(z)]
        symbols.append(symbol)
        atonums.append(int(atonum.split(".")[0]))
    # data in internal coordinates
    lzmat = None
    if pos_zmat is not None:
       pos_zmat += 2
       lzmat = []
       for idx,line in enumerate(lines[pos_zmat:]):
           line = line.strip()
           if line == "": break
           symbol, atk, atj, ati, dist, angle, dihedral = line.split()
           if symbol == "-": symbol = "XX"
           atk, atj, ati = int(atk)-1, int(atj)-1, int(ati)-1
           dist, angle, dihedral = float(dist), float(angle), float(dihedral)
           if   idx == 0: connections, values = ()            , ()
           elif idx == 1: connections, values = (atk,)        , (dist,)
           elif idx == 2: connections, values = (atk,atj)     , (dist,angle)
           else         : connections, values = (atk,atj,ati) , (dist,angle,dihedral)
           data = (symbol,connections,values)
           lzmat += [data]
    # Return
    return xcc, lzmat, atonums, ch, mtp, E, calclevel.strip()
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
#-----------------------------------------------#
def gen_zmatrix_string(lzmat,zmatvals):
    string = ""
    for idx,zmatline in enumerate(lzmat):
        symbol, connections, keys = zmatline
        # Dummy atoms
        if symbol.upper() in "XX,X": symbol = "DA"
        values = [zmatvals.get(key,key) for key in keys]
        # complete tuples
        connections = [at+1 for at in connections]
        while len(connections) != 3: connections = tuple(list(connections)+[0])
        while len(values     ) != 3: values      = tuple(list(values     )+[0.0])
        # generate line
        string += "%2s  "%symbol
        for at in connections:
            try   : string += "%3i  "%at
            except: string += "%3s  "%at
        for value in values:
            try   : string += "%13.7f  "%value
            except: string += "%13s  "%value
        string += "\n"
    return string
#==========================================================#





#==========================================================#
def read_orca(orca_out):
    xcc, lzmat, atonums, ch, mtp, E, calclevel = read_orcaout(orca_out)
    # Looking for files with gcc and Fcc
    orca_engrad = ".".join(orca_out.split(".")[:-1]) + ".engrad"
    orca_hess   = ".".join(orca_out.split(".")[:-1]) + ".hess"
    # Read them
    gcc, dummy = read_orcaengrad(orca_engrad)
    Fcc        = read_orcahess(orca_hess)
    # Get masses
    masses     = atonums2masses(atonums)
    # return data
    return xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel
#==========================================================#
