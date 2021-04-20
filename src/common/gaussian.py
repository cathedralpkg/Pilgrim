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
| Sub-module :  gaussian           |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Interface for the Electronic Structure calculation
using GAUSSIAN
'''

#=============================================#
import sys
import os
import time
import subprocess
#---------------------------------------------#
import common.Exceptions as Exc
from   common.fncs       import xyz
from   common.fncs       import clean_lines
from   common.fncs       import flatten_llist
from   common.fncs       import atonums2masses
from   common.fncs       import symbol_and_atonum
from   common.fncs       import symbols_and_atonums
from   common.fncs       import extract_string
from   common.fncs       import get_atonums
from   common.physcons   import ANGSTROM
from   common.files      import read_file
from   common.pgs        import get_pgs
#=============================================#


KEYGAU  = "GauExe"
KEYFCHK = "GauFchk"
#==========================================================#
#                  TO BE MODIFIED BY USER                  #
#==========================================================#
#EXE  = "/home/programs/G09_64david/g09/g09"               #
#FCHK = "/home/programs/G09_64david/g09/formchk"           #
#----------------------------------------------------------#
# in .bashrc:                                              #
#  export GauExe="/home/programs/G09_64david/g09/g09"      #
#  export GauFchk="/home/programs/G09_64david/g09/formchk" #
#==========================================================#



#=======================================================#
# FUNCTION FOR READING THE FCHK FILE OF GAUSSIAN PROG.  #
#=======================================================#
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
#=======================================================#


#=======================================================#
# FUNCTIONS FOR READING THE LOG FILE OF GAUSSIAN PROG.  #
#=======================================================#
def split_gaulog_into_gaublocks(filename):
    # Key to split the system into blocks
    str_end='Normal termination'
    # Read log/out file and join lines as string
    with open(filename,'r') as asdf: lines = asdf.readlines()
    text  = "".join(lines)
    # Divide by str_end (last element has to be excluded)
    blocks = text.split(str_end)[:-1]
    # For some reason, sometimes Gaussian prints a block
    # without information. In these cases, the block consists
    # of a few lines. Here, we exclude that cases
    #print [len(block.split("\n")) for block in blocks]
    blocks = [block for block in blocks if len(block.split("\n")) > 300]
    # Remove the lines list and the whole text
    del lines, text
    # No normal termination?
    return blocks
#-------------------------------------------------------#
def get_data_from_maintext(mtext):
    key_geom1 = 'Z-Matrix orientation'
    key_geom2 = 'Input orientation'
    key_geom3 = 'Standard orientation'
    key_force = "     Forces ("
    key_zmat  = 'Final structure in terms of initial Z-matrix:'
    key_end   = "------------------"
    key_1bar  = "\\"
    key_oniom = "ONIOM: extrapolated energy"

    # (a) Find cartesian coordinates
    if   key_geom1 in mtext: geom = mtext.split(key_geom1)[-1]
    elif key_geom2 in mtext: geom = mtext.split(key_geom2)[-1]
    elif key_geom3 in mtext: geom = mtext.split(key_geom3)[-1]
    else                   : geom,xcc = None, None
    if geom is not None:
       # convert to list of lines and get the lines associated to geometry
       geom = "\n".join(geom.split("\n")[5:])
       idx  = geom.find(key_end)
       geom = geom[:idx].strip()
       # convert to list of floats
       geom = [line.split() for line in geom.split("\n") if line.strip() != ""]
       xcc  = [[float(x),float(y),float(z)] for (_,atnum,_,x,y,z) in geom]
       xcc  = flatten_llist(xcc)
    # (b) Find forces --> gradient
    if key_force in mtext:
       force = mtext.split(key_force)[-1]
       # convert to list of lines and get the lines associated to forces
       force = "\n".join(force.split("\n")[3:])
       idx = force.find(key_end)
       force = force[:idx]
       # convert to list of floats
       force = [line.split() for line in force.split("\n") if line.strip() != ""]
       gcc  = [[-float(gx),-float(gy),-float(gz)] for (_,atnum,gx,gy,gz) in force]
       gcc  = flatten_llist(gcc)
    else: gcc = None
    # (c) Find z-matrix
    if  key_zmat in mtext:
        lines_zmat = mtext.split(key_zmat)[-1].strip().split("\n")
        for idx,line in enumerate(lines_zmat):
            line = line.strip()
            if line == "" or key_1bar in line:
               idx1 = idx
               break
        zmat = [line.strip() for line in lines_zmat[:idx1] if "Variables:" not in line]
    else: zmat = None
    # (d) ONIOM energy?
    E_ONIOM = None
    for line in mtext.split("\n")[::-1]:
        if key_oniom in line:
           E_ONIOM = float(line.split()[-1])
           break
    # Convert xcc to bohr
    if xcc is not None: xcc = [xi/ANGSTROM for xi in xcc]
    # Return data
    return xcc, gcc, zmat, E_ONIOM
#-------------------------------------------------------#
def get_data_from_archive(summary):
    # Keywords to look for
    key_hag  = "#"
    key_1bar ='\\'
    key_2bar ='\\\\'
    key_ver  ='Version'
    key_en1  ='State='
    key_en2  ='RMSD='
    key_imag ='NImag='
    # logical variables
    Lzmat = False
    Lxyz  = False
    Lhess = False
    # (a) the command line
    idx1 = summary.find(key_hag,0)+len(key_hag)
    idx2 = summary.find(key_2bar,idx1)
    commands = summary[idx1:idx2]
    # (b) the comment line
    idx3 = summary.find(key_2bar,idx2)+len(key_2bar)
    idx4 = summary.find(key_2bar,idx3)
    comment = summary[idx3:idx4]
    # (c) charge and multiplicity
    idx5 = summary.find(key_2bar,idx4)+len(key_2bar)
    idx6 = summary.find(key_1bar,idx5)
    ch,mtp = [int(value) for value in summary[idx5:idx6].split(",")]
    # (d) z-matrix or Cartesian coordinates
    idx7 = summary.find(key_ver,idx6)
    geom = [string for string in summary[idx6+len(key_1bar):idx7].split(key_1bar) if string != ""]
    if len(geom[0]) <= 4:
       Lzmat   = True
       zmat    = list(geom)
       symbols = [line.split(",")[0] for line in zmat if "=" not in line]
       xcc     = None
    else:
       Lxyz    = True
       zmat    = None
       symbols = [line.split(",")[0]   for line in geom]
       # sometimes, this line has 5 elements instead of 4
       # for this reason, coordinates are extracted with [-3:]
       # instead of [1:]
       xyz = [line.split(",")[-3:] for line in geom]
       xyz = [[float(x),float(y),float(z)] for (x,y,z) in xyz]
       xcc = flatten_llist(xyz)
    # (e) Energy and other info
    idx8a = summary.find(key_ver,idx7)
    idx8b = summary.find(key_en1,idx7)
    idx8  = max(idx8a,idx8b)
    idx9  = summary.find(key_1bar,idx8)+len(key_1bar)
    idx10 = summary.find(key_en2,idx9)
    str_energies = summary[idx9:idx10].replace(key_1bar," ")
    energies = str_energies.split()
    energies = [line.split("=") for line in energies]
    # remove S**2 (for open-shell)
    energies = [(float(energy),level.strip()) for level,energy in energies if not level.strip().startswith("S2")]
    # (f) Hessian matrix
    Lhess = key_imag in summary
    if Lhess:
       idx11=summary.find(key_imag,0)+len(key_imag)
       idx12=summary.find(key_2bar,idx11)
       num_imag = int(summary[idx11:idx12])
       idx12 += len(key_2bar)
       idx13=summary.find(key_2bar,idx12)
       # low-triangle hessian
       Fcc = [float(value) for value in summary[idx12:idx13].split(",")]
    else:
       num_imag = -1
       Fcc  = None
    # (g) Gradient may appera after hessian matrix
    if Lhess:
       idx13 += len(key_2bar)
       idx14=summary.find(key_2bar,idx13)
       gcc = [float(value) for value in summary[idx13:idx14].split(",")]
    else:
       gcc = None
    # Convert xcc to bohr
    if xcc is not None: xcc = [xi/ANGSTROM for xi in xcc]
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,num_imag,zmat
#-------------------------------------------------------#
def get_data_from_gaublock(gaublock):
    '''
    gaublock is a string 
    gaublock contains the info from begining till Normal termination
    '''
    # Divide the block into the summary part and the rest
    key_start='GINC'
    key_end  ='@'
    mtext    = gaublock.split(key_start)[0]
    summary  = gaublock.split(key_start)[1].split(key_end)[0]
    # Remove the initial blank space of each line in summary
    # Also remove the line breaks
    summary  = "".join([line.strip() for line in summary.split("\n")])
    # Data in summary (aka archive)
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,num_imag,zmat = get_data_from_archive(summary)
    # Data in main text (excluding archive)
    xcc_mt, gcc_mt, zmat_mt, E_oniom = get_data_from_maintext(mtext)
    # If data not in archive --> get from main text
    if xcc  is None and  xcc_mt is not None: xcc  =  xcc_mt
    if gcc  is None and  gcc_mt is not None: gcc  =  gcc_mt
    if zmat is None and zmat_mt is not None: zmat = zmat_mt
    # Return data
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,E_oniom,num_imag,zmat
#-------------------------------------------------------#
def read_gaussian_log(filename,target_level=None):
    if not os.path.exists(filename): return
    # split lines into blocks (in case of Link1)
    blocks = split_gaulog_into_gaublocks(filename)
    # Get info of each block
    data = [get_data_from_gaublock(block) for block in blocks]
    # There is nothing to return
    if data == []: return [None]*12
    # Localize data with hessian matrix
    IDX = -1
    for idx,data_i in enumerate(data):
        Fcc =  data_i[7]
        if Fcc is not None:
           IDX = idx
           break
    # Return the best set of data (the last with the hessian or the last block)
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energies,E_oniom,num_imag,zmat = data[IDX]
    # If user does not ask for level, send one of lowest energy
    if target_level is None:
       energies.sort()
       energy,level = energies[0]
    else:
       IDX = None
       exception = Exc.LevelNotFound()
       exception._var = target_level
       for idx,(energy,level) in enumerate(energies):
           if level.lower() == target_level.lower():
               IDX = idx
               break
       if IDX is None: raise exception
       energy, level = energies[IDX]
    # oniom?
    if E_oniom is not None:
       energy = E_oniom
       level  = "ONIOM"
    # Return data
    return commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energy,num_imag,zmat,level
#-------------------------------------------------------#
def gen_zmatrix_string(lzmat,zmatvals,constants=[]):
    string = ""
    all_keys = []
    for idx,zmatline in enumerate(lzmat):
        symbol, connections, keys = zmatline
        # Dummy atoms
        if symbol.upper() in "XX,X,DA": symbol = "X"
        # generate line
        string += "%2s  "%symbol
        for (at,key) in zip(connections,keys):
            string += "%3i  %-7s"%(at+1,key)
        string += "\n"
        # avoid duplicates
        all_keys += [key for key in keys if key not in all_keys]
    # Exclude that keys that are, indeed, numbers
    the_keys = []
    for key in all_keys:
        try   : key = float(key)
        except: the_keys.append(key)
    # Write variables and constants
    string += "Variables:\n"
    for key in the_keys:
        if key.startswith("-") or key in constants: continue
        try   : string += "%-7s %.5f\n"%(key,zmatvals[key])
        except: pass
    if len(constants) != 0:
       string += "Constants:\n"
       for key in constants:
           try   : string += "%-7s %.5f\n"%(key,zmatvals[key])
           except: pass
    return string
#-------------------------------------------------------#
def convert_zmat(lines):
    '''
    basically, a modification of read_xyz_zmat (common.files)
    '''
    lines_values  = [line.replace("="," ") for line in lines if "="     in line]
    lines_zmatrix = [line.replace(","," ") for line in lines if "=" not in line]
    # symbols from zmat
    symbols = []
    atonums = []
    lzmat   = []
    negkeys = []
    for idx,line in enumerate(lines_zmatrix):
        line = line.split()
        # Expected number of columns in this line
        if   idx == 0: expected_cols = 1
        elif idx == 1: expected_cols = 3
        elif idx == 2: expected_cols = 5
        else         : expected_cols = 7
        # Get symbol
        symbol,atonum = symbol_and_atonum(line[0])
        # Get other data
        connections = tuple([int(at_i)-1 for at_i  in line[1:expected_cols:2]])
        keys        = tuple([  key_i     for key_i in line[2:expected_cols:2]])
        # add keys with negative sign
        negkeys += [key_i for key_i in keys if key_i.startswith("-")]
        # save data
        symbols.append(symbol)
        atonums.append(atonum)
        lzmat.append( (symbol,connections,keys) )
    # Get diccionary with values
    zmatvals = {line.split()[0]:float(line.split()[1]) for line in lines_values}
    # Generate another dict
    zmatatoms = {}
    for at1,(symbol,connections,keys) in enumerate(lzmat):
        at2,at3,at4,k12,k123,k1234 = [None for dummy in range(6)]
        if   len(connections) == 1: at2,k12  = connections[0], keys[0]
        elif len(connections) == 2: at2,at3,k12,k123 = connections[0:2]+keys[0:2]
        elif len(connections) == 3: at2,at3,at4,k12,k123,k1234 = connections[0:3]+keys[0:3]
        else: continue
        if k12   is not None: zmatatoms[k12]   = (at1,at2)
        if k123  is not None: zmatatoms[k123]  = (at1,at2,at3)
        if k1234 is not None: zmatatoms[k1234] = (at1,at2,at3,at4)
    # any keyword with negative value
    for key_i in negkeys:
        key = key_i[1:]
        if key in zmatvals.keys(): zmatvals[key_i] = -zmatvals[key]
    # Return data
    return (lzmat,zmatvals,zmatatoms), symbols
#=======================================================#


#=======================================================#
def zmat_from_loginp(log):
    '''
    read Z-matrix from the input lines in the log file
    '''
    with open(log,'r') as asdf: lines = "".join(asdf.readlines())
    if " Symbolic Z-matrix:" not in lines: return None
    lines = lines.split("Symbolic Z-matrix:")[1].split("NAtoms= ")[0].strip()
    zmat = []
    sep = ","
    for line in lines.split("\n"):
        if "Charge" in line: continue
        if ":" in line:
            sep = "="
            continue
        line = line.strip()
        line = sep.join(line.split())
        zmat.append(line)
    return zmat
#=======================================================#




#=======================================================#
def get_fccards_string(gcc,Fcc):
    n1 = len(gcc)
    n2 = len(Fcc)
    assert n2 == n1*(n1+1)/2
    string = " +0.00000000\n"
    # cartesian forces (-gradient)
    for idx in range(0,len(gcc),6):
        line = "".join(["%+12.8f"%(-gcc_i) for gcc_i in gcc[idx:idx+6]])
        string += line+"\n"
    # hessian matrix
    for idx in range(0,len(Fcc),6):
        line = "".join(["%+12.8f"%Fcc_i for Fcc_i in Fcc[idx:idx+6]])
        string += line+"\n"
    return string
#=======================================================#


#=======================================================#
# reading method for Pilgrim                            #
#-------------------------------------------------------#
def read_gauout(filename):
    # read gaussian file
    data_gaulog = read_gaussian_log(filename)
    # split data
    ch      = data_gaulog[2]
    mtp     = data_gaulog[3]
    symbols = data_gaulog[4]
    xcc     = data_gaulog[5]
    gcc     = data_gaulog[6]
    Fcc     = data_gaulog[7]
    V0      = data_gaulog[8]
    level   = data_gaulog[11]
    # symbols and atomic numbers
    symbols,atonums = symbols_and_atonums(symbols)
    # atomic mass
    atomasses = atonums2masses(atonums)
    # return data
    return xcc, atonums, ch, mtp, V0, gcc, Fcc, atomasses, level
#-------------------------------------------------------#
def read_gauout_old(gauout):
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
       for line in extract_string(lines,key1,key2).split("\n")[3:-3]:
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
    strings = extract_string(lines,key1,key2,accumulate=True)
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
        atonums = get_atonums(symbols)
        masses  = atonums2masses(atonums)
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
#=======================================================#



   

#==========================================================#
def set_EXE():
    global EXE
    # Defined in this file
    if 'EXE' in globals():
        return
    # Try to export it from bashrc
    elif KEYGAU in os.environ:
        EXE = os.environ[KEYGAU]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def set_FCHK():
    global FCHK
    # Defined in this file
    if 'FCHK' in globals():
        return
    # Try to export it from bashrc
    elif KEYFCHK in os.environ:
        FCHK = os.environ[KEYFCHK]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    if not os.path.exists(EXE): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def check_FCHK():
    if not os.path.exists(FCHK): raise Exc.ExeNotFound(Exception)
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
    # Execute Gaussian
    command = "%s <%s 1>%s 2>%s"%(EXE,ifile,ofile,err)
    # Try up to NN times
    NN = 5
    for ii in range(NN):
        try   : status = os.system(command)
        except: status = 0
        # interrupted by ctrl+c (status = 2?)
        if status == 2: raise KeyboardInterrupt
        # Get log status
        logstatus = log_status(ofile)
        # Act according log status
        # (a) Normal termination
        if   logstatus == 1       : break
        # (b) empty file, no file or open-new-file error
        elif logstatus in (-1,0,2):
             time.sleep(1)
             continue
        # (c) other type of error
        else: raise exception
    return status
#----------------------------------------------------------#
def log_status(ofile):
    '''
    -1 --> file does not exists
     0 --> files exists but it is empty
     1 --> normal termination
     2 --> open-new-file
     3 --> other...
    '''
    if not os.path.exists(ofile)      : return -1
    with open(ofile,'r') as asdf: lines = asdf.readlines()
    if lines == []                    : return  0
    lline = lines[-1].lower()
    if   "normal termination" in lline: return  1
    elif "open-new-file"      in lline: return  2
    else                              : return  3
#----------------------------------------------------------#
def normal_termination(ofile):
    lines = read_file(ofile)
    if len(lines) == 0: return False
    lastline = lines[-1].lower()
    if "normal termination" in lastline: return True
    else                               : return False
#----------------------------------------------------------#
def genfchk(chk,fchk,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       chk   = folder + chk
       fchk  = folder + fchk
       err   = folder + err
    # fchk tool?
    set_FCHK()
    check_FCHK()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = chk
    # Execute fchk tool
    command = "%s %s %s 1>%s 2>&1"%(FCHK,chk,fchk,err)
    try   : status = os.system(command)
    except: raise exception
    return status
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
    ifile  = wname + ".gjf"
    ofile  = wname + ".log"
    chk    = wname + ".chk"
    fchk   = wname + ".fchk"
    err    = wname + ".err"
    return wname, ifile, ofile, chk, fchk, err
#==========================================================#


