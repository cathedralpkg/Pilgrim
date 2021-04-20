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
| Sub-module :  pilgau             |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

InTerFace betweem Pilgrim and Gaussian
'''

#---------------------------------------------------------------#
import os                            
import time
import common.Exceptions   as Exc    
import common.fncs         as fncs   
import common.physcons     as pc     
import common.gaussian     as ITF    
#---------------------------------------------------------------#


#===============================================================#
key1 = "[Pilgrim_geometry]"
key2 = "[Pilgrim_name]" 
key3 = "[Pilgrim_gradhess]"
#---------------------------------------------------------------#
def pilgrim_template(ch=0,mtp=1,nproc=1,mem=1):
    '''
    Input for Gaussian software
    This calculation IS a single point calculation!
    '''
    level   = "hf/sto-3g"
   #level   = "mpwb95/6-31+G(d,p) IOp(3/76=0560004400) int=ultrafine"

    string  = "%%nproc=%i   \n"%nproc
    string += "%%mem=%iGB   \n"%mem
    string += "%%chk=%s.chk \n"%key2
    string += "#p %s        \n"%level
    string += "scf=verytight\n"
    string += "NoSymm       \n"
    string += "%s           \n"%key3
    string += "\n"
    string += "Input file for MEP calculation\n"
    string += "\n"
    string += "%i %i        \n"%(ch,mtp)
    string += "%s           \n"%key1
    string += "\n"
    return string
#---------------------------------------------------------------#
def pilgrim_templateHL(ch=0,mtp=1,nproc=1,mem=1):
    '''
    Input for Gaussian software
    This calculation IS a single point calculation!
    '''
    level   = "hf/6-31G"

    string  = "%%nproc=%i   \n"%nproc
    string += "%%mem=%iGB   \n"%mem
    string += "%%chk=%s.chk \n"%key2
    string += "#p %s        \n"%level
    string += "scf=verytight\n"
    string += "NoSymm       \n"
    string += "\n"
    string += "Input file for MEP calculation\n"
    string += "\n"
    string += "%i %i        \n"%(ch,mtp)
    string += "%s           \n"%key1
    string += "\n"
    return string
#---------------------------------------------------------------#
def pilgrim_spc(xcc,symbols,bhessian,mname,eargs):
    # initialize
    folder = None
    clean  = False
    oniomH = []
    oniomM = []
    oniomL = []
    frozen = ([],[])
    # expand extra-args
    if   len(eargs) == 1:
        spc_template = eargs[0]
    elif len(eargs) == 2:
        spc_template, folder = eargs
    elif len(eargs) == 3:
        spc_template, folder, clean = eargs
    elif len(eargs) == 5:
        spc_template, folder, clean, frozen, oniom_layers = eargs
        oniomH,oniomM,oniomL = oniom_layers
    # no template??
    if spc_template is None: raise Exc.NoTemplateGiven(Exception)
    # create folder
    if folder is not None:
        if not os.path.exists(folder): os.mkdir(folder)
    # Calculate gradient&hessian or just gradient
    if bhessian: inkey3 = "freq=noraman "
    else       : inkey3 = "force        "
    # names of files
    wname, ifile, ofile, chk, fchk, err = ITF.iofiles(mname,folder)
    # Input
    string_ifile = ""
    for line in spc_template.split("\n"):
        if key1 in line:
           line = ""
           for idx,symbol in enumerate(symbols):
               x,y,z  = fncs.xyz(xcc,idx)
               if   idx in oniomH: layer = " H"
               elif idx in oniomM: layer = " M"
               elif idx in oniomL: layer = " L"
               else              : layer = "  "
               linedata = (symbol,x*pc.ANGSTROM,y*pc.ANGSTROM,z*pc.ANGSTROM,layer)
               line += "%2s  %+11.6f  %+11.6f  %+11.6f %s\n"%linedata
           # add frozen atoms
           frozen_symbols, frozen_xcc = frozen
           for idx,symbol in enumerate(frozen_symbols):
               x,y,z  = fncs.xyz(frozen_xcc,idx)
               # oniom layer?
               at = len(symbols)+idx
               if   at  in oniomH: layer = " H"
               elif at  in oniomM: layer = " M"
               elif at  in oniomL: layer = " L"
               else              : layer = "  "
               linedata = (symbol,x*pc.ANGSTROM,y*pc.ANGSTROM,z*pc.ANGSTROM,layer)
               line += "%2s  -1 %+11.6f  %+11.6f  %+11.6f %s\n"%linedata
        if key2 in line:
           pos  = line.find(key2)
           line = line[0:pos] + wname + line[pos+len(key2):]
        if key3 in line:
           pos  = line.find(key3)
           line = line[0:pos] + inkey3 + line[pos+len(key3):]
        # Add \n to line
        if not line.endswith("\n"): line += "\n"
        string_ifile += line
    with open(ifile,'w') as asdf: asdf.write(string_ifile)
    # Send calculation
    status = ITF.execute(ifile,ofile,err)
    # Check output
    exception = Exc.CalcFails(Exception)
    exception._var = ofile
    if not ITF.normal_termination(ofile): raise exception
    # (a) read log file
    xcc, atonums, ch, mtp, E, gcc, Fcc, masses, clevel = ITF.read_gauout_old(ofile)
    # (b) chk exists! read data from fchk
    if os.path.exists(chk):
       try:
          # Generate fchk
          status = ITF.genfchk(chk,fchk,err)
          # Read data
          xcc, atonums, ch, mtp, E, gcc, Fcc, masses, clevel = ITF.read_fchk(fchk)
       except: pass
    # Remove files
    if clean:
       files = os.listdir(folder)
       files = [fff for fff in files if     fff.startswith(name)   ]
       files = [fff for fff in files if not fff.endswith(".gjf")   ]
       files = [fff for fff in files if not fff.endswith(".log")   ]
       files = [fff for fff in files if not fff.endswith(".out")   ]
       files = [fff for fff in files if not fff.endswith(".chk")   ]
       files = [fff for fff in files if not fff.endswith(".fchk")  ]
       for fff in files: os.remove(folder+fff)
    return xcc, atonums, ch, mtp,  E, gcc, Fcc, masses
#===============================================================#




