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
| Sub-module :  strings            |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#---------------------------------------------------------------#
import os
import getpass
import itertools
import sys
import time
import numpy             as np
#---------------------------------------------------------------#
import modpilgrim.names      as PN
import modpilgrim.steepdesc  as sd
#---------------------------------------------------------------#
import common.fncs       as fncs
import common.internal   as intl
import common.partfns    as partfns
#---------------------------------------------------------------#
from   common.physcons   import AMU
from   common.physcons   import KB, METER, JOULE
from   common.physcons   import KCALMOL
from   common.physcons   import ML
from   common.physcons   import PRE0, VOL0
from   common.physcons   import SECOND
from   common.criteria   import EPS_AMU, EPS_MEPS
#---------------------------------------------------------------#
from   modpilgrim.fit2anarc  import activation1
from   modpilgrim.fit2anarc  import activation2
from   modpilgrim.fit2anarc  import activation3
from   modpilgrim.fit2anarc  import activation4
from   modpilgrim.fit2anarc  import activation5
#---------------------------------------------------------------#

# methods
KEYS_TC = "tst,tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")
KEYS_X  = ("tst,"+\
          "mststzct,mststsct,mscvt,mscvtzct,mscvtsct,"+\
          "mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct").split(",")

# Nice string
KEYNICE = {}
KEYNICE["tst"     ] = "MS-TST"
KEYNICE["mststzct"] = "MS-TST/ZCT"
KEYNICE["mststsct"] = "MS-TST/SCT"
KEYNICE["mscvt"   ] = "MS-CVT"
KEYNICE["mscvtzct"] = "MS-CVT/ZCT"
KEYNICE["mscvtsct"] = "MS-CVT/SCT"
KEYNICE["mptstzct"] = "MP-TST/ZCT"
KEYNICE["mptstsct"] = "MP-TST/SCT"
KEYNICE["mpcvt"   ] = "MP-CVT"
KEYNICE["mpcvtzct"] = "MP-CVT/ZCT"
KEYNICE["mpcvtsct"] = "MP-CVT/SCT"



PROGNAME  = ""
PROGNAME += "\n"
#PROGNAME += "  __________.___.____     __________________.___   _____    \n"
#PROGNAME += "  \______   \   |    |   /  _____/\______   \   | /     \   \n"
#PROGNAME += "   |     ___/   |    |  /   \  ___ |       _/   |/  \ /  \  \n"
#PROGNAME += "   |    |   |   |    |__\    \_\  \|    |   \   /    Y    \ \n"
#PROGNAME += "   |____|   |___|_______ \______  /|____|_  /___\____|__  / \n"
#PROGNAME += "                        \/      \/        \/            \/  \n"
PROGNAME += "  __________.___.____     __________________.___  _     _    \n"
PROGNAME += "  \______   \   |    |   /  _____/\______   \   |/ \___/ \   \n"
PROGNAME += "   |     ___/   |    |  /   \  ___ |       _/   |         \  \n"
PROGNAME += "   |    |   |   |    |__\    \_\  \|    |   \   |  /\_/\   \ \n"
PROGNAME += "   |____|   |___|_______ \______  /|____|_  /___|_/     \  / \n"
PROGNAME += "                        \/      \/        \/             \/  \n"
PROGNAME += "\n"


VERSION   = "2021.2 (2021-04-20)"

PROGHEAD  = " -------------------------------------------------------------\n"
PROGHEAD += "  Program version: Pilgrim v%s\n"%VERSION
PROGHEAD += " -------------------------------------------------------------\n"
PROGHEAD += "                                                            \n"
PROGHEAD += "          A Thermal Rate Constant Calculator and            \n"
PROGHEAD += "              Kinetic Monte Carlo Simulator                 \n"
PROGHEAD += "                                                            \n"
PROGHEAD += " -------------------------------------------------------------\n"

AUTHORINFO = '''
 -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
 |                                                           |
 | Main authors:                                             |
 |     * Ferro-Costas, David       (1)                       |
 |     * Fernandez-Ramos, Antonio  (1)                       |
 |                                                           |
 | In collaboration with:                                    |
 |     * Truhlar, Donald G.        (2)                       |
 |                                                           |
 | (1) Centro Singular de Investigacion en Quimica Bioloxica |
 |     e Materiais Moleculares (CIQUS), Universidade de      |
 |     Santiago de Compostela, Galicia, Spain                |
 |                                                           |
 | (2) Department of Chemistry and Supercomputer Institute,  |
 |     University of Minnesota, Minneapolis, Minnesota       |
 |                                                           |
 -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'''


STRING_END_TXT  = "                                                      Date: %s \n"
STRING_END_TXT += "                                                      Elapsed time: %.1f %s\n"
STRING_END_TXT += "                                                     -----------------------------\n"


#===============================================================#
def exe_info():
    #------------#
    cdate   = time.strftime("%Y-%m-%d")
    ctime   = time.strftime("%H:%M:%S")
    user    = getpass.getuser()
    host    = os.uname()[1]
    pwd     = os.getcwd()
    vinfo   = "%i.%i.%i"%(sys.version_info[0:3])
    #------------#
    string  = "  Current date (YY-MM-DD)   : %s\n"%cdate
    string += "  Current time (HH:MM:SS)   : %s\n"%ctime
    string += "  Python interpreter version: %s\n"%vinfo
    #string += "  User Name / Host Name     : %s in %s\n"%(user,host)
    #string += "  Current directory         : %s\n"%pwd
    return string
#---------------------------------------------------------------#
def init_txt():
    string  = PROGHEAD
    string += exe_info()
    string += " -------------------------------------------------------------\n"
    string += "\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def end_txt(elapsed_time,units):
    # string date and elapsed time
    cdate   = time.strftime("%Y-%m-%d")
    ctime   = time.strftime("%H:%M:%S")
    etime = "%.1f %s"%(elapsed_time,units)
    # length
    ml    = max(len(cdate),len(ctime),len(etime))
    # texts and division
    text1 = "| Current date: %%%is |"%ml%cdate
    text2 = "| Current time: %%%is |"%ml%ctime
    text3 = "| Elapsed time: %%%is |"%ml%etime
    divis = "-"*len(text1)
    # The string
    string  = ""
    string += " "*60 + divis + "\n"
    string += " "*60 + text1 + "\n"
    string += " "*60 + text2 + "\n"
    string += " "*60 + text3 + "\n"
    string += " "*60 + divis + "\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --pfn                          #
#===============================================================#
def title_pfn(ctc,pof):
    title = "Analysis of STRUC: %s"%ctc
    # the string
    string  = ""
    string += "-"*(len(title)+2)+"\n"
    string += " %s \n"%title
    string += "-"*(len(title)+2)+"\n"
    string += "\n"
    string += "    Pilgrim output file: %s\n"%pof
    return string
#---------------------------------------------------------------#
def getstring_pfn1(ltemp,pf_tot,pf_PIB,pf_RR,pf_HO,pf_ele):
    cols = ["T (K)","Qtr","Qrot","Qvib","Qel","Qtot"]
    string  = ""
    string += "-".join(["-"*12     for col in cols])+"\n"
    string += "|".join([fncs.fill_string(col,12) for col in cols])+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    for idx,T in enumerate(ltemp):
        row = pf_PIB[idx], pf_RR[idx], pf_HO[idx], pf_ele[idx], pf_tot[idx]
        row = [" %10.2f "%T]+[" %10.3E "%value for value in row]
        string += "|".join(row)+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    string += "  Qtr : translational pfn per volume unit (in au)\n"
    string += "  Qrot: rotational pfn (rigid-rotor)\n"
    string += "  Qvib: vibrational pfn (harmonic-oscillator) relative to V1\n"
    string += "  Qel : electronic pfn\n"
    string += "  Qtot: total pfn per unit volume (in au)\n"
    string += "  \n"
    string += "  Both Qrot and Qtot include rotational symmetry number\n"
    return string
#---------------------------------------------------------------#
def spfn_igibbs(ltemp,gibbs1cm3,gibbs1bar):
    string  = "-------------------------------------------\n"
    string += "  T (K)  |   V = 1 cm^3   |   V = kbT/p0   \n"
    string += "-------------------------------------------\n"
    for idx,T in enumerate(ltemp):
        g1cm3 = gibbs1cm3[idx]
        g1bar = gibbs1bar[idx]
        string += " %7.2f | %14.8f | %14.8f \n"%(T,g1cm3,g1bar)
    string += "-------------------------------------------\n"
    string += "  V : volume per molecule\n"
    string += "  p0: 1bar\n"
    return string
#---------------------------------------------------------------#
def getstring_pfn2(ltemp,pf_tot,gibbs1cm3,gibbs1bar):
    cols = ["T (K)","QMS_HO","GFE [V = 1 cm^3]","GFE [V = kbT/p0]"]
    string  = ""
    string += "-".join(["-"*18     for col in cols])+"\n"
    string += "|".join([fncs.fill_string(col,18) for col in cols])+"\n"
    string += "-".join(["-"*18     for col in cols])+"\n"
    for idx,T in enumerate(ltemp):
        row = [ "    %10.2f    "%T,"    %10.3E    "%pf_tot[idx],"  %14.8f  "%gibbs1cm3[idx],"  %14.8f  "%gibbs1bar[idx] ]
        string += "|".join(row)+"\n"
    string += "-".join(["-"*18     for col in cols])+"\n"
    string += "   QMS_HO is calculated with regard to min(V1)\n"
    return string
#---------------------------------------------------------------#
def getstring_ctc(sptype,molecules,itcs,V0,V1):
    nconfs = sum([weight for name,weight in itcs])
    fscal  = float(molecules[0]._fscal)
    string  = "Number of conformers: %i\n"%nconfs
    string += "\n"
    string += "   V0 = electronic energy\n"
    string += "   V1 = electronic energy + zero-point energy (ZPE)\n"
    string += "\n"
    string += "   ZPE is calculated using scaled frequencies\n"
    string += "   Frequency scale factor: %.5f\n"%fscal
    string += "\n"
    string += "   min(V0) = %.8f hartree\n"%V0
    string += "   min(V1) = %.8f hartree\n"%V1
    string += "\n"
    string += "   Relative energies (in kcal/mol):\n"
    # table of ctc
    lls      = [4,10,10,7,10,6,5]
    cols     = ["name","V0-min(V0)","V1-min(V1)","ZPE","mass (amu)","weight","PGS"]
    if sptype == 1:
       cols.append("imag.freq.")
       lls.append(10)
    head     = "|".join( [fncs.fill_string(col,ll+2) for ll,col in zip(lls,cols)] )
    division = "-"*len(head)
    # begin table
    string += "   "+division+"\n"
    string += "   "+head+"\n"
    string += "   "+division+"\n"
    refmass = None
    bool_wrongmass = False
    for molecule,itc in zip(molecules,itcs):
        name,weight = itc
        weight = "%i"%weight
        V0i  = "%7.2f"%((molecule._V0  -V0)*KCALMOL)
        V1i  = "%7.2f"%((molecule._ccV1-V1)*KCALMOL)
        zpe  = "%7.2f"%((molecule._ccV1 - molecule._V0)*KCALMOL)
        mass = "%7.2f"%( molecule._mass*AMU)
        pg   = molecule._pgroup
        row  = [name, V0i, V1i, zpe, mass, weight, pg]
        if sptype == 1:
           row.append( "%7.2f"%fncs.afreq2cm(molecule._ccfreqs[0]) )
        line = "|".join( [fncs.fill_string(val,ll+2) for ll,val in zip(lls,row)] )
        string += "   "+line+"\n"
        if refmass is None: refmass = molecule._mass
        diff = abs(refmass-molecule._mass)*AMU
        if diff > EPS_AMU: bool_wrongmass = True
    string += "   "+division+"\n"
    string += "   weight: equals 2 if the structure has a conformational enantiomer,\n"
    string += "           equals 1 otherwise\n"
    string += "   PGS   : point group of symmetry\n"
    if bool_wrongmass:
       string += "\n"
       string += "   ERROR: Mass differs from one conformer to another!\n"
    return string
#---------------------------------------------------------------#
def string_contributions(itcs,ltemp,dchi):
    # get contributions
    nn = 8
    string = ""
    for idx in range(0,len(itcs),nn):
        targets = itcs[idx:idx+nn]
        line = [" T (K) "]+["%5s"%itc for itc,weight in targets]
        head = " "+" | ".join(line)+" "
        divi = "-"*len(head)
        string += divi+"\n"
        string += head+"\n"
        string += divi+"\n"
        for idx,T in enumerate(ltemp):
            line    =  ["%7.2f"%T]+["%5.3f"%dchi[itc][idx] for itc,weight in targets]
            string += " "+" | ".join(line)+" "+"\n"
        string += divi+"\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for MEP                #
#===============================================================#
def smep_title(target,pathvars,pof):
    # title
    title = "Minimum Energy Path for %s"%target
    divis = "="*len(title)
    # begin string
    string  = divis+"\n"
    string += title+"\n"
    string += divis+"\n"
    string += "\n"
    string += "  Pilgrim output file: %s\n"%pof
    string += "\n"
    # reference energies
    if pathvars._reaction is not None:
       string += "  This transition state IS INVOLVED in reaction: '%s'\n"%pathvars._reaction
    else:
       string += "  This transition state IS NOT INVOLVED in any reaction\n"
    if pathvars._beyondmep:
       ctc, itc = PN.name2data(target)
       # E ref
       Eref = pathvars._eref
       string += "     Eref: %.6f hartree\n"%Eref
      ## E0
      #E0   = pathvars._e0
      #if E0 is not None:
      #   E0rel= (E0-Eref)*KCALMOL
      #   string += "     E0  : %.6f hartree ==> %.4f kcal/mol\n"%(E0,E0rel)
    else:
       string += "  WARNING: unable to get Eref from any reaction...\n"
    return string
#---------------------------------------------------------------#
def smep_init(target,software,PARALLEL,var_first,var_sdbw,var_sdfw):
    # split name into info
    ctc, itc = PN.name2data(target)
    # variables first step
    ds,mu,d3,idir = var_first
    string  = "Variables for first step\n"
    string += "   ds        %.4f bohr\n"%ds
    string += "   mu        %.4f amu\n"%(mu*AMU)
    if d3 in [None,False,"no"]: string += "   cubic     no\n"
    else                      : string += "   cubic     %s \n"%d3
    string += "   idir      %s %s\n"%idir
    string += "\n"
    # variables MEP
    method,mu,ds,sbw,hsteps,epse,epsg = var_sdbw
    method,mu,ds,sfw,hsteps,epse,epsg = var_sdfw
    string += "Variables for steepest descent\n"
    string += "   method    %s\n"%method
    string += "   mu        %.4f amu\n"%(mu*AMU)
    string += "   ds        %.4f bohr\n"%ds
    string += "   sbw       %+.4f bohr\n"%sbw
    string += "   sfw       %+.4f bohr\n"%sfw
    string += "   hsteps    %i\n"%hsteps
    string += "   epse      %.2E hartree\n"%epse
    string += "   epsg      %.2E hartree/bohr\n"%epsg
    string += "\n"
    # software and parallization
    string += "Other variables:\n"
    string += "   software  %s (for single-point calcs)\n"%software
    if PARALLEL:
       string += "   paral     yes (both sides of MEP will be calculated in unison)\n"
    else:
       string += "   paral     no  (both sides of MEP will NOT be calculated in unison)\n"
    return string
#---------------------------------------------------------------#
def smep_ff(f1,f2,f3,rstfile,xyzfile):
    string  = "Folders of interest:\n"
    string += "   temporal for calcs: %s\n"%f1
    string += "   folder with rst files : %s\n"%f2
    string += "   folder with xyz files : %s\n"%f3
    string += "\n"
    string += "Files of interest:\n"
    string += "   rst file: %s\n"%rstfile
    string += "   xyz file: %s\n"%xyzfile
    return string
#---------------------------------------------------------------#
def smep_rst(rstfile,drst):
    lbw, lfw, sbw, sfw, Ebw, Efw = sd.rstlimits(drst)
    if not os.path.exists(rstfile): return ""
    string  = "Restart file (rst) found!\n"
    if sbw is None:
        string += "   no data in it :(\n"
    else:
        ml = max(len(lbw),len(lfw))
        lbw = "%%-%is"%ml%lbw
        lfw = "%%-%is"%ml%lfw
        string += "   first point (%s): %+7.3f bohr (%+11.6f hartree)\n"%(lbw,sbw,Ebw)
        string += "   last  point (%s): %+7.3f bohr (%+11.6f hartree)\n"%(lfw,sfw,Efw)
    return string
#---------------------------------------------------------------#
def smep_ts(ts):
    ''' ts is an instance of Molecule'''
    string  = "\n"
    string += "Information about transition structure:\n"
    string += "\n"
    string += "   File with info    : %s\n"%ts._gts
    string += ts.info_string(3)
    return string
#---------------------------------------------------------------#
def smep_oniom(oniomlayers,natoms,software):
    oniomh, oniomm, onioml = oniomlayers
    num_oniom = len(set(oniomh+oniomm+onioml))
    if num_oniom == 0: return ""
    string  = "ONIOM layers defined!\n"
    string += "\n"
    string += "   number of atoms in layer H     : %i\n"%len(oniomh)
    string += "   number of atoms in layer M     : %i\n"%len(oniomm)
    string += "   number of atoms in layer L     : %i\n"%len(onioml)
    string += "   total number of atoms in layers: %i\n"%num_oniom
    if num_oniom != natoms:
        string += "   ERROR: number of atoms is %i\n"%natoms
    if software != "gaussian":
        string += "   ERROR: ONIOM layers is only implemented in Gaussian!\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def smep_first(symbols,xms,v0,v1,oniomlayers):
    oniomh, oniomm, onioml = oniomlayers
    natoms = fncs.howmanyatoms(v0)
    string = "mass-scaled coordinates (in bohr):\n"
    for at,symbol in enumerate(symbols):
        vx,vy,vz = fncs.xyz(xms,at)
        if   at in oniomh: soniom = " [ONIOM LAYER: H]"
        elif at in oniomm: soniom = " [ONIOM LAYER: M]"
        elif at in onioml: soniom = " [ONIOM LAYER: L]"
        else             : soniom = ""
        string += "  %2s   %+10.6f  %+10.6f  %+10.6f %s\n"%(symbol,vx,vy,vz,soniom)
    string += "\n"
    string += "mass-scaled v0 vector:\n"
    for at,symbol in enumerate(symbols):
        vx,vy,vz = fncs.xyz(v0,at)
        string += "  %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,vx,vy,vz)
    if v1 is not None:
       string += "\n"
       string += "v1 vector (mass-scaled):\n"
       for at,symbol in enumerate(symbols):
           vx,vy,vz = fncs.xyz(v1,at)
           string += "  %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,vx,vy,vz)
    return string
#---------------------------------------------------------------#
def smep_table(drst,Eref):
    string  = "Reference energy (Eref) set to: %.6f hartree\n"%Eref
    string += "\n"
    string += "     -------------------------------------------------------------------\n"
    string += "      s (bohr) |   E (hartree)   | E-Eref (hartree) | E-Eref (kcal/mol) \n"
    string += "     -------------------------------------------------------------------\n"
    points = sd.sorted_points(drst,hess=True)
    for point in points:
        s_i, E_i = drst[point][0:2]
        dE_i      = (E_i-Eref)
        dE_i_kcal = (E_i-Eref)*KCALMOL
        string += "      %+8.4f | %+15.7f | %+16.7f |    %+12.3f   \n"%(s_i,E_i,dE_i,dE_i_kcal)
    string += "     -------------------------------------------------------------------\n"
    return string
#---------------------------------------------------------------#
def smep_tableDLEVEL(drst,tdleveldata,Eref):
    points, xx, yyll, yyhl = tdleveldata

    string  = "Reference energy (Eref) set to: %.6f hartree\n"%Eref
    string += "\n"
    string += "     --------------------------------------------------------------------------------------\n"
    string += "      s (bohr) | E_LL  (hartree) || E_HL  (hartree) | E-Eref (hartree) | E-Eref (kcal/mol) \n"
    string += "     --------------------------------------------------------------------------------------\n"
    for idx,point in enumerate(points):
        if drst[point][4] is None: continue
        s_i, Ell_i, Ehl_i = xx[idx],yyll[idx],yyhl[idx]
        dE_i      = (Ehl_i-Eref)
        dE_i_kcal = (Ehl_i-Eref)*KCALMOL
        string += "      %+8.4f | %+15.7f || %+15.7f | %+16.7f |    %+12.3f   \n"%(s_i,Ell_i,Ehl_i,dE_i,dE_i_kcal)
    string += "     --------------------------------------------------------------------------------------\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for adiab pot          #
#===============================================================#
def sadipot_ics(ics,inter,icsbw=None,icsfw=None):
    nepl = 4
    string = ""
    if inter == "no": return string
    for IDX,the_ics in enumerate([ics,icsbw,icsfw]):
        if the_ics is None: continue
        if IDX == 0: string += "- Internal coordinates were found!\n"
        if IDX == 1: string += "- Internal coordinates for BACKWARD direction were found!\n"
        if IDX == 2: string += "- Internal coordinates for FORWARD  direction were found!\n"
        string += "\n"
        ics_st,ics_ab,ics_lb,ics_it,ics_pt = intl.unmerge_ics(the_ics)
        for tics,lics in [("st",ics_st),("ab",ics_ab),("lb",ics_lb),("it",ics_it),("pt",ics_pt)]:
            nics = len(lics)
            if   nics == 0: continue
            elif nics == 1:
                if tics == "st": string += "  1 stretch:\n"
                if tics == "ab": string += "  1 bent bond angle:\n"
                if tics == "lb": string += "  1 linear bond angle:\n"
                if tics == "it": string += "  1 improper torsion:\n"
                if tics == "pt": string += "  1 proper torsion:\n"
            else:
                if tics == "st": string += "  %i stretches:\n"%len(lics)
                if tics == "ab": string += "  %i bent bond angles:\n"%len(lics)
                if tics == "lb": string += "  %i linear bond angles:\n"%len(lics)
                if tics == "it": string += "  %i improper torsions:\n"%len(lics)
                if tics == "pt": string += "  %i proper torsions:\n"%len(lics)
            for idx in range(0,len(lics),nepl):
                string += " "*4+" ".join(["%-12s"%intl.ic2string((tics,ic)) \
                                          for ic in lics[idx:idx+nepl]])+"\n"
            string += "\n"
    return string
#---------------------------------------------------------------#
def sadipot_table(ls, lV1, sAG, VAG, lV0, lcc_zpe, lic_zpe, Eref):
    string  =  "-----------------------------------------------------------------------\n"
    string +=  " s (bohr) |  V_MEP  | ZPE(cc) | ZPE(ic) |    VaG    || Eref + VaG (au) \n"
    string +=  "-----------------------------------------------------------------------\n"
    for idx,s_i in enumerate(ls):
        s_i   = "%+7.3f"%s_i
        V0    = "%+7.3f"%(lV0[idx]*KCALMOL)
        V1    = "%+9.3f"%(lV1[idx]*KCALMOL)
        V1_au = "%.8f"%(lV1[idx]+Eref)
        cczpe = "%7.3f"%(lcc_zpe[idx] *KCALMOL)
        arrow = "   "
        try:
          iczpe = "%7.3f"%(lic_zpe[idx]*KCALMOL)
          #if float(cczpe)-float(iczpe) > 0.1: arrow = "<--"
        except:
          iczpe = "   -   "
        string +=  " %7s  | %7s | %7s | %7s | %9s || %15s  %s\n"%(s_i,V0,cczpe,iczpe,V1,V1_au,arrow)
    string +=  "-----------------------------------------------------------------------\n"
    # Maximum
    sAG = "%+7.3f"%sAG
    VAG_au = "%.8f"%(VAG+Eref)
    VAG = "%+9.3f"%(VAG*KCALMOL)
    string +=  " %7s  |    Maximum of VaG (VAG)     | %9s || %15s\n"%(sAG,VAG,VAG_au)
    string +=  "-----------------------------------------------------------------------\n"
    string +=  "  V_MEP, ZPE and VaG in kcal/mol\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def sadipot_freqs(ls,lcc_frqs,lic_frqs):
    nepl = 5
    # Generate matrix with all data
    matrix_cc = fncs.ll2matrix(lcc_frqs,varcomp=None,pos=0)
    matrix_cc = [list(xxx) for xxx in zip(*matrix_cc)]
    nrows, ncols = len(matrix_cc), len(matrix_cc[0])
    # same for internal
    matrix_ic = fncs.ll2matrix(lic_frqs,varcomp=None,pos=0)
    if matrix_ic is not None: matrix_ic = [list(xxx) for xxx in zip(*matrix_ic)]
    # Print frequencies
    string = ""
    for idxi in range(0,ncols,nepl):
        svalues = " idx | "+"|".join(["   %+8.4f   "%ss for ss in ls[idxi:idxi+nepl]])
        string +=  "-"*len(svalues)+"\n"
        string +=          svalues +"\n"
        string +=  "-"*len(svalues)+"\n"
        # add all frequencies
        for row in range(nrows):
            string += " %3i | "%row
            for col in range(idxi,min(idxi+nepl,ncols)):
                fcc = matrix_cc[row][col]
                if fcc is not None: string +=  "%5.0f "%fncs.afreq2cm(fcc)
                else              : string +=  "  -   "
                # internal
                if matrix_ic is None: fic = None
                else                : fic = matrix_ic[row][col]
                if fic is not None: string += "(%5.0f)"%fncs.afreq2cm(fic)
                else              : string += "(  -  )"
                if col+1 == min(idxi+nepl,ncols): string += "\n"
                else                            : string += " |"
    return string
#---------------------------------------------------------------#
def sadipot_checks(ok1,ok2,ok3,useics,svalues_with_imag):
    string  = ""
    string += "- cc is used to indicate Cartesian coordinates\n"
    string += "- ic is used to indicate Internal  coordinates\n"
    string += "\n"
    # print which ZPE is used to calculate Vadi
    if useics == "yes": string += "- ZPE(ic) will be used to calculate Vadi\n"
    else              : string += "- ZPE(cc) will be used to calculate Vadi\n"
    string += "\n"
    #string = "Comparing Cartesian coordinates (cc) vs internal coordinates (ic) frequencies...\n"
    if useics == "yes": ftype = "ic"
    else              : ftype = "cc"
    # info about what failed
    if not ok1      : string += "- WARNING! Different number of cc and ic frequencies at some point(s)!\n"
    if not ok2      : string += "- WARNING! Transition structure cc and ic frequencies differ!\n"
    if not ok3      : string += "- WARNING! ZPE(cc)-ZPE(ic) > 0.1 kcal/mol at some point(s)!\n"
    # imaginary freqs
    if len(svalues_with_imag) !=0:
       string += "- WARNING! There are imaginary %s frequencies along the MEP!!\n"%ftype
       string += "           Values of s (bohr) with imaginary freq(s):\n"
       nn = 5
       for idx in range(0,len(svalues_with_imag),nn):
           string += "           "+" ".join(["%+9.4f"%svalue for svalue in svalues_with_imag[idx:idx+nn]])+"\n"
    else:
       string += "- Fine! There are no imaginary %s frequencies along the MEP!\n"%ftype
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for CVT                #
#===============================================================#
def scvt_gibbs(svals,temps,mgibbs,pathvars,gibbsTS):
    #---------------------------------------
    def get_table_cvtgibbs(svals,temps,mgibbs):
        nrows = len(svals)
        nepl = 5
        matrix_string = ""
        for idx in range(0,len(temps),nepl):
            colnames = ["  s value  "]+[" %7.2f K "%T for T in temps[idx:idx+nepl]]
            thead    = " "+"|".join(colnames)+" "
            tdiv     = "-"*len(thead)
            matrix_string += tdiv +"\n"
            matrix_string += thead+"\n"
            matrix_string += tdiv +"\n"
            for row in range(nrows):
                values  = [" %+9.4f "%svals[row]]
                values += [" %9.3f "%(gibbs*KCALMOL) for gibbs in mgibbs[row][idx:idx+nepl]]
                line    = " "+"|".join(values)+" "
                matrix_string += line+"\n"
            matrix_string += tdiv +"\n\n"
        return matrix_string
    #---------------------------------------
    nrows, ncols = mgibbs.shape
    string = ""
    if gibbsTS is not None:
      #string += "Matrix of gibbs free energy with regards to the saddle point structure\n"
       string += "Gibbs free energy profile [with regard to the saddle point structure]\n"
       string += "    - shape: %i x %i\n"%(nrows,ncols)
       string += "    - units: kcal/mol\n"
      #string += "    - sigma_rot was set to 1 along MEP for the\n"
      #string += "      calculation of Gibbs free energy profile\n"
       string += "\n"
       string += get_table_cvtgibbs(svals,temps,mgibbs)
    # gibbs matrix with regards reactants
    if pathvars._GibbsR is not None:
       DGRTS = [gTS-gR for gR,gTS in zip(pathvars._GibbsR,gibbsTS)]
       for col in range(ncols):
           for row in range(nrows):
               mgibbs[row][col] += DGRTS[col]
      #string += "Matrix of gibbs free energy with regards to '%s' reactant(s)\n"%pathvars._reaction
       string += "Gibbs free energy profile [with regards to '%s' reactant(s)]\n"%pathvars._reaction
       string += "    - shape: %i x %i\n"%(nrows,ncols)
       string += "    - units: kcal/mol\n"
      #string += "    - sigma_rot was set to 1 along MEP for the\n"
      #string += "      calculation of Gibbs free energy profile\n"
       string += "\n"
       string += get_table_cvtgibbs(svals,temps,mgibbs)
    return string
#---------------------------------------------------------------#
def scvt_coefs(lcvt_s, lcvt_gamma, temps):
    string  = "CVT variational coefficient: \n"
    string += "    - variation in Gibbs free energy (DGFE) in kcal/mol\n"
    string += "           DDGFE = -R T ln(Gamma^CVT) \n"
    string += "\n"
    string += "    -----------------------------------------------\n"
    string += "      T (K)  |  s_CVT  |  Gamma^CVT  |   DDGFE     \n"
    string += "    -----------------------------------------------\n"
    for idx,T in enumerate(temps):
        s_i, gamma_i = lcvt_s[idx], lcvt_gamma[idx]
        gibbs        = - KB * T * np.log(gamma_i) * KCALMOL
        string += "     %7.2f | %+7.4f | %11.4E | %9.4f \n"%(T, s_i,gamma_i,gibbs)
    string += "    -----------------------------------------------\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for SCT                #
#===============================================================#
def ssct_init(E0,VadiSpl,pathvars,v1mode):
    sbw, Vbw = VadiSpl.get_alpha()
    sts, Vts = VadiSpl.get_saddle()
    sfw, Vfw = VadiSpl.get_omega()
    sAG, VAG = VadiSpl.get_max()
    V1R      = pathvars._V1R
    V1P      = pathvars._V1P
    string  = "Energy limits for tunnelling [E0,VAG]:\n"
    string += "\n"
    string += "  ---------------------------------------\n"
    string += "            | s (bohr) | Vadi (kcal/mol) \n"
    string += "  ---------------------------------------\n"
    if V1R is None:
       string += "  reactants |          |    --    \n"
    else:
       V1R = V1R - pathvars._eref
       string += "  reactants |          | %+8.3f \n"%(V1R*KCALMOL)
    string += "  bw        | %+8.4f | %+8.3f \n"%(sbw,Vbw*KCALMOL)
    string += "  saddle    | %+8.4f | %+8.3f \n"%(sts,Vts*KCALMOL)
    string += "  fw        | %+8.4f | %+8.3f \n"%(sfw,Vfw*KCALMOL)
    if V1P is None:
       string += "  products  |          |    --    \n"
    else:
       V1P = V1P - pathvars._eref
       string += "  products  |          | %+8.3f \n"%(V1P*KCALMOL)
    string += "  =======================================\n"
    string += "  E0        |          | %+8.3f \n"%(E0*KCALMOL)
    string += "  VAG       | %+8.4f | %+8.3f \n"%(sAG,VAG*KCALMOL)
    string += "  =======================================\n"
    string += "\n"
    if type(pathvars._e0) == float:
       string += "  E0 defined by user in '%s'\n"%PN.IFILE3
       string += "\n"
    string += "v1mode  : %s\n"%v1mode
    string += "muintrpl: %s %i\n"%(pathvars._muintrpl)
    return string
#---------------------------------------------------------------#
def ssct_mueff(svals,VadiSpl,lkappa,ltbar,lmueff,mu,toignore=[]):
    string  = "Summary table\n"
    string += "  - Progress along the path (s) in bohr\n"
    string += "  - Vibrationally adiabatic potential (VaG) in kcal/mol\n"
    string += "  - kappa (curvature) in bohr^-1\n"
    string += "  - Turning point (turnpoint) in bohr\n"
    string += "  - Effective mass (mueff) in a.u.\n"
    string += "  ----------------------------------------------------------------------\n"
    string += "      s    |    VaG    |   kappa   | turnpoint |   mueff   |  mueff/mu  \n"
    string += "  ----------------------------------------------------------------------\n"
    for idx in range(len(svals)):
         s       = "%+7.3f"%svals[idx]
         Vadi    = "%+9.3f"%(VadiSpl(svals[idx])*KCALMOL)
         kappa   = "%+9.2E"%lkappa[idx]
         tbar    = "%9.5f"%ltbar[idx]
         mu_eff  = "%9.4f"%(lmueff[idx])
         ratiomu = "%7.4f"%(lmueff[idx]/mu)
         string += "   %s | %s | %s | %s | %s | %s  "%(s,Vadi,kappa,tbar,mu_eff,ratiomu)
         if idx in toignore: string += "**"
         string += "\n"
    string += "  ----------------------------------------------------------------------\n"
    string += "   ** kappa, turnpoint and mueff were interpolated\n"
    return string
#---------------------------------------------------------------#
def ssct_E0VAG(E0,VAG):
    string  = "Transmission probabilities will be calculated between E0 and VAG:\n"
    string += "   E0  = %8.4f kcal/mol\n"%( E0 *KCALMOL )
    string += "   VAG = %8.4f kcal/mol\n"%( VAG*KCALMOL )
    string += "\n"
    return string
#---------------------------------------------------------------#
def ssct_E0_above_VAG(E0,VAG):
    string  = "WARNING! E0 > VAG:\n"
    string += "   E0  = %7.3f kcal/mol\n"%(E0 *KCALMOL)
    string += "   VAG = %7.3f kcal/mol\n"%(VAG*KCALMOL)
    string += "   --> ZCT/SCT probabilities will not be calculated!\n"
    return string
#---------------------------------------------------------------#
def ssct_probs(E_list,probs_ZCT,probs_SCT,rpoints_SCT,sbw,sfw):
    string = ""
    # Get P(E) for all these energies
    head1 = " E [kcal/mol]  |  P^ZCT(E)  |  P^SCT(E)  |  Classical turning "
    head2 = "               |            |            |  points [bohr]     "
    divi = "-"*len(head1)
    string += "  "+divi+"\n"
    string += "  "+head1+"\n"
    string += "  "+head2+"\n"
    string += "  "+divi+"\n"
    write_message = False
    for idx in range(len(E_list)):
        E           = "%10.4f"%(E_list[idx]*KCALMOL)
        pE_ZCT      = "%10.3e"%(probs_ZCT[idx])
        pE_SCT      = "%10.3e"%(probs_SCT[idx])
        rpoints     = rpoints_SCT[idx]
        # string return points
        rp_string = ""
        for idx,(start,end) in enumerate(rpoints):
            rp_string += "[%+.3f,%+.3f]U"%(start,end)
            if (idx+1)%2==0 and idx+1 != len(rpoints): rp_string+= "\n"+" "*46
        # are probabilities converged?
        try:
            sbw1  = rpoints[ 0][0]
            sfw1  = rpoints[-1][1]
            bool1 = abs(sbw1-sbw) < EPS_MEPS
            bool2 = abs(sfw1-sfw) < EPS_MEPS
            if bool1 and bool2: endline = "B"
            elif bool1        : endline = "L"
            elif bool2        : endline = "B"
            else              : endline = ""
        except: endline = ""
        if endline != "": write_message = True
        # add to string
        string += "     %s  | %s | %s |  %s  %s\n"%(E,pE_ZCT,pE_SCT,rp_string[:-1],endline)
    string += "  "+divi+"\n"
    string += "   Number of tunnelling energies: %i\n"%len(E_list)

    if write_message:
       string += "\n"
       string += "   WARNING! Some tunnelling probabilities are not converged:\n"
       string += "     * 'L' --> at the left-side of VaG(s).\n"
       string += "     * 'R' --> at the right-side of VaG(s).\n"
       string += "     * 'B' --> at both sides of VaG(s).\n"
       string += "   This fact may lack of importance if the temperature\n"
       string += "   is high enough so these probabilities do not play\n"
       string += "   any role.\n"
       #string += "   When the MEP is not complete, the value of VaG(s)\n"
       #string += "   may be greater than E.\n"
    return string
#---------------------------------------------------------------#
def ssct_diffs(lE_SCT,diffs_SCT):
    maxdiff = 10.0
    nepl    = 5
    bigdiff = [E for E,diff in zip(lE_SCT,diffs_SCT) if diff > maxdiff]
    string = ""
    if len(bigdiff) > 0:
       string += "Action integral (theta) is calculated by:\n"
       string += "   * Gaussian quadrature (theta_g)\n"
       string += "   * trapezoidal rule    (theta_t)\n"
       string += "and theta = min(theta_g,theta_t)\n"
       string += "\n"
       string += "WARNING: differences > %i%% for the following energies:\n"%maxdiff
       for idx in range(0,len(bigdiff),nepl):
           string += "         "
           string += " ".join(["%8.4f"%(E*KCALMOL) for E in bigdiff[idx:idx+nepl]])
           string += "\n"
      #string += "   * theta = theta_t in those cases\n"
       string += "\n"
    return string
#---------------------------------------------------------------#
def ssct_qrc(pathvars):
    string  = "Quantum reaction coordinate keyword (qrc) activated!\n"
    string += "\n"

    # Any error?
    qrccase = pathvars._qrccase
    error1 = "No reaction related to the transition state!"
    error2 = "The reaction is not unimolecular!"
    error3 = "Reactant is not defined!"
    error4 = "The reactant gts file was not found"
    error5 = "Unable to get product(s) energy!"
    if qrccase == 1: string += "   * ERROR: %s\n"%error1; return string, False
    if qrccase == 2: string += "   * ERROR: %s\n"%error2; return string, False
    if qrccase == 3: string += "   * ERROR: %s\n"%error3; return string, False
    if qrccase == 4: string += "   * ERROR: %s\n"%error4; return string, False
    if qrccase == 5: string += "   * ERROR: %s\n"%error5; return string, False

    # Print reaction and reactant information
    reaction   = pathvars._reaction
    reactioneq = pathvars._reactioneq
    qrcname    = pathvars._qrcname
    mode       = pathvars._qrc[0]+1
    afreq_cm   = fncs.afreq2cm(pathvars._qrcafreq)
    numst      = pathvars._qrc[1]

    string += "   * reaction name & eq  : %s (%s)\n"%(reaction,reactioneq)
    if not pathvars._exorgic:
        string += "     WARNING!\n"
        string += "     QRC must be used in exoergic direction!\n"
        string += "     Reaction will be reversed for QRC calculation!\n"
    string += "   * reactant name       : %s\n"%qrcname
    string += "   * reactant vib. freq. : %i (%.2f cm^-1)\n"%(mode,afreq_cm)
    string += "   * number of states    : %i\n"%numst
    string += "   * contribution to Kappa^SAG from E0 to VAG will\n"
    string += "     be obtained from discrete set of energies\n"
    string += "\n"
    string += "   * calculating transmission probabilities...\n"
    return string, True
#---------------------------------------------------------------#
def ssct_kappa(temps,KAPPA,lIi,RTE,E0,bqrc,case="sct"):
    string  = "%s transmission coefficient:\n"%(case.upper())
    string += "----------------------------------------------------\n"
    string += "   T (K)   |  %%I1   |  %%I2   |  Kappa^%3s  |   RTE   \n"%(case.upper())
    string += "----------------------------------------------------\n"
    for idx,T in enumerate(temps):
        I1, I2, I3 = lIi[idx]
        kappa, rte = KAPPA[idx], RTE[idx]*KCALMOL
        I1 *= 100./kappa
        I2 *= 100./kappa
        I3 *= 100./kappa
        string += " %9.2f | %6.2f | %6.2f | %11s | %7.3f "%(T,I1,I2+I3,fncs.eformat(kappa,3),rte)
        if bqrc[idx]             : string += "++ "
        if (rte-E0*KCALMOL) < 1.5: string += "** "
        string += "\n"
    string += "----------------------------------------------------\n"
    string += " RTE: Representative Tunnelling Energy (in kcal/mol)           \n"
    string += " %I1: contribution of tunnelling\n"
    string += " %I2: contribution of non-classical reflection\n"
    string += " ++ : indicates that QRC was used at the given temperature\n"
    string += " ** : indicates that RTE is close to E0 (less than 1.5 kcal/mol)\n"
    return string
#---------------------------------------------------------------#
def ssct_onlykappa(temps,kappa_ZCT,kappa_SCT):
    string  = "Transmission coefficients:\n"
    string += "---------------------------------------\n"
    string += "   T (K)   |  Kappa^ZCT  |  Kappa^SCT  \n"
    string += "---------------------------------------\n"
    for idx,T in enumerate(temps):
        zct, sct = kappa_ZCT[idx], kappa_SCT[idx]
        string += " %9.2f | %11s | %11s \n"%(T,fncs.eformat(zct,4),fncs.eformat(sct,4))
       #string += " %9.2f | %11.4E | %11.4E \n"%(T,zct,sct)
    string += "---------------------------------------\n"
    return string
#---------------------------------------------------------------#
def ssct_convergence(convlist_sct,convlist_lim,scterr,converged=False):
    string  = "Convergence of Kappa^SCT for the lowest temperature \n"
    string += "-------------------------------------------------\n"
    string += "\n"
    string += "     scterr = %.2f %%\n"%scterr
    string += "\n"
    string += "     --------------------------------------------------------------------\n"
    string += "       step   |   sbw   |   sfw   |  Kappa^SCT  |  diff(%)  | converged? \n"
    string += "     --------------------------------------------------------------------\n"
    nn = len(convlist_sct)
    for idx, SCT in enumerate(convlist_sct):
        sbw, sfw = convlist_lim[idx]
        if idx == 0:
            diff  = "    -    "
            sconv = "    -    "
        else:
            SCT_pr = convlist_sct[idx-1]
            diff = "%9.2f"%float(100*abs(SCT-SCT_pr)/SCT)
            if idx+1 != nn: sconv = "    NO    "
            elif converged: sconv = "    YES   "
            else          : sconv = "    NO    "
        string += "          %3i | %+7.3f | %+7.3f | %11s | %s | %s \n"%(idx,sbw,sfw,fncs.eformat(SCT,4),diff,sconv)
       #string += "          %3i | %+7.3f | %+7.3f | %+11.4E | %s | %s \n"%(idx,sbw,sfw,SCT,diff,sconv)
    string += "     --------------------------------------------------------------------\n"
    string += "     When data files for products are not given,\n"
    string += "     diff must be smaller than scterr twice in a row\n"
    # add square around for visualization
    with_square  = "="*82+"\n"
    for line in string.split("\n"):
        with_square += "|| " + "%-77s"%line + "||\n"
    with_square += "="*82+"\n"
    string = with_square
    # return string
    return string
#===============================================================#

#===============================================================#
#                STRINGS USED IN --path: for CAG                #
#===============================================================#
def scag_table(ltemp, dE_cagtst, cagtst, dE_cagcvt, cagcvt):
    string  = "    --------------------------------------------------------------------\n"
    string += "      T (K)  |  VAG-VaG(s=0)  |  TST/CAG  || VAG-VaG(s_CVT) |  CVT/CAG  \n"
    string += "    --------------------------------------------------------------------\n"
    for idx,T in enumerate(ltemp):
        dE1 = "%14.4f"%(dE_cagtst[idx]*KCALMOL)
        cf1 = "%9.3E"%(cagtst[idx])
        if dE_cagcvt is None or cagcvt is None:
           dE2 = "      -       "
           cf2 =   "    -    "
        else:
           dE2 = "%14.4f"%(dE_cagcvt[idx]*KCALMOL)
           cf2 = "%9.3E"%(cagcvt[idx])
        string += "     %7.2f | %s | %s || %s | %s\n"%(T,dE1,cf1,dE2,cf2)
    string += "    --------------------------------------------------------------------\n"
    string += "      Energy differences in kcal/mol\n"
    return string
#===============================================================#


#===============================================================#
#                                                               #
#===============================================================#
def spath_allcoefs(ltemp,dcoefs):
    keys  = ["CVT"  ,"ZCT"  ,"SCT"  ,"CAGTST" ,"CAGCVT" ]
    cols1 = ["Gamma","Kappa","Kappa"," Kappa "," Kappa "]
    cols2 = [" CVT "," ZCT "," SCT ","TST/CAG","CVT/CAG"]
    # print coefs again
    table = ""
    thead1 = ["  T (K) "]+[fncs.fill_string(col,11) for col in cols1]
    thead2 = ["        "]+[fncs.fill_string(col,11) for col in cols2]
    thead1 = " "+" | ".join(thead1)+" "
    thead2 = " "+" | ".join(thead2)+" "
    tdivi = "-"*len(thead1)
    table += tdivi+"\n"
    table += thead1+"\n"
    table += thead2+"\n"
    table += tdivi+"\n"
    for idx,T in enumerate(ltemp):
        tline = ["%8.2f"%T]
        for coef in keys:
            coef = coef.lower()
            if coef in dcoefs.keys():
                value = "%11s"%fncs.eformat(dcoefs[coef][idx],4)
               #value = "%11.4E"%dcoefs[coef][idx]
            else:
                value = "     -     "
            tline.append(value)
        tline = " "+" | ".join(tline)+" "
        table += tline+"\n"
    table += tdivi+"\n"
    # return string
    string  = "SUMMARY OF CALCULATED COEFFICIENTS:\n"
    string += "\n"
    for line in table.split("\n"): string += "    "+line+"\n"
    string += "\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --rcons                        #
#===============================================================#
def srcons_init(rcname,pof,Rs,TS,Ps):
    # Title
    title    = "| Reaction to analyze: %s |"%rcname
    division = "-"*len(title)
    # the string
    string  = division+"\n"
    string += title   +"\n"
    string += division+"\n"
    string += "\n"
    string += "  Pilgrim output file: %s\n"%pof
    string += "\n"
    string += "  Chemical equation: "
    string += " + ".join(Rs)
    string += " --> "
    if TS is not None: string += TS
    string += " --> "
    string += " + ".join(Ps)
    string += "\n\n"
    # Identify reactants, transition state and products
    string += "     reactant(s)      ==> %s\n"%(" + ".join(Rs))
    if TS is not None: string += "     transition state ==> %s\n"%TS
    string += "     product(s)       ==> %s\n"%(" + ".join(Ps))
    string += "\n"
    if len(Rs) == 2 and Rs[0] == Rs[1]:
        string += "  The two reactants are equal.\n"
        string += "  Forward rate constants will be multiplied by 2\n"
        string += "\n"
    if len(Ps) == 2 and Ps[0] == Ps[1]:
        string += "  The two products are equal.\n"
        string += "  Backward rate constants will be multiplied by 2\n"
        string += "\n"
    return string
#---------------------------------------------------------------#
def srcons_conservation(chemreac):
    string  = "Conservation of charge and mass:\n"
    string += "\n"
    string += "   ----------------------------------------\n"
    string += "                     | charge | mass (amu) \n"
    string += "   ----------------------------------------\n"
    if chemreac._nR != 0:
       string += "    reactant(s)      |  %3i   | %10.3f \n"%(chemreac._chR,chemreac._massR*AMU)
    if chemreac._ts is not None:
       string += "    transition state |  %3i   | %10.3f \n"%(chemreac._chTS,chemreac._massTS*AMU)
    if chemreac._nP != 0:
       string += "    product(s)       |  %3i   | %10.3f \n"%(chemreac._chP,chemreac._massP*AMU)
    string += "   ----------------------------------------\n"
    if chemreac._problem:
       string += "    charge and/or mass NOT conserved!!\n"
    return string
#---------------------------------------------------------------#
def weights_as_dict(dctc):
    dweights = {}
    for ctc,cluster in dctc.items():
        for itc,weight in cluster._itcs:
            dweights[ (ctc,itc) ] = weight
    return dweights
#---------------------------------------------------------------#
def srcons_relenergies(chemreac):

    string  = "Relative energies (kcal/mol):\n"
    string += "\n"
    string += "   "+"V0(i) is the electronic energy of the i-th conformer\n"
    string += "   "+"V1(i) = V0(i)+ZPE(i)\n"
    string += "   "+"ZPE(i) is the harmonic oscillator ZPE of the i-th conformer\n"
    string += "\n"
    string += "   "+"min{V0(i)} of reactants ==> V0 = %.8f hartree\n"%chemreac._V0R
    string += "   "+"min{V1(i)} of reactants ==> V1 = %.8f hartree\n"%chemreac._V1R
    string += "\n"

    dweights = weights_as_dict(chemreac._dctc)

    # table lines
    tlines = []
    for targets in [chemreac._reacts,[chemreac._ts],chemreac._prods]:
        if len(targets) == 0: continue
        if None in targets: continue
        # generate combinations
        combinations = []
        for target in targets:
            # the itcs
            ctc,itcs,ms = chemreac._itcs[target]
            combinations.append( [PN.struckey(ctc,itc) for itc,weight in itcs] )
        while len(combinations) != 1:
              combined = list(itertools.product(combinations[0],combinations[1]))
              combinations =  [[ele1+"+"+ele2 for ele1,ele2 in combined]] + combinations[2:]
        combinations = combinations[0]
        # now print
        tlines_i = []
        for combination in combinations:
            V0_comb = 0.0
            V1_comb = 0.0
            tot_weight = 1
            for species in combination.split("+"):
                ctc, itc = PN.name2data(species)
                V0, V1 = chemreac.return_V0V1(ctc,itc)
                V0_comb += V0
                V1_comb += V1
                tot_weight *= dweights[(ctc,itc)]
            V0_comb -= chemreac._V0R
            V1_comb -= chemreac._V1R
            row  = (combination,V0_comb*KCALMOL,V1_comb*KCALMOL,tot_weight)
            tlines_i.append( (V0_comb,len(combination),row) )
        tlines_i.sort()
        tlines += [(ll,row) for (V0,ll,row) in tlines_i]
        # add None for printing
        tlines.append( (0,None) )

    # add to string
    ml    = max([ll for ll,line in tlines])
    col1  = ("%%-%is"%ml)%"SP"
    thead = " %s | V0(i)-V0 | V1(i)-V1 | weight "%col1
    tdivi = "-"*len(thead)
    string += "   "+tdivi+"\n"
    string += "   "+thead+"\n"
    string += "   "+tdivi+"\n"
    for ll,row in tlines:
        if row is None: string +=  "   "+tdivi+"\n"
        else          : string += ("   "+" %%-%is | %%8.2f | %%8.2f | %%2i "%ml)%row+"\n"
    string += "   "+"SP: stationary point\n"
    return string
#---------------------------------------------------------------#
#def srcons_conformers(data,ml=4):
#    # head of table and division line
#    head = (" %%%is | conformer | weight"%ml)%"Name"
#    division = "-"*len(head)
#    # listing conformers and weights
#    string  = "Conformational flexibility:\n"
#    string += "\n"
#    string += "   "+division+"\n"
#    string += "   "+head+"\n"
#    string += "   "+division+"\n"
#    for target,itcs in data:
#        for idx,(itc,weight) in enumerate(itcs):
#            if idx == 0: string += "   "+(" %%%is "%ml)%target+ "|    %3s    |  %2i  \n"%(itc,weight)
#            else       : string += "   "+(" %%%is "%ml)%""    + "|    %3s    |  %2i  \n"%(itc,weight)
#        string += "   "+division+"\n"
#    return string
#---------------------------------------------------------------#
def srcons_anhtable(chemreac):
    iblank = "    "
    if len(chemreac._anhctcs) == 0: return ""

    string  = "-------------------------\n"
    string += " TORSIONAL ANHARMONICITY \n"
    string += "-------------------------\n"
    string += "\n"
    string += "    Anharmonic factors found for:\n"
    for ctc in chemreac._anhctcs: string += iblank+"* %s\n"%ctc
    string += "\n"
    string += "    Anharmonic corrections will be calculated for:\n"
    string += "    * for equilibrium constant    : ANHC(Keq)\n"
    string += "    * for rate constant (forward ): ANHC(k,fw)\n"
    string += "    * for rate constant (backward): ANHC(k,bw)\n"
    string += "\n"
    thead = "  T  (K)  |  ANHC(Keq)  |  ANHC(k,fw)  |  ANHC(k,bw)  "
    tdivi = "-"*len(thead)
    string += iblank+tdivi+"\n"
    string += iblank+thead+"\n"
    string += iblank+tdivi+"\n"
    for idx,T in enumerate(chemreac._ltemp):
        col1 = "     -    "
        col2 = "     -    "
        col3 = "     -    "
        if chemreac._nR*chemreac._nP != 0: col1 = "%11.3E"%chemreac._ANHKeq[idx]
        if chemreac._ts is not None:
           if chemreac._nR != 0: col2 = "%11.3E"%chemreac._ANHkfw[idx]
           if chemreac._nP != 0: col3 = "%11.3E"%chemreac._ANHkbw[idx]
        tline = " %8.2f | %s | %s | %s "%(T,col1,col2,col3)
        string += iblank+tline+"\n"
    string += iblank+tdivi+"\n"
    return string
#---------------------------------------------------------------#
def srcons_keq(chemreac):
    if chemreac._nR * chemreac._nP == 0: return ""
    dn = chemreac._nR-chemreac._nP

    # calculate Gibbs
    string  = "----------------------\n"
    string += " EQUILIBRIUM CONSTANT \n"
    string += "----------------------\n"
    string += "\n"
    # select equilibrium constant
    KEQ = list(chemreac._Keq)
    # begin table
    string += "  - Keq  : the equilibrium constant\n"
    string += "  - GFER : the Gibbs free energy of reaction (kcal/mol)\n"
    string += "  - R2P  : from reactant(s) to product(s)\n"
    string += "  - P2R  : from product(s) to reactant(s)\n"
    string += "\n"
    string += "    Keq(P2R)  = 1/Keq(R2P)\n"
    string += "    GFER(P2R) = - GFER(R2P)\n"
    string += "\n"
    if len(chemreac._anhctcs) != 0:
       string += "  - The values listed below include torsional anharmonicity\n"
       string += "    (see TORSIONAL ANHARMONICITY section)\n"
    else:
       string += "  - The values listed below do not include torsional anharmonicity\n"
    if chemreac._wfw*chemreac._wbw != 1:
       string += "  - Constants are corrected by: %i/%i\n"%(chemreac._wfw,chemreac._wbw)
    string += "\n"

    table_1cm3  = "   --------------------------------------------------\n"
    table_1cm3 += "            |       for V=1cm^3 per molecule         \n"
    table_1cm3 += "   --------------------------------------------------\n"
    table_1cm3 += "     T (K)  |  Keq (R2P)  |  Keq (P2R)  | GFER (R2P) \n"
    table_1cm3 += "   --------------------------------------------------\n"

    table_1bar  = "\n"
    table_1bar += "   --------------------------------------------------\n"
    table_1bar += "            |  for V=kB*T/p0 per molecule, p0=1bar   \n"
    table_1bar += "   --------------------------------------------------\n"
    table_1bar += "     T (K)  |  Keq (R2P)  |  Keq (P2R)  | GFER (R2P) \n"
    table_1bar += "   --------------------------------------------------\n"

    for idx,T in enumerate(chemreac._ltemp):
        from1cm3to1bar = (KB*T/VOL0/PRE0)**(chemreac._nP-chemreac._nR)
        keq_1cm3 = KEQ[idx]
        keq_1bar = keq_1cm3 * from1cm3to1bar
        gibbs_1cm3 = partfns.Kc2GFE([T],[keq_1cm3])[0] * KCALMOL
        gibbs_1bar = partfns.Kc2GFE([T],[keq_1bar])[0] * KCALMOL
        # string format
        table_1cm3 += "    %7.2f | %+11.3E | %+11.3E | %10.3f \n"%(T,keq_1cm3,1./keq_1cm3,gibbs_1cm3)
        table_1bar += "    %7.2f | %+11.3E | %+11.3E | %10.3f \n"%(T,keq_1bar,1./keq_1bar,gibbs_1bar)
    table_1cm3 += "   --------------------------------------------------\n"
    table_1bar += "   --------------------------------------------------\n"
    string += table_1cm3
    string += table_1bar
    string += "\n"
    return string
#---------------------------------------------------------------#
def srcons_indivtcoefs(chemreac):
    if chemreac._ts is None: return ""
    ANY = False
    # data to print
    ltemp       = chemreac._ltemp
    ts          = chemreac._ts
    ctc,itcs,ms = chemreac._itcs[chemreac._ts]
    oneconf = (len(itcs) == 1)

    string  = "---------------------------------\n"
    string += "CVT/SCT TRANSMISSION COEFFICIENTS\n"
    string += "---------------------------------\n"
    string += "\n"
    string += "    * Gamma_CVT = k^CVT / k^TST\n"
    string += "    * kappa_CAG/CVT\n"
    string += "    * kappa_SCT\n"
    string += "    * kappa_CVT/SCT = kappa_CAG/CVT * kappa_SCT \n"
    string += "    * gamma_CVT/SCT = Gamma_CVT * kappa_CVT/SCT\n"
    string += "\n"

    table  = ""
    thead1 = " T  (K) ,   Gamma   ,   kappa   ,   kappa   ,   kappa   ,   gamma   ".split(",")
    thead2 = "        ,    CVT    ,  CAG/CVT  ,    SCT    ,  CVT/SCT  ,  CVT/SCT  ".split(",")
    if not oneconf: thead1 = thead1[0:1]+["Conf"]+thead1[1:]
    if not oneconf: thead2 = thead2[0:1]+["    "]+thead2[1:]
    thead1 = " "+" | ".join(thead1)+" "
    thead2 = " "+" | ".join(thead2)+" "
    tdivi = "-"*len(thead1)
    table += tdivi+"\n"
    table += thead1+"\n"
    table += thead2+"\n"
    table += tdivi+"\n"
    for idxT,T in enumerate(chemreac._ltemp):
        # for each itc
        for idxITC,(itc,weights) in enumerate(itcs):
            tsname = PN.struckey(ctc,itc)
            # fill columns
            if idxITC == 0: tline = ["%8.2f"%T]
            else          : tline = [" "*8]
            if not oneconf: tline += ["%4s"%itc]
            # Transmission coefficients
            cvt    = chemreac._dall.get("cvt"   ,{}).get(tsname,None)
            sct    = chemreac._dall.get("sct"   ,{}).get(tsname,None)
            cagcvt = chemreac._dall.get("cagcvt",{}).get(tsname,None)
            if None in (cvt,sct,cagcvt): return ""
            cvtsct = cagcvt[idxT] * sct[idxT]
            total = cvt[idxT] * cvtsct
            tline += ["%11.3E"%cvt[idxT]]
            tline += ["%11.3E"%cagcvt[idxT]]

            tline += ["%11s"%fncs.eformat(sct[idxT],3)]
            tline += ["%11s"%fncs.eformat(cvtsct,3)]
            tline += ["%11s"%fncs.eformat(total,3)]
           #tline += ["%11.3E"%sct[idxT]]
           #tline += ["%11.3E"%cvtsct]
           #tline += ["%11.3E"%total]

            table += " "+" | ".join(tline)+" \n"
        # division when changing temperature
        if not oneconf: table += tdivi+"\n"
    if oneconf: table += tdivi+"\n"
    for line in table.split("\n"): string += "    "+line+"\n"
    return string
#---------------------------------------------------------------#
def srcons_avertcoefs(chemreac):
    if chemreac._ts is None: return ""
    ANY = False

    string  = "-------------------------------\n"
    string += "TOTAL TRANSMISSION COEFFICIENTS\n"
    string += "-------------------------------\n"
    string += "\n"
    string += "  The averaged transmission coefficient\n"
    string += "  for a given method (X) is:\n"
    string += "\n"
    string += "       <gamma>^X = k^X / k^TST   \n"
    string += "\n"
    string += "  where\n"
    string += "\n"
    string += "       k^TST : rate constant calculated with MS-TST\n"
    string += "       k^X   : rate constant calculated with method X\n"
    string += "\n"
    string += "  It can be also expressed as: \n"
    string += "\n"
    string += "       <gamma>^X = \sum_j chi_j^TST gamma_j^X  \n"
    string += "\n"
    string += "  with\n"
    string += "\n"
    string += "       chi_j^TST = w_j * (Q^{RR-HO}_j / Q^{MS-HO}) * exp(-U_j/kB/T)\n"
    string += "\n"
    string += "  where\n"
    string += "\n"
    string += "    gamma_j^X  : the transmission coefficient associated\n"
    string += "                 to the j-th transition state conformer\n"
    string += "    chi_j^TST  : the contribution of the j-th conformer\n"
    string += "                 to the MS-TST rate constant\n"
    string += "    w_j        : weight of j-th conformer (1 or 2)\n"
    string += "    Q^{RR-HO}_j: rigid-rotor harmonic-oscillator partition function\n"
    string += "    Q^{MS-HO}  : multi-structural harmonic-oscillator partition function\n"
    string += "    U_j        : relative energy with regard to the most stable conformer\n"
    string += "                 (considering the ZPE)\n"
    string += "\n"

    # (a) INDIVIDUAL + AVERAGED
    keys  = "tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")
    cols  = "TST/ZCT,TST/SCT,CVT,CVT/ZCT,CVT/SCT".split(",")
    stringA, ANYA = table_transcoeffs(chemreac,keys,cols)

    if ANYA:
       for line in stringA.split("\n"): string += "    "+line+"\n"

    if ANYA: return string
    else   : return ""
#---------------------------------------------------------------#
def table_transcoeffs(chemreac,keys,cols):
    # data to print
    ltemp       = chemreac._ltemp
    ts          = chemreac._ts
    ctc,itcs,ms = chemreac._itcs[chemreac._ts]
    oneconf = (len(itcs) == 1)
    table = ""
    thead1 = [" T  (K) "]+[fncs.fill_string("gamma",11) for col in cols]
    thead2 = ["        "]+[fncs.fill_string(col    ,11) for col in cols]
    if not oneconf: thead1 = thead1[0:1]+["Conf"]+thead1[1:]
    if not oneconf: thead2 = thead2[0:1]+["    "]+thead2[1:]
    thead1 = " "+" | ".join(thead1)+" "
    thead2 = " "+" | ".join(thead2)+" "
    tdivi = "-"*len(thead1)
    table += tdivi+"\n"
    table += thead1+"\n"
    table += thead2+"\n"
    table += tdivi+"\n"
    ANY = False
    for idxT,T in enumerate(chemreac._ltemp):
        # averaged
        tline  = ["%8.2f"%T]
        if not oneconf: tline += ["all "]
        for X in keys:
            values = chemreac._dtcoef[X]["averaged"]
            if values is None: value = "     -     "
            else:
                value = "%11s"%fncs.eformat(values[idxT],3)
               #value = "%11.3E"%values[idxT]
                ANY = True
            tline.append(value)
        table += " "+" | ".join(tline)+" \n"
        if oneconf: continue
        # for each itc
        for idxITC,(itc,weights) in enumerate(itcs):
            tline  = [" "*8]
            tline += ["%-4s"%itc]
            for X in keys:
                values = chemreac._dtcoef[X][itc]
                if values is None: value = "     -     "
                else:
                   value = "%11s"%fncs.eformat(values[idxT],3)
                  #value = "%11.3E"%values[idxT]
                   ANY   = True
                tline.append(value)
            table += " "+" | ".join(tline)+" \n"
        # division when changing temperature
        table += tdivi+"\n"
    if oneconf: table += tdivi+"\n"
    return table, ANY
#---------------------------------------------------------------#
def srcons_tscontributions(chemreac):
    if chemreac._ts is None: return ""
    ltemp       = chemreac._ltemp
    ts          = chemreac._ts
    ctc,itcs,ms = chemreac._itcs[chemreac._ts]
    if not ms: return ""
    if len(itcs) == 1: return ""

    keys  = "tst,tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")
    cols  = "TST,TST/ZCT,TST/SCT,CVT,CVT/ZCT,CVT/SCT".split(",")
    cols  = {k:c for k,c in zip(keys,cols)}
    keys  = clean_keys_chis(keys,chemreac._tschi)

    string  = "----------------------------------\n"
    string += "TRANSITION STRUCTURE CONTRIBUTIONS\n"
    string += "----------------------------------\n"
    string += "\n"
    if keys == ["tst"]:
       string += "  The contribution of the j-th transition state conformer\n"
       string += "  to the MS-TST rate constant is calculated as:\n"
       string += "\n"
       string += "     chi_j^TST = w_j * (Q^{RR-HO}_j / Q^{MS-HO}) * exp(-U_j/kB/T)\n"
       string += "\n"
       string += "  where\n"
       string += "\n"
       string += "    w_j        : weight of j-th conformer (1 or 2)\n"
       string += "    Q^{RR-HO}_j: rigid-rotor harmonic-oscillator partition function\n"
       string += "    Q^{MS-HO}  : multi-structural harmonic-oscillator partition function\n"
       string += "    U_j        : relative energy with regard to the most stable conformer\n"
       string += "                 (considering the ZPE)\n"
       string += "\n"
    else:
       string += "  The contribution of the j-th transition state conformer\n"
       string += "  to the X rate constant (chi_j^X) is calculated as:\n"
       string += "\n"
       string += "    chi_j^X = (gamma_j^X/<gamma>^X) * chi_j^TST\n"
       string += "\n"
       string += "  where:\n"
       string += "\n"
       string += "    gamma_j^X  : transmission coefficient of the j-th\n"
       string += "                 conformer for method X\n"
       string += "    <gamma>^X  : averaged transmission coefficient\n"
       string += "                 for method X\n"
       string += "    chi_j^TST  : contribution of the j-th conformer to\n"
       string += "                 the TST rate constant\n"
       string += "\n"


    table = ""
    thead1 = [" T  (K) ","Conf"]+[fncs.fill_string("chi_j",7) for key in keys]
    thead2 = ["        ","    "]+[fncs.fill_string(cols[key],7) for key in keys]
    thead1 = " "+" | ".join(thead1)+" "
    thead2 = " "+" | ".join(thead2)+" "
    tdivi = "-"*len(thead1)
    table += tdivi+"\n"
    table += thead1+"\n"
    table += thead2+"\n"
    table += tdivi+"\n"
    for idxT,T in enumerate(chemreac._ltemp):
        # for each itc
        for idxITC,(itc,weights) in enumerate(itcs):
            if idxITC == 0: tline = ["%8.2f"%T,"%4s"%itc]
            else          : tline = [" "*8    ,"%4s"%itc]
            for X in keys:
                values = chemreac._tschi[X][itc]
                if values is None: value = "   -   "
                else             : value = "%7.5f"%values[idxT]
                tline.append(value)
            table += " "+" | ".join(tline)+" \n"
        # division when changing temperature
        table += tdivi+"\n"
    for line in table.split("\n"): string += "    "+line+"\n"
    return string
#---------------------------------------------------------------#
def srcons_rateconstants(chemreac,direction="fw"):
    if chemreac._ts is None: return ""

    ctc,itcs,ms = chemreac._itcs[chemreac._ts]

    if direction == "fw": nR,weight,tt = chemreac._nR, chemreac._wfw, "FORWARD "
    if direction == "bw": nR,weight,tt = chemreac._nP, chemreac._wbw, "BACKWARD"

    if nR == 0: return ""

    string  = "-----------------------\n"
    string += "%s RATE CONSTANTS\n"%tt
    string += "-----------------------\n"
    string += "\n"

    # start table
    if   nR == 1: string += "    - units: sec^-1\n"
    elif nR == 2: string += "    - units: cm^3/molecule/s\n"
    else        : string += "    - units: (cm^3)^%i/molecule^%i/s\n"%(nR-1,nR-1)
    if len(chemreac._anhctcs) != 0:
       string += "    - the rate constants listed below include torsional anharmonicity\n"
       string += "      (see TORSIONAL ANHARMONICITY section)\n"
    else:
       string += "    - the rate constants listed below do not include torsional anharmonicity\n"
    string += "\n"
    if weight != 1.0:
       string += "    - Weight for the rate constant is applied: %.3f\n"%weight
       string += "\n"
    # tables with multi-structural & multi-path
    tables = ""
    if direction == "fw": dks = chemreac._kfw
    if direction == "bw": dks = chemreac._kbw
    tables += "\n"
    keys1 = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct".split(",")
    keys2 = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")
    keys1 = clean_keys_rcons(keys1,dks)
    keys2 = clean_keys_rcons(keys2,dks)
    # (a) table with multi-structural
    if ms and keys1 != keys2: tables += table_rcons(chemreac,dks,keys1,nR)
    # (b) table with multi-path
    tables += table_rcons(chemreac,dks,keys2,nR)
    for line in tables.split("\n"): string += "    "+line+"\n"
    return string
#----------------------------------------------------------#
def clean_keys_chis(keys,dtschi):
    # determine which constants will be plotted
    KEYS = []
    for X in keys:
        save = False
        if dtschi.get(X,None) in (None,{}): continue
        for itc,chis in dtschi[X].items():
            if chis is not None: save = True
        if save: KEYS.append(X)
    return KEYS
#----------------------------------------------------------#
def clean_keys_rcons(keys,dks):
    # determine which constants will be plotted
    KEYS = []
    for X in keys:
        if dks.get(X,None) is None: continue
        KEYS.append(X)
    return KEYS
#----------------------------------------------------------#
def table_rcons(chemreac,dks,keys,nR):
    ctc,itcs,ms = chemreac._itcs[chemreac._ts]
    human_units = ML**(nR-1.0) / SECOND
    oneconf = (len(itcs) == 1)

    thead = ["  T (K)  "]+["%11s"%KEYNICE[X] for X in keys]
    if not oneconf: thead = thead[0:1]+["  TS  "]+thead[1:]
    thead = " | ".join(thead)+" "
    division   = len(thead)*"-"
    string  = division  +"\n"
    string += thead+"\n"
    string += division  +"\n"

    for idx,T in enumerate(chemreac._ltemp):
        values = ["%9.2f"%T]
        if not oneconf: values += ["total "]
        # (a) THE TOTAL ONE
        for X in keys:
            try   : values += ["%11s"%fncs.eformat(dks[X][idx]*human_units,3)]
           #try   : values += ["%11.3E"%(dks[X][idx]*human_units)]
            except: values += ["     -     "]
        string += " | ".join(values)+"\n"
        # (b) MORE THAN ONE CONFORMER
        if oneconf: continue
        for itc,weights in itcs:
            values = [" "*9,"%-6s"%itc]
            for X in keys:
                if   X.startswith("ms" ): key2 = "tst"
                elif X.startswith("tst"): key2 = "tst"
                else                    : key2 = X[2:]
                try:
                   chi_i = chemreac._tschi[key2][itc][idx]
                   k_i = dks[X][idx] * chi_i
                   values += ["%11s"%fncs.eformat(k_i*human_units,3)]
                  #values += ["%11.3E"%(k_i*human_units)]
                except: values += ["     -     "]
            string += " | ".join(values)+"\n"
        string += division+"\n"
    if oneconf: string += division+"\n"
    string += "\n"
    return string
#----------------------------------------------------------#
def srcons_actgibbs(chemreac,direction="fw"):
    if chemreac._ts is None: return

    ctc,itcs,ms = chemreac._itcs[chemreac._ts]

    if direction == "fw": nR,weight,tt = chemreac._nR, chemreac._wfw, "FORWARD "
    if direction == "bw": nR,weight,tt = chemreac._nP, chemreac._wbw, "BACKWARD"

    if nR == 0: return ""

    string  = "-------------------------------------------------\n"
    string += "%s GIBBS FREE ENERGIES OF ACTIVATION (GFEA)\n"%tt
    string += "-------------------------------------------------\n"
    string += "\n"

    # table with multi-structural & multi-path
    if direction == "fw": dks = chemreac._kfw
    if direction == "bw": dks = chemreac._kbw
    keys1 = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct".split(",")
    keys2 = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")
    keys1 = clean_keys_rcons(keys1,dks)
    keys2 = clean_keys_rcons(keys2,dks)
    # (a) table with multi-structural
    if 1-nR == 0: tables_1cm3  = "    * GFEA = -R T ln(h k / kB T)\n"
    else        : tables_1cm3  = "    * GFEA = -R T ln(h k V^%i / kB T)\n"%(1-nR)
    tables_1cm3 += "    - units: kcal/mol\n"
    tables_1cm3 += "    - reference volume: 1 cm^3 per molecule\n"
    tables_1cm3 += "\n"
    if 1-nR == 0: tables_1bar  = "    * GFEA = -R T ln(h k / kB T)\n"
    else        : tables_1bar  = "    * GFEA = -R T ln(h k V^%i / kB T)\n"%(1-nR)
    tables_1bar += "    - units: kcal/mol\n"
    tables_1bar += "    - reference volume: kB*T/p0 per molecule, with p0 = 1 bar\n"
    if len(chemreac._anhctcs) != 0:
       tables_1cm3 += "    - the free energies listed below include torsional anharmonicity\n"
       tables_1bar += "    - the free energies listed below include torsional anharmonicity\n"
    else:
       tables_1cm3 += "    - the free energies listed below do not include torsional anharmonicity\n"
       tables_1bar += "    - the free energies listed below do not include torsional anharmonicity\n"
    tables_1bar += "\n"
    if ms and keys1 != keys2:
       table_1cm3, table_1bar = table_gibbs(chemreac,dks,keys1,nR)
       tables_1cm3 += table_1cm3
       tables_1bar += table_1bar
    # (b) table with multi-path
    table_1cm3, table_1bar = table_gibbs(chemreac,dks,keys2,nR)
    tables_1cm3 += table_1cm3
    tables_1bar += table_1bar
    # start table
    for line in tables_1cm3.split("\n"): string += "    "+line+"\n"
    for line in tables_1bar.split("\n"): string += "    "+line+"\n"

    return string
#----------------------------------------------------------#
def table_gibbs(chemreac,dks,keys,nR):
    table_head = " | ".join( ["  T (K)  "]+["%11s"%KEYNICE[X] for X in keys] )
    division   = len(table_head)*"-"
    string_1cm3  = division  +"\n"
    string_1cm3 += table_head+"\n"
    string_1cm3 += division  +"\n"
    string_1bar  = str(string_1cm3)
    for idx,T in enumerate(chemreac._ltemp):
        values_1cm3 = ["%9.2f"%T]
        values_1bar = ["%9.2f"%T]
        for X in keys:
            try:
               rate = dks[X][idx]
               gfe_1cm3  = partfns.rate2GFE([T],[rate],nR)[0]
               gfe_1bar  = gfe_1cm3 + KB*T *(1-nR) * np.log(VOL0*PRE0/KB/T)
               gfe_1cm3 *= KCALMOL
               gfe_1bar *= KCALMOL
               if abs(gfe_1cm3) > 0.01: values_1cm3 += ["%11.3f"%(gfe_1cm3)]
               else                   : values_1cm3 += ["%11.3E"%(gfe_1cm3)]
               if abs(gfe_1bar) > 0.01: values_1bar += ["%11.3f"%(gfe_1bar)]
               else                   : values_1bar += ["%11.3E"%(gfe_1bar)]
            except:
               values_1cm3 += ["     -     "]
               values_1bar += ["     -     "]
        string_1cm3 += " | ".join(values_1cm3)+"\n"
        string_1bar += " | ".join(values_1bar)+"\n"
    string_1cm3 += division+"\n"
    string_1cm3 += "\n"
    string_1bar += division+"\n"
    string_1bar += "\n"
    return string_1cm3,string_1bar
#----------------------------------------------------------#
def srcons_rconsX2(chemreac):
    string  = "--------------------------\n"
    string += "|  Reactants = Products  |\n"
    string += "--------------------------\n"
    string += "\n"
    string += "  * Reactants and products are identical.  \n"
    string += "  * Rate constants must be multiplied by 2.\n"
    string += "\n"
    ctc, itcs, ms = chemreac._itcs[chemreac._ts]
    rcons = {X:2*k for X,k in chemreac._kfw.items() if k is not None}
    keys1 = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct".split(",")
    keys2 = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")
    keys1 = clean_keys_rcons(keys1,rcons)
    keys2 = clean_keys_rcons(keys2,rcons)
    if ms and keys1 != keys2:
       lines = table_rcons(chemreac,rcons,keys1,chemreac._nR)
       for line in lines.split("\n"): string += "   "+line+"\n"
    lines = table_rcons(chemreac,rcons,keys2,chemreac._nR)
    for line in lines.split("\n"): string += "   "+line+"\n"
    return string
#===============================================================#


#===============================================================#
#                                                               #
#===============================================================#
def sfit_rconsX(ltemp,rcons,lX,Rs,Ps):
    #print_string(PS.sfit_rconsX(ltemp,dall["rcons"],Xbw,Ps,Rs),6)
    nR = len(Rs)
    if True in [X.split(".")[2] == "fw" for X in lX]: direction = "FORWARD"
    else                                            : direction = "BACKWARD"
    # print table of forward rate constants
    string = "    * Chemical equation (%s): '%s --> %s'\n"%(direction," + ".join(Rs)," + ".join(Ps))
    human_units = ML**(nR-1.0) / SECOND
    if   nR == 1: sunits = "in sec^-1"
    elif nR == 2: sunits = "in cm^3 / molecule / s"
    else        : sunits = "in (cm^3)^%i / molecule^%i / s"%(nR-1,nR-1)
    string += "\n"
    string += "      units rate constant: %s\n"%sunits
    keys1   = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")
    keys2   = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct".split(",")
    string += table_rcons(ltemp,rcons,keys1,nR)
    try: string += table_rcons(ltemp,rcons,keys2,nR)
    except: pass
    return string
#----------------------------------------------------------#
def sfit_rcons(ltemp,rcons,nR):
    human_units = ML**(nR-1.0) / SECOND
    # get name of correction coefficients
    rctypes = list(rcons.keys())
    sorted1 = [rctype for rctype in TRANSCOEFFS if rctype     in rctypes    ]
    sorted2 = [rctype for rctype in rctypes     if rctype not in TRANSCOEFFS]
    rctypes = sorted1 +  sorted(sorted2)
    
    string = "Rate Constants "
    if nR == 1  : string += "in sec^-1\n"
    else        : string += "in (cm^3)^%i/(molecule^%i sec)\n"%(nR-1,nR-1)
    # table head and table division
    keys    = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct".split(",")
    string += table_rcons(ltemp,rcons,keys,nR)
    keys    = "tst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct".split(",")
    string += table_rcons(ltemp,rcons,keys,nR)
    return string
#----------------------------------------------------------#
def sfit_reaction(rcname):
    # Generate strings
    title = "Reaction name: %s"%rcname
    string  = "-"*len(title)+"\n"
    string +=         title +"\n"
    string += "-"*len(title)+"\n"
    string += "\n"
    return string
#----------------------------------------------------------#
def sfit_anafit(ltemp,dfit,key):
    X, rcname, direc = key.split(".")
    cols   = [1,2,3,4,5]
    rows   = ["A","B","n","Tr","T0"]
    # print tables
    rcstr  = fncs.fill_string(KEYNICE[X],18)
    string = ""
    string += "----------------------------------------------------------------------------------\n"
    string += "                                %-18s                                \n"%rcstr
    string += "----------------------------------------------------------------------------------\n"
    string += "                                FITTING PARAMETERS                                \n"
    string += "----------------------------------------------------------------------------------\n"
    string += " Parameters |     (1)     |     (2)     |     (3)     |     (4)     |     (5)     \n"
    string += "----------------------------------------------------------------------------------\n"
    # add parameters
    for idx,row in enumerate(rows):
        string += "   %-6s   |"%row
        for col in cols:
            params,r2 = dfit.get(col,(None,None))
            try:
                if col in [4,5]: params = params[0:3]+[params[4],params[3]]
                param = "%11.3E"%params[idx]
            except:
                param = " "*11
            string += " %s "%param + "|"
        string = string[:-1]+"\n"
    string += "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
    # add 1-r^2
    string += "   1-r^2    |"
    for col in cols:
        params,r2 = dfit.get(col,(None,None))
        try   : r2 = "%7.2E"%(1-r2)
        except: r2 = " "*7
        string += " %11s |"%r2
    string = string[:-1]+"\n"
    # activation energies
    string += "----------------------------------------------------------------------------------\n"
    string += "                          ACTIVATION ENERGIES (KCAL/MOL)                          \n"
    string += "----------------------------------------------------------------------------------\n"
    string += "    T(K)    |     (1)     |     (2)     |     (3)     |     (4)     |     (5)     \n"
    string += "----------------------------------------------------------------------------------\n"
    # fitting parameters
    params_1 = dfit.get(1,[None])[0]
    params_2 = dfit.get(2,[None])[0]
    params_3 = dfit.get(3,[None])[0]
    params_4 = dfit.get(4,[None])[0]
    params_5 = dfit.get(5,[None])[0]
    for idxT,T in enumerate(ltemp):
        Ea_1 = "     -     "
        Ea_2 = "     -     "
        Ea_3 = "     -     "
        Ea_4 = "     -     "
        Ea_5 = "     -     "
        if params_1 is not None:
           A, B = params_1
           Ea_1 = KCALMOL*activation1(T,A,B)
        if params_2 is not None:
           A, B, n = params_2
           Ea_2 = KCALMOL*activation2(T,A,B,n)
        if params_3 is not None:
           A, B, n, Tr = params_3
           Ea_3 = KCALMOL*activation3(T,A,B,n,Tr)
        if params_4 is not None:
           A, B, n, T0, Tr = params_4
           Ea_4 = KCALMOL*activation4(T,A,B,n,T0,Tr)
        if params_5 is not None:
           A, B, n, T0, Tr = params_5
           Ea_5 = KCALMOL*activation5(T,A,B,n,T0,Tr)
        # correct format of value
        args = [T]
        for Ea_i in [Ea_1,Ea_2,Ea_3,Ea_4,Ea_5]:
            if type(Ea_i) == str: args.append( Ea_i )
            elif Ea_i < 0.010   : args.append( "%11.3E"%Ea_i )
            else                : args.append( "%11.3f"%Ea_i )
        string += " %10.2f | %s | %s | %s | %s | %s \n"%tuple(args)
    string += "----------------------------------------------------------------------------------\n"
    return string
#----------------------------------------------------------#
def sfit_fitting(rcname,fitting,ltemp):
    string = ""
    for direc in ["fw","bw"]:
        if len(fitting[direc].keys()) == 0: continue
        string += "- - - - - - - - -\n"
        if direc == "fw": string += "FORWARD  REACTION\n"
        if direc == "bw": string += "BACKWARD REACTION\n"
        string += "- - - - - - - - -\n"
        string += "\n"
        # methods
        for X in KEYS_X:
            if X not in fitting[direc].keys(): continue
            key, dfit, ks_human = fitting[direc][X]
            if dfit is None: continue
            string += sfit_anafit(ltemp,dfit,key)
            string += "\n"
        string += "\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --kies                         #
#===============================================================#
def skies_summary(dchem,dall):

    sinit = "   "
    string  = " Available reactions \n"
    string += "---------------------\n"
    string += "\n"
    string += sinit+"Nmethods: the number of methods the rate constant was calculated with\n"
    string += "\n"

    # table 
    ml = max([len(rcname) for rcname in dchem.keys()]+[8])
    ml2 = ml+4
    row_format = " %%%is | %%3s | %%%is | %%5s "%(ml,ml2)

    titles =["reaction","dir",fncs.fill_string("keyword",ml2),"Nmethods"]
    head = row_format%tuple(titles)
    divi = "-"*len(head)

    string += sinit+divi+"\n"
    string += sinit+head+"\n"
    string += sinit+divi+"\n"
    for rcname in dchem.keys():
        for direc in "fw,bw".split(","):
            nmeths = 0
            for X in KEYS_X:
                key1 = "%s.%s.%s"%(X,rcname,direc)
                if dall["rcons"].get(key1,None) is not None: nmeths += 1
            key_to_intro = fncs.fill_string(rcname+"."+direc,ml2)
            row = [fncs.fill_string(rcname,ml),direc,key_to_intro,nmeths]
            string += sinit+row_format%tuple(row)+"\n"
        string += sinit+divi+"\n"
    string += "\n"
    string += sinit+"Reactions are selected as in the 'keyword' column ($reaction.$dir)\n"
    string += "\n"

    return string
#---------------------------------------------------------------#
def skies_definitions():
    string  = " Calculation of Kinetic Isotopic Effects (KIEs) \n"
    string += "------------------------------------------------\n"
    string += "\n"
    string += "  Total KIE (kie_tot) is given by:\n"
    string += "\n"
    string += "     kie_tot = kie_tr * kie_rv * kie_vtun * kie_tor \n"
    string += "\n"
    string += "  where its contributions are:\n"
    string += "     * translational            --> kie_tr\n"
    string += "     * rovibrational            --> kie_rv\n"
    string += "     * variational & tunnelling --> kie_vtun\n"
    string += "     * anharmonic               --> kie_tor\n"
    string += "\n"
    string += "  With more than one transition state conformers\n"
    string += "  the following magnitudes are also calculated:\n"
    string += "\n"
    string += "    * kie_j  --> KIE calculated considering only j-th conformer\n"
    string += "    * kie#_j --> weighted individual KIE for confomer j-th\n"
    string += "    * P_j,D  --> ratio between individual and total rate constants\n"
    string += "                 of the isotopically substituted species\n"
    string += "    * P_j,H  --> ratio between the weighted individual KIE and the\n"
    string += "                 total KIE\n"
    string += "\n"
    string += "  The following relationships hold:\n"
    string += "\n"
    string += "      kie_tot = sum_j( P_j,D * kie_j   )\n"
    string += "      kie_tot = sum_j( P_j,H * kie_tot )\n"
    string += "      kie_tot = sum_j(         kie#_j  )\n"
    string += "      kie#_j  = P_j,D * kie_j\n"
    string += "      P_j,H   = kie#_j / kie_tot\n"
    string += "\n"
    string += "  Reference:\n"
    string += "       L. Simon-Carballido, T. V. Alves,\n"
    string += "       A. Dybala-Defratyka and A. Fernandez-Ramos,\n"
    string += "       J. Phys. Chem. B (2016), 120, 1911-1918\n"
    string += "\n"
    return string
#--------------------------------------------------#
def skies_basiccontris(ltemp,kie_tr,kie_rv,kie_tor):
    string  = "Contributions to the total KIE (method-independent):\n"
    string += "\n"
    string += "    ---------------------------------------------\n"
    string += "      T (K)  |  kie_tr   |  kie_rv   |  kie_tor  \n"
    string += "    ---------------------------------------------\n"
    for idx,T in enumerate(ltemp):
        line_data  = ["%7.2f"%T]
        line_data += ["%9.3f"%kie_tr[idx]]
        line_data += ["%9.3f"%kie_rv[idx]]
        line_data += ["%9.3f"%kie_tor[idx]]
        string += "     %7s | %9s | %9s | %9s \n"%tuple(line_data)
    string += "    ---------------------------------------------\n"
    return string
#---------------------------------------------------------------#
def skies_vtun_tot(ltemp,kies_vtun,kies_tot):
    string  = "Total KIE and vtun contribution (method-dependent):\n"
    string += "\n"
    string += "    ------------------------------------------------\n"
    string += "      T (K)  |   method   |  kie_vtun  |  kie_tot   \n"
    string += "    ------------------------------------------------\n"
    for idxT,T in enumerate(ltemp):
        for idxX,X in enumerate(KEYS_X):
            if X not in kies_vtun.keys(): continue
            kie_vtun    = kies_vtun[X][idxT]
            kie_tot     = kies_tot[X][idxT]
            if idxX == 0: line_data  = ["%7.2f"%T]
            else        : line_data  = ["       "]
            line_data += ["%-10s"%KEYNICE[X]]
            if   X == "tst"      : colvtun = "     -    "
            elif kie_vtun < 1000 : colvtun = "%10.3f"%kie_vtun
            else                 : colvtun = "%10.3E"%kie_vtun
            line_data += [colvtun]
            line_data += ["%10.3f"%kie_tot]
            string += "     %s | %s | %s | %s \n"%tuple(line_data)
        string += "    ------------------------------------------------\n"
    return string
#---------------------------------------------------------------#
def skies_indkies(ltemp,ikies_rv,ikies_vtun,pjh,pjd,kies_j,kiest_j):
    string = ""
    row_format = " %7s | %4s | %6s | %10s | %7s | %7s | %7s | %7s "

    for idxX,X in enumerate(KEYS_X):
        string += "Contribution of each transition structure to KIE using %s:\n"%KEYNICE[X]
        string += "\n"
        string += "    -------------------------------------------------------------------------------\n"
        string += "      T (K)  |  SP  |  rv,j  |   vtun,j   |  P_j,D  |  P_j,H  |  kie_j  |  kie#_j  \n"
        string += "    -------------------------------------------------------------------------------\n"
        for idxT,T in enumerate(ltemp):
            for idxITC,itc in enumerate(sorted(ikies_rv.keys())):
                # data for row
                if X.startswith("ms"):
                    XG = X[2:]
                    col_rv    = ikies_rv[itc][idxT]
                    col_vt    = ikies_vtun[X][itc][idxT]
                    col_pD    = pjd[X][itc][idxT]
                    col_pH    = pjh[X][itc][idxT]
                    col_kiej  = kies_j[X][itc][idxT]
                    col_kietj = kiest_j[X][itc][idxT]
                else:
                    col_rv    = ikies_rv[itc][idxT]
                    col_vt    = ikies_vtun[X][itc][idxT]
                    col_pD    = pjd[X][itc][idxT]
                    col_pH    = pjh[X][itc][idxT]
                    col_kiej  = kies_j[X][itc][idxT]
                    col_kietj = kiest_j[X][itc][idxT]

                if idxITC == 0: line_data  = ["%7.2f"%T]
                else          : line_data  = ["       "]
                line_data += ["%4s"%itc]
                line_data += ["%6.3f"%col_rv]
                if X == "tst"     : col_vt    = "    -     "
                elif col_vt < 1000: col_vt    = "%10.3f"%col_vt
                else              : col_vt    = "%10.3E"%col_vt
                line_data += [col_vt]
                line_data += ["%7.3f"%col_pD]
                line_data += ["%7.3f"%col_pH]
                line_data += ["%7.3f"%col_kiej]
                line_data += ["%7.3f"%col_kietj]
                string += "     %s | %s | %s | %s | %s | %s | %s | %s \n"%tuple(line_data)
            string += "    -------------------------------------------------------------------------------\n"
        string += "\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --kmc                          #
#===============================================================#
def skmc_init(ipops,POP0,excess,rcs,psteps,volume,timeunits,dof):
    # some lengths for nice format
    try   : l1 = max([len(key) for key       in rcs.keys()])
    except: l1 = 4
    try   : l2 = max([len(v1)  for (v1,v2)   in rcs.values()])
    except: l2 = 4
    try   : l3 = max([len(key) for key,val   in ipops.items() if val != 0.0])
    except: l3 = 4
    # generate string
    string  = ""
    string += "---------------------\n"
    string += "KMC basic information\n"
    string += "---------------------\n"
    string += "  Writing step set to    : %s\n"%psteps
    string += "  Simulation volume      : %.3E mL\n"%(volume*ML)
    string += "  Time units             : %s\n"%(timeunits)
    string += "\n"
    string += "  Initial populations, pop(i;t=0):\n"
    for species,ipop in sorted(ipops.items()):
        if ipop == 0.0: continue
        species = ("%%-%is"%l3)%species
        string += "      pop(%s;t=0) = %.3E particles"%(species,ipop)
        if species in excess: string += " (excess)\n"
        else                : string += "\n"
    string += "\n"
    string += "  Population of limiting reactant (POP0): %.3E particles\n"%POP0
    string += "\n"
    string += "  Rate constants:\n"
    for rname,(ktype,weight,coefs) in sorted(rcs.items()):
        if rname.endswith(".both"): rname = rname[:-5]
        rname = ("%%-%is"%l1)%rname
        ktype = ("%%-%is"%l2)%ktype
        string += "      %s for %s"%(ktype,rname)
        if weight != 1: string += " (x%i)"%weight
        string += "\n"
    string += "---------------------\n"
    string += "\n"

    return string
#---------------------------------------------------------------#
def skmc_processes(T,processes):
    string  = ""
    string += "Processes to be considered:\n"
    ml1 = max([len(" + ".join(Rs)) for Rs,Ps,k in processes])
    ml2 = max([len(" + ".join(Ps)) for Rs,Ps,k in processes])
    for Rs,Ps,k in processes:
        nR = len(Rs)
        if nR == 1:
           k *= (1.0/SECOND)
           units = "sec^-1"
        else:
            k *= (ML**(nR-1)) / SECOND
            units = "(cm^3)^%i/(molecule^%i sec)"%(nR-1,nR-1)
        str_Rs = ("%%-%is"%ml1)%" + ".join(Rs)
        str_Ps = ("%%-%is"%ml2)%" + ".join(Ps)
        string += "    %s --> %s (%s %s)\n"%(str_Rs,str_Ps,fncs.eformat(k,3),units)
       #string += "    %s --> %s (%.3E %s)\n"%(str_Rs,str_Ps,k,units)
    return string
#---------------------------------------------------------------#
def skmc_results(xvalues,yvalues,timeunits):
    stime    = "time (%s)"%timeunits
    nrpl     = 6
    ib       = "    "

    len1 = max([len(string) for string in list(yvalues.keys())+[stime]])
    molecules = sorted(yvalues.keys())

    # xvalues to units
    if timeunits == "fs" : xvalues = [t_i*SECOND*1e15  for t_i in xvalues]
    if timeunits == "ps" : xvalues = [t_i*SECOND*1e12  for t_i in xvalues]
    if timeunits == "mcs": xvalues = [t_i*SECOND*1e6   for t_i in xvalues]
    if timeunits == "ms" : xvalues = [t_i*SECOND*1e3   for t_i in xvalues]
    if timeunits == "s"  : xvalues = [t_i*SECOND       for t_i in xvalues]
    if timeunits == "min": xvalues = [t_i*SECOND/60.   for t_i in xvalues]
    if timeunits == "hr" : xvalues = [t_i*SECOND/3600. for t_i in xvalues]
    string  = ""
    string += "Evolution of population(s) with time, pop_i(t):\n"
    string += "\n"
    for idx1 in range(0,len(xvalues),nrpl):
        idx2 = idx1+nrpl
        table_head  = (" %%-%is | "%len1)%stime
        table_head += " | ".join( ["%10s"%fncs.eformat(tt,2) for tt in xvalues[idx1:idx2]])+" "
       #table_head += " | ".join( ["%.2E"%tt for tt in xvalues[idx1:idx2]])
        string += ib+table_head+"\n"
        string += ib+"-"*len(table_head)+"\n"
        for molecule in molecules:
            pops = yvalues[molecule][idx1:idx2]
            molecule = ("%%-%is"%len1)%molecule
            string += ib+" %s | "%molecule + " | ".join( ["%10s"%fncs.eformat(pop,2) for pop in pops])+" \n"
           #string += ib+" %s | "%molecule + " | ".join( ["%.2E"%pop for pop in pops])+"\n"
        string += "\n"
    return string
#---------------------------------------------------------------#
def skmc_finaltimes(ltemp,ftimes,timeunits):
    if timeunits == "fs" : ftimes = [t_i*SECOND*1e15  for t_i in ftimes]
    if timeunits == "ps" : ftimes = [t_i*SECOND*1e12  for t_i in ftimes]
    if timeunits == "mcs": ftimes = [t_i*SECOND*1e6   for t_i in ftimes]
    if timeunits == "ms" : ftimes = [t_i*SECOND*1e3   for t_i in ftimes]
    if timeunits == "s"  : ftimes = [t_i*SECOND       for t_i in ftimes]
    if timeunits == "min": ftimes = [t_i*SECOND/60.   for t_i in ftimes]
    if timeunits == "hr" : ftimes = [t_i*SECOND/3600. for t_i in ftimes]
    string = ""
    string += "------------------------\n"
    string += " Simulation time in %3s \n"%timeunits
    string += "------------------------\n"
    string += "\n"
    string += "     ---------------------\n"
    string += "       T (K)  | sim. time \n"
    string += "     ---------------------\n"
    for T,ftime in zip(ltemp,ftimes):
        string += "      %7.2f | %9s \n"%(T,fncs.eformat(ftime,2))
       #string += "      %7.2f | %9.2E \n"%(T,ftime)
    string += "     ---------------------\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def skmc_finalratio(ltemp,species,fratios):
    NEPT    = 5 # number of elements per table
    ml      = max([len(xx) for xx in species]+[7])
    string = ""
    string += "----------------------------\n"
    string += " Final ratios (pop(i)/POP0) \n"
    string += "----------------------------\n"
    string += "\n"
    for idx in range(0,len(species),NEPT):
        head     = [" T (K) "]+[("%%%is"%ml)%xx for xx in species[idx:idx+NEPT]]
        head     = " "+" | ".join(head)+" "
        division = "-"*len(head)
        string += "     "+division+"\n"
        string += "     "+head    +"\n"
        string += "     "+division+"\n"
        for temp in ltemp:
            key1 = "%7.2f"%temp
            line = [key1]
            for species_i in species[idx:idx+NEPT]:
                value = fratios[key1][species_i]
                if value < 10: line.append( ("%%%i.3f"%ml)%value )
                else         : line.append( ("%%%i.1E"%ml)%value )
            line = " "+" | ".join(line)+" "
            string += "     "+line+"\n"
        string += "     "+division+"\n\n"
    return string
#===============================================================#

