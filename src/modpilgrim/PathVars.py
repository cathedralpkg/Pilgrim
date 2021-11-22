'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.5
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
| Sub-module :  PathVars           |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#---------------------------------------------------------------#
import os
#---------------------------------------------------------------#
import common.Exceptions as Exc
#---------------------------------------------------------------#
import modpilgrim.names      as PN
import modpilgrim.pilrw       as RW
import modpilgrim.steepdesc  as sd
#---------------------------------------------------------------#
from   common.physcons   import AMU, KCALMOL, HBAR
from   common.fncs       import afreq2cm
from   common.fncs       import afreq2zpe
from   common.fncs       import cm2afreq
from   common.fncs       import do_parallel
from   common.Molecule   import Molecule
from   common.criteria   import EPS_MEPINCR, EPS_KCALMOL
#---------------------------------------------------------------#


#===============================================================#
def var_info(key):
    '''
    description of variables, if needed
    '''
    info = {}
    info["fwdir"   ] = "ic & action     ; ic increases(++)/decreases(--) in the fw direction   "
    info["sbw"     ] = "float           ; minimum value of s in bohr                           "
    info["sfw"     ] = "float           ; maximum value of s in bohr                           "
    info["ds"      ] = "float           ; step size in bohr                                    "
    info["hsteps"  ] = "integer         ; number of steps for Hessian update                   "
    info["cubic"   ] = "yes/no/float    ; use cubic first step                                 "
    info["mtype"   ] = "pm/es           ; MEP-following algorithm; pm: Page-McIver, es: Euler) "
    info["mu"      ] = "float           ; scaled mass in amu                                   "
    info["eref"    ] = "auto/float      ; ref. energy for MEP; float value or auto (in hartree)"
    info["paral"   ] = "yes/no          ; get both sides of MEP at the same time (check manual)"
    info["epse"    ] = "float           ; conv. criterium for energy (hartree)                 "
    info["epsg"    ] = "float           ; conv. criterium for ms gradient norm (hartree/bohr)  "
    info["useics"  ] = "yes/no          ; projects Hessian to internal coordinates             "
    info["cvt"     ] = "yes/no          ; calculates CVT variational coefficient               "
    info["sct"     ] = "yes/no          ; calculates SCT transmission coefficient              "
    info["muintrpl"] = "lineal/cubic int; type of intrpl of mueff at the TS structure          "
    info["v1mode"  ] = "grad/hess       ; indicates how to calculate the Bmf terms             "
    info["scterr"  ] = "float           ; criterium for kappa_SCT variation (in percentage)    "
    info["e0"      ] = "auto/float      ; lowest E for SCT; float value or auto (in Hartree)   "
    info["onioml"  ] = "list of integers; atoms in the L layer in ONIOM (only GAUSSIAN)        "
    info["oniomm"  ] = "list of integers; atoms in the L layer in ONIOM (only GAUSSIAN)        "
    info["oniomh"  ] = "list of integers; atoms in the L layer in ONIOM (only GAUSSIAN)        "
    return info.get(key,"")
#---------------------------------------------------------------#
def list_of_atoms(string):
    latoms = []
    for atoms in string.split():
        # a single value
        if "-" not in atoms:
           latoms += [int(atoms)-1]
        # a range of values
        else:
            at1,at2 = atoms.split("-")
            latoms += range(int(at1)-1,int(at2))
    return list(set(latoms))
#---------------------------------------------------------------#
def string_of_atoms(latoms):
    string = ""
    # sort and correct to start at position 1
    latoms.sort()
    latoms = [at+1 for at in latoms]
    # generate string
    ii     = 0
    at0    = latoms[0]
    string = "%i"%at0
    for idx,at in enumerate(latoms):
        if idx == 0: continue
        if at == at0+1:
           ii   = 1
           if idx+1 == len(latoms): string += "-%s"%at
        elif ii == 0:
           string += " %s"%at
        else:
           string += "-%s %s"%(latoms[idx-1],at)
           ii = 0
        at0  = at
    return string
#---------------------------------------------------------------#
def get_reaction_energies(TS,dchem,dof):
    '''
    * checks the reactions which involved the target transition
      state (TS) in order to extract V0 and V1 for reactants
      and products 
    '''
    # initialize variables
    reaction = None
    Eref     = None
    V0R      = None
    V0P      = None
    V1R      = None
    V1P      = None
    GibbsR   = None
    # select reaction
    ctc, itc = PN.name2data(TS)
    for rname in dchem.keys():
        Rs, ts, Ps = dchem[rname]
        ctc2, itc2 = PN.name2data(ts)
        if ctc == ctc2:
           reaction = rname
           if itc == itc2: break
    # no reaction?
    if reaction is None: return reaction, V0R, V0P, V1R, V1P, GibbsR
    # read dof
    dall = RW.read_alldata(dof)[0]
    # get energy from reaction
    Rs = dchem[reaction][0]
    Ps = dchem[reaction][2]
    # reactants
    if len(Rs) != 0:
        V0R, V1R = 0.0, 0.0
        for R in Rs:
            ctc, itc = PN.name2data(R)
            if itc is None: key = PN.struckey(ctc,"msho")
            else          : key = R
            data = dall["pfn"].get(key,None)
            if data is None:
               V0R, V1R = None, None
               break
            V0, V1, pfns = data
            V0R += V0
            V1R += V1
            # Gibbs energy
            if key in dall["gibbs1cm3"].keys():
               gibbs = dall["gibbs1cm3"][key]
               if GibbsR is None: GibbsR = gibbs
               else             : GibbsR = [ gi+gj for gi,gj in zip(GibbsR,gibbs)]
    # products
    if len(Ps) != 0:
        V0P, V1P = 0.0, 0.0
        for P in Ps:
            ctc, itc = PN.name2data(P)
            if itc is None: key = PN.struckey(ctc,"msho")
            else          : key = P
            data = dall["pfn"].get(key,None)
            if data is None:
               V0P, V1P = None, None
               break
            V0, V1, pfns = data
            V0P += V0
            V1P += V1
    # Output
    return reaction, V0R, V0P, V1R, V1P, GibbsR
#---------------------------------------------------------------#
class PathVars():
      '''
      class that manages the variables of the path
      '''
      def __init__(self,pathtype="mep"):

         # type of path
         self._pathtype = pathtype

         # path variables (via input) - basic
         self._sbw     = "default"
         self._sfw     = "default"
         self._ds      = "default"
         self._hsteps  = "default"

         # path variables (via input) - advanced
         self._fwdir   = "default"
         self._mtype   = "default"
         self._mu      = "default"
         self._cubic   = "default"
         self._eref    = "default"
         self._paral   = "default"
         self._epse    = "default"
         self._epsg    = "default"

         # path variables (via input) - highly advanced
         self._cvt     = "default"
         self._sct     = "default"
         self._useics  = "default"
         self._v1mode  = "default"
         self._scterr  = "default"
         self._lowfq   = "default"
         self._e0      = "default"
         self._muintrpl= "default"
         self._qrc     = "default"
         self._onioml  = "default"
         self._oniomm  = "default"
         self._oniomh  = "default"

         # other variables of interest
         self._masses   = None
         self._V0R      = None
         self._V0P      = None
         self._dlevel   = None
         self._ics      = None
         self._icsbw    = None
         self._icsfw    = None
         self._cagtst   = "no"
         self._cagcvt   = "no"
         self._freqscal = 1.0
         self._convbw   = False
         self._convfw   = False
         self._keeptmp  = False

         # specific variables for itc
         self._spec    = {}

         self._reaction  = None
         self._reactioneq= None
         self._exorgic   = True
         self._V1R       = None
         self._V1P       = None
         self._GibbsR    = None
         self._beyondmep = True
         self._qrclE     = []
         self._qrcafreq  = None
         self._qrcauto   = True
         self._qrcname   = None
         self._qrccase   = 0

      def converged_in_bw(self,sbw):
          self._convbw = True
          self._sbw    = sbw

      def converged_in_fw(self,sfw):
          self._convfw = True
          self._sfw    = sfw

      def set_masses(self,masses):
          self._masses = masses

      def setvar(self,var,value):
          '''
          save a string at each given variable - except for highly advanced
          requires a latter setup
          '''
          value = value.lower()

          # split variable in case of specific
          var = var.split(".")
          if   len(var) == 1: var, itc = var[0], "*"
          elif len(var) == 2: var, itc = var[0], var[1]
          else              : return

          # save anything in specific?
          if itc != "*":
             self._spec[itc] = self._spec.get(itc,[]) + [(var,value)]
             return

          # path variables (via input) - basic
          if   var == "sbw"     : self._sbw     = str(-abs(float(value)))
          elif var == "sfw"     : self._sfw     = str(+abs(float(value)))
          elif var == "ds"      : self._ds      = str(+abs(float(value)))
          elif var == "hsteps"  : self._hsteps  = str(+abs(  int(value)))

          # path variables (via input) - advanced
          elif var == "bwdir" :
               if "++" in value : self._fwdir   = value.replace("++","--")
               else             : self._fwdir   = value.replace("--","++")
          elif var == "fwdir"   : self._fwdir   = value
          elif var == "mtype"   : self._mtype   = value
          elif var == "mu"      : self._mu      = value
          elif var == "cubic"   : self._cubic   = value
          elif var == "eref"    : self._eref    = value
          elif var == "paral"   : self._paral   = value
          elif var == "epse"    : self._epse    = value
          elif var == "epsg"    : self._epsg    = value
          elif var == "cvt"     : self._cvt     = value
          elif var == "sct"     : self._sct     = value
          elif var == "useics"  : self._useics  = value
          elif var == "v1mode"  : self._v1mode  = value
          elif var == "scterr"  : self._scterr  = value
          elif var == "e0"      : self._e0      = value
          elif var == "keeptmp" : self._keeptmp = True

          # path variables (via input) - highly advanced
          elif var == "muintrpl":
             try:
               mode, num = value.split()
               mode = mode.lower()
               num  = min(int(num),2)
               if mode not in ["linear","cubic"]: mode = "linear"
               self._muintrpl = (mode,num)
             except:
               self._muintrpl = "default"
          elif var == "qrc":
               nnn = len(value.split())
               mode,num,auto = "1","1000","auto"
               if nnn == 1: mode          = value.split()[0]
               if nnn == 2: mode,num      = value.split()
               if nnn == 3: mode,num,auto = value.split()
               if auto == "always": self._qrcauto = False
               else               : self._qrcauto = True
               self._qrc = (int(mode)-1,int(num))
          elif var == "onioml":
               self._onioml = list_of_atoms(value)
          elif var == "oniomm":
               self._oniomm = list_of_atoms(value)
          elif var == "oniomh":
               self._oniomh = list_of_atoms(value)
          elif var == "lowfq":
             if self._lowfq == "default": self._lowfq = {}
             keys  = self._lowfq.keys()
             value = value.split()
             if len(value) == 2:
                mode, cm = value
                direc    = "+-"
             elif len(value) == 3:
               mode, cm, direc = value
             else: return
             mode = int(mode)
             if direc != "++" and direc != "--": direc = "+-"
             # save as double dict
             if "+" in direc:
                self._lowfq["fw"] = self._lowfq.get("fw",{})
                self._lowfq["fw"][mode] = cm2afreq(float(cm))
             if "-" in direc:
                self._lowfq["bw"] = self._lowfq.get("bw",{})
                self._lowfq["bw"][mode] = cm2afreq(float(cm))

      def apply_specific(self,itc="*"):
          if itc not in self._spec.keys(): return
          for var,value in self._spec[itc]: self.setvar(var,value)

      def setup1(self):
          '''
          convert user data to value
          '''
          for var in [self._fwdir,self._sbw,self._sfw,self._ds,self._hsteps,\
                      self._cvt,self._sct,self._mtype,self._mu,self._cubic,\
                      self._eref,self._paral,self._epse,self._epsg,\
                      self._useics,self._v1mode,self._scterr,self._e0]:

              if type(var) != str: continue
              if var in ["default","auto"]: continue

              try:
                 # path variables (via input) - basic
                 if var is self._sbw   : self._sbw    = float(self._sbw)
                 if var is self._sfw   : self._sfw    = float(self._sfw)
                 if var is self._ds    : self._ds     = float(self._ds)
                 if var is self._hsteps: self._hsteps = int(self._hsteps)
                 # path variables (via input) - advanced
                 if var is self._fwdir : self._fwdir  = tuple(self._fwdir.split())
                 if var is self._cvt   : self._cvt    = str(self._cvt)
                 if var is self._sct   : self._sct    = str(self._sct)
                 if var is self._mtype : self._mtype  = str(self._mtype)
                 if var is self._mu    : self._mu     = float(self._mu)/AMU
                 if var is self._cubic : self._cubic  = str(self._cubic)
                 if var is self._eref  : self._eref   = float(self._eref)
                 if var is self._paral : self._paral  = str(self._paral)
                 if var is self._epse  : self._epse   = float(self._epse)
                 if var is self._epsg  : self._epsg   = float(self._epsg)
                 if var is self._useics: self._useics = str(self._useics)
                 if var is self._v1mode: self._v1mode = str(self._v1mode)
                 if var is self._scterr: self._scterr = float(self._scterr)
                 if var is self._e0    : self._e0     = float(self._e0 )
              except:
                 # path variables (via input) - basic
                 if var is self._sbw   : self._sbw    = "default"
                 if var is self._sfw   : self._sfw    = "default"
                 if var is self._ds    : self._ds     = "default"
                 if var is self._hsteps: self._hsteps = "default"
                 # path variables (via input) - advanced
                 if var is self._fwdir : self._fwdir  = "default"
                 if var is self._cvt   : self._cvt    = "default"
                 if var is self._sct   : self._sct    = "default"
                 if var is self._mtype : self._mtype  = "default"
                 if var is self._mu    : self._mu     = "default"
                 if var is self._cubic : self._cubic  = "default"
                 if var is self._eref  : self._eref   = "default"
                 if var is self._paral : self._paral  = "default"
                 if var is self._epse  : self._epse   = "default"
                 if var is self._epsg  : self._epsg   = "default"
                 if var is self._useics: self._useics = "default"
                 if var is self._v1mode: self._v1mode = "default"
                 if var is self._scterr: self._scterr = "default"
                 if var is self._e0    : self._e0     = "default"

      def setup2(self):
          '''
          setup of basic vars (those in string4inp)
          '''
          # path variables (via input) - basic
          if self._sbw    == "default": self._sbw    = -0.50
          if self._sfw    == "default": self._sfw    = +0.50
          if self._ds     == "default": self._ds     =  0.01
          if self._hsteps == "default": self._hsteps = 10

      def setup3(self):
          '''
          apply default value to those not defined!
          '''
          # use ics but no ics?
          if self._useics == "yes" and (self._ics is None or len(self._ics) == 0):
             raise Exc.NoICS(Exception)

          # path variables (via input) - advanced
          if self._fwdir    == "default"     : self._fwdir    = None # ("1-2","++")
          if self._cvt      == "default"     : self._cvt      = "yes"
          if self._sct      == "default"     : self._sct      = "yes"
          if self._mtype    == "default"     : self._mtype    = "pm"
          if self._mu       == "default"     : self._mu       = 1.0/AMU
          if self._cubic    == "default"     : self._cubic    = "no"
          if self._eref in ["default","auto"]: self._eref     = None
          if self._paral    == "default"     : self._paral    = "no"
          if self._epse     == "default"     : self._epse     = 1e-8
          if self._epsg     == "default"     : self._epsg     = 1e-4
          if self._useics   == "default"     : self._useics   = "yes"
          if self._v1mode   == "default"     : self._v1mode   = "grad"
          if self._scterr   == "default"     : self._scterr   = 0.0
          if self._e0 in ["default","auto"]  : self._e0       = None
          if self._lowfq    == "default"     : self._lowfq    = {}
          if self._muintrpl == "default"     : self._muintrpl = ("linear",0)
          if self._qrc      == "default"     : self._qrc      = None
          if self._onioml   == "default"     : self._onioml   = []
          if self._oniomm   == "default"     : self._oniomm   = []
          if self._oniomh   == "default"     : self._oniomh   = []

          # setup parallel
          self._paral = do_parallel(self._paral)
          
          # setup cag
          self.set_cag()

          # useics = yes but there are no ics?
          if self._useics == "yes" and (self._ics is None or len(self._ics) == 0):
             self._useics = "no"


      def string4inp(self,target):
          # setup of basic data
          self.setup1()
          self.setup2()
          # the string
          string  = "%s (%s):\n"%(target,self._pathtype)
          string += "   sbw    = %+7.4f     sfw    = %+7.4f\n"%(self._sbw,self._sfw)
          string += "   ds     = %7.5f     hsteps = %7i\n"%(self._ds,self._hsteps)
          if self._paral  != "default": string += "   paral  = %s\n"%self._paral
          if self._scterr != "default": string += "   scterr = %.3f\n"%self._scterr
          # return data
          return string

      def string4pif(self,target):
          # setup of basic data
          self.setup1()
          self.setup2()
          # string for MEP
          string  = "start_mep %s     \n"%target
          # path variables (via input) - basic
          string += "  sbw      %-+8.4f \n"%self._sbw
          string += "  sfw      %-+8.4f \n"%self._sfw
          string += "  ds       %-8.5f  \n"%self._ds
          string += "  hsteps   %i     \n"%self._hsteps

          # path variables (via input) - advanced
          if self._fwdir    not in ["default","auto"]: string += "  fwdir    %s  %s \n"%self._fwdir
          if self._cvt      not in ["default","auto"]: string += "  cvt      %-11s  \n"%self._cvt
          if self._sct      not in ["default","auto"]: string += "  sct      %-11s  \n"%self._sct
          if self._mtype    not in ["default","auto"]: string += "  mtype    %-11s\n"%self._mtype
          if self._mu       not in ["default","auto"]: string += "  mu       %.3f\n"%(self._mu*AMU)
          if self._cubic    not in ["default","auto"]: string += "  cubic    %s\n"%self._cubic
          if self._eref     not in ["default","auto"]: string += "  eref     %.6f\n"%self._eref
          if self._paral    not in ["default","auto"]: string += "  paral    %-11s\n"%self._paral
          if self._epse     not in ["default","auto"]: string += "  epse     %.2e\n"%self._epse
          if self._epsg     not in ["default","auto"]: string += "  epsg     %.2e\n"%self._epsg
          if self._useics   not in ["default","auto"]: string += "  useics   %-11s\n"%self._useics
          if self._v1mode   not in ["default","auto"]: string += "  v1mode   %-11s\n"%self._v1mode
          if self._scterr   not in ["default","auto"]: string += "  scterr   %.3f\n"%self._scterr
          if self._e0       not in ["default","auto"]: string += "  e0       %.6f\n"%self._e0

          if self._muintrpl not in ["default","auto"]: string += "  muintrpl %s  %i\n"%self._muintrpl
          if self._qrc      not in ["default","auto"]: string += "  qrc      %i  %i  %s\n"%(self._qrc[0]+1,self._qrc[1],"auto" if self._qrcauto else "always")
          if self._lowfq    not in ["default","auto"]:
             for direc in self._lowfq.keys():
                 for mode,freq in self._lowfq[direc].items():
                     if direc == "fw": string +="  lowfq    %i %.2f ++\n"%(mode,afreq2cm(freq))
                     if direc == "bw": string +="  lowfq    %i %.2f --\n"%(mode,afreq2cm(freq))

          if self._onioml   not in ["default","auto"]:
             string += "  onioml   %s\n"%string_of_atoms(self._onioml)
          if self._oniomm   not in ["default","auto"]:
             string += "  oniomm   %s\n"%string_of_atoms(self._oniomm)
          if self._oniomh   not in ["default","auto"]:
             string += "  oniomh   %s\n"%string_of_atoms(self._oniomh)
          if self._keeptmp: string += "  keeptmp      \n"
          # add everything in specific
          if self._spec != {}: string += "  # specific\n"
          for itc in self._spec.keys():
              for var,value in self._spec[itc]:
                  var = "%s.%3s"%(var,itc)
                  string += "  %-16s  %s\n"%(var,value)

          string += "end_mep\n"
          return string

      def set_cag(self):
          if self._sct == "yes":
             self._cagtst = "yes"
             if self._cvt == "yes":
                self._cagcvt = "yes"

      def set_eref_from_reaction(self,tsname,dchem,dof):
          thebool = self._eref in [None,"auto","default"]
          ctc, itc = PN.name2data(tsname)
          # Get data from reactions
          rname, V0R, V0P, V1R, V1P, GibbsR = get_reaction_energies(tsname,dchem,dof)
          self._V1R = V1R
          self._V1P = V1P
          if None not in (self._V1R,self._V1P) and self._V1P > self._V1R: self._exorgic = False
          # save GibbsR for CVT gibbs
          self._GibbsR = GibbsR
          # go case by case
          if thebool and V0R is not None: self._eref = V0R
          # save reaction name and check if beyond mep
          self._reaction = rname
          if self._eref in [None,"auto","default"]: self._beyondmep = False
          try:
              Rs,TS,Ps = dchem[self._reaction]
              self._reactioneq = '%s --> %s --> %s'%("+".join(Rs),TS,"+".join(Ps))
          except: self._reactioneq = None
          return rname


      def get_layers(self):
          return (self._oniomh,self._oniomm,self._onioml)

      def isONIOMok(self,natoms,software):
          num_oniom = len(set(self._oniomh+self._oniomm+self._onioml))
          if   num_oniom == 0         : return True
          elif num_oniom != natoms    : return False
          elif software  != "gaussian": return False
          else                        : return True

      def prepare_qrc(self,dchem,dctc,dimasses):
          '''
          also modifies self._qrccase:
              self._qrccase == 0: everything is ok / qrc is not activated
              self._qrccase == 1: no reaction is associated to the TS
              self._qrccase == 2: reaction is not unimolecular
              self._qrccase == 3: ctc for reactant not defined
              self._qrccase == 4: gts file for reactant not found
              self._qrccase == 5: unable to get energy of products
          '''

          # Will qrc be used?
          if self._qrc      is None: return
          if self._reaction is None: self._qrccase = 1; return
          # assert unimolecular
          reactants = dchem[self._reaction][0]
          products  = dchem[self._reaction][2]
          if self._exorgic:
             if len(reactants) != 1: self._qrccase = 2; return
             self._qrcname = reactants[0]
          else:
             if len(products ) != 1: self._qrccase = 2; return
             self._qrcname = products[0]
          if self._V1P is None: self._qrccase = 5; return
          #---------------------------------#
          # Now, generate Molecule instance #
          #---------------------------------#
          ctc, itc = PN.name2data(self._qrcname)
          if ctc not in dctc.keys(): self._qrccase = 3; return
          cluster = dctc[ctc]
          if itc is None: itc = cluster._itcs[0][0]
          if   itc in cluster._diso.keys(): imods = cluster._diso[itc]
          elif "*" in cluster._diso.keys(): imods = cluster._diso["*"]
          else                            : imods = None
          gtsfile = dctc[ctc].gtsfile(itc)
          if not os.path.exists(gtsfile): self._qrccase = 4; return
          # Generate Molecule instance
          molecule = Molecule()
          molecule.set_from_gts(gtsfile)
          # apply fscal
          molecule.setvar(fscal=cluster._fscal)
          # apply isotopic masses
          if imods is not None:
             molecule.apply_imods(imods,dimasses)
          # calculate frequencies
          molecule.setup()
          #------------------------------#
          # Prepare list of QRC energies #
          #------------------------------#
          mode , nE      = self._qrc
          self._qrcafreq = molecule._ccfreqs[mode]
          self._qrclE    = [n*HBAR*self._qrcafreq for n in range(nE)]

      def sct_convergence(self):
          if self._sct == "no"   : return False
          if self._scterr is None: return False
          if self._scterr <= 0.0 : return False
          return True

      def increase_svals(self,V0bw,V0fw,V1bw,V1fw):
          ds       = self._ds
          hsteps   = self._hsteps
          dV0      = (V0fw-V0bw)*KCALMOL
          # Decide where to increase
          if abs(dV0) < EPS_KCALMOL     : increase = "+-"
          elif self._V1R > V1fw         : increase = "-"
          elif self._V1P is None        : increase = "+-"
          elif dV0 > +abs(EPS_MEPINCR)  : increase = "+"
          elif dV0 < -abs(EPS_MEPINCR)  : increase = "-"
          else                          : increase = "+-"
          # Modify MEP limits
          if "-" in increase and not self._convbw: self._sbw -= ds*hsteps
          if "+" in increase and not self._convfw: self._sfw += ds*hsteps

      def set_ics(self,ics,icsbw=None,icsfw=None):
          self._ics   = ics
          self._icsbw = icsbw
          self._icsfw = icsfw

      def get_ics(self): return self._ics

      def tuple_rst(self):
          d3 = sd.cubic2float(self._cubic)
          return (self._mtype,self._mu,self._ds,self._hsteps,d3)

      def tuple_first(self):
          d3 = sd.cubic2float(self._cubic)
          return (self._ds,self._mu,d3,self._fwdir)

      def tuple_sdbw(self):
          return (self._mtype,self._mu,self._ds,self._sbw,\
                  self._hsteps,self._epse,self._epsg)

      def tuple_sdfw(self):
          return (self._mtype,self._mu,self._ds,self._sfw,\
                  self._hsteps,self._epse,self._epsg)
#===============================================================#





