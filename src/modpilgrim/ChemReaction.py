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
| Sub-module :  ChemReaction       |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#==============================================#
import numpy as np
#----------------------------------------------#
import common.partfns   as partfns
import common.fncs      as fncs
from   common.Molecule  import Molecule
from   common.physcons  import AMU, KB, VOL0
#----------------------------------------------#
import modpilgrim.names as PN
#==============================================#


RCONS  = "tst,mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct"
RCONS +=    ",mststzct,mststsct,mscvt,mscvtzct,mscvtsct"
RCONS  = RCONS.split(",")


#==============================================#
class ChemReaction():

      #========================================#
      def __init__(self,rcname,reacts,ts,prods,ltemp,dctc):

          # save input data
          self._rcname    = rcname
          self._ltemp     = ltemp
          self._beta      = 1.0/(np.array(ltemp)*KB)
          self._reacts    = reacts
          self._ts        = ts
          self._prods     = prods

          #---------------------------------------#
          # Some basic operations with input data #
          #---------------------------------------#

          remove        = []
          self._itcs    = {}
          self._gtsfile = {} # name of gts files
          self._anhfile = {} # name of anh files
          self._disos   = {} # isotopic modifications
          self._weight  = {}
          for target in self._reacts+[self._ts]+self._prods:
              ctc,itcs,ms = None,None,None
              # get itcs
              if target is not None:
                 ctc, itc = PN.name2data(target)
                 # list of itcs
                 if   itc is not None   : itcs = [(itc,1)]
                 elif ctc in dctc.keys(): itcs = list(dctc[ctc]._itcs)
                 else                   : itcs = None
                 # multi-structural?
                 ms = itc is None # (True also with 1 conformer if itc is None)
              # itcs NOT found
              if itcs is None:
                 remove.append(target)
                 continue
              # save data
              self._anhfile[ctc] = dctc[ctc]._anh
              self._disos[ctc]   = dctc[ctc]._diso
              self._itcs[target] = (ctc,itcs,ms)
              for itc,weight in itcs:
                  self._weight [(ctc,itc)] = weight
                  self._gtsfile[(ctc,itc)] = dctc[ctc].gtsfile(itc)

          # clean up
          self._reacts = [target for target in self._reacts if target not in remove]
          self._prods  = [target for target in self._prods  if target not in remove]
          if self._ts in remove: self._ts = None

          # number of reactants/product and correction weight
          self._nR,self._wfw = len(self._reacts),1
          self._nP,self._wbw = len(self._prods ),1
          if self._nR == 2 and self._reacts[0] == self._reacts[1]: self._wfw = 2
          if self._nP == 2 and self._prods[0]  == self._prods[1] : self._wbw = 2


          #------------------------------#
          # initialize rest of variables #
          #------------------------------#
          self._dall      = {}
          self._massch    = {}

          self._massR     = 0.0
          self._massTS    = 0.0
          self._massP     = 0.0
          self._chR       = 0
          self._chTS      = 0
          self._chP       = 0

          self._V0R       = 0.0
          self._V1R       = 0.0
          self._V0P       = 0.0
          self._V1P       = 0.0
          self._V0TS      = 0.0
          self._V1TS      = 0.0

          self._QtR       = np.array([1.0 for t in self._ltemp])
          self._QtP       = np.array([1.0 for t in self._ltemp])
          self._QtTS      = np.array([1.0 for t in self._ltemp])
          self._chi0      = {} # contribution of most stable conformer [considering V0]

          self._tsitc0    = None # for ctc of ts [considering V0]
          self._tschi     = {X:{} for X in "tst,tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")}
          self._dtcoef    = {X:{} for X in "tstzct,tstsct,cvt,cvtzct,cvtsct".split(",")}

          self._anhctcs   = set([])
          self._FMSQH     = np.array([1.0 for t in self._ltemp]) # SS considering V0
          self._ANHR      = np.array([1.0 for t in self._ltemp])
          self._ANHP      = np.array([1.0 for t in self._ltemp])
          self._ANHTS     = np.array([1.0 for t in self._ltemp])
          self._ANHKeq    = np.array([1.0 for t in self._ltemp])
          self._ANHkfw    = np.array([1.0 for t in self._ltemp])
          self._ANHkbw    = np.array([1.0 for t in self._ltemp])

          self._Keq       = np.array([1.0 for t in self._ltemp])
          
          self._kfw       = {X:None for X in RCONS}
          self._kbw       = {X:None for X in RCONS}

          self._problem   = False
      #========================================#

      #========================================#
      def external_data(self,dall):
          self._dall     = dall
      #========================================#

      #========================================#
      def get_masses_charges(self,dimasses={}):
          for ctc,itcs,ms in self._itcs.values():
              # select first itc
              itc,weight = itcs[0]
              # generate instance of Molecule
              molecule = Molecule()
              # read gts
              gtsfile  = self._gtsfile[(ctc,itc)]
              molecule.set_from_gts(gtsfile)
              # isotopic modification
              diso = self._disos[ctc]
              if   itc in diso.keys(): imod = diso[itc]
              elif "*" in diso.keys(): imod = diso["*"]
              else                   : imod = None
              molecule.apply_imods(imod,dimasses)
              # data to save: mass & charge
              mass   = float(molecule._mass)
              charge = int(molecule._ch)
              # save data
              self._massch[ctc] = (mass,charge)
      #========================================#

      #========================================#
      def return_V0V1(self,ctc,itc):
          key = PN.struckey(ctc,itc)
          V0,V1,pfns = self._dall["pfn"][key]
          return V0, V1
      #========================================#

      #========================================#
      def check_conservation(self):
          
          # reactants
          if self._nR != 0:
             self._massR,self._chR = 0.0, 0
             for target in self._reacts:
                 ctc,itcs,ms  = self._itcs[target]
                 self._massR += self._massch[ctc][0]
                 self._chR   += self._massch[ctc][1]

          # transition state
          if self._ts is not None:
             self._massP,self._chP = 0.0, 0
             ctc,itcs,ms  = self._itcs[self._ts]
             self._massTS = self._massch[ctc][0]
             self._chTS   = self._massch[ctc][1]

          # products
          if self._nP != 0:
             for target in self._prods:
                 ctc,itcs,ms  = self._itcs[target]
                 self._massP += self._massch[ctc][0]
                 self._chP   += self._massch[ctc][1]

          # Now check
          eps_mass = 0.001/AMU
          if self._nR*self._nP != 0:
             if self._chR != self._chP: self._problem = True
             if abs(self._massR-self._massP) > eps_mass: self._problem = True
          if self._ts is not None and self._nR != 0:
             if self._chR != self._chTS: self._problem = True
             if abs(self._massR-self._massTS) > eps_mass: self._problem = True
          if self._ts is not None and self._nP != 0:
             if self._chP != self._chTS: self._problem = True
             if abs(self._massP-self._massTS) > eps_mass: self._problem = True
      #========================================#

      #========================================#
      def calculate_chi0(self,target):
          ctc,itcs,ms = self._itcs[target]
          # Get total partition function
          key        = PN.struckey(ctc,"msho")
          V0,V1,PFN  = self._dall["pfn"][key ]
          # Get contribution of each conformer
          itc0,V0_0,V1_0,PFN_0 = None,float("inf"),float("inf"),None
          for itc,weight_i in itcs:
              # individual partition function
              key_i = PN.struckey(ctc,itc)
              V0_i,V1_i,PFN_i = self._dall["pfn"][key_i]
              if V0_i < V0_0: itc0,V0_0,V1_0,PFN_0 = itc,V0_i,V1_i,PFN_i
          # contribution of most stable conformer (wo considering weight)
          dE        = (V1_0-V1)
          exp_arg   = [-dE/KB/T for T in self._ltemp]
          ratio_pfn = [pfi/pftot for pfi,pftot in zip(PFN_0,PFN)]
          chi_i     = [aa*fncs.exp128(bb) for aa,bb in zip(ratio_pfn,exp_arg)]
          self._chi0[target] = (itc0,np.array(chi_i))
      #========================================#


      #========================================#
      def return_pfns(self,target):
          ctc,itcs,ms = self._itcs[target]
          # key for partition functions
          if ms: key = PN.struckey(ctc,"msho")
          else : key = PN.struckey(ctc,itcs[0][0])
          # Energies and partition functions from data
          V0,V1,pfns = self._dall["pfn"][key]
          # anharmonicity (if commented in pif.struc --> not include it)
          anh = None
          if ms and self._anhfile[ctc] is not None:
             anh = self._dall["anh"].get(ctc,None)
          # Save data
          if anh is None: anh = np.array([1.0 for t in self._ltemp])
          else          : self._anhctcs.add(ctc)
          return V0, V1, pfns, anh
      #========================================#

      #========================================#
      def obtain_pfns(self):
          # (a) reactants
          for target in self._reacts:
              self.calculate_chi0(target)
              V0, V1, pfns, anh = self.return_pfns(target)
              self._V0R    += V0
              self._V1R    += V1
              self._QtR    *= pfns
              self._ANHR   *= anh
              self._FMSQH  *= self._chi0[target][1]

          # (b) products
          for target in self._prods:
              V0, V1, pfns, anh = self.return_pfns(target)
              self._V0P    += V0
              self._V1P    += V1
              self._QtP    *= pfns
              self._ANHP   *= anh

          # (c) transition state
          if self._ts is None: return
          # (c.1) get partition functions
          self.calculate_chi0(self._ts)
          self._V0TS, self._V1TS, self._QtTS, self._ANHTS = self.return_pfns(self._ts)
          # (c.2) Conversion betweem SS-QH and MS-QH rate constants
          self._FMSQH /= self._chi0[self._ts][1]
          # (c.2) get individual contributions of TS and itc of min(V0)
          ctc, itcs, ms = self._itcs[self._ts]
          if ms:
             V0,V1,PFN = self._dall["pfn"][PN.struckey(ctc,"msho")]
             # save most stable conformer
             self._tsitc0 = self._chi0[self._ts][0]
             # calculate chi for each conformer of the ts
             for itc,weight in itcs:
                 # see pfn
                 V0i,V1i,PFNi = self._dall["pfn"][PN.struckey(ctc,itc)]
                 dE        = (V1i-V1)
                 exp_arg   = [-dE/KB/T for T in self._ltemp]
                 ratio_pfn = [weight*pfi/pftot for pfi,pftot in zip(PFNi,PFN)]
                 chi_i     = [aa*fncs.exp128(bb) for aa,bb in zip(ratio_pfn,exp_arg)]
                 # save data
                 self._tschi["tst"][itc] = np.array(chi_i)
          else :
              itc,weight = itcs[0]
              self._tschi["tst"][itc] = np.array( [1.0 for T in self._ltemp] )
              self._tsitc0 = None
      #========================================#

      #========================================#
      def obtain_pfnparts(self,which):
          if which == "reactants": targets, V1 = self._reacts, self._V1R
          if which == "products" : targets, V1 = self._prods , self._V1P
          if which == "ts"       : targets, V1 = [self._ts]  , self._V1TS

          Qtr = np.array([1.0 for T in self._ltemp])
          Qrv = np.array([1.0 for T in self._ltemp])
          Qel = np.array([1.0 for T in self._ltemp])
          for target in targets:
              # get contributions
              tra, rovib, ele = self._dall["ctr"][target]
              # Apply data
              Qtr *= np.array(tra)
              Qrv *= np.array(rovib)
              Qel *= np.array(ele)
          return Qtr, Qrv, Qel, V1
      #========================================#

      #========================================#
      def calculate_anharmonicity(self):
          self._ANHKeq = self._ANHP  / self._ANHR
          self._ANHkfw = self._ANHTS / self._ANHR
          self._ANHkbw = self._ANHTS / self._ANHP
      #========================================#

      #========================================#
      def calculate_eqconstant(self):
          # calculate eq constant for v=1cm^3/molecule
          self._Keq = partfns.Qs2Kc(self._ltemp,self._QtR,self._QtP,self._V1R,self._V1P,self._nR,self._nP)
          self._Keq = float(self._wfw)/float(self._wbw) * np.array(self._Keq)
          # include anharmonicity
          self._Keq = self._Keq * self._ANHKeq
      #========================================#

      #========================================#
      def calculate_transcoeffs(self):
          if self._ts is None: return
          ctc, itcs, ms = self._itcs[self._ts]
          # Get individual transmission coefficients
          for itc,weight in itcs:
              tsname = PN.struckey(ctc,itc)
              cvt    = self._dall.get("cvt"   ,{}).get(tsname,None)
              zct    = self._dall.get("zct"   ,{}).get(tsname,None)
              sct    = self._dall.get("sct"   ,{}).get(tsname,None)
              cagtst = self._dall.get("cagtst",{}).get(tsname,None)
              cagcvt = self._dall.get("cagcvt",{}).get(tsname,None)
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
              self._dtcoef["tstzct"][itc] = tc_tstzct
              self._dtcoef["tstsct"][itc] = tc_tstsct
              self._dtcoef["cvt"   ][itc] = tc_cvt
              self._dtcoef["cvtzct"][itc] = tc_cvtzct
              self._dtcoef["cvtsct"][itc] = tc_cvtsct

          # Calculate contributions to other methods
          for X in self._dtcoef.keys():
              # calculate numerator of contribution
              for itc,weight in itcs:
                  chi_i = self._tschi["tst"][itc]
                  tc_i  = self._dtcoef[X][itc]
                  if tc_i is not None: self._tschi[X][itc] = chi_i*tc_i
                  else               : self._tschi[X][itc] = None
              # averaged transmission coefficient
              chis_X = [chi_i_X for chi_i_X in self._tschi[X].values()]
              try   : self._dtcoef[X]["averaged"] = sum(chis_X)
              except: self._dtcoef[X]["averaged"] = None
              # calculate contribution
              averaged = self._dtcoef[X]["averaged"] 
              for itc,weight in itcs:
                  chi_i = self._tschi["tst"][itc]
                  tc_i  = self._dtcoef[X][itc]
                  if averaged is None or tc_i is None: self._tschi[X][itc] = None
                  else                               : self._tschi[X][itc] /= averaged
      #========================================#

      #========================================#
      def calculate_rateconstants(self):
          if self._ts is None: return
          ctc, itcs, ms = self._itcs[self._ts]

          # (A) MS-TST rate constant
          Kpseudo = partfns.Qs2Kc(self._ltemp,self._QtR,self._QtTS,self._V1R,self._V1TS)
          kTST    = partfns.Kc2rate(self._ltemp,Kpseudo)
          kTST    = float(self._wfw) * np.array(kTST)

  #       # (B) SS-TST rate constant
  #       kSSQH   = np.array(kTST) / self._FMSQH
  #       self._kfw["sstst"] = kSSQH

          # (C) MS-TST rate constant with anharmonicity
          kTST    = self._ANHkfw * np.array(kTST)
          self._kfw["tst"] = kTST

          # (D) Multi-path
          tc_tstzct = self._dtcoef["tstzct"]["averaged"]
          tc_tstsct = self._dtcoef["tstsct"]["averaged"]
          tc_cvt    = self._dtcoef["cvt"   ]["averaged"]
          tc_cvtzct = self._dtcoef["cvtzct"]["averaged"]
          tc_cvtsct = self._dtcoef["cvtsct"]["averaged"]

          if tc_tstzct is not None: self._kfw["mptstzct"] = tc_tstzct * kTST
          if tc_tstsct is not None: self._kfw["mptstsct"] = tc_tstsct * kTST
          if tc_cvt    is not None: self._kfw["mpcvt"   ] = tc_cvt    * kTST
          if tc_cvtzct is not None: self._kfw["mpcvtzct"] = tc_cvtzct * kTST
          if tc_cvtsct is not None: self._kfw["mpcvtsct"] = tc_cvtsct * kTST

          
          # (E) Multi-structural
          if ms:
             tc_tstzct = self._dtcoef["tstzct"][self._tsitc0]
             tc_tstsct = self._dtcoef["tstsct"][self._tsitc0]
             tc_cvt    = self._dtcoef["cvt"   ][self._tsitc0]
             tc_cvtzct = self._dtcoef["cvtzct"][self._tsitc0]
             tc_cvtsct = self._dtcoef["cvtsct"][self._tsitc0]

             if tc_tstzct is not None: self._kfw["mststzct"] = tc_tstzct * kTST
             if tc_tstsct is not None: self._kfw["mststsct"] = tc_tstsct * kTST
             if tc_cvt    is not None: self._kfw["mscvt"   ] = tc_cvt    * kTST
             if tc_cvtzct is not None: self._kfw["mscvtzct"] = tc_cvtzct * kTST
             if tc_cvtsct is not None: self._kfw["mscvtsct"] = tc_cvtsct * kTST

  #       # (F) REST OF SINGLE STRUCTURE
  #       tc_tstzct = self._dtcoef["tstzct"][self._tsitc0]
  #       tc_tstsct = self._dtcoef["tstsct"][self._tsitc0]
  #       tc_cvt    = self._dtcoef["cvt"   ][self._tsitc0]
  #       tc_cvtzct = self._dtcoef["cvtzct"][self._tsitc0]
  #       tc_cvtsct = self._dtcoef["cvtsct"][self._tsitc0]

  #       if tc_tstzct is not None: self._kfw["sststzct"] = tc_tstzct * kSSQH
  #       if tc_tstsct is not None: self._kfw["sststsct"] = tc_tstsct * kSSQH
  #       if tc_cvt    is not None: self._kfw["sscvt"   ] = tc_cvt    * kSSQH
  #       if tc_cvtzct is not None: self._kfw["sscvtzct"] = tc_cvtzct * kSSQH
  #       if tc_cvtsct is not None: self._kfw["sscvtsct"] = tc_cvtsct * kSSQH

          # (G) BACKWARDS
          if self._nP != 0:
             for X in self._kfw.keys():
                 kfw = self._kfw[X]
                 if kfw is None: continue
                 self._kbw[X] = VOL0**(self._nP-self._nR) *kfw / self._Keq
      #========================================#

#==============================================#

