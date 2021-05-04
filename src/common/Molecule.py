'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.3
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
| Sub-module :  Molecule           |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the Molecule class
'''

#=============================================#
import os
import numpy              as np
#---------------------------------------------#
import common.fncs        as     fncs
import common.partfns     as     pf
import common.internal    as     intl
import common.Exceptions  as     Exc
from   common.criteria    import EPS_IC
from   common.dicts       import dpt_im
from   common.files       import read_gtsfile
from   common.files       import write_gtsfile
from   common.files       import write_xyz, write_molden
from   common.pgs         import get_pgs
from   common.physcons    import AMU
from   common.physcons    import KCALMOL
from   common.physcons    import EV
from   common.physcons    import ANGSTROM
from   common.physcons    import H2CM
from   common.gaussian    import read_fchk
from   common.gaussian    import read_gauout
#=============================================#

class Molecule():

      # Initialization method
      def __init__(self,label=None):

          self._label    = label

          # Unidimensional
          self._mform    = "-"
          self._mu       = None
          self._ch       = None
          self._mtp      = None
          self._V0       = None
          self._pgroup   = None
          self._rotsigma = None
          self._natoms   = None
          self._nel      = None # number of electrons
          self._rtype    = None
          self._linear   = None

          # Multi-dimensional
          self._atnums   = None
          self._symbols  = None
          self._masses   = None
          self._les      = None # list of electronic states
          self._itensor  = None
          self._imoms    = None
          self._rotTs    = None

          # Arrays of importance
          self._xcc      = None
          self._gcc      = None
          self._Fcc      = None
          self._xms      = None
          self._gms      = None
          self._Fms      = None

          # related to frequencies
          self._fscal    = 1.0
          self._nvdof    = None
          self._cczpe    = None
          self._ccfreqs  = None
          self._ccFevals = None
          self._ccFevecs = None
          self._iczpe    = None
          self._icfreqs  = None
          self._icFevals = None
          self._icFevecs = None

          # other stuff for very particular occasion
          self._gts      = None

      def __str__(self): return self._mform

      def setvar(self,xcc=None,gcc=None,Fcc=None,\
                      atonums=None,symbols=None,masses=None,\
                      ch=None,mtp=None, V0=None,\
                      pgroup=None,rotsigma=None,\
                      fscal=None,les=None):

          if xcc      is not None: self._xcc      = xcc
          if gcc      is not None: self._gcc      = gcc
          if Fcc      is not None: self._Fcc      = Fcc

          if atonums  is not None: self._atnums   = atonums
          if symbols  is not None: self._symbols  = symbols
          if masses   is not None: self._masses   = masses

          if ch       is not None: self._ch       = int(ch)
          if mtp      is not None: self._mtp      = int(mtp)
          if V0       is not None: self._V0       = V0

          if pgroup   is not None: self._pgroup   = pgroup
          if rotsigma is not None: self._rotsigma = rotsigma

          if fscal    is not None: self._fscal    = fscal
          if les      is not None: self._les      = les

      def genderivates(self):
          self._mform   = fncs.get_molformula(self._symbols)
          self._natoms  = len(self._atnums)
          self._mass    = sum(self._masses)
          self._nel     = sum(self._atnums)-self._ch
          if self._les is None: self._les = [ (self._mtp,0.0) ]

      def prepare(self):
          # check atnums
          if self._atnums is not None and type(self._atnums[0]) == str:
             self._symbols = list(self._atnums)
          # check symbols
          if self._symbols is not None and type(self._symbols[0]) == int:
             self._atnums  = list(self._symbols)
          # Get both atnums and symbols if None
          if self._atnums  is None: self._atnums  = fncs.symbols2atonums(self._symbols)
          if self._symbols is None: self._symbols = fncs.atonums2symbols(self._atnums)
          # check masses
          if self._masses is None:
             self._masses = fncs.atonums2masses(self._atnums)
          # derivated magnitudes
          self.genderivates()
          # check Fcc
          if self._Fcc not in (None,[]) and len(self._Fcc) != 3*self._natoms:
             self._Fcc = fncs.lowt2matrix(self._Fcc)

      def calc_pgroup(self,force=False):
          calculate = False
          if force                 : calculate = True
          if self._pgroup   is None: calculate = True
          if self._rotsigma is None: calculate = True
          if calculate: self._pgroup,self._rotsigma = get_pgs(self._atnums,self._masses,self._xcc)

      def remove_frozen(self):
          frozen = fncs.detect_frozen(self._Fcc,self._natoms)
          if len(frozen) == 0: return [],[]
          # coordinates and symbols of frozen moiety
          bN   = [at in frozen for at in range(self._natoms)]
          b3N  = [at in frozen for at in range(self._natoms) for ii in range(3)]
          frozen_xcc     = np.array(self._xcc)[b3N]
          frozen_symbols = np.array(self._symbols)[bN]
          # now system is just the flexible moiety
          bN   = [at not in frozen for at in range(self._natoms)]
          b3N  = [at not in frozen for at in range(self._natoms) for ii in range(3)]
          self._xcc      = np.array(self._xcc)[b3N].tolist()
          self._symbols  = np.array(self._symbols)[bN].tolist()
          self._atnums   = np.array(self._atnums)[bN].tolist()
          self._masses   = np.array(self._masses)[bN].tolist()
          self._pgroup   = None
          self._rotsigma = None
          # Gradient and hessian
          if self._gcc is not None and len(self._gcc) != 0:
             self._gcc     = np.array(self._gcc)[b3N].tolist()
          if self._Fcc is not None and len(self._Fcc) != 0:
             n3 = self._natoms*3
             self._Fcc = [[self._Fcc[idx1][idx2] for idx1 in range(n3) if b3N[idx1]]\
                                                 for idx2 in range(n3) if b3N[idx2]]
          # set origin for frozen moiety
          com = fncs.get_com(self._xcc,self._masses)
          frozen_xcc = fncs.set_origin(frozen_xcc,com)
          # prepare system
          self.prepare()
          return frozen_xcc, frozen_symbols

      def mod_masses(self,masses):
          self._masses = list(masses)
          self._mass   = sum(self._masses)
          # re-calculate point group
          self.calc_pgroup(force=True)

      def apply_imods(self,imods,imasses):
          '''
          example: imods   = ["H2(4,5)","C13(all_C)"]
                   imasses = {"H2":2.0141/AMU, "C13":13.0034/AMU}
          '''
          if imods is None: return

          for imod in imods:
              isymbol = imod.split("(")[0]
              if   isymbol in imasses.keys(): imass = imasses[isymbol]
              elif isymbol in  dpt_im.keys(): imass =  dpt_im[isymbol]
              else:
                 exception = Exc.WrongInIsomass
                 exception._var = isymbol
                 raise exception
              atoms   = imod.split("(")[1].split(")")[0]
              if "all_" in atoms:
                 atype = atoms.split("all_")[1].strip()
                 for idx,symbol in enumerate(self._symbols):
                     if symbol == atype: self._masses[idx] = imass
              else:
                 list_of_atoms = []
                 for atom in atoms.split(","):
                     if "-" in atom:
                        at1,atn = atom.split("-")
                        list_of_atoms += range(int(at1),int(atn)+1)
                     else: list_of_atoms.append(int(atom))
                 list_of_atoms = sorted(list(set(list_of_atoms)))
                 for idx in list_of_atoms: self._masses[idx-1] = imass
          # re-calculate total mass and point group
          self.mod_masses(self._masses)

      def setup(self,mu=1.0/AMU,projgrad=False):
          self._mu = mu
          # derivated magnitudes (again, in case sth was modified)
          # for example, when set from gts and masses are added latter
          self.genderivates()
          # shift to center of mass and reorientate molecule
          idata = (self._xcc,self._gcc,self._Fcc,self._masses)
          self._xcc, self._gcc, self._Fcc = fncs.center_and_orient(*idata)
          # symmetry
          self.calc_pgroup(force=False)
          # Generate mass-scaled arrays
          self._xms = fncs.cc2ms_x(self._xcc,self._masses,self._mu)
          self._gms = fncs.cc2ms_g(self._gcc,self._masses,self._mu)
          self._Fms = fncs.cc2ms_F(self._Fcc,self._masses,self._mu)
          #-------------#
          # Atomic case #
          #-------------#
          if self._natoms == 1:
              self._nvdof    = 0
              self._linear   = False
             #self._xms      = list(self._xcc)
             #self._gms      = list(self._gcc)
             #self._Fms      = list(self._Fcc)
              self._ccfreqs  = []
              self._ccFevals = []
              self._ccFevecs = []
          #----------------#
          # Molecular case #
          #----------------#
          else:
             # Calculate inertia
             self._itensor = fncs.get_itensor_matrix(self._xcc,self._masses)
             self._imoms, self._rotTs, self._rtype, self._linear = \
                     fncs.get_itensor_evals(self._itensor)
             # Vibrational degrees of freedom
             if self._linear: self._nvdof = 3*self._natoms - 5
             else           : self._nvdof = 3*self._natoms - 6
             # calculate frequencies
             if self._Fcc is None        : return
             if len(self._Fcc) == 0      : return
             if self._ccfreqs is not None: return

             v0   = self._gms if projgrad else None
             data = fncs.calc_ccfreqs(self._Fcc,self._masses,self._xcc,self._mu,v0=v0)
             self._ccfreqs, self._ccFevals, self._ccFevecs = data
             # Scale frequencies
             self._ccfreqs = fncs.scale_freqs(self._ccfreqs,self._fscal)

      def get_imag_main_dir(self):
          ic, fwsign = intl.ics_idir(self._xcc,self._symbols,\
                       self._masses,self._ccfreqs,self._ccFevecs)
          return ic, fwsign

      def icfreqs(self,ics,bool_pg=False):
          #----------------#
          # Molecular case #
          #----------------#
          if self._natoms != 1:
             ituple = (self._Fcc,self._masses,self._xcc,self._gcc,ics,bool_pg)
             self._icfreqs, self._icFevals, self._icFevecs = intl.calc_icfreqs(*ituple)
          #-------------#
          # Atomic case #
          #-------------#
          else:
             self._icfreqs  = []
             self._icFevals = []
             self._icFevecs = []
          # scale frequencies
          self._icfreqs = [freq*self._fscal for freq in self._icfreqs]

      def ana_freqs(self,case="cc"):
          if case == "cc":
             # Keep record of imaginary frequencies
             if self._ccFevecs is not None:
                self._ccimag = [ (frq,self._ccFevecs[idx]) for idx,frq in enumerate(self._ccfreqs)\
                                 if frq < 0.0]
             else:
                self._ccimag = [ (frq,None)                for idx,frq in enumerate(self._ccfreqs)\
                                if frq < 0.0]
             # Calculate zpe
             self._cczpes = [fncs.afreq2zpe(frq) for frq in self._ccfreqs]
             self._cczpe  = sum(self._cczpes)
             self._ccV1   = self._V0 + self._cczpe
          if case == "ic":
             # Keep record of imaginary frequencies
             if self._icFevecs is not None:
                self._icimag = [ (frq,self._icFevecs[idx]) for idx,frq in enumerate(self._icfreqs)\
                                 if frq < 0.0]
             else:
                self._icimag = [ (frq,None)                for idx,frq in enumerate(self._icfreqs)\
                                if frq < 0.0]
             # Calculate zpe
             self._iczpes = [fncs.afreq2zpe(frq) for frq in self._icfreqs]
             self._iczpe  = sum(self._iczpes)
             self._icV1   = self._V0 + self._iczpe

      def clean_freqs(self,case="cc"):
          # select case
          if case == "cc": freqs = self._ccfreqs
          else           : freqs = self._icfreqs
          # keep track of those to save
          keep = []
          for idx,freq in enumerate(freqs):
              if abs(fncs.afreq2cm(freq)) < EPS_IC: continue
              keep.append(idx)
          # keep only those > EPS_IC
          if case == "cc":
             self._ccfreqs  = [self._ccfreqs[idx]  for idx in keep]
             if self._ccFevals is not None:
                self._ccFevals = [self._ccFevals[idx] for idx in keep]
             if self._ccFevecs is not None:
                self._ccFevecs = [self._ccFevecs[idx] for idx in keep]
          if case == "ic":
             self._icfreqs  = [self._icfreqs[idx]  for idx in keep]
             if self._icFevals is not None:
                self._icFevals = [self._icFevals[idx] for idx in keep]
             if self._icFevecs is not None:
                self._icFevecs = [self._icFevecs[idx] for idx in keep]

      def deal_lowfq(self,lowfq={},case="cc"):
          # for Cartesian Coordinates
          if   case == "cc":
             # frequencies were not projected along MEP
             if   self._nvdof - len(self._ccfreqs) == 0:
                for idx,newfreq in lowfq.items():
                    self._ccfreqs[idx] = max(self._ccfreqs[idx],newfreq)
             # frequencies were projected along MEP
             elif self._nvdof - len(self._ccfreqs) == 1:
                for idx,newfreq in lowfq.items():
                    self._ccfreqs[idx-1] = max(self._ccfreqs[idx-1],newfreq)
          # for Internal Coordinates
          elif case == "ic":
             # frequencies were not projected along MEP
             if   self._nvdof - len(self._icfreqs) == 0:
                for idx,newfreq in lowfq.items():
                    self._icfreqs[idx] = max(self._icfreqs[idx],newfreq)
             # frequencies were projected along MEP
             elif self._nvdof - len(self._icfreqs) == 1:
                for idx,newfreq in lowfq.items():
                    self._icfreqs[idx-1] = max(self._icfreqs[idx-1],newfreq)

      def calc_pfns(self,temps,case="cc",fmode=0,imag=1E10):
          '''
          fmode = -1 or 0 (0 is default)
          '''
          # Calculate translational partition function (per unit volume)
          ph_tra = np.array([pf.pf_partinbox(self._mass,T) for T in temps])
          # Calculate rotational partition function (Rigid-Rotor)
          if self._natoms > 1:
             pf_rot = np.array([pf.pf_rigidrotor(self._imoms,T,self._rotsigma) for T in temps])
          else:
             pf_rot = np.array([1.0 for T in temps])
          # Calculate vibrational partition function (Harmonic-Oscillator)
          if self._nvdof != 0:
             # remove freq if required
             nf     = self._nvdof + fmode
             if case == "cc": afreqs = list(self._ccfreqs)
             if case == "ic": afreqs = list(self._icfreqs)
             while len(afreqs) > nf: afreqs = afreqs[1:]
             # Calculate vib pfn
             pf_vib = np.array([pf.pf_harmosc(afreqs,T,imag=imag) for T in temps])
          else:
             pf_vib = np.array([1.0 for T in temps])
          # Calculate electronic partition function
          pf_ele = np.array([pf.pf_electr(self._les,T) for T in temps])
          # Total partition function
          qtot = ph_tra * pf_rot * pf_vib * pf_ele
          if case == "cc": return qtot, self._ccV1, (ph_tra,pf_rot,pf_vib,pf_ele)
          if case == "ic": return qtot, self._icV1, (ph_tra,pf_rot,pf_vib,pf_ele)

      def info_string(self,ib=0):
          root_mass = sum(fncs.symbols2masses(self._symbols))
          string  = "Molecular formula     : %s\n"%self._mform
          string += "Number of atoms       : %i\n"%self._natoms
          string += "Number of electrons   : %i\n"%self._nel
          string += "Vibrational DOFs      : %i\n"%self._nvdof
          string += "Charge                : %i\n"%self._ch
          string += "Multiplicity          : %i\n"%self._mtp
          string += "Electronic energy (V0): %.8f hartree\n"%self._V0
          string += "Total mass [root]     : %.4f amu\n"%(root_mass *AMU)
          string += "Total mass            : %.4f amu\n"%(self._mass*AMU)
          if self._pgroup   is not None: string += "Point group symmetry  : %s\n"%(self._pgroup)
          if self._rotsigma is not None: string += "Rotational sym num    : %i\n"%(self._rotsigma)
          string += "Cartesian coordinates (Angstrom):\n"
          for at,symbol in enumerate(self._symbols):
              mass   = self._masses[at]*AMU
              x,y,z  = fncs.xyz(self._xcc,at)
              x *= ANGSTROM
              y *= ANGSTROM
              z *= ANGSTROM
              string += "  %2s   %+10.6f  %+10.6f  %+10.6f  [%7.3f amu]\n"%(symbol,x,y,z,mass)

          try:
              str2  = "Moments and product of inertia (au):\n"
              if len(self._imoms) == 1:
                 str2 += "        %+10.3E\n"%self._imoms[0]
              if len(self._imoms) == 3:
                 prodinert = self._imoms[0]*self._imoms[1]*self._imoms[2]
                 dataline = (self._imoms[0],self._imoms[1],self._imoms[2],prodinert)
                 str2 += "        %+10.3E  %+10.3E  %+10.3E  [%10.3E]\n"%dataline
              string += str2
          except: pass

          try:
              str2  = "Vibrational frequencies [1/cm] (scaled by %.3f):\n"%self._fscal
              for idx in range(0,len(self._ccfreqs),6):
                  str2 += "  %s\n"%("  ".join("%8.2f"%fncs.afreq2cm(freq) \
                                      for freq in self._ccfreqs[idx:idx+6]))
              if len(self._ccfreqs) != 0: string += str2
          except: pass

          try:
              str2  = "Vibrational zero-point energies [kcal/mol]:\n"
              for idx in range(0,len(self._cczpes),6):
                  str2 += "  %s\n"%("  ".join("%8.2f"%(zpe*KCALMOL) \
                                      for zpe in self._cczpes[idx:idx+6]))
              zpe_au   = self._cczpe
              zpe_kcal = self._cczpe * KCALMOL
              zpe_eV   = self._cczpe * EV
              zpe_cm   = self._cczpe * H2CM
              str2 += "Vibrational zero-point energy: %+14.8f hartree  = \n"%zpe_au
              str2 += "                               %+14.2f kcal/mol = \n"%zpe_kcal
              str2 += "                               %+14.2f eV       = \n"%zpe_eV
              str2 += "                               %+14.2f cm^-1 \n"%zpe_cm
              str2 += "V0 + zero-point energy (V1)  : %+14.8f hartree\n"%self._ccV1
              if self._cczpe != 0.0: string += str2
          except: pass

          # add blank spaces
          string = "\n".join([" "*ib+line for line in string.split("\n")])
          return string

      #=======================================#
      # Calculation of geometric parameters   #
      #=======================================#
      def dihedral(self,at1,at2,at3,at4):
          x1 = self._xcc[3*at1:3*at1+3]
          x2 = self._xcc[3*at2:3*at2+3]
          x3 = self._xcc[3*at3:3*at3+3]
          x4 = self._xcc[3*at4:3*at4+3]
          return fncs.dihedral(x1,x2,x3,x4)
      #=======================================#

      #=======================================#
      # Set variables from external files     #
      #=======================================#
      def set_from_gts(self,gtsfile):
          if not os.path.exists(gtsfile): return
          # read file
          self._gts = gtsfile
          xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma,freq_list = read_gtsfile(self._gts)
          # set variables
          self.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
          self.setvar(atonums=atonums,masses=masses)
          self.setvar(ch=ch,mtp=mtp,V0=E,pgroup=pgroup,rotsigma=rotsigma)
          # Prepare system
          self.prepare()
          # only for developers: freq list
          if freq_list is not None and len(freq_list) != 0: self._ccfreqs = freq_list

      def set_from_fchk(self,fchk):
          if not os.path.exists(fchk): return
          # read file
          xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel = read_fchk(fchk)
          # set variables
          self.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
          self.setvar(atonums=atonums,masses=masses)
          self.setvar(ch=ch,mtp=mtp,V0=E)
          # Prepare system
          self.prepare()

      def set_from_gauout(self,gauout):
          if not os.path.exists(gauout): return
          # read file
          xcc, atonums, ch, mtp, E, gcc, Fcc, masses, calclevel = read_gauout(gauout)
          # set variables
          self.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
          self.setvar(atonums=atonums,masses=masses)
          self.setvar(ch=ch,mtp=mtp,V0=E)
          # Prepare system
          self.prepare()

      #=======================================#
      # Generation of different kind of files #
      #=======================================#
      def genfile_xyz(self,filename):
          try   : write_xyz(filename,self._xcc,self._symbols)
          except: return 0
          return 1
      #---------------------------------------#
      def genfile_molden(self,filename):
          try   : write_molden(filename,self._xcc,self._symbols,self._ccfreqs,self._ccFevecs)
          except: return 0
          return 1
      #---------------------------------------#
      def genfile_gts(self,filename,level=""):
          write_gtsfile(self._xcc,self._atnums,self._ch,self._mtp,\
                  self._V0,self._pgroup,self._rotsigma,self._gcc,\
                  self._Fcc,filename,level=level)

          try   : write_gtsfile(self._xcc,self._atnums,self._ch,self._mtp,\
                  self._V0,self._pgroup,self._rotsigma,self._gcc,\
                  self._Fcc,filename,level=level)
          except: return 0
          return 1
      #=======================================#

