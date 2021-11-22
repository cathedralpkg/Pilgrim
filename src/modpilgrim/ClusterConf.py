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
| Sub-module :  ClusterConf        |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the ClusterConf class
'''

#====================================================#
import os
import modpilgrim.names  as     PN
import common.Exceptions as     Exc
from   common.files      import findingts_energy_pg
from   common.internal   import unmerge_ics
from   common.internal   import string2ic
from   common.physcons   import KCALMOL
from   common.Molecule   import Molecule
from   common.fncs       import numimag
#====================================================#

class ClusterConf():

      def __init__(self,name):
          self._ctc      = name
          self._root     = name
          self._mformu   = None
          self._itcs     = []
          self._ch       = None
          self._mtp      = None
          self._es       = []
          self._type     = None
          self._fscal    = 1.0
          self._dics     = {}
          self._dicsbw   = {}
          self._dicsfw   = {}
          self._diso     = {}
          self._anh      = None
          self._lV0      = []
          self._lpg      = []
          self._V0       = None
          self._itcminV0 = None
          # list of molecule
          self._molecules = []

      def gtsfile(self,itc):
          return PN.get_gts(self._root,itc)

      def gtsfiles(self):
          return sorted([self.gtsfile(itc) for itc,weight in self._itcs])

      def setvar(self,key,val,vtype="int"):
          variables  = "ctc,root,mformu,itcs,ch,mtp,es,type,fscal,dics,"
          variables += "dicsbw,dicsfw,diso,anh,lV0,lpg,V0,itcminV0,molecules"
          if key not in variables.split(","): return
          if   vtype == "str": update_string = "self._%s='%s'"%(key,val)
          elif vtype == "int": update_string = "self._%s=int(%s)"%(key,val)
          else               : update_string = "self._%s=%s"%(key,val)
          exec(update_string)

      def set_from_piflines(self,lines):
          if type(lines) == str: lines = lines.split("\n")
          for line in lines:
              # common data
              if   line.startswith("root "     ): self._root   = line.split()[1]
              elif line.startswith("mformu "   ): self._mformu = line.split()[1]
              elif line.startswith("ch "       ): self._ch     = int(line.split()[1])
              elif line.startswith("mtp "      ): self._mtp    = int(line.split()[1])
              elif line.startswith("type "     ): self._type   = int(line.split()[1])
              elif line.startswith("freqscal " ): self._fscal  = float(line.split()[1])
              elif line.startswith("anharfile "): self._anh = line.split()[1]
              elif line.startswith("elestate " ):
                 deg, relE = line.split()[1:]
                 self._es += [(int(deg),float(relE))]
              # internal coordinates
              elif line.startswith("icsbw"):
                 key = line.split()[0]
                 if "." in key: itc = key.split(".")[1]
                 else         : itc = "*"
                 ics = [string2ic(ic) for ic in line.split()[1:]]
                 self._dicsbw[itc] = self._dicsbw.get(itc,[])+ics
              elif line.startswith("icsfw"):
                 key = line.split()[0]
                 if "." in key: itc = key.split(".")[1]
                 else         : itc = "*"
                 ics = [string2ic(ic) for ic in line.split()[1:]]
                 self._dicsfw[itc] = self._dicsfw.get(itc,[])+ics
              elif line.startswith("ics"):
                 key = line.split()[0]
                 if "." in key: itc = key.split(".")[1]
                 else         : itc = "*"
                 ics = [string2ic(ic) for ic in line.split()[1:]]
                 self._dics[itc] = self._dics.get(itc,[])+ics
              # isotopic modification
              elif line.startswith("iso"):
                 key = line.split()[0]
                 if "." in key: itc = key.split(".")[1]
                 else         : itc = "*"
                 imods = line.split()[1:]
                 self._diso[itc] = self._diso.get(itc,[])+imods
              # individual itc
              elif line.startswith("itc") or line.startswith("conformer"):
                 itc_name   = line.split()[1]
                 itc_weight = int(line.split()[3])
                 if itc_weight > 2: itc_weight = 1
                 if itc_weight < 1: itc_weight = 1
                 itc_tuple  = (itc_name,itc_weight)
                 self._itcs.append( itc_tuple )
          # sort itcs
          self._itcs.sort()
   
      def get_weight(self,itc):
          if type(itc) == int: itc = "%003i"%itc
          for itc_i,weight_i in self._itcs:
              if itc == itc_i: return weight_i
          return 0

      def add_itc(self,itc,V0,pg,weight=1):
          if type(itc) == int: itc = "%003i"%itc
          self._lV0.append(V0)
          self._lpg.append(pg)
          self._itcs.append( (itc,weight) )

      def remove_itc(self,itc):
          itcs = [itc for (itc,weight) in self._itcs]
          if itc not in itcs: return
          # index of this itc
          idx = itcs.index(itc)
          # remove it
          self._itcs.pop(idx)
          if self._lV0 != []: self._lV0.pop(idx)
          if self._lV0 != []: self._lpg.pop(idx)
          # remove only if exists
          if self._molecules != []: self._molecules.pop(idx)
          if itc in self._dics    : self._dics.pop(itc)
          if itc in self._dicsbw  : self._dicsbw.pop(itc)
          if itc in self._dicsfw  : self._dicsfw.pop(itc)


     #def set_from_gtsfiles(self,gtsfiles):
     #    '''
     #    returns status, string
     #    status = -1 ==> inconsistences between conformers
     #    status =  0 ==> one or several gts files do not exist
     #    status =  1 ==> everything is ok
     #    string is "" except for status = -1
     #    '''
     #    self._itcs = []
     #    self._lpg  = []
     #    lch        = []
     #    lmtp       = []
     #    limag      = []
     #    lmformu    = []
     #    self._lV0  = []
     #    if len(gtsfiles) == 0: return 0, None
     #    # save lists
     #    bool_imagfreqsOK = True
     #    for gts in gtsfiles:
     #        ctc, itc, ext = gts.split("/")[-1].split(".")
     #        self._root = ctc
     #        if not os.path.exists(gts): return 0, ""
     #        molecule = Molecule()
     #        molecule.set_from_gts(gts)
     #        molecule.setup()
     #        # save molecule
     #        self._molecules.append(molecule)
     #        # complete CTC data
     #        self._lV0.append(float(molecule._V0))
     #        self._lpg.append(str(molecule._pgroup))
     #        self._itcs.append( (itc,1) )
     #        # save some special data of this gts
     #        lch.append( int(molecule._ch) )
     #        lmtp.append( int(molecule._mtp) )
     #        limag.append( int(numimag(molecule._ccfreqs)) )
     #        lmformu.append( str(molecule._mform) )
     #        if int(numimag(molecule._ccfreqs)) not in [0,1]: bool_imagfreqsOK = False
     #    # check
     #    len1 = len(list(set(lch    )))
     #    len2 = len(list(set(lmtp   )))
     #    len3 = len(list(set(limag  )))
     #    len4 = len(list(set(lmformu)))
     #    if len1*len2*len3*len4 != 1 or not bool_imagfreqsOK:
     #       # table head and division
     #       ml1 = max([len(name) for name in lmformu ]+[7])
     #       ml2 = max([len(name) for name in gtsfiles]+[10])
     #       line_format = " %%-%is | %%-%is | %%-%is | %%-%is | %%-%is | %%-%is "%(5,6,3,7,ml1,ml2)
     #       thead = line_format%("itc","charge","mtp","n.imag.","m.form.","gts file")
     #       tdivi = "-"*len(thead)
     #       # start string
     #       string  = "Main properties of each conformer in '%s'\n"%self._ctc
     #       string += "   "+tdivi+"\n"
     #       string += "   "+thead+"\n"
     #       string += "   "+tdivi+"\n"
     #       for idx,(itc,weight) in enumerate(self._itcs):
     #           col1 = itc
     #           col2 = "%i"%lch[idx]
     #           col3 = "%i"%lmtp[idx]
     #           col4 = "%i"%limag[idx]
     #           col5 = lmformu[idx]
     #           col6 = gtsfiles[idx]
     #           ldata = (col1,col2,col3,col4,col5,col6)
     #           string += "   "+line_format%ldata+"\n"
     #       string += "   "+tdivi+"\n"
     #       return -1,string
     #    else:
     #       self._mformu = lmformu[0]
     #       self._ch     = lch[0]
     #       self._mtp    = lmtp[0]
     #       self._type   = limag[0]
     #       self._V0     = min(self._lV0)
     #       self._es     = [(self._mtp,0.0)]
     #       return 1,""

      def find_minV0(self):
          if self._lV0 == []: return
          for idx,(itc,weight) in enumerate(self._itcs):
              V0 = self._lV0[idx]
              # if first, set V0 and itcminV0
              if idx == 0:
                 self._V0       = V0
                 self._itcminV0 = itc
              # compare
              if V0 >= self._V0: continue
              self._V0       = V0
              self._itcminV0 = itc

      def get_minV0(self):
          self._lV0 = []
          self._lpg = []
          min_V0    = float("inf")
          min_itc   = None
          for idx,(itc,weight) in enumerate(self._itcs):
              gts = self.gtsfile(itc)
              if not os.path.exists(gts):
                 self._lV0 = []
                 break
                #exception = Exc.NoGTSfile(Exception)
                #exception._var = gts
                #raise exception
              # read energy and point group
              V0,pg = findingts_energy_pg(gts)
              # save data
              self._lV0.append( V0 )
              self._lpg.append( pg )
              # conformer of minimum V0?
              if V0 < min_V0: min_V0, min_itc = V0, itc
          if min_itc is not None:
             self._V0       = min_V0
             self._itcminV0 = min_itc

      def read_gtsfiles(self,dimasses={}):
          # initialize variables
          self._lpg       = []
          self._lV0       = []
          self._molecules = []
          # get list of gts files
          gtsfiles = self.gtsfiles()
          if len(gtsfiles) == 0: return
          # read and save
          for idx,gts in enumerate(gtsfiles):
              # no gts file
              if not os.path.exists(gts):
                 self._molecules.append(None)
                 self._lpg.append(None)
                 self._lV0.append(None)
                 continue
              # isotopic modification?
              itc,weight = self._itcs[idx]
              if   itc in self._diso.keys(): imods = self._diso[itc]
              elif "*" in self._diso.keys(): imods = self._diso["*"]
              else                         : imods = {}
              # create molecule
              molecule = Molecule(label=itc)
              molecule.set_from_gts(gts)
              molecule.apply_imods(imods,dimasses)
              molecule.setup()
              # save molecule
              self._molecules.append(molecule)
              # complete CTC data
              self._lV0.append(float(molecule._V0))
              self._lpg.append(str(molecule._pgroup))

      def get_piflines(self):
          string = ""
          # get min(V0)
          if self._V0 is None: self.get_minV0()
          # write itcs in CTC
          string += "root      %s\n"%self._root
          string += "# conformers\n"
          for idx,(itc,weight) in enumerate(self._itcs):
              if weight > 2: weight = 1
              if weight < 1: weight = 1
              if self._V0 is not None:
                 V0 = (self._lV0[idx] - self._V0)*KCALMOL
                 pg = self._lpg[idx]
                 string += "conformer %-4s * %i # %6.3f kcal/mol, %s\n"%(itc,weight,V0,pg)
              else:
                 string += "conformer %-4s * %i \n"%(itc,weight)
          # anharmonicity file
          if self._anh is not None:
             string += "# anharmonicity\n"
             string += "anharfile %s\n"%self._anh
          elif os.path.exists(PN.ANHDIR):
             anhfiles = []
             for anhfile in os.listdir(PN.ANHDIR):
                 if anhfile.startswith(self._ctc+"."):
                    anhfiles.append(anhfile)
             if len(anhfiles) == 1:
                string += "anharfile %s\n"%anhfiles[0]
             else:
                for anhfile in anhfiles: string += "#anharfile %s\n"%anhfile
          # basic info
          string += "# basic data\n"
          string += "mformu    %s\n"%self._mformu
          string += "ch        %s\n"%self._ch
          string += "mtp       %s\n"%self._mtp
          string += "type      %s\n"%self._type
          string += "freqscal  %.3f\n"%self._fscal
          for deg,relE in self._es:
              string += "elestate  %i  %.10E\n"%(deg,relE)
          # isotopic modifications
          if self._diso != {}:
             string += "# isotopic modifications\n"
             for itc,imods in sorted(self._diso.items()):
                 if itc == "*": string += "iso       %s\n"%("  ".join(imods))
                 else         : string += "iso.%3s   %s\n"%(itc,"  ".join(imods))
          # basic info
          for case,dics in enumerate([self._dics,self._dicsbw,self._dicsfw]):
              if dics == {}: continue
              if case == 0: string += "# internal coordinates\n"
              if case == 1: string += "# internal coordinates (backward)\n"
              if case == 2: string += "# internal coordinates (forward)\n"
              for ics_case,ctc_ics in sorted(dics.items()):
                  ics_st,ics_ab,ics_lb,ics_it,ics_pt = unmerge_ics(ctc_ics)
                  if case == 0: ss = "ics"
                  if case == 1: ss = "icsbw"
                  if case == 2: ss = "icsfw"
                  if ics_case == "*": ss += "       "
                  else              : ss += ".%3s   "%ics_case
                  if case == 0: ss += "  "
                  for idx in range(0,len(ics_st),5):
                      string += ss+ "  ".join( ["-".join(["%i"%(at+1) for at in atoms])\
                                                for atoms in ics_st[idx:idx+5]]) + "\n"
                  for idx in range(0,len(ics_ab),5):
                      string += ss+ "  ".join( ["-".join(["%i"%(at+1) for at in atoms])\
                                                for atoms in ics_ab[idx:idx+5]]) + "\n"
                  for idx in range(0,len(ics_lb),5):
                      string += ss+ "  ".join( ["=".join(["%i"%(at+1) for at in atoms])\
                                                for atoms in ics_lb[idx:idx+5]]) + "\n"
                  for idx in range(0,len(ics_it),5):
                      string += ss+ "  ".join( ["_".join(["%i"%(at+1) for at in atoms])\
                                                for atoms in ics_it[idx:idx+5]]) + "\n"
                  for idx in range(0,len(ics_pt),5):
                      string += ss+ "  ".join( ["-".join(["%i"%(at+1) for at in atoms])\
                                               for atoms in ics_pt[idx:idx+5]]) + "\n"
          return string

