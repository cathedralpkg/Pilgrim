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
| Sub-module :  optGATHER          |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#===============================================================#
import os
import sys
#---------------------------------------------------------------#
import modpilgrim.names as PN
import modpilgrim.pilrw  as RW
#---------------------------------------------------------------#
from   modpilgrim.ClusterConf       import ClusterConf
#---------------------------------------------------------------#
import common.Exceptions as     Exc
#---------------------------------------------------------------#
from   common.dicts      import dpt_im
from   common.orca       import read_orca
from   common.gaussian   import read_fchk
from   common.gaussian   import read_gauout     as read_gauout_v1
from   common.gaussian   import read_gauout_old as read_gauout_v2
from   common.files      import read_gtsfile
from   common.files      import write_gtsfile
from   common.files      import write_molden
from   common.fncs       import numimag
from   common.fncs       import symbols_and_atonums
from   common.fncs       import clean_dummies
from   common.fncs       import is_string_valid
from   common.fncs       import fill_string
from   common.fncs       import print_string
from   common.Molecule   import Molecule
from   common.pgs        import get_pgs
#===============================================================#


#===============================================================#
def gtsname(gtsfile,case="full"):
    # name without path
    gtsfile_name = gtsfile.split("/")[-1]
    gtsfile_full = PN.DIR1+gtsfile.split("/")[-1]
    if   case == "full": return gtsfile_full
    elif case == "name": return gtsfile_name
    else               : return None
#---------------------------------------------------------------#
def known_files(files):
    allowed_ext = ["log","out","fchk","gts"]
    files = [ifile for ifile in files if ifile.split(".")[-1].lower() in allowed_ext]
    return files
#---------------------------------------------------------------#
def get_gtsfiles_from_dir(ctc):
    return sorted([PN.DIR1+gts for gts in os.listdir(PN.DIR1) \
                   if gts.startswith(ctc+".") and gts.endswith(".gts")])
#---------------------------------------------------------------#
def userfile_to_gtsfile(filename,gtsfile):
    '''
    Read a file given by the user and generate gts file
    if possible
    '''
    gtsfile = gtsname(gtsfile,"full")
    read_methods = []
    read_methods.append(read_gtsfile)    #(a) a gts file
    read_methods.append(read_fchk   )    #(b) a fchk file
    read_methods.append(read_gauout_v1 ) #(c.1) a Gaussian output file
    read_methods.append(read_gauout_v2 ) #(c.2) a Gaussian output file
    read_methods.append(read_orca   )    #(d) an Orca output file
    for read_method in read_methods:
        try:
          xcc,atonums,ch,mtp,E,gcc,Fcc,masses,clevel = read_method(filename)[0:9]
          # symbols and atonums
          symbols,atonums = symbols_and_atonums(atonums)
          # remove dummy atoms
          symbols_wo,xcc = clean_dummies(symbols,xcc=xcc)
          masses = clean_dummies(symbols,masses=masses)[1]
          if gcc is not None: gcc = clean_dummies(symbols,gcc=gcc)[1]
          if Fcc is not None: Fcc = clean_dummies(symbols,Fcc=Fcc)[1]
          # update atonums to remove dummies
          symbols,atonums = symbols_and_atonums(symbols_wo)
          # clevel
          if read_method == read_gtsfile: clevel = ""
          # in case no data in masses
          if masses is None or len(masses) == 0 or sum(masses) == 0.0:
             masses = atonums2masses(atonums)
          # Some checking
          if len(atonums) != 1 and Fcc in ([],None):
             #raise Exc.FccNotFound
             return -1, None
          # Generate Molecule instance
          molecule = Molecule()
          molecule.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
          molecule.setvar(atonums=atonums,masses=masses)
          molecule.setvar(ch=ch,mtp=mtp,V0=E)
          molecule.prepare()
          # Deal with frozen atoms
          ffrozen = gtsfile+".frozen"
          frozen_xcc, frozen_symbols = molecule.remove_frozen()
          RW.write_frozen(ffrozen,frozen_xcc,frozen_symbols)
          # write gts
          molecule.calc_pgroup(force=True)
          molecule.genfile_gts(gtsfile,level=clevel)
          return 1, E
        except Exc.FileType:
          continue
        except:
          continue
    return 0, None
#---------------------------------------------------------------#
def deal_with_tracking():
    dtrack,(fname,status) = RW.read_track()
    if status == 1: print("  - File '%s' exists and is not empty\n"%fname)
    topop = []
    ml1 = max([len(k) for k in dtrack.keys()  ]+[1])
    ml2 = max([len(v) for v in dtrack.values()]+[1])
    for outfile,gtsfile in sorted(dtrack.items()):
        gtsfile_full = gtsname(gtsfile,"full")
        gtsfile_name = gtsname(gtsfile,"name")
        exist1 = os.path.exists(PN.UFOLDER+outfile)
        exist2 = os.path.exists(gtsfile_full)
        if not exist1 and exist2:
           print("       * File '%s' does not exists, but '%s' does"%(outfile,gtsfile_name))
           answer = input("         remove gts file '%s' (y/N)? "%(gtsfile_full))
           if answer.strip().lower() in ["y","yes"]:
              os.remove(gtsfile_full)
        elif exist1 and not exist2:
           print("       * File '%s' exists, but '%s' does not..."%(outfile,gtsfile_name))
           answer = input("         Should Pilgrim create the gts file (y/N)? ")
           if answer.strip().lower() in ["y","yes"]:
              status,E = userfile_to_gtsfile(PN.UFOLDER+outfile,gtsfile_full)
              if   status ==  1: print("         Created!")
              else             : print("         Creation failed!")
        elif not exist1 and not exist2:
           print("       * Neither '%s' nor '%s' exist! Removing from '%s'"%(outfile,gtsfile_name,fname))
           topop.append(outfile)
        else:
           print("       * %s ==> %s"%(PN.UFOLDER+"%%-%is"%ml1%outfile,PN.DIR1+"%%-%is"%ml2%gtsfile_name))
    print("")
    # remove lines in tracking
    for key in topop:
        if key in dtrack.keys(): dtrack.pop(key)
    # return dictionary
    return dtrack
#---------------------------------------------------------------#
def deal_with_pif():
    (dctc,dimasses),(fname,status) = RW.read_ctc()
    if status == 1 and dctc !={}:
       print("  - File '%s' exists and is not empty"%fname)
       ls_struc(dctc)
       print("       Checking existency of gts files...")
       print("")
       notexist  = {}
       notlisted = {}
       for ctc in sorted(dctc.keys()):
           gtsfilesA = [gtsname(gts,"name") for gts in dctc[ctc].gtsfiles()]
           gtsfilesB = [gtsname(gts,"name") for gts in get_gtsfiles_from_dir(dctc[ctc]._root)]
           set_notexist  = set(gtsfilesA).difference(gtsfilesB)
           set_notlisted = set(gtsfilesB).difference(gtsfilesA)
           if len(set_notexist)  != 0: notexist[ctc]  = set_notexist
           if len(set_notlisted) != 0: notlisted[ctc] = set_notlisted
       if notexist != {}:
          print("          * Some files should exist according to '%s', but they do not."%PN.IFILE1)
          print("            Actions:")
          print("              (0) do nothing [default]")
          print("              (1) remove from %s"%PN.IFILE1)
          print("            Files:")
          for ctc,gtslist in notexist.items():
              for gts in gtslist:
                  answer = input("            - %s (%s); action? "%(gts,ctc)).strip().lower()
                  if answer == "1":
                     target_itc = gts.split(".")[1]
                     itcs = dctc[ctc]._itcs
                     itcs = [(itc,weight) for (itc,weight) in itcs if itc != target_itc]
                     dctc[ctc]._itcs = itcs
                     print("              removed!")
       else:
          print("          * All gts files listed in '%s' exist in '%s'"%(PN.IFILE1,PN.DIR1))
       print("")
       if notlisted != {}:
          print("          * Some gts files are not listed in '%s' but they exist."%PN.IFILE1)
          print("            Actions:")
          print("              (0) do nothing [default]")
          print("              (1) add to %s"%PN.IFILE1)
          print("              (2) remove file")
          print("            Files:")
          for ctc,gtslist in notlisted.items():
              for gts in gtslist:
                  answer = input("              - %s ; action? "%gts).strip().lower()
                  if answer == "1":
                     itc = gts.split("/")[-1].split(".")[1]
                     dctc[ctc]._itcs.append( (itc,1) )
                     print("                added!")
                  if answer == "2":
                     os.remove(gtsname(gts,"full"))
                     print("                removed!")
       else: 
          print("          * All gts files in '%s' are listed in '%s'"%(PN.DIR1,PN.IFILE1))
       print("")
    # return data
    return dctc, dimasses
#---------------------------------------------------------------#
def dctc_from_DIR1(gtsfiles):
    # list of gts files
    print("  - File '%s' does not exist but '%s' is not empty!"%(PN.IFILE1,PN.DIR1))
    print("")
    # classify gts files according to ctc
    dctc = {}
    for gtsfile in gtsfiles:
        gts_name = gtsname(gtsfile,"name")
        gts_full = gtsname(gtsfile,"full")
        ctc, itc = gts_name.split(".")[0:2]
        dctc[ctc] = dctc.get(ctc,[]) + [gts_full]
    # generate real dctc
    print("     number of structures: %i"%len(dctc.keys()))
    print("")
    ml = max([len(key) for key in dctc.keys()]+[1])
    for ctc,gtsfiles in dctc.items():
        print("       * %%-%is (num gts files: %%2i)"%ml%(ctc,len(gtsfiles)))
        CTC = ClusterConf(ctc)
        status, string = CTC.set_from_gtsfiles(gtsfiles)
        if status == -1:
           print("         INCONSISTENCES!")
           print_string(string,9)
           answer = input("         Remove gts files (y/N)?")
           if answer.strip().lower() in ["y","yes"]:
              for gts in gtsfiles: os.remove(gts)
           else: raise Exc.ABORTED
        dctc[ctc] = CTC
    print("")
    return dctc
#---------------------------------------------------------------#
def convert_udata_to_gts(ff,dtrack,dctc):
    global wrong1
    global wrong2
    global wrong5
    inconsistence = None
    print(   "         |--> %s"%ff)
    iblank = "         |    "
    # Is ff a folder or a file?
    if os.path.isdir(PN.UFOLDER+ff):
        # Get CTC name & list of files
        ctc = ff[:-1]
        files = sorted([ff+filename for filename in os.listdir(PN.UFOLDER+ff)])
    else:
        # Get CTC name & list of files
        ctc = ff.split(".")[0]
        files = [ff]

    # valid name??
    if not is_string_valid(ctc,extra="_"):
        print(iblank+"invalid name!")
        print(iblank)
        wrong1   = True
        return dtrack, dctc, inconsistence

    # only files with known extension
    files = known_files(files)
    if len(files) == 0:
       print(iblank+"no valid file(s) or empty folder")
       print(iblank)
       wrong5 = True
       return dtrack, dctc, inconsistence

    in_dtrack = [ff in dtrack.keys() for ff in files].count(False) == 0
    if ctc in dctc.keys() and in_dtrack:
       print(iblank+"already registered")
       print(iblank)
       return dtrack, dctc, inconsistence

    # Convert to gts (provisional name)
    count = 1
    gtslist = []
    for ifile in files:
        # already in dtrack
        if ifile in dtrack.keys():
           print(iblank+"|--> %s (already in %s)"%(ifile.split("/")[-1],PN.IFILE0))
           continue
        # provisional gts name
        while True:
              #prov_gts = PN.DIR1+"%s.gts_%i"%(ctc,count)
              prov_gts = PN.DIR1+"%s.prov%i.gts"%(ctc,count)
              if not os.path.exists(prov_gts): break
              count +=1
        # read file and create gts
        status,E = userfile_to_gtsfile(PN.UFOLDER+ifile,prov_gts)
        if status == 0:
           print(iblank+"|--> %s (sth failed!)"%(ifile.split("/")[-1]))
           wrong2 = True
        elif status == -1:
           print(iblank+"|--> %s (Force Constant Matrix NOT FOUND!)"%(ifile.split("/")[-1]))
           wrong2 = True
        else: 
           print(iblank+"|--> %s (new)"%(ifile.split("/")[-1]))
        if E is not None: gtslist.append( (E,ifile,prov_gts) )
    print(iblank)

    # sort by energy (only those created)
    if len(gtslist) == 0: return dtrack, dctc, inconsistence
    gtslist.sort()

    # Rename gts files and update dtrack
    info, idx = [], 1
    generated_gtsfiles = []
    for E,ifile,prov_gts in gtslist:
        while True:
              gtsfile = ctc+".%003i.gts"%(idx)
              idx += 1
              if gtsfile in dtrack.values(): continue
              if not os.path.exists(PN.DIR1+gtsfile): break 
        os.rename(prov_gts,PN.DIR1+gtsfile)
        if os.path.exists(prov_gts+".frozen"):
           os.rename(prov_gts+".frozen",PN.DIR1+gtsfile+".frozen")
        info.append( (E,ifile,gtsfile) )
        dtrack[ifile] = gtsfile
        generated_gtsfiles.append(gtsfile)

    # consider all gts files and prepare CTC
    CTC,status,string,gtsfiles = gen_and_check_ctc(ctc,dtrack)
    # save ctc in dctc
    if status == 1: dctc[ctc] = CTC
    # save inconsistence
    if status == -1: inconsistence = (string,generated_gtsfiles)
    return dtrack, dctc, inconsistence
#---------------------------------------------------------------#
def gts2molecule(gtsfile):
    # name of file (without path to folder)
    gtsfile_name = gtsname(gtsfile,"name")
    gtsfile_full = gtsname(gtsfile,"full")
    # Get ctc and itc
    ctc,itc,ext = gtsfile_name.split(".")
    # read gts and prepare Molecule
    if not os.path.exists(gtsfile_full): return None,(ctc,itc)
    molecule = Molecule()
    molecule.set_from_gts(gtsfile_full)
    # setup
    molecule.setup()
    # return data
    return molecule,(ctc,itc)
#---------------------------------------------------------------#
def molecule2molden(molecule,ctc,itc):
    # gn t name for molden file
    molden  = PN.get_gtsmolden(ctc,itc)
   ## does molden file exist?
   #if os.path.exists(molden): return
    # create molden file
    idata = (molden,molecule._xcc,molecule._symbols,molecule._ccfreqs,molecule._ccFevecs)
    write_molden(*idata)
#---------------------------------------------------------------#
def gen_and_check_ctc(ctc,dtrack,gtsfiles=None):
    #--------------------#
    # Read all gts files #
    #--------------------#
    # if gtsfiles not given, get them from folder
    if gtsfiles is None: gtsfiles = get_gtsfiles_from_dir(ctc)
    # initialize
    created_molden = []
    # generate cluster
    CTC = ClusterConf(ctc)
    status, string = CTC.set_from_gtsfiles(gtsfiles)
    for idx,molecule in enumerate(CTC._molecules):
        itc,weigh = CTC._itcs[idx]
        molecule2molden(molecule,ctc,itc)
    return CTC, status, string, gtsfiles
#---------------------------------------------------------------#
def ls_struc(dctc):
    '''
    extra option to show what's in .ctc file
    to show: ls ctc
    '''
    ib = "       "
    htable = ["species name","m.form.","num.ifreqs.","ch","mtp","num.confs.","iso.mod."]
    # some numbers and list for nice formatting
    ml1 = max([len(ctc)         for ctc in dctc.keys()  ]+[len(htable[0])])
    ml2 = max([len(CTC._mformu) for CTC in dctc.values()]+[len(htable[1])])
    ll  = [ml1,ml2,len(htable[2]),3,3,len(htable[5]),len(htable[6])]
    # fill strings in htable
    htable = [fill_string(string,ml) for ml,string in zip(ll,htable)]
    # start loop
    htable  = " "+" | ".join(htable)+" "
    divis   = "-"*len(htable)
    stable  = "\n"
    stable += ib+divis+"\n"
    stable += ib+htable+"\n"
    stable += ib+divis+"\n"
    # separate minima from saddle point and from other (xx)
    ctc_root_min = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 0 and dctc[ctc]._diso=={}])
    ctc_root_ts  = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 1 and dctc[ctc]._diso=={}])
    ctc_isoX_min = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 0 and dctc[ctc]._diso!={}])
    ctc_isoX_ts  = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 1 and dctc[ctc]._diso!={}])
    ctc_xx  = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type not in [0,1]])
    
    numSP0, numSP0all = 0, 0
    numSP1, numSP1all = 0, 0
    numSPX, numSPXall = 0, 0
    for idx,ctcs in enumerate([ctc_root_min+ctc_isoX_min,ctc_root_ts+ctc_isoX_ts,ctc_xx]):
        if ctcs == []: continue
        for ctc in ctcs:
            mformu = dctc[ctc]._mformu
            itcs   = dctc[ctc]._itcs
            ch     = dctc[ctc]._ch
            mtp    = dctc[ctc]._mtp
            sptype = dctc[ctc]._type
            diso   = dctc[ctc]._diso
            if diso == {}:
               isostr = "none"
            else:
               try   : isostr = " ".join(diso["*"])
               except: isostr = "yes"
            # count conformers
            num1  = len(itcs)
            num2  = int(sum([weight for itc,weight in itcs]))
            # variable to string
            ncnfs = "%i (%i)"%(num2,num1)
            ch    = "%i"%ch
            mtp   = "%i"%mtp
            nif   = "%i"%sptype
            # fill strings
            ltable = ["%%-%is"%ml1%ctc,"%%-%is"%ml2%mformu,nif,ch,mtp,ncnfs,isostr]
            ltable = [fill_string(string,ml) for ml,string in zip(ll,ltable)]
            # add to table
            ltable  = " "+" | ".join(ltable)+" "
            stable += ib+ltable+"\n"
            if   sptype == 0: numSP0 += num1; numSP0all += num2
            elif sptype == 1: numSP1 += num1; numSP1all += num2
            else            : numSPX += num1; numSPXall += num2
        stable += ib+divis+"\n"
    stable += ib+"  * num(minimum) = %i (%i)\n"%(numSP0all,numSP0)
    stable += ib+"  * num(saddle ) = %i (%i)\n"%(numSP1all,numSP1)
    print(stable)
#===============================================================#



#===============================================================#
def main(idata,status):

    global wrong1
    global wrong2
    global wrong3
    global wrong4
    global wrong5
    wrong1      = False
    wrong2      = False
    wrong3      = False
    wrong4      = False
    wrong5      = False

    #-------------------#
    # Deal with folders #
    #-------------------#
    print("  - Folders of interest:")
    print("")
    print("       * %-9s ==> folder of user's input data"%PN.UFOLDER)
    print("       * %-9s ==> folder with gts files"%PN.DIR1)
    print("       * %-9s ==> folder with molden files"%PN.DIR5)
    print("")
    print("  - Existence of folders of interest:")
    print("")
    if os.path.exists(PN.UFOLDER): print("       * %-9s exists"%PN.UFOLDER)
    else                         : print("       * %-9s does not exist"%PN.UFOLDER)
    if os.path.exists(PN.DIR1   ): print("       * %-9s exists"%PN.DIR1)
    else                         : print("       * %-9s does not exist"%PN.DIR1)
    if os.path.exists(PN.DIR5   ): print("       * %-9s exists"%PN.DIR5)
    else                         : print("       * %-9s does not exist"%PN.DIR5)
    print("")

    # Data in folders
    if os.path.exists(PN.UFOLDER):
       files   = [ff     for ff in os.listdir(PN.UFOLDER) if not os.path.isdir(PN.UFOLDER+ff)]
       folders = [ff+"/" for ff in os.listdir(PN.UFOLDER) if     os.path.isdir(PN.UFOLDER+ff)]
       files   = known_files(files) # only those with known extensions
    else:
       files   = []
       folders = []

    if os.path.exists(PN.DIR1):
       gtsfiles = [gts for gts in os.listdir(PN.DIR1) if gts.endswith(".gts")]
    else:
       gtsfiles = []


    #--------------------------------#
    # Basic cases to abort execution #
    #--------------------------------#

    # NO UFOLDER, NO DIR1
    if not os.path.exists(PN.UFOLDER) and not os.path.exists(PN.DIR1):
       print("  - Neither %s nor %s exist!"%(PN.UFOLDER,PN.DIR1))
       raise Exc.ABORTED
    # NO UFOLDER, EMPTY DIR1
    if not os.path.exists(PN.UFOLDER) and len(gtsfiles) == 0:
       print("  - %s does not exist and %s contains no data!"%(PN.UFOLDER,PN.DIR1))
       print("")
       raise Exc.ABORTED
    # EMPTY UFOLDER, NO DIR1
    if len(files+folders) == 0 and not os.path.exists(PN.DIR1):
       print("  - %s contains no data and %s does not exist!"%(PN.UFOLDER,PN.DIR1))
       print("")
       raise Exc.ABORTED
    # EMPTY UFOLDER, EMPTY DIR1
    if len(files+folders) == 0 and len(gtsfiles) == 0:
      print("  - %s contains no data and %s contains no data either!"%(PN.UFOLDER,PN.DIR1))
      print("")
      raise Exc.ABORTED


    #----------------#
    # Create folders #
    #----------------#
    if not os.path.exists(PN.DIR1   ): os.mkdir(PN.DIR1   )
    if not os.path.exists(PN.DIR5   ): os.mkdir(PN.DIR5   )


    #--------------------#
    # DEAL WITH tracking #
    #--------------------#
    dtrack = deal_with_tracking()
    # All removed except UDATA and tracking --> gtsfiles need to be listed again
    if os.path.exists(PN.DIR1):
       gtsfiles = [gts for gts in os.listdir(PN.DIR1) if gts.endswith(".gts")]

    #---------------------------------#
    # DEAL WITH pif.struc or 1-GTS/   #
    #---------------------------------#
    if os.path.exists(PN.IFILE1):
       dctc, dimasses = deal_with_pif()
    elif len(gtsfiles) != 0:
       dctc     = dctc_from_DIR1(gtsfiles)
       dimasses = {}
    else:
       dctc, dimasses = {},  {}
    # isotopic masses
    if dimasses == {}: dimasses = dpt_im

    #-----------------------#
    # User data is analyzed #
    #-----------------------#
    if files+folders != []:
       print("  - Going through user's input data (%s) to create gts files"%PN.UFOLDER)
       print("")
       print("      number of folders = %i "%(len(folders)))
       print("      number of files   = %i (outside folders)"%(len(files)))
       print("")

       print("      %s"%(PN.UFOLDER))
       inconsistences = []
       for ff in sorted(folders)+sorted(files):
           dtrack,dctc,inconsistence = convert_udata_to_gts(ff,dtrack,dctc)
           if inconsistence is not None: inconsistences.append(inconsistence)
       print("")
    else:
       print("  - No data inside '%s'"%PN.UFOLDER)
       print("")
       inconsistences = []


    #-------------------------------#
    # Conformational inconsistences #
    #-------------------------------#
    if len(inconsistences) != 0:
       print("  - There are some inconsistences!")
       print("")
       wrong4 = True
       dtrack2 = {v:k for k,v in dtrack.items()}
       for string,gtsfiles in inconsistences:
           print_string(string,7)
           for gts in gtsfiles:
               gts_name = gtsname(gts,"name")
               gts_full = gtsname(gts,"full")
               ifile    = dtrack2[gts_name]
               print("          removing '%s' (from %s)..."%(gts_full,dtrack2[gts_name]))
               dtrack.pop(ifile)
               os.remove(gts_full)
           print("")
       print("")


    #------------------------#
    # Final checking of dctc #
    #------------------------#
    print("  - Checking final data for %s"%PN.IFILE1)
    ls_struc(dctc)
    for ctc,CTC in sorted(dctc.items()):
        gtsfiles = CTC.gtsfiles()
        CTC, status, string, gtsfiles = gen_and_check_ctc(ctc,dtrack,gtsfiles)
        if status == -1:
           wrong4 = True
           print_string(string,7)
           print("          Inconsistence found!")
           answer = input("          Remove gts files (y/N)? ")
           if answer.strip().lower() == "y":
               for gts in gtsfiles:
                   print("             removing %s"%gts)
                   os.remove(gts)
           print("")
        if status == 0:
           print("       ERROR: Unable to find the following gts file(s):")
           print("")
           wrong3 = True
           for gts in gtsfiles:
               if os.path.exists(gts): continue
               print("         * %s (%s)"%(gts,ctc))
           print("")


    #----------------------------------#
    # (Re)Write tracking and pif files #
    #----------------------------------#
    if dtrack != {}:
       print("  - Writing/Updating information in '%s'"%PN.IFILE0)
       RW.write_track(dtrack)
    if dctc != {}:
       print("  - Writing/Updating information in '%s'"%PN.IFILE1)
       RW.write_ctc(dctc,dimasses)


    #----------------#
    # Print warnings #
    #----------------#
    if wrong1 or wrong2 or wrong3 or wrong4 or wrong5:
       print("")
       print("  - WARNING!! It seems that something did not go well")
       if wrong1: print("       * invalid name for some file/folder(s)!")
       if wrong2: print("       * reading process failed for a/some file(s) in %s!"%PN.UFOLDER)
       if wrong3: print("       * gts file(s) not found!")
       if wrong4: print("       * inconsistences in conformational cluster(s)!")
       if wrong5: print("       * some folder(s) in %s is(are) empty!"%PN.UFOLDER)
    else:
       print("")
       print("  - Everything seems to be OK")
#===============================================================#




