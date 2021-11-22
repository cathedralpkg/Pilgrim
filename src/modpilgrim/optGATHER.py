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
| Sub-module :  optGATHER          |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#===============================================================#
import os
#---------------------------------------------------------------#
import common.Exceptions as     Exc
import common.fncs       as     fncs
#---------------------------------------------------------------#
from   common.dicts      import dpt_im
from   common.files      import write_gtsfile
from   common.files      import write_molden
from   common.Molecule   import Molecule
from   common.pgs        import get_pgs
#---------------------------------------------------------------#
from   common.files      import read_gtsfile
from   common.gaussian   import read_fchk
from   common.gaussian   import read_gauout     as read_gauout_v1
from   common.gaussian   import read_gauout_old as read_gauout_v2
from   common.orca       import read_orca
#---------------------------------------------------------------#
READ_METHODS = [read_gtsfile,read_fchk,read_gauout_v1,read_gauout_v2,read_orca]
EXTENSIONS   = ["gts","fchk","log","out"]
#---------------------------------------------------------------#
import modpilgrim.names       as     PN
import modpilgrim.pilrw       as     RW
from   modpilgrim.ClusterConf import ClusterConf
#===============================================================#

# WRONG: 1 (invalid name)
#        2 (folder without valid files)
#        3 (reading failed)
#        4 (Fcc not found)
#        5 (inconsistences CTC)
#        6 (unexpected number of imag. freqs)
WRONG = {i+1:False for i in range(6)}



#===============================================================#
def get_ufiles_ufolders():
    '''
    lists files in UDATA/ folder;
    returns a list with the folders [ufolder],
            a list with the files that are not inside a folder [ufiles],
            and a list with all files (i.e. files inside each folder + ufiles) [all_ufiles]
    + 
    '''
    # list files and folders
    ufiles   = [ff     for ff in os.listdir(PN.UFOLDER) if not os.path.isdir(PN.UFOLDER+ff)]
    ufolders = [ff+"/" for ff in os.listdir(PN.UFOLDER) if     os.path.isdir(PN.UFOLDER+ff)]
    # only files with correct extension
    ufiles   = [ifile for ifile in ufiles if ifile.split(".")[-1].lower() in EXTENSIONS]
    # sort
    ufiles.sort()
    ufolders.sort()
    # all the files (including those in folders)
    all_ufiles = list(ufiles)
    for ufolder in ufolders:
        all_ufiles += sorted([ufolder+esfile for esfile in os.listdir(PN.UFOLDER+ufolder) \
                              if  esfile.split(".")[-1].lower() in EXTENSIONS])
    return ufiles,ufolders,all_ufiles
#---------------------------------------------------------------#
def remove_definitive(lgts):
    '''
    Given a list of gts files, remove each file and the corresponding files associated
    (i.e. the molden file and the frozen file, in case any of them exists)
    * molden files can be where the gts files are or they can be at 5-MOLDEN
    * frozen files can only be at the same directory as the gts file
    '''
    toremove  = [PN.DIR1+gts                           for gts in lgts if gts is not None]
    toremove += [PN.DIR1+gts+".frozen"                 for gts in lgts if gts is not None]
    toremove += [PN.DIR1+gts+".molden"                 for gts in lgts if gts is not None]
    toremove += [PN.DIR5+gts.replace(".gts",".molden") for gts in lgts if gts is not None]
    for fname in toremove:
        if os.path.exists(fname): os.remove(fname)
#---------------------------------------------------------------#
def deal_removed_ufiles(all_ufiles,dtrack,dctc):
    '''
    compares files in UDATA with those listed in tracking.
    If any file in tracking is not in UDATA/, this function is in charge of
    remove the corresponding associated information (its line in tracking,
    gts file, the conformer in pif.struc)
    '''
    # any file removed by user?
    removed = []
    for file_in_tracking in dtrack.keys():
        if file_in_tracking not in all_ufiles:
            removed.append(file_in_tracking)
    # deal with removed files
    if len(removed) != 0:
       print("   --> The following files were removed from %s by user:"%PN.UFOLDER)
       print("")
       ii = len(str(len(removed)))
       for idx,fname in enumerate(removed):
           idx = "%%%ii"%ii%idx
           print("       [%s] %s"%(idx,fname))
       print("")
       print("       To continue, the corresponding gts files (in %s) must be removed"%PN.DIR1)
       answer = input("       remove (y/N)? ").strip().lower()
       print("")
       # if no permission to remove, abort execution
       if answer not in ("yes","y"): raise Exc.ABORTED
       # deal with one file at a time
       for idx,fname in enumerate(removed):
           # remove from dtrack
           gts = dtrack.pop(fname)
           # ctc and itc from gts name
           ctc,itc,ext= gts.split(".")
           # remove from dctc
           if ctc in dctc.keys():
              dctc[ctc].remove_itc(itc)
              # update referene energy for this CTC
              dctc[ctc].find_minV0()
           # remove gts
           idx = "%%%ii"%ii%idx
           print("       [%s] %s"%(idx,PN.DIR1+gts))
           remove_definitive([gts])
       print("")
    # return new dicts
    return dtrack,dctc
#---------------------------------------------------------------#
def deal_nonlisted_gts(dtrack):
    '''
    remove all gts files (and associated files) that exist
    but that are not listed in tracking
    '''
    # gts files
    gtsfiles = []
    if os.path.exists(PN.DIR1):
       gtsfiles = [gts for gts in os.listdir(PN.DIR1) if gts.endswith(".gts")]
    # check which ones are non-listed (nl) in dtrack
    nl_gtsfiles = []
    for gts in gtsfiles:
        if gts not in dtrack.values(): nl_gtsfiles.append(gts)
    # remove files
    if len(nl_gtsfiles) != 0:
       print("   --> The following files are not listed in %s"%PN.IFILE0)
       print("")
       for gts in nl_gtsfiles: print("       * %s"%gts)
       print("")
       print("       To continue, these files must be removed")
       answer = input("       remove (y/N)? ").strip().lower()
       print("")
       if answer in ("yes","y"): remove_definitive(nl_gtsfiles)
       else: raise Exc.ABORTED
#---------------------------------------------------------------#
def cleanup_dctc(dctc):
    '''
    removes conformers in dctc whose gts file does not exist
    '''
    for ctc in dctc.keys():
        remove = []
        for itc,weight in dctc[ctc]._itcs:
            # the corresponding gts file
            gts = "%s.%s.gts"%(ctc,itc)
            if not os.path.exists(PN.DIR1+gts): remove.append(itc)
        # remove conformers
        for itc in remove: dctc[ctc].remove_itc(itc)
        # redefine minimum energy using conformers we kept
        if remove != []  : dctc[ctc].find_minV0()
    # remove CTCs without conformers
    dctc = {ctc:CTC for ctc,CTC in dctc.items() if len(CTC._itcs) != 0}
    return dctc
#===============================================================#


#===============================================================#
def temporal_gtsname(ctc,count=0):
    '''
    gives a temporary name for the gts file in such a way
    that the gts file does not exist yet
    '''
    while True:
      prov_gts = "%s.prov%i.gts"%(ctc,count)
      count   +=1
      if os.path.exists(PN.DIR1+prov_gts): continue
      return prov_gts,count
#---------------------------------------------------------------#
def print_ctc_info(CTC,iblank="      "):
    '''
    print main information of the given CTC
    '''
    string  = "ctc name         : %s\n"%CTC._ctc
    string += "charge           : %i\n"%CTC._ch
    string += "spin multiplicity: %i\n"%CTC._mtp
    string += "molecular formula: %s\n"%CTC._mformu
    string += "num. imag. freqs.: %i\n"%CTC._type
    string += "num. conformers  : %i\n"%len(CTC._itcs)
    for line in string.split("\n"): print(iblank+line)
#---------------------------------------------------------------#
def generate_gts_file(esfile,gtsfile,read_method):
    '''
    Generate gts file (and molden & frozen files)
    from given electronic structure (ES) file;
    The ES file is read using read_method
    '''
    # Extra files that (maybe) will be generated
    file_frozen = gtsfile+".frozen"       # only if any frozen atom
    file_molden = gtsfile+".molden"       # molden file
    # read file
    xcc,atonums,ch,mtp,E,gcc,Fcc,masses,clevel = read_method(PN.UFOLDER+esfile)[0:9]
    # clevel is not really clevel when using read_gtsfile as read_method
    if read_method == read_gtsfile: clevel = ""
    # symbols and atonums
    symbols,atonums = fncs.symbols_and_atonums(atonums)
    # is masses available?
    if masses is None or len(masses) == 0 or sum(masses) == 0.0:
       masses = atonums2masses(atonums)
    # is Fcc available?
    if len(atonums) != 1 and (Fcc is None or len(Fcc) == 0):
       status   = -1
       cache    = None
       molec    = None
    else:
       # Generate Molecule instance
       molec = Molecule()
       molec.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
       molec.setvar(atonums=atonums,masses=masses)
       molec.setvar(ch=ch,mtp=mtp,V0=E)
       molec.prepare()
       # Deal with frozen atoms [atom i is frozen if Fij=0 forall j)
       frozen_xcc, frozen_symbols = molec.remove_frozen()
       # Calculate point group (must be done after remove frozen)
       molec.calc_pgroup(force=True)
       # Deal with Hessian
       molec.setup()
       # write gts file
       molec.genfile_gts(PN.DIR1+gtsfile,level=clevel)
       status = 1
       # write frozen (if there are frozen atoms)
       RW.write_frozen(PN.DIR1+file_frozen,frozen_xcc,frozen_symbols)
       # write molden file
       idata = (PN.DIR1+file_molden,molec._xcc,molec._symbols,molec._ccfreqs,molec._ccFevecs)
       write_molden(*idata)
       # save to cache
       nimag  = int(fncs.numimag(molec._ccfreqs))
       pgroup = str(molec._pgroup)
       mform  = str(molec._mform)
       cache = [esfile,gtsfile,E,ch,mtp,nimag,mform,pgroup]
    # delete variable, just in case
    del xcc, gcc, Fcc
    del symbols, atonums, masses
    del molec
    # return
    return status, cache
#---------------------------------------------------------------#
def userfile_to_gtsfile(esfile,gtsfile):
    '''
    Read the electronic structure (ES) file given by
    the user and generate the corresponding gts file;
    This function try with different methods to read the ES file

    returns
        * status: -1 (Fcc not found), 0 (not read) or 1 (all ok)
        * cache : a tuple with data if status == 1, else None
    '''
    for read_method in READ_METHODS:
        try                : return generate_gts_file(esfile,gtsfile,read_method)
        except Exc.FileType: continue
        except             : continue
    return 0, None
#---------------------------------------------------------------#
def print_table_ufiles(cache):
    '''
    print table with basic info about each electronic structure file
    associated to a given CTC;
    useful to see if there is any incompatibility between conformers
    (e.g. conformers with different charge, or different spin multiplicity)
    '''
    # table head and division
    ml1 = max([len(cache_i[6]) for cache_i in cache]+[ 7])
    ml2 = 21
    line_format = " %%-%is | %%-%is | %%-%is | %%-%is | %%-%is "%(6,3,7,ml1,ml2)
    thead = line_format%("charge","mtp","n.imag.","m.form.","user file")
    tdivi = "-"*len(thead)
    # start string
    string  = "\n"
    string += tdivi+"\n"
    string += thead+"\n"
    string += tdivi+"\n"
    for idx,cache_i in enumerate(cache):
        col1 = "%i"%cache_i[3] # charge
        col2 = "%i"%cache_i[4] # spin multiplicity
        col3 = "%i"%cache_i[5] # number of imaginary frequencies
        col4 = cache_i[6]      # molecular formula
        col5 = cache_i[0]      # name of ES file
        # format col5
        col5 = col5.split("/")[-1]
        if len(col5) > ml2: col5 = col5[:(ml2-3)//2]+"..."+col5[-(ml2-3)//2:]
        # add to table
        ldata = (col1,col2,col3,col4,col5)
        string += line_format%ldata+"\n"
    string += tdivi+"\n"
    # print string in prompt
    for line in string.split("\n"): print("    "+line)
#---------------------------------------------------------------#
def rename_files(ctc,prov_gts,iitc=0,gts=None):
    '''
    rename provisional gts files (also molden and frozen files)
    to definitive name
    '''
    assert iitc >= 0
    # decide new name for gts file
    if gts is None:
       itc = int(iitc)
       while True:
           itc += 1
           new_gts = "%s.%003i.gts"%(ctc,itc)
           if not os.path.exists(PN.DIR1+new_gts): break
    else:
       new_gts = gts
       itc     = int(gts.split(".")[-2])
    # assert file does not exists
    if os.path.exists(PN.DIR1+new_gts):
       print("    * File '%s' already exists!"%(PN.DIR1+new_gts))
       raise Exc.ABORTED
    # provisional files
    prov_frozen = prov_gts+".frozen"       # only if any frozen atom
    prov_molden = prov_gts+".molden"       # molden file
    # other files
    new_frozen = "%s.%003i.gts.frozen"%(ctc,itc)
    new_molden = "%s.%003i.molden"%(ctc,itc)
    # Rename files
    if os.path.exists(PN.DIR1+prov_gts   ): os.rename(PN.DIR1+prov_gts   ,PN.DIR1+new_gts   )
    if os.path.exists(PN.DIR1+prov_frozen): os.rename(PN.DIR1+prov_frozen,PN.DIR1+new_frozen)
    if os.path.exists(PN.DIR1+prov_molden): os.rename(PN.DIR1+prov_molden,PN.DIR5+new_molden)
    return new_gts, itc
#---------------------------------------------------------------#
def is_provisional(gts):
    return "prov" in gts.split(".")[-2]
#---------------------------------------------------------------#
def remove_provisionals(lgts):
    '''
    uses remove_definitive to remove gts files BUT this function
    only passes those that are provisional
    '''
    prov_gts = [gts for gts in lgts if is_provisional(gts)]
    remove_definitive(prov_gts)
#---------------------------------------------------------------#
def print_warnings():
    '''
    print warnings according to values of WRONG global variable
    '''
    #----------------#
    # Print warnings #
    #----------------#
    if True in WRONG.values():
       print("")
       print("  - WARNING!! It seems that something did not go well")
       if WRONG[1]: print("    * invalid name for some file/folder(s)!")
       if WRONG[2]: print("    * some folder(s) in %s does(do) NOT contain valid files!"%PN.UFOLDER)
       if WRONG[3]: print("    * reading process failed for a/some file(s) in %s!"%PN.UFOLDER)
       if WRONG[4]: print("    * file(s) without Hessian matrix!")
       if WRONG[5]: print("    * inconsistences in conformational cluster(s)!")
       if WRONG[6]: print("    * unexpected number of imaginary frequencies!")
    else:
       print("")
       print("  - Everything seems to be OK! :)")
#===============================================================#


#===============================================================#
def deal_with_file(esfile,dtrack,dctc):
    '''
    deal with an electronic structure (es) file inside UDATA/
    that is not inside a subfolder (i.e. a system with a single conformer)
    '''
    # Return data
    print("  - Reading file: %s"%(PN.UFOLDER+esfile))
    # Get CTC
    ctc = ".".join(esfile.split(".")[:-1])
    # Check CTC name
    if not fncs.is_string_valid(ctc,extra="_"):
       print("    '%s' is an invalid name for a CTC!\n"%ctc)
       WRONG[1] = True
       return dtrack,dctc

    # previously assignated?
    if esfile in dtrack:
       gtsfile = dtrack[esfile]
       print("    * gts file already assignated in %s! (%s)"%(PN.IFILE0,PN.DIR1+gtsfile))
       if os.path.exists(PN.DIR1+gtsfile):
          print("      - gts file already exists! Skipping generation\n")
          return dtrack,dctc
    else: gtsfile = None

    # GTS generation
    prov_gts, __   = temporal_gtsname(ctc)
    status  ,cache = userfile_to_gtsfile(esfile,prov_gts)
    remove = False
    if status ==  0:
       print("    * [reading failed]")
       remove,WRONG[3] = True,True
    elif status == -1:
       print("    * [Fcc NOT found] ")
       remove,WRONG[4] = True,True
    elif cache[5] not in (0,1):
       print("    * Incorrect number of imaginary frequencies in this CTC!")
       print("      n.imag. = %i"%cache[5])
       remove,WRONG[6] = True,True
    print()
    # Remove generated files
    if remove:
       remove_provisionals([prov_gts])
       return dtrack,dctc
    # rename files
    new_gts,itc = rename_files(ctc,prov_gts,0,gtsfile)
    # Create CTC instance
    CTC = ClusterConf(ctc)
    V0,ch,mtp,nifreq,mformu,pg = cache[2:]
    CTC.setvar("ch"    ,ch         ,"int")
    CTC.setvar("mtp"   ,mtp        ,"int")
    CTC.setvar("type"  ,nifreq     ,"int")
    CTC.setvar("mformu",mformu     ,"str")
    CTC._es = [(mtp,0.0)]
    CTC.add_itc(itc,V0,pg)
    # update dictionaries
    dtrack[esfile] = new_gts
    dctc[ctc] = CTC
    print_ctc_info(CTC)
    return dtrack,dctc
#---------------------------------------------------------------#
def deal_with_folder(folder,dtrack,dctc):
    '''
    deal with an electronic structure (es) file of a subfolder inside UDATA/
    '''

    print("  - Reading files in: %s"%(PN.UFOLDER+folder))

    # Get CTC
    ctc = folder[:-1]

    # Check CTC name
    if not fncs.is_string_valid(ctc,extra="_"):
       print("    '%s' is an invalid name for a CTC!\n"%ctc)
       WRONG[1] = True
       return dtrack,dctc

    # Files for each conformer of the CTC
    esfiles = sorted([folder+esfile for esfile in os.listdir(PN.UFOLDER+folder) \
                                    if  esfile.split(".")[-1].lower() in EXTENSIONS])
    if len(esfiles) == 0:
       print("    empty folder...\n")
       WRONG[2] = True
       return dtrack,dctc

    # initialize CTC
    CTC      = ClusterConf(ctc)
    ctc_vars = None

    # GTS generation
    count  = 0
    cache  = []
    remove = False
    for idx,esfile in enumerate(esfiles):

        # previously assignated?
        gtsfile = dtrack.get(esfile,None)
        if gtsfile is not None and os.path.exists(PN.DIR1+gtsfile):
           # read file
           print("    %s [gts exists: %s]"%(esfile,PN.DIR1+gtsfile))
           molecule = Molecule()
           molecule.set_from_gts(PN.DIR1+gtsfile)
           molecule.setup()
           status = 2
           # add data to cache
           V0      = molecule._V0
           ch      = molecule._ch
           mtp     = molecule._mtp
           nimag   = int(fncs.numimag(molecule._ccfreqs))
           mform   = molecule._mform
           pgroup  = molecule._pgroup
           cache_i = [esfile,gtsfile,V0,ch,mtp,nimag,mform,pgroup]
        else:
           prov_gts,count = temporal_gtsname(ctc,count)
           status,cache_i = userfile_to_gtsfile(esfile,prov_gts)
        # save cache_i
        cache.append(cache_i)

        # use first conformer to update ctc_vars
        if ctc_vars is None: ctc_vars = cache_i[3:7] 

        # Check reading and consistence within CTC
        if status ==  0:
           print("    %s [reading failed]"%esfile)
           print("    -- STOPPED --")
           WRONG[3],remove = True,True
        elif status == -1:
           print("    %s [Fcc NOT found] "%esfile)
           print("    -- STOPPED --")
           WRONG[4],remove = True,True
        elif ctc_vars != cache_i[3:7]:
           print("    %s [inconsistence found] "%esfile)
           print("    -- STOPPED --")
           WRONG[5],remove = True,True
           print_table_ufiles(cache)
        elif ctc_vars[2] not in (0,1):
           print("    %s [unexpected number of imag. freqs.] "%esfile)
           print("    -- STOPPED --")
           WRONG[6],remove = True,True
           print_table_ufiles(cache)
        elif status == 1:
           print("    %s [ok]            "%esfile)
        # sth wrong so we have to remove?
        if remove: break
    print()

    # Remove files if sth wrong
    if remove:
       remove_provisionals([cache_i[1] for cache_i in cache])
       return dtrack, dctc

    # sort cache by energy
    cache.sort(key = lambda x: x[2])
    # already in dctc
    if ctc in dctc: CTC0 = dctc.pop(ctc)
    else          : CTC0 = None
    # Update CTC instance
    CTC.setvar("ch"    ,ctc_vars[0],"int")
    CTC.setvar("mtp"   ,ctc_vars[1],"int")
    CTC.setvar("type"  ,ctc_vars[2],"int")
    CTC.setvar("mformu",ctc_vars[3],"str")
    CTC._es = [(ctc_vars[1],0.0)]
    # Rename gtsfiles
    itc = 0
    for idx,cache_i in enumerate(cache):
        prov_gts = cache_i[1]
        if is_provisional(prov_gts):
           new_gts,itc = rename_files(ctc,prov_gts,itc)
           itc_i = itc
           # update dtrack
           dtrack[cache_i[0]] = new_gts
           # update gtsfile in cache
           cache[idx][1] = new_gts
        else:
           itc_i = prov_gts.split(".")[-2]
        # weight
        weight = 1
        if CTC0 is not None: weight = max(CTC0.get_weight(itc_i),weight)
        # update CTC
        CTC.add_itc(itc_i,cache_i[2],cache_i[7],weight)
    CTC.find_minV0()
    # update with info of CTC0
    if CTC0 is not None:
       CTC._es     = CTC0._es
       CTC._fscal  = CTC0._fscal
       CTC._dics   = CTC0._dics
       CTC._dicsfw = CTC0._dicsfw
       CTC._dicsbw = CTC0._dicsbw
       CTC._diso   = CTC0._diso
       CTC._anh    = CTC0._anh

    # update dctc
    dctc[ctc] = CTC
    # Print info of this CTC
    print_ctc_info(CTC)
    # Delete cache
    del cache
    # Return data
    return dtrack,dctc
#===============================================================#


#===============================================================#
from   common.fncs       import fill_string

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
def main(_,__):

    # Assert folder with electronic-structure files exists!
    if not os.path.exists(PN.UFOLDER):
       print("   - %s: folder NOT found!"%PN.UFOLDER)
       print("")
       return

    # list user files
    ufiles,ufolders,all_ufiles = get_ufiles_ufolders()

    # Read tracking & pif.struc
    dtrack          = RW.read_track()[0]
    (dctc,dimasses) = RW.read_ctc()[0]

    # dictionary with isotopic masses
    if dimasses == {}: dimasses = dpt_im

    # deal with removed files (if any) - this may happen when executed more than once
    dtrack,dctc = deal_removed_ufiles(all_ufiles,dtrack,dctc)

    # remove non-listed gts files
    deal_nonlisted_gts(dtrack)

    # clean-up dctc
    dctc = cleanup_dctc(dctc)

    # Folder creation
    for folder in [PN.DIR1,PN.DIR5]:
        if not os.path.exists(folder): os.mkdir(folder)

    # Deal with CTCs defined through single files
    for esfile in ufiles:
        dtrack,dctc = deal_with_file(esfile,dtrack,dctc)

    # Deal with CTCs defined through folders
    for folder in ufolders:
        dtrack,dctc = deal_with_folder(folder,dtrack,dctc)

    # print summary
    print("  - SUMMARY")
    ls_struc(dctc)

    # (re)write file: tracking
    if dtrack != {}:
       print("  - Writing/Updating information in '%s'"%PN.IFILE0)
       RW.write_track(dtrack)

    # (re)write file: pif.struc
    if dctc != {}:
       print("  - Writing/Updating information in '%s'"%PN.IFILE1)
       RW.write_ctc(dctc,dimasses)

    print_warnings()
#===============================================================#


