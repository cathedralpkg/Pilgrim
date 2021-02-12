#!/usr/bin/python3
'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.1
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
| Program    :  pilgrim            |
| Last Update:  2020/04/18 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''


#===========================================================#
# FIRSTLY: check modules!!!                                 #
#-----------------------------------------------------------#
import sys                                                  #
if sys.version_info.major < 3: exit("Python 3 required!!")  #
#-----------------------------------------------------------#
import modpilgrim.checkmods   as     checkmods              #
checkmods.checkmods()                                       #
#-----------------------------------------------------------#
import os                                                   #
import time                                                 #
#-----------------------------------------------------------#
from modpilgrim.strings import PROGNAME, VERSION            #
from modpilgrim.strings import PROGHEAD, AUTHORINFO         #
from modpilgrim.strings import exe_info                     #
#-----------------------------------------------------------#
import common.Exceptions       as     Exc                   #
from   common.fncs             import classify_args         #
from   common.fncs             import print_string          #
from   common.fncs             import time2human            #
#-----------------------------------------------------------#
import modpilgrim.names       as PN                         #
#-----------------------------------------------------------#
import modpilgrim.optGATHER      as gather                  #
import modpilgrim.optINPUT       as inpmenu                 #
import modpilgrim.optICS         as ics                     #
import modpilgrim.optPFN         as pfn                     #
import modpilgrim.optPATH        as path                    #
import modpilgrim.optRCONS       as rcons                   #
import modpilgrim.optKMC         as kmc                     #
import modpilgrim.optFIT         as fit                     #
import modpilgrim.optPLOT        as plot                    #
import modpilgrim.optHLCALC      as hlcalc                  #
import modpilgrim.optKIES        as kies                    #
#import modpilgrim.optHARVEST     as harvest                 #
import modpilgrim.optSUMMARY     as summary                 #
from   modpilgrim.exceptions     import deal_with_exception #
from   modpilgrim.diverse        import get_input_data      #
from   modpilgrim.diverse        import dlevel_to_files     #
#-----------------------------------------------------------#
from   modpilgrim.helps          import HELP_main           #
from   modpilgrim.helps          import HELP_gather         #
from   modpilgrim.helps          import HELP_input          #
from   modpilgrim.helps          import HELP_ics            #
from   modpilgrim.helps          import HELP_pfn            #
from   modpilgrim.helps          import HELP_path           #
from   modpilgrim.helps          import HELP_hlcalc         #
from   modpilgrim.helps          import HELP_rcons          #
from   modpilgrim.helps          import HELP_fit            #
from   modpilgrim.helps          import HELP_kmc            #
from   modpilgrim.helps          import HELP_summary        #
#from   modpilgrim.helps          import HELP_harvest       #
from   modpilgrim.helps          import HELP_plot           #
from   modpilgrim.helps          import HELP_kies           #
#===========================================================#

OPTIONS1 = "pfn,path,hlcalc,rcons,kmc,fit,plot,kies,summary"
OPTIONS2 = "ls,gather,input,ics"
OPTIONS  = OPTIONS1+","+OPTIONS2

#==========================================================#
#                    Checking arguments                    #
#==========================================================#
def check_ics(dargs):
    key = "ics"
    if key in dargs.keys():
        if   len(dargs[key]) == 0:
             dargs[key].append("2")
             dargs[key].append("*")
        elif len(dargs[key]) == 1:
             dargs[key].append("*")
        elif len(dargs[key]) == 2:
             pass
        else:
            print("Maximum number of variables for --ics is two!!")
            exit()
    return dargs
#----------------------------------------------------------#
def check_sft(dargs):
    key = "software"
    if key in dargs.keys():
        try   : return dargs[key][0]
        except: return None
    else      : return "gaussian"
#----------------------------------------------------------#
def check_dlevel(dargs):
    key = "dlevel"
    if key in dargs.keys(): return True
    else                  : return False
#----------------------------------------------------------#
def check_boolms(dargs):
    key = "ts0"
    if key in dargs.keys(): return True
    else                  : return False
#==========================================================#



#==========================================================#

def main():
    cdate = time.strftime("%Y-%m-%d")
    ctime = time.strftime("%H:%M:%S")
    t1 = time.time()

    # Read user arguments
    user_args = sys.argv[1:]
    if len(user_args) == 0:
       print(PROGNAME)
       print(PROGHEAD)
       print(AUTHORINFO)
       print_string(HELP_main,2)
       return
    dargs = classify_args(user_args)

    # User asked for version
    if "version" in dargs.keys() or "-v" in user_args:
       print("Current version: %s"%VERSION)
       return

    # User asked for help
    elif "help"    in dargs.keys():
       print(PROGNAME)
       print(PROGHEAD)
       print(AUTHORINFO)
       if len(dargs["help"]) == 0   : print_string(HELP_main   ,2)
       if "gather"  in dargs["help"]: print_string(HELP_gather ,1)
       if "input"   in dargs["help"]: print_string(HELP_input  ,1)
       if "ics"     in dargs["help"]: print_string(HELP_ics    ,1)
       if "pfn"     in dargs["help"]: print_string(HELP_pfn    ,1)
       if "path"    in dargs["help"]: print_string(HELP_path   ,1)
       if "hlcalc"  in dargs["help"]: print_string(HELP_hlcalc ,1)
       if "rcons"   in dargs["help"]: print_string(HELP_rcons  ,1)
       if "fit"     in dargs["help"]: print_string(HELP_fit    ,1)
       if "kmc"     in dargs["help"]: print_string(HELP_kmc    ,1)
       if "plot"    in dargs["help"]: print_string(HELP_plot   ,1)
       if "kies"    in dargs["help"]: print_string(HELP_kies   ,1)
      #if "harvest" in dargs["help"]: print_string(HELP_harvest,1)
       if "summary" in dargs["help"]: print_string(HELP_summary,1)
       return

    # Print logo and current date
    string_head  = PROGNAME+"\n"+PROGHEAD+"\n"+AUTHORINFO+"\n"+"\n"
    string_head += " -----------------------------------------------------------\n"
    string_head += exe_info()
    string_head += " -----------------------------------------------------------\n\n"
    print(string_head)

    # check some arguments
    dargs    = check_ics(dargs)
    software = check_sft(dargs)
    dlevel   = check_dlevel(dargs)
    boolms   = check_boolms(dargs)

    # --software used properly
    if software is None:
       print("  Software option has to be followed by an argument!")
       print("")
       exit()

    # Act according user argument
    IN_OPTS = False
    for option in OPTIONS.split(","):
        if option not in dargs.keys(): continue
        IN_OPTS = True
        # read input files and print table
        idata, status, string = get_input_data()
        print_string(string,nbs=3)
        # data files according to case
        datafiles,string = dlevel_to_files(dlevel)
        if option in OPTIONS1.split(","): print_string(string,nbs=3)
        case = (datafiles,dlevel,software)
        # execute pilgrim with option
        print(" ===========================================================")
        print(" ||  EXECUTING PILGRIM WITH --%-6s                      ||"%option)
        print(" ===========================================================")
        print("")
        if   option == "ls"     :  gather.ls_struc(idata[0][0])
        elif option == "gather" :  gather.main(idata,status                      )
        elif option == "input"  : inpmenu.main(idata,status                      )
        elif option == "ics"    :     ics.main(idata,status,*dargs["ics"   ]     )
        elif option == "pfn"    :     pfn.main(idata,status,case,dargs["pfn"   ])
        elif option == "path"   :    path.main(idata,status,case,dargs["path"  ],boolms)
        elif option == "hlcalc" :  hlcalc.main(idata,status,case,dargs["hlcalc"])
        elif option == "rcons"  :   rcons.main(idata,status,case,dargs["rcons" ])
        elif option == "kmc"    :     kmc.main(idata,status,case,dargs["kmc"   ])
        elif option == "kies"   :    kies.main(idata,status,case                )
        elif option == "fit"    :     fit.main(idata,status,case,dargs["fit"   ])
        elif option == "plot"   :    plot.main(             case,dargs["plot"  ])
       #elif option == "harvest": harvest.main(idata,status,case,dargs["harvest"])
        elif option == "summary": summary.main(idata,status,case,dargs["summary"])
        print("")
        print(" ===========================================================")
        break

    # User asked for nothing. Print help!
    if not IN_OPTS: print_string(HELP_main,2); return


    # Print elapsed time
    t2 = time.time()
    timeline = "| Total Elapsed Time: %5.1f %5s |"%time2human(t2-t1,"secs")
    print(" "*25+"-"*len(timeline))
    print(" "*25+        timeline)
    print(" "*25+"-"*len(timeline))
    print("")
#==========================================================#


#==========================================================#
if __name__ == '__main__':
   try: main()
   except Exception as exception: deal_with_exception(exception)
#==========================================================#


