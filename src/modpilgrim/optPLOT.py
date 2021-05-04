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
| Module     :  modpilgrim         |
| Sub-module :  optPLOT            |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#===============================================================#
import os
import modpilgrim.names as PN
#---------------------------------------------------------------#
from   modpilgrim.plotting import read_plotfile
#===============================================================#

SMENU = '''    -------------------------------------------------
    Select indices of data to plot:
      * to select all, just type 'all'
      * to select a range, use 'idx1-idx2' syntax
      * examples:
        >> all
        >> 1 3 5
        >> 1-3 5 8
    -------------------------------------------------
    Reduce list of plots:
      * to list plots with a given keyword in
        their name, use:
           >> select keyword
        For example, with:
           >> select mueff
        plots related to effective masses are listed
      * to exclude plots with a given keyword, use:
           >> exclude keyword
      * to go back to the original list, type:
           >> reset          
    -------------------------------------------------
    To exit use 'exit()', 'end()', 'exit', or 'end'
    -------------------------------------------------'''


TITLESIZE =15
LEGENDSIZE=14
XLABELSIZE=14
YLABELSIZE=14
XTICKSIZE =14
YTICKSIZE =14
#===============================================================#
def main(case,pltcase="show"):

    (dof,hlf,plotfile),dlevel,software = case

    if "pdf" in pltcase: pltcase = "pdf"
    else               : pltcase = "show"

    plotfile = PN.get_plf(dlevel)
    print("    Data will be read from '%s'"%plotfile)
    print("")

    # file exists?
    if not os.path.exists(plotfile):
       print("    ERROR! File does not exist!")
       print("")
       return
    # read file
    plotdata = read_plotfile(plotfile)
    allplots = sorted(plotdata.keys())
    # file contains data?
    if len(allplots) == 0:
       print("    ERROR! File contains no data!")
       print("")
       return

    # print which data is available
    current = list(allplots)
    while True:
       try: 
          print("    Available data:")
          for idx,mlabel in enumerate(current):
              print("       (%02i) %s"%(idx+1,mlabel))
          print("")
          print(SMENU)
          print("")

          # ask user what to plot
          answer = input("       >> ").strip()
          print("")
          # exit?
          if   answer in ["end","end()","exit","exit()",".."]: break
          # select?
          elif "select" in answer:
              keys       = answer.split("select")[1].split()
              newcurrent = []
              for plot in current:
                  nfalse = [(key in plot) for key in keys].count(False)
                  if nfalse == 0: newcurrent.append(plot)
              current = list(newcurrent)
              continue
          # exclude?
          elif "exclude" in answer:
              keys       = answer.split("exclude")[1].split()
              newcurrent = []
              for plot in current:
                  ntrue = [(key in plot) for key in keys].count(True)
                  if ntrue == 0: newcurrent.append(plot)
              current = list(newcurrent)
              continue
          # reset to all keys
          elif answer.split()[0] == "reset":
             current = list(allplots)
             continue
          # extract targets
          targets = answer.split()
          indices = []
          for target in targets:
              if "-" in target:
                  idx1, idx2 = target.split("-")
                  indices += range(int(idx1),int(idx2)+1)
              elif target == "all":
                  indices += range(1,len(current)+1)
              else:
                  indices.append(int(target))
          # sort indices
          indices.sort()

          # ask for name for pdf file
          if pltcase == "pdf":
             print("    Name for the pdf file (no extension needed)? ")
             pdffile = input("       >> ").strip()
             if not pdffile.endswith(".pdf"): pdffile += ".pdf"
             # does the pdf file exist?
             if os.path.exists(pdffile):
                 print("       file '%s' already exists..."%pdffile)
                 answer = input("       continue (y/N)? ").strip().lower()
                 if answer not in ["yes","y","si","s"]: return
             print("")

          # Import modules
          print("    Importing matplotlib.pyplot")
          try:
             import matplotlib.pyplot as plt
          except:
             import pylab             as plt

          if pltcase == "pdf":
             print("    Importing PdfPages from matplotlib.backends.backend_pdf")
             from   matplotlib.backends.backend_pdf import PdfPages
             thepdf = PdfPages(pdffile,keep_empty=False)
          print("")

          # plot data
          if pltcase == "pdf" : print("    Adding plots to '%s' file"%pdffile)
          if pltcase == "show": print("    Showing plots")
          for idx in indices:
              mlabel = current[idx-1]
              ptype, xydata, title, xlabel, ylabel = plotdata[mlabel]
              print("      %s..."%mlabel)
              plt.xlabel(xlabel,fontsize=XLABELSIZE)
              plt.ylabel(ylabel,fontsize=YLABELSIZE)
              plt.xticks(fontsize=XTICKSIZE)
              plt.yticks(fontsize=YTICKSIZE)
              plt.title(title,fontsize=TITLESIZE)
              for xx,yy,pformat,plabel in xydata:
                  case1 = (pformat == "") and (plabel == "")
                  case2 = (pformat != "") and (plabel == "")
                  case3 = (pformat == "") and (plabel != "")
                  case4 = (pformat != "") and (plabel != "")
                  if   xx is None: plt.axhline(y=yy, color='r', linestyle='--')
                  elif yy is None: plt.axvline(x=xx, color='r', linestyle='--')
                  elif case1 and ptype == "plot": plt.plot(xx,yy)
                  elif case2 and ptype == "plot": plt.plot(xx,yy,pformat)
                  elif case3 and ptype == "plot": plt.plot(xx,yy,label=plabel)
                  elif case4 and ptype == "plot": plt.plot(xx,yy,pformat,label=plabel)
                  elif           ptype == "bar" : plt.bar(range(len(xx)),yy,tick_label=xx)
              if len(xydata) > 1: plt.legend(loc="best",framealpha=0.5,fontsize=LEGENDSIZE)
              if pltcase == "show": plt.show()
              if pltcase == "pdf": thepdf.savefig(bbox_inches='tight')
              plt.close()
          if pltcase == "pdf": thepdf.close()
          print("")
       except KeyboardInterrupt:
          return
#===============================================================#



