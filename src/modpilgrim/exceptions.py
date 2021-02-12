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
| Module     :  modpilgrim         |
| Sub-module :  exceptions         |
| Last Update:  2020/04/18 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

import traceback
import modpilgrim.names  as PN
import common.Exceptions as Exc

#===============================================================#
def deal_with_exception(exception):
    # initialize error string
    serr = "\n"
    # add nice heading
    serr += "================================\n"
    serr += "|| Pilgrim ended unexpectedly ||\n"
    serr += "================================\n"
    serr += "\n"
    # type of exception
    exctype = type(exception)
    # add to string according to error
    if   exctype == Exc.NoTemplateGiven:
         serr += "  * Unable to find template for the electronic structure calculation!\n"
    elif exctype == Exc.ReadProblem:
         serr += "  * Problem reading line in '%s':\n"%exception._file
         serr += "    %s\n"%exception._var
    elif exctype == Exc.ICFails:
         serr += "  * Problems with internal coordinates...\n"
    elif exctype == Exc.RstReadProblem:
        serr += "  * Problems reading rst file: %s\n"%exception._var
    elif exctype == Exc.RstDiffVar:
         serr += "  * Problems with data in '%s'\n"%(exception._rst)
         serr += "    Variable '%s' differs from that in '%s'\n"%(exception._var,PN.IFILE3)
    elif exctype == Exc.RstDiffPoint:
         l_i,s_i = exception._var
         serr += "  * Problems with data in rst file!\n"
         serr += "    Geometry %s (%.4f bohr) does not agree with generated one"%(l_i,s_i)
    elif exctype == Exc.RstDiffTS:
         serr += "  * Problems with data in rst file!\n"
         serr += "    Geometry of TS differs!"
    elif exctype == Exc.ExeNotDef:
         serr += "  * Variable(s) for software is(are) not defined in your .bashrc file\n"
    elif exctype == Exc.ExeNotFound:
         serr += "  * Executable file for software was not found\n"
    elif exctype == Exc.NoDLEVELpfn:
         serr += "  * It seems that the partition functions were not calculated at dual-level\n"
    elif exctype == Exc.NoDLEVELdata:
         serr += "  * It seems that high-level energies are missing\n"
    elif exctype == Exc.DLEVELsthWrong:
         serr += "  * Something went wrong when applying dual-level energies...\n"
    elif exctype == Exc.CalcFails:
         serr += "  * The electronic structure calculation failed (%s)!\n"%(exception._var)
    elif exctype == Exc.NoTemps:
         serr += "  * Unable to find working temperatures!\n"
    elif exctype == Exc.NoICS:
         serr += "  * User asked for the adiabatic potential in internal coordinates.\n"
         serr += "    However, these coordinates are not defined in '%s'\n"%PN.IFILE1
    elif exctype == Exc.IncompData:
         serr += "  * Inconsistences with the list of temperatures: '%s'\n"%exception._var
    elif exctype == Exc.NoRateCons:
         rctype, rcname = exception._var
         serr += "  WARNING: Unable to find '%s' rate constant for '%s'\n"%(rctype,rcname)
    elif exctype == Exc.WrongInIsomass:
         serr += "  * Isotopic mass '%s' not defined in '%s'\n"%(exception._var,PN.IFILE1)
    elif exctype == Exc.NoData:
         serr += "  * It seems that the partition functions are missing from the data files\n"
         serr += "  * Did you run Pilgrim with the --pfn option?\n"
    elif exctype == Exc.LostConformer:
         serr += "  * It seems that the partition functions are missing from the data files\n"
         serr += "    for conformer: %s\n"%exception._var
    elif exctype == Exc.TSnotFOUND:
         serr += "  WARNING: transition state not found in rst file!\n"
    elif exctype == Exc.OnlyMEP:
         serr += "  WARNING: This transition state is not defined in a reaction\n"
         serr += "           and/or the reaction reactant is incorrectly defined \n"
         serr += "           and/or the --pfn option was not run.\n"
         serr += "           Therefore, TST correction coefficients will not be calculated!\n"
    elif exctype == Exc.FileIsNotGTS:
         serr += "  WARNING: it seems that %s is not a gts file\n"%exception._var
    elif exctype == Exc.WrongONIOMlayers:
         serr += "  * Problems in the definition of ONIOM layers!\n"
    elif exctype == Exc.ErrorQRC:
         serr += "  * Problems with QRC!\n"
    elif exctype == Exc.ABORTED:
         serr += "  * Execution aborted!\n"
    else:
         serr += "  * Unexpected error: %s\n"%str(exception)
         serr += "\n"
         message = str(traceback.format_exc())
         for line in message.split("\n"):
             if line == "": continue
             serr += "    "+line+"\n"
    serr += "\n"
    # print lines of exceptions
    for line in serr.split("\n"): print(" "+line)
#----------------------------------------------------------#

