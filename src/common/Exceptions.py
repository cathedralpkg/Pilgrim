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
| Sub-module :  Exceptions         |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains classes used to
raise exceptions
'''

#=============================================#
class END(Exception)             : pass
class ABORTED(Exception)         : pass
class FileType(Exception)        : pass
class ReadProblem(Exception)     : pass
class UnknownSoft(Exception)     : pass
class ExeNotDef(Exception)       : pass
class ExeNotFound(Exception)     : pass
class CalcFails(Exception)       : pass
class RstReadProblem(Exception)  : pass
class RstDiffTS(Exception)       : pass
class RstDiffPoint(Exception)    : pass
class NoDLEVELdata(Exception)    : pass
class RstDiffVar(Exception)      : pass
class SDpmNoHess(Exception)      : pass
class SDdiffds(Exception)        : pass
class NoICS(Exception)           : pass
class NoTemps(Exception)         : pass
class NoData(Exception)          : pass
class IncompData(Exception)      : pass
class WrongVar(Exception)        : pass
class WrongInIsomass(Exception)  : pass
class NoDLEVELpfn(Exception)     : pass
class DLEVELsthWrong(Exception)  : pass
class ICFails(Exception)         : pass
class NoTemplateGiven(Exception) : pass
class NoReacMol(Exception)       : pass
class LostConformer(Exception)   : pass
class DiffMassConf(Exception)    : pass
class LackOfGts(Exception)       : pass
class NoRateCons(Exception)      : pass
class OnlyMEP(Exception)         : pass
class TSnotFOUND(Exception)      : pass
class NoGTSfile(Exception)       : pass
class FileIsNotGTS(Exception)    : pass
class FccNotFound(Exception)     : pass
class LevelNotFound(Exception)   : pass
class WrongDimension(Exception)  : pass
class WrongONIOMlayers(Exception): pass
class UnableGenRandAng(Exception): pass
class ErrorHConstraint(Exception): pass
class ErrorSConstraint(Exception): pass
class ErrorQRC(Exception)        : pass
class ErrorTorsion1(Exception)   : pass
class ErrorTorsionN(Exception)   : pass
class ErrorTorsionRepe(Exception): pass
#=============================================#

