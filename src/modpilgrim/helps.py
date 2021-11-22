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
| Sub-module :  helps              |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

import modpilgrim.names as PN

#===============================================================#
# help string for main function                                 #
#===============================================================#
HELP_main  = '''

Hi! I am PILGRIM, a program that calculates thermal
rate constants and simulates reaction mechanisms

   * Setting the system up:

       $ pilgrim --gather

   * Input files preparation:

       $ pilgrim --input

   * Running Pilgrim:

       $ pilgrim --option [target[.idx]] [--dlevel] [--software (ESSO)]

       Main options:
           --pfn   | Calculates the partition function at the 
                   |    stationary points
           --path  | Calculates the MEP & VTST coefficients
           --rcons | Calculates the thermal rate constants
           --kmc   | Performs a kinetic Monte Carlo simlulation

       Options for special actions:
           --ts0      | Selects the TS conformer of smallest energy
                      |    It goes together with --path
           --dlevel   | Dual-level calculations and interpolation 
                      |    via ISPE
           --software | Selects ESSO ('gaussian' (default) or 'orca')
                      |    It goes together with --path OR --hlcalc

       Additional options and tools
           --ics     | Generates internal coordinates
           --hlcalc  | Performs high-level calculations for the
                     |    subsequent dual-level interpolation
           --fit     | Fits the thermal rate constants to prescribed
                     |    analytical expressions
           --plot    | Generates diverse plots
           --summary | Summarizes output files

       Information options:
           --help [option]| Displays help messages for a given option,
                          |    for instance $pilgrim --help input
           --ls           | Lists all the CTC-blocks defined in %s
           --version      | Displays the version of the program

       Abbreviations and further explanations:

           ESSO     : Electronic structure software

           CTC-block: A section of %s contained between the keywords:

                          start_ctc label 
                          end_ctc

                      'label' --> a name given to the CTC-block
                                  We'll use the generic label 'ctcsp' 
                                  if the structures inside the CTC-block
                                  is a stationary point (minima or TSs). 
                                  When we want to address only TS structures 
                                  the generic label 'ctcts' will be used
                                  Individual structures inside the CTC-block
                                  will be addressed as 'ctcsp.idx' OR
                                  'ctcts.idx'
           'chemname' : a generic label given to a chemical reaction. 


           List of the targets allowed by the different options. The table 
           also shows the case when 'target' is omitted (default).
         
               --option  |      default      |  allowed 'target'
               =================================================
               --pfn     |    all ctcsp      |   'ctcsp'
               --path    |    all ctcts      |   'ctcts[.idx]'
               --rcons   |    all chemname   |   'chemname'
               --kmc     |      ---          |      ---

               --ics     |    all ctcts      |   'ctcts[.idx]'
               --hlcalc  |    all ctcsp      |   'ctcsp[.idx]
               --fit     |    all chemname   |   'chemname'
               --plot    |      ---          |      ---
               --summary |      ---          | pfn, path, rcons


Let's find together how your system behaves!
Go ahead!
'''%(PN.IFILE1,PN.IFILE1)
#===============================================================#



#===============================================================#
# help string for --gather                                      #
#===============================================================#
HELP_gather  = '''
==================================================================
|               Information for option:  --gather                |
==================================================================
 _____________________                                        

 Command line:
 $ pilgrim.py --gather
 _____________________

   INPUT FILES:
       The electronic structure files (ESFILs).

   OUTPUT FILES:
       pif.struc

   FORMAT OF THE ESFILs:
       This option gathers the ESFILs inside %s and converts them 
       to the Pilgrim internal format (gts).
       Currently, Pilgrim is able to read/convert ESFILs
       with the following extensions:


	 =======================================================
	  Program      |  ESFIL extension for each structure (1)
	 =======================================================
	  GAUSSIAN (2) |        .fchk OR  .out  OR .log
	 ------------------------------------------------------
	  ORCA (3)     |        .out  AND .hess
	 ------------------------------------------------------
	  Q2DTor       |        .gts
	  Pilgrim      |        .gts
	  =====================================================
          (1) Other extensions not listed in the Table will be
              ignored by Pilgrim.
          (2) For a given structure use .fchk OR .out OR .log.
          (3) For a given structure the .out and the .hess 
              files should have the same name.


       Notice that ESFILs with other extension will be ignored
       by Pilgrim!

       We highlight that each of these ESFILs must contain: 
	  - The optimized structure
	  - The Hessian matrix

       The generated gts files will be stored inside %s and
       an the '%s' input file will be generated.  Check it
       carefully!

   WHERE TO PLACE THE ESFILs?
   
       Example: To study the reaction:

	    EtOH + H  --> TSc --> EtO + H2

       EtOH (ethanol) has three conformers (one anti and two gauge).
       The two gauche structures are conformational enantiomers and
       therefore only one of them should be included. 
       In summary, include only distinguishable conformations, and from
       those include only one structure if there is another structure
       that is a conformational enantiomer of the former.
       In this case the ESFILs are called (the name of the
       files is arbitrary): 
	    t-et.fchk
	    g-et.fchk

       H is the hydrogen atom with a ESFIL called:
	    H.out
    
       TSc is a transition state with two conformations (one anti and two
       gauche of which only one should be included). The ESFILs
       are called:
	    t-tsc.fchk
	    g-tsc.fchk

       EtO is a radical with only one distinguishable conformation. 
       The ESFILs is called:
	    prod.fchk

       The molecule of H2 has a ESFIL called:
            H2.fchk

       The files should be placed as indicated below:

	   %s/
	     |---> EtOH/
		   |--> t-et.fchk
		   |--> g-et.fchk
	     |---> H.fchk
	     |---> TSc
		   |--> t-tsc.fchk
		   |--> g-tsc.fchk
	     |---> EtO/
		   |--> prod.fchk
             |---> H2.fchk

       Notice that compounds with a single structure can hang directly
       from UDATA/

       The labels of the species are: 

	     EtOH --> EtOH
	     H    --> H
	     TS   --> TSc
	     EtO  --> EtO
             H2   --> H2

           NOTE:  We use the generic label 'spname' to name each 
           of the species defined by the user (minima and TSs).
       
	   These spnames will be used by the option --gather to create 
           four blocks inside pif.struc. Each block contains a cluster of
           torsional conformers (CTC) of each species, and it has a label. 
           This label has the generic name of 'ctcsp' in the case of a 
           statinary point (minima or TS) and 'ctcts' in the particular 
           case of a TS.

              start_ctc EtOH
                ...
              end_ctc
              start_ctc H
                ...
              end_ctc
              start_ctc EtO
                ...
              end_ctc
              start_ctc H2
                ...
              end_ctc
              start_ctc TSc
                ...
              end_ctc

           In this case the 'ctcsp' names are: 
              EtOH, H, EtO, H2 and TSc
           The 'ctcsts name is: TSc

==================================================================
'''%(PN.UFOLDER,PN.DIR1,PN.IFILE1,PN.UFOLDER)
#===============================================================#



#===============================================================#
# help string for --ics                                         #
#===============================================================#
HELP_ics  = '''
==================================================================
|                 Information for option:  --ics                 |
==================================================================

___________________________________

 Command line:
 $ pilgrim.py --ics [mode [target[.idx]]]
___________________________________


   INPUT FILES:
       The gts files.

   OUTPUT FILES:
       Adds internal coordinates to the TS structures in pif.struc

     This option generates (or at least tries to generate) a set of
     internal coordinates for each transition state.  It admits up
     to two arguments (mode and target, see above):

       'mode':
       --------------
	 1: for standard generation of a redundant set
	 2: for the generation of a smaller redundant set (default)
	 3: for the generation of a non-redundant set
	-1: to check a previously generated set
	-2: to delete the current set

       'target':
       ----------------
	 * 'target = ctcts[.idx]' (use --help for explanation)
	   Example:
	     For a TS labelled 'tstry' with three conformers and gts files:
	     tstry1.001.gts, tstry.002.gts and tstry.003.gts
		 $ pilgrim --ics 2 tstry.001 
		     - It uses mode 2 for structure 001 of 'tstry'
		 $ pilgrim --ics 2 tstry
		     - It uses modes 2 for structures 001, 002 and 003 
		       of 'tstry'
		 
       If 'target' is omitted (default; all 'ctcts' inside pif.struc 
       will be considered) 

==================================================================
'''
#===============================================================#



#===============================================================#
# help string for --pfn                                         #
#===============================================================#
HELP_pfn = '''
==================================================================
|                 Information for option:  --pfn                 |
==================================================================

 ______________________________________

 Command line:
 $ pilgrim.py --pfn [target] [--dlevel]
 ______________________________________

   PIF INPUT FILES:
       pif.struc
       pif.temp

   OUTPUT FILES (in %s):
       pfn.target.slevel.txt
     With --dlevel 
       pfn.target.dlevel.txt

       'target':
       ----------------
	 * 'target = ctcsp' (use --help for explanation)
       (Default) If 'target' is omitted all 'ctcsp'
	   inside pif.struc will be considered
   
     If high-level energies are available (see option --hlcalc),
     the --dlevel option can be included in the command line.
     In this case, the program uses the high-level energies in the
     calculation of the partition functions.

==================================================================
'''%(PN.DIR3)
#===============================================================#



#===============================================================#
# help string for --path                                        #
#===============================================================#
HELP_path = '''
==================================================================
|                 Information for option: --path                 |
==================================================================

 _______________________________________________________________

 Command line:
 $ pilgrim.py --path [target[.idx]] [--software ESSO] [--dlevel]
 _______________________________________________________________


   PIF INPUT FILES:
       pif.struc
       pif.temp
       pif.path
       pif.calcs

   RESTART FILE (if available):
       %starget.idx.rst

   OUTPUT FILES (in %s):
       path.target.idx.slevel.txt
     With --dlevel 
       path.target.idx.dlevel.txt

     This command line allows the calculation of the minimum energy path
     (MEP) associated with a given transition state structure. 

     Pilgrim also calculates the CVT, SCT and CAG coefficients.

     'target':
     ----------------
       * 'target = ctcts[.idx]' (use --help for explanation)
           index [.idx] of a single TS inside the CTC-block
       (Default) If 'target' is omitted, all 'ctcts'
	   inside pif.struc will be considered


     --software ESSO
	 Selects the electronic structure software (ESSO);
	 at this moment ESSO can be 'gaussian' (default) or 'orca'

     --dlevel
	 Corrects the MEP with high-level energies previously
         calculated with the --hlcalc option; the correction is
	 carried out by the Interpolated Single Point Energies (ISPE)
	 algorithm.


     Example:
       For a TS with 'ctcts' name equal to 'TSc' and with two
       conformers, the command line:

	   $ pilgrim --path TSc.001 

	       - It calculates the MEP for structure 001 of 'TSc'

	   $ pilgrim --ics 2 TSc

	       - It calculates the MEP for structures 001 and 002 
		 of 'TSc'. The calculations are performed in 
                 sequential form. First tstry.001, and then tstry.002
		 However, the user can also send all jobs at 
                 the same time:

		   $ pilgrim --path TSc.001 

		   $ pilgrim --path TSc.002 

   If good computational resources are available, we highly recommend to
   execute this option specifically for each transition state,
   sending each execution to a queue or to the background with
   the nohup command.

==================================================================
'''%(PN.DIR4,PN.DIR3)
#===============================================================#



#===============================================================#
# help string for --hlcalc                                      #
#===============================================================#
HELP_hlcalc = '''
==================================================================
|                 Information for option: --hlcalc               |
==================================================================

 ________________________________________________

 Command line:
 $ pilgrim.py --hlcalc [target[.idx]] [--software ESSO]
 ________________________________________________

   PIF INPUT FILES:
       pif.struc
       pif.path
       pif.calcs

   OUTPUT FILES:
       %shighlevel.txt

     This option carries out the high-level calculations of the
     stationary points of the system, as well as the selected MEP
     points, as indicated in %s.

     'target':
     ----------------
     * 'target = ctcsp[.idx]' (use --help for explanation)
     (Default) If 'target' is omitted all 'ctcsp'
        inside pif.struc will be considered

     --software ESSO
	 Selects the electronic structure software (ESSO);
	 at this moment ESSO can be 'gaussian' (default) or 'orca'

==================================================================
'''%(PN.DIR2,PN.IFILE7)
#===============================================================#



#===============================================================#
# help string for --rcons                                       #
#===============================================================#
HELP_rcons = '''
==================================================================
|                Information for option:  --rcons                |
==================================================================

 ________________________________________

 Command line:
 $ pilgrim.py --rcons [target] [--dlevel]
 ________________________________________

   PIF INPUT FILES:
       pif.struc
       pif.temp
       pif.path

   OUTPUT FILES (in %s):
       rcons.target.slevel.txt
     With --dlevel 
       rcons.target.dlevel.txt

     With this option, the rate constant of the selected chemical 
     reaction is calculated.

     'target':
     ----------------
       * 'target = chemname' (use --help for explanation)
       (Default) If 'target' is omitted all 'chemname'
          inside %s will be considered

     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.

==================================================================
'''%(PN.DIR3,PN.IFILE5)
#===============================================================#



#===============================================================#
# help string for --fit                                         #
#===============================================================#
HELP_fit = '''
==================================================================
|                 Information for option:  --fit                 |
==================================================================

_______________________________________

 Command line:
 $ pilgrim.py --fit [target] [--dlevel]
_______________________________________

     'target':
     ----------------
       * 'target = chemname' (use --help for explanation)
       (Default) If 'target' is omitted all 'chemname'
          inside %s will be considered

     Fits previously calculated rate constants of the available chemical 
     reactions to analytical expressions.

     At the end of the execution, a "Fitting sumamry" is printed.
     These lines can be used in the %s file to use analytical rate 
     expressions instead of calculated values.
     
     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.

==================================================================
'''%(PN.IFILE5,PN.IFILE6)
#===============================================================#



#===============================================================#
# help string for --kmc                                         #
#===============================================================#
HELP_kmc  = '''
==================================================================
|                 Information for option:  --kmc                 |
==================================================================

_____________________________

 Command line:
 $ pilgrim.py --kmc [kmcname] [--dlevel]
______________________________

   INPUT FILES:
       pif.temp
       pif.kmc

   OUTPUT FILES (in %s):
       kmc.kmcname.slevel.txt
     With --dlevel 
       kmc.kmcname.dlevel.txt

     With this option, a kinetic Monte Carlo simulation is
     performed. 

     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.


==================================================================
'''%PN.DIR3
#===============================================================#

#===============================================================#
# help string for --kies                                        #
#===============================================================#
HELP_kies  = '''
==================================================================
|                 Information for option:  --kies                |
==================================================================

_____________________________

 Command line:
 $ pilgrim.py --kies
______________________________

   INPUT FILES:
       pif.temp
       pif.struc
       pif.chem

     With this option, Kinetic Isotopic Effects (KIEs) are
     calculated.
     KIEs are split into different contributions.
     The program asks for the root and the isotopic reactions.

     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.


==================================================================
'''
#===============================================================#


#===============================================================#
# help string for --plot                                        #
#===============================================================#
HELP_plot  = '''
==================================================================
|               Information for option:  --plot                  |
==================================================================

 _________________________________________

 Command line:
 $ pilgrim.py --plot [show/pdf] [--dlevel]
 _________________________________________

 With this option, the data stored in the plot-file inside
 folder %s is plotted.

 The option may be accompanied by one of the next arguments:
   * 'show': the plots pop up in the screen (default)
   * 'pdf ': generates a pdf file with the plots
 If you are using Pilgrim in remote (via ssh) this option
 may be quite slow.

 Once executed, a list of plots is shown. Type the numbers
 of the plots to be selected.

 Depending on the system, the list of plots may be huge.
 In such cases, the following commands are useful:

    * select

      Lists the plots whose name contains the
      word 'keyword' just type:
      >> select keyword

      For example:
      >> select sct_int
      will only list the plots related to the SCT integrand

    * exclude

      Excludes the plots whose name contains the
      word 'keyword' just type:
      >> exclude keyword

    * reset

      To revert the effect of the select and exclude commands,
      type:
      >> reset

     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.


==================================================================
'''%(PN.DIR6)
#===============================================================#



#===============================================================#
# help string for --summary                                     #
#===============================================================#
HELP_summary = '''
==================================================================
|                 Information for option:  --summary             |
==================================================================

_____________________________

 Command line:
 $ pilgrim.py --summary [pfn/path/rcons [argument(s)]] [--dlevel]
______________________________

 With this option, Pilgrim reads the pfn/path/rcons output
 files and summarizes the most important data into tables.

 * Summary of files from --pfn execution:

   $ pilgrim.py --summary pfn ctcsp[.itc] [temperature(s)]

   > replace ctcsp with the selected system name
   > temperature(s) is(are) optional

   > information:
     - weight of each conformer
     - V0 for each conformer
     - V1 for each conformer
     - Partition functions for each conformer
     - Partition functions for the whole system
     - Gibbs free energy for each conformer
     - Gibbs free energy for the whole system

 * Summary of files from --path execution:

   $ pilgrim.py --summary path ctcts[.itc] [temperature(s)]

   > replace ctcts with the selected transition state
   > temperature(s) is(are) optional

   > information:
     - imaginary frequency
     - E0
     - maximum of VaG (VAG)
     - position of VAG (sAG)
     - MEP limits
     - SCT probability at E0
     - s_CVT & Gamma^CVT
     - kappa^SCT
     - gamma^CVT/SCT
     - Representative tunnelling energy for SCT
     - If a specific itc is selected, V_MEP, VaG and mueff
       along the MEP is also shown

 * Summary of files from --rcons execution:

   $ pilgrim.py --summary rcons labelrc [rcname(s)] [temperature(s)]

   > replace labelrc with the selected rate constant
   > reaction name(s) is(are) optional
   > temperature(s) is(are) optional

   > information:
     - barrier height (V0)
     - barrier height with zero-point energy correction (V1)
     - Gibbs free energies
     - rate constants

==================================================================
'''
#===============================================================#

#===============================================================#
# help string for --harvest                                     #
#===============================================================#
HELP_harvest = '''
==================================================================
|                 Information for option:  --harvest             |
==================================================================

_____________________________

 Command line:
 $ pilgrim.py --harvest [rcname] [--dlevel]
______________________________

   INPUT FILES:
       pif.chem

   OUTPUT FILES (in %s):
       summary.rcname.slevel.txt
     With --dlevel 
       summary.rcname.dlevel.txt

     With this option, 

     --dlevel
         Uses the dual-level data obtained by the ISPE algorithm.


==================================================================
'''
#===============================================================#

#===============================================================#
# help string for --input                                       #
#===============================================================#
MENU  = '''
 There are several variables ($var) and commands ($cmd) available
 in this interactive menu. The command line should have the 
 following syntax:

     > $cmd $var [$values]

 where the square brackets indicate that the $cmd-$var combination 
 may require the specification of values ($values)
 For more information, use the 'help' command on each
 variable inside the interactive menu.

 List of commands ($cmd) and variables ($var):

   ----------------------------------------------------------------
   $cmd\$var | struc | isomass | temp | chem | path | kmc | dlevel 
   ----------------------------------------------------------------
   help      | x     | x       |  x   |  x   |  x   |  x  |  x     
   ls        | x     | x       |  x   |  x   |  x   |  x  |  x
   add       |       | x+      |  x+  |  x+  |  x+  |  x+ |  x     
   mod       | x+    |         |      |      |  x+  |  x+ |        
   rm        | x+    | x       |  x+  |  x+  |  x+  |  x  |  x     
   ----------------------------------------------------------------
    x: the combination $cmd $var is available
    +: the combination $cmd $var requires $values

 Information about variables ($var):

   -----------------------------------------------------------------
    $var    | addresses...          | which contains...
   -----------------------------------------------------------------
    struc   | pif.struc             | structures & isot. masses
    isomass | pif.struc             | structures & isot. masses
    temp    | pif.temp              | temperatures
    chem    | pif.chem              | reactions     
    path    | pif.path & pif.calcs  | MEP parameters          
    kmc     | pif.kmc               | variables in the KMC
    dlevel  | pif.dlevel            | structures for high-level
   -----------------------------------------------------------------

 To go back to an upper level in the menu or to exit, use
 one of the next strings: end / .. / exit
'''
#---------------------------------------------------------------#
HELP_input  = '''
==================================================================
|                 Information for option:  --input               |
==================================================================

 Execution:
 $ pilgrim.py --input

 This option initializes an interactive menu which can be used
 to generate the Pilgrim input files (pifs).
 %s
==================================================================
'''%MENU
#---------------------------------------------------------------#
HELP_input_temp = '''
* 'add' : adds temperatures

  > add temp 100 200 300 400 550 780

  An equispaced number of temperatures can be included
  with 'range':

  > add temp range(100,200,20)

  which is equivalent to:

  > add temp 100 120 140 160 180 200

* 'rm' : removes temperatures

  > rm  temp 100 300 500 800

  To remove all the temperatures, use:

  > rm temp all

  To remove temperatures in an interval, use:

  > rm temp from $Ti to $Tf

  with $Ti and $Tf being the initial and final
  temperatures. For example:

  > rm temp from 100 to 500

* 'ls' : lists the previously defined temperatures 
  > ls temp
'''
#---------------------------------------------------------------#
HELP_input_chem = '''
* 'add' : creates a new reaction

  > add  $chemname : $reaction_equation

  where:
    - $chemname is just a name to identify the chemical reaction
    - $reaction_equation follows the structure:

          reactants --> transition state --> products

      Use '+' if there is more than one reactant (or product).
      If an element in the reaction is not defined in %s
      it will be ignored in the calculation of the rate
      constant, BUT it will be considered in the KMC simulation.

  Valid examples:
     > add myreaction1 : R1 + R2 --> TS --> P1 + P2
     > add myreaction2 : nh3_pyram --> nh3_planar --> nh3_pyram
     > add myreaction3 : R + R --> TS --> P1 + P2
     > add myreaction4 : R --> TS --> P1 + P2

* 'rm': removes a reaction

  > rm  chem  $chemname

* 'ls': lists all the available chemical reactions

  > ls  chem

'''%(PN.IFILE1)
#---------------------------------------------------------------#
HELP_input_path = '''
* 'add' : creates a new path:

  >  add  path  $targets

  where '$targets = ctcts.[idx]' OR
        '$target = *' all ctcts will be considered

    After adding a path the user enters a submenu, where
    some path variables can be modified. The syntax is:

    >> $pathvar = $value

    The following variables can be modified using this submenu:
        sbw, sfw, hsteps, ds, paral*, sctmns*, scterr*

    Variables with the * are not shown unless modified

* 'mod' : modifies a defined path

  > mod path $targets

* 'rm' : removes a path

  > rm path $target

* 'ls' : lists the previously defined paths

  > ls path
'''
#---------------------------------------------------------------#
HELP_input_kmc = '''
* 'add' : Initializes the KMC variables

  > add kmc kmcname

  This option initializes the KMC variables ($kmcvar), 
  whose values can be modified using the command line:

  >> $kmcvar = $value

  $kmcvar can be:
     psteps       | prints the number of molecules of each species 
                    after a given number of steps
     volume       | the volume of the reaction vessel (in cm**3)
     timeunits    | time evolution units: fs, ps, mcs, ms, s, min, hr
     pop0(ctcsp)  | initial number of particles of the 'ctcsp' species

  For example:

  >> volume = 2.00

  RATE CONSTANTS SPECIFICATION (inside 'add'):

    Command line for the rate constants calculated by Pilgrim:

      >> k($chemname[.$direc])[*$wgtreac] = $chemtype

      * $chemname is the name of the reaction,
  
      * $wgtreac = 1 by default; The user can set it up to 2 if there
         is a another chiral channel to the one is being considered

      * $direc is the direction of the reaction, which can be:

         - 'fw' for the forward  process (i.e. reactants --> products)

         - 'bw' for the backward process (i.e. products  --> reactants)

         - omitted; in such case, the syntax is just:

           >> k($chemname)[*$wgtreac] = $chemtype

            . bimolecular react. : only the 'fw' process is considered.
            . unimolecular react. : 'bw' and 'fw' processes are considered.

      * $chemtype is the type of rate constant.

         Valid values for the multi-structural rate constants are:
            tst    mscvt    mststzct    mscvtzct    mststsct    mscvtsct 

         Valid values for the multi-path are constants are:
                   mpcvt    mptstzct    mpcvtzct    mptstsct    mpcvtsct 

    Command line for the rate constants specified by analytical expressions:

      >> k($chemname.$direc) = analyticX $coefs

      X=1,2,3,4 and the coefficients ($coefs) have to included in the 
      following order:

        k($reaction_name.$direc) = analytic1 $A $B
        k($reaction_name.$direc) = analytic2 $A $B $n
        k($reaction_name.$direc) = analytic3 $A $B $n $Tr
        k($reaction_name.$direc) = analytic4 $A $B $n $Tr $T0

* 'mod' : modifies a defined KMC mechanism

  > mod kmc $kmcname

* 'rm': removes a KMC mechanism

  > rm  kmc  $kmcname

* 'ls' : list the kmc variables and reactions
  > ls kmc
'''
#---------------------------------------------------------------#
HELP_input_dlevel = '''
* 'add' : generates the high-level input files using the 
          command line:

  > add dlevel

  and the '%s' and '%s' files will be generated/modified.

* 'ls' : lists the structures which will be included in the
  high-level calculations

  >> ls dlevel

* 'rm' : removes a structure from the high-level calculation

  >> rm dlevel $target

NOTICE THAT:

  * %s should be modified to define the high-level method !!

  * Edit the file '%s' to include the MEP points 
    for which high-level calculation will be carried out

'''%(PN.IFILE4,PN.IFILE7,PN.IFILE4,PN.IFILE7)
#---------------------------------------------------------------#
HELP_input_struc = '''

* 'mod' : modifies a given structure using the command line:

  > mod struc $ctcsp

  After this command, the menu enters in a second environment
  which displays the variables to be modified. 
  The variables $var inside a $ctcsp can be modified 
  using the syntax:

  >> $var = $val

    The variables are:

    * root: points to the spname files defined by the user

    * freqscal: the frequency scale factor

    * weight idx : the weight of conformer with index idx.. 
        By default, all weights are set to 1. 
        However,  the user should switch it to 2 for conformers 
        with C1 symmetry that have another conformational enantiomer.
        To modify the weight of  a given conformer with index idx use:

           >> weight(idx) = $value

        To modify several weights at the same time, list them
        in parenthesis. For example:

           >> weight(001,003,004) = 2

        When  '-' is used, a range of conformers can be modified 
        at the same time. Thus, the command:

           >> weight(001-003) = 2

        is equivalent to:

           >> weight(001,002,003) = 2

        To modify the weight of all conformers having C1 point
        group symmetry, just type:

           >> weight(all) = 2

        ONLY weights of structures having C1 point group symmetry 
        can be modified!

    * anharfile: points towards an anharmonicity output file 
        of the MSTor or Q2DTor programs. Those files should
        be inside %s

    * iso: adds isotopic modifications:

        >> iso = $isomod(s)

        where each $isomod presents the following structure:
          
          $isomod(s) = $IMASS($atidxs),

        $IMASS is the isotopic mass.
        $atidxs is the numbering of the atoms with the 
        new $IMASS isotopic mass. 
        For example:

          >> iso C13(1,4,5)

        indicates that the mass of atoms 1, 4 and 5 is the one of C13.

        'all_X' : selects all the atoms of a given species
                  X being the symbol of an atom. For example, 
                  to consider deuterium in all hydrogen atoms, use:

                  >> iso D(all_H)

        To know how to define the $IMASS labels, check the command 
        isomass in the help menu:

          > help isomass 

  CREATING a NEW 'CTC-BLOCK' WITH A ISOTOPIC MODIFICATION

    For that, use the syntax:

    > mod struc $ctcsp

       >> copywith $isomod as $new_ctcsp

    For example, to copy the molecular hydrogen (H2) as D2 you can use:

    >  mod struc H2

       >> copywith D(all_H) as D2

  MODIFICATION OF THE VARIABLES OF SEVERAL 'CTC-BLOCKs' AT THE SAME TIME

    This can be done by using:

    > mod struc $ctcsp1 $ctcsp2 ...

    or if all are selected:

    > mod struc *

    Note that only the freqscal variable and the weights of the
    C1 symmetry conformers can be modified for several/all the targets
    at the same time.
  
* 'ls' : Lists the CTC-blocks

  > ls struc

* 'rm' : removes a CTC-block

  > rm struc $ctcsp

'''%PN.ANHDIR
#---------------------------------------------------------------#
HELP_input_isomass = '''
* 'add' : adds an isotopic mass

  >  add isomass label mass

  For example:

  > add isomass Deu  2.0141

* 'ls' : lists the isotopic masses

  > ls isomass

* 'rm' : removes a defined isotopic mass

  > rm isomass $target

  where '$target = label'

'''
#===============================================================#




INFO_pifstruc = '''
#----------------------------------------------------------------#
# Info about variables                                           #
#                                                                #
#     * root     : points towards the spname file of the root    #
#                  species defined by the user                   #
#     * conformer idx * weight : the index of the conformer and  # 
#                                its weight                      #
#                                it may content, as comment,     #
#                                the relative energy (kcal/mol)  #
#                                and the point group symmetry    #
#     * anharfile: anharmonicity file (Q2DTor/MSTor)             #
#     * mformu   : the molecular formula                         #
#     * ch       : the charge of the system                      #
#     * mtp      : the multiplicity of the system                #
#     * type     : 0 for minima, 1 for 1st order saddle points,  #
#                  et cetera                                     #
#     * freqcal  : scale factor for the vibrational frequencies  #
#     * elestate : each line like this adds a electronic state   #
#                  for the calculation of the eleectronic        #
#                  partition function; this keyword is followed  #
#                  by the multiplicity and the relative energy   #
#                  of the state (with regards to the             #
#                  ground-state). It is given in hartree.        #
#     * ics      : each line with this keyword contains one (or  #
#                  more) internal coordinates.                   #
#                  Stretching, angular bendings and proper       #
#                  torsions are indicated with '-':              #
#                      example:   1-2-3                          #
#                  linear bendings with '=':                     #
#                      example:   1=2=3                          #
#                  and improper torsions with '_' OR '-':        #
#                      example:   1_2_3_4                        #
#     * iso      : the isotopic modification                     #
#                                                                #
#----------------------------------------------------------------#\n'''

INFO_pifchem = '''
#-----------------------------------------------#
# To define a reaction, just use the next
# syntax:
#      chemname : chemreac
# where:
#     * 'chemname' is the label of the chemical reaction
#     * 'chemreac' is the chemical reaction
#
# The chemical reaction contains reactant(s),
# transition state and product(s) of the
# reaction, separated by arrows (-->).
# If there are more than one reactant/product,
# use the '+' symbol.
#
# Examples of valid lines:
#     myreac1 : R     --> TS --> P
#     myreac2 : R1+R2 --> TS --> P
#     myreac3 : R     --> TS --> P1+P2
#     myreac4 : R1+R2 --> TS --> P1+P2
#
#-----------------------------------------------#\n'''

