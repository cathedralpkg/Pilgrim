# level: Hartree-Fock(GTOs)  STO-3G
# Atomic number and non-scaled cartesian coordinates [bohr]
start_cc
   006   +5.81660000E-02  +1.05693900E+00  +1.03973000E-01
   006   -2.42485700E+00  -4.80307000E-01  -2.10440000E-02
   001   -4.03947000E+00  +7.76145000E-01  +1.30240000E-01
   001   -2.54890700E+00  -1.49820800E+00  -1.79788100E+00
   001   -2.51675300E+00  -1.83936000E+00  +1.51512200E+00
   001   +8.22780000E-02  +2.44360700E+00  -1.42825200E+00
   001   +1.22680000E-01  +2.10783900E+00  +1.89051700E+00
   008   +2.27485400E+00  -4.70803000E-01  -2.11984000E-01
   001   +2.20302200E+00  -1.67375000E+00  +1.22200000E+00
end_cc

# Charge, multiplicity, energy [hartree],
# point group and rotational symmetry number
start_basic
   charge        0
   multiplicity  1
   energy       -152.13306611    # Total energy in hartree
   pointgroup    C1              # Point group
   rotsigma      1               # Rotational sigma
end_basic

# Non-scaled cartesian gradient [hartree/bohr]
start_grad
   +1.19519400E-06  +5.17739000E-06  -1.17474200E-06
   -9.85763000E-07  -2.76614900E-06  +3.40069200E-06
   +4.26228200E-06  -2.09590500E-06  -1.09197000E-07
   -8.24612000E-07  -1.40618000E-06  -1.55387100E-06
   -3.19713500E-06  -1.27076600E-06  +3.63207000E-07
   +2.25244900E-06  +3.49554300E-06  -2.12324000E-06
   +3.68197800E-06  -3.08165700E-06  -1.29912000E-07
   -5.48588800E-06  +1.55555300E-06  +5.53600000E-09
   -8.98503000E-07  +3.92172000E-07  +1.32152700E-06
end_grad

# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]
# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...
start_hess
   +7.06759572E-01  -3.40070536E-02  +8.10890235E-01  -1.31472961E-02  -3.69564374E-02
   +8.78772545E-01  -2.38422963E-01  -9.23470775E-02  -1.24441543E-02  +7.51431675E-01
   -9.00818067E-02  -1.79853660E-01  -3.29655892E-03  -6.91874080E-02  +8.23628095E-01
   -6.80377748E-03  -4.30580649E-03  -1.27679475E-01  -6.94888013E-03  +6.53515763E-04
   +8.60828956E-01  -2.28823607E-02  +2.55972795E-02  +2.77792252E-03  -3.12383152E-01
   +1.90939910E-01  +2.22177428E-02  +3.23295254E-01  -2.16602064E-02  +1.87018676E-02
   +1.91417236E-03  +1.87446582E-01  -2.23850959E-01  -1.72978666E-02  -1.98560557E-01
   +2.26817992E-01  -2.05566448E-03  +1.75173400E-03  -5.51877183E-04  +2.24056320E-02
   -1.76034848E-02  -8.16584610E-02  -2.38737383E-02  +1.86724455E-02  +7.10567068E-02
   +2.36412107E-04  -2.03918030E-02  -3.72977540E-02  -7.71832888E-02  -8.05925867E-03
   -1.97713528E-02  +6.67880519E-03  +1.52944363E-02  +3.43516154E-02  +7.17904025E-02
   +1.25616020E-03  -1.04856946E-02  -1.98396119E-02  -7.72189384E-03  -1.73098806E-01
   -1.65018440E-01  -7.21223300E-03  -8.95007400E-03  -2.31869003E-02  +1.28879479E-02
   +1.71094570E-01  +5.51767036E-04  +3.46861418E-03  +4.91667985E-03  -1.72142210E-02
   -1.63460999E-01  -3.71813762E-01  -3.38250447E-03  +5.16577515E-03  +3.19716724E-03
   +1.94104614E-02  +1.76039320E-01  +3.81012047E-01  +7.05089758E-04  -2.72538499E-02
   +3.29657434E-02  -7.75112346E-02  -8.65427227E-03  +1.30883530E-02  +6.13435254E-03
   +2.20098973E-02  -3.09158674E-02  -3.51604365E-04  +2.02530863E-03  -1.32542757E-03
   +7.28355721E-02  +1.81551181E-03  -1.35787399E-02  +1.68786849E-02  -6.66992656E-03
   -2.44710809E-01  +1.88374583E-01  -7.25847896E-03  -1.22724348E-02  +2.05427437E-02
   +9.10191181E-04  +2.01561092E-02  -2.19167926E-02  +1.23584030E-02  +2.48322692E-01
   -6.65396029E-04  -6.49629258E-03  +6.65079116E-03  +9.91011934E-03  +1.87352165E-01
   -2.96095244E-01  +2.02805443E-03  -8.68387324E-03  +6.59438151E-03  +1.11996930E-03
   +3.05643001E-02  -1.91984774E-02  -1.12747239E-02  -2.01978848E-01  +3.00326829E-01
   -7.37512553E-02  -9.11503167E-03  +9.42597256E-03  +3.48162583E-05  -2.59894697E-02
   +3.08005973E-02  +2.04148918E-03  +9.28192942E-04  -7.12359050E-05  +1.72335383E-03
   +1.29659579E-03  -1.91793683E-04  -6.46737131E-03  -3.56516818E-03  -1.95518710E-03
   +8.03404842E-02  -5.59355233E-03  -2.39552604E-01  +1.80225379E-01  +5.54409727E-05
   -1.17791268E-02  +1.73311348E-02  +5.30427912E-04  +9.41487829E-04  -4.04381323E-04
   +1.65318121E-03  +1.15701144E-03  +1.17758644E-04  -3.06709157E-03  -6.14704750E-04
   -6.17723047E-04  -8.11145532E-04  +2.50942425E-01  +5.51427856E-03  +1.81239953E-01
   -2.93201908E-01  -1.88373728E-03  -6.52597820E-03  +6.69830610E-03  +1.88838545E-04
   -6.97795973E-05  +1.31755805E-04  +3.11597246E-04  -6.25858833E-05  +9.20967112E-04
   -1.69708288E-03  -9.43414513E-04  +7.76415772E-04  -3.56285542E-03  -1.95955330E-01
   +2.92259989E-01  -7.21722793E-02  -8.45373383E-03  -1.50371527E-02  -2.17545650E-04
   -1.99313278E-02  -3.61236258E-02  +2.29877255E-03  +1.00366222E-03  +2.86785578E-04
   -6.65872504E-03  -3.99209754E-03  +1.58145993E-03  +1.61517297E-03  +1.19795101E-03
   +4.47522708E-04  +1.79426178E-03  +3.49388109E-05  +2.01756699E-03  +8.20519113E-02
   -1.27973201E-02  -1.68170540E-01  -1.49218905E-01  -2.59594151E-03  -8.77363310E-03
   -2.18312997E-02  +4.19423140E-04  +6.05058494E-04  +4.89449765E-04  -3.70665452E-03
   -8.57921591E-04  +6.08649974E-04  +1.84891486E-03  +1.31298779E-03  -9.34431504E-05
   +1.69424802E-03  +1.91613334E-02  +2.86780193E-02  +3.88426250E-03  +1.70027554E-01
   -1.64602835E-02  -1.47189452E-01  -3.51027023E-01  +2.63438321E-03  +3.87755513E-03
   +3.72111936E-03  -1.89276543E-04  +1.65727582E-04  +2.17183551E-04  +1.08338758E-03
   +5.31899574E-04  +1.03396933E-03  +2.93464300E-05  +3.97736493E-05  +9.12499374E-04
   -2.65225673E-03  -2.08463453E-02  -1.94757868E-02  +1.03120569E-02  +1.66109019E-01
   +3.59641663E-01  -2.91979475E-01  +1.19011164E-01  +8.11289277E-02  -4.81532334E-02
   +3.02791364E-02  +5.00344160E-03  -5.62331911E-03  -6.40898017E-03  -3.05703013E-04
   +4.03758546E-03  +1.63841613E-03  +5.88241900E-04  +3.32811053E-03  +1.56164482E-03
   -3.77867032E-05  +1.97539204E-03  +2.33476338E-03  -3.58642724E-03  -1.09801943E-02
   +1.35548196E-02  +5.30445863E-03  +4.02512764E-01  +1.73060988E-01  -1.93609533E-01
   -2.49401824E-02  -6.66944284E-03  +1.49316710E-02  +3.32552340E-03  -3.70137877E-03
   -1.60202466E-03  -1.59716200E-05  +1.54700068E-03  +1.07168663E-03  +2.20556858E-04
   +1.67261518E-05  +8.13734891E-04  -2.66345693E-04  +2.90185688E-02  -1.80297015E-02
   -4.06444950E-03  +2.84280214E-02  -1.48918679E-02  -3.43574882E-03  -1.92896183E-01
   +4.77772007E-01  +1.40182159E-02  +2.88480833E-02  -1.20231978E-01  +2.89424887E-03
   -1.08199685E-04  +5.90320759E-03  -7.55104545E-05  +2.91151223E-04  +8.70545900E-04
   +8.95302358E-04  +8.54801418E-04  +1.07002260E-04  -5.11848548E-04  -8.85864269E-04
   -1.37133391E-04  -3.30387594E-02  +2.08005258E-02  +1.02796553E-02  +3.59662463E-02
   -2.37250279E-02  +2.39599754E-03  -1.00469683E-01  -2.96557992E-01  +4.87901387E-01
   -8.50900943E-03  +4.69365786E-02  -4.83816606E-02  +2.39878141E-03  +6.87393735E-04
   -1.46173619E-03  +4.39802282E-04  -5.24964785E-05  +1.78352466E-04  -2.72631391E-04
   -1.77998700E-04  -1.88805067E-05  -2.87888142E-04  -3.50213751E-04  +4.28083434E-04
   -7.69086704E-03  +4.86411786E-03  +2.69711605E-03  +2.27252512E-03  -2.30272484E-03
   -6.07832594E-05  -5.50958902E-02  -2.87736019E-02  +8.03230791E-02  +6.67414961E-02
   -1.19876504E-02  -2.43582940E-02  +3.52072027E-02  -2.31409403E-03  +3.49969836E-03
   -1.22812024E-03  -7.56489135E-04  -3.89489574E-04  -2.45395532E-04  -1.33805204E-04
   -8.58726794E-05  -2.44822072E-04  +7.18137766E-04  +5.72857104E-04  +2.21579033E-04
   +6.53979112E-03  -2.22322424E-03  -2.29877909E-03  -2.16732936E-03  +1.58767853E-03
   +7.46880231E-04  +3.09282833E-02  -2.66440524E-01  +2.70520564E-01  -2.08376702E-02
   +2.87838543E-01  +1.90310503E-02  -2.03868691E-02  +2.30258323E-03  +6.40978590E-04
   -8.87246321E-04  +9.69631294E-05  +3.10904023E-04  -1.61479180E-04  +1.41875314E-04
   -1.02587509E-04  +1.21300282E-04  -1.70665437E-04  -3.57615223E-04  -1.08955604E-04
   +1.67156876E-04  +1.24559007E-03  -6.51240725E-04  +1.61148134E-03  +5.55425483E-04
   -1.01823638E-03  +2.58406173E-03  +1.23850045E-02  +3.25753831E-01  -3.87031586E-01
   -3.37016148E-02  -3.02672899E-01  +3.80283102E-01
end_hess

