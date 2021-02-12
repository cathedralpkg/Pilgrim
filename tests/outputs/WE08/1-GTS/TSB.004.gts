# level: Hartree-Fock(GTOs)  STO-3G
# Atomic number and non-scaled cartesian coordinates [bohr]
start_cc
   006   +2.02657000E+00  -1.22741200E+00  -4.00000000E-06
   006   +1.91598000E-01  +1.02004600E+00  -1.00000000E-06
   008   -2.34780900E+00  +5.27870000E-02  +1.50000000E-05
   001   +4.29886400E+00  -3.39917000E-01  -2.00000000E-06
   001   +1.87544700E+00  -2.37854800E+00  +1.69159800E+00
   001   +1.87544800E+00  -2.37854500E+00  -1.69160700E+00
   001   +5.40736000E-01  +2.19834500E+00  +1.67024100E+00
   001   +5.40719000E-01  +2.19833500E+00  -1.67025200E+00
   001   -3.41421600E+00  +1.59285800E+00  +1.30000000E-05
   001   +5.93477700E+00  +3.03976000E-01  -1.00000000E-06
end_cc

# Charge, multiplicity, energy [hartree],
# point group and rotational symmetry number
start_basic
   charge        0
   multiplicity  2
   energy       -152.56374374    # Total energy in hartree
   pointgroup    Cs              # Point group
   rotsigma      1               # Rotational sigma
end_basic

# Non-scaled cartesian gradient [hartree/bohr]
start_grad
   -2.68831365E-04  +1.49635706E-04  +2.28353000E-07
   +4.60583440E-05  -1.33305295E-04  +6.35117000E-07
   -1.28337087E-04  +1.45438612E-04  -2.74040000E-08
   +2.12185008E-04  +1.12547592E-04  -8.45700000E-09
   -2.83878870E-05  -5.10098750E-05  -3.62297180E-05
   -2.83065670E-05  -5.12855870E-05  +3.58758020E-05
   +9.82122500E-06  +1.04325960E-05  -7.95906200E-06
   +9.86513900E-06  +1.06074090E-05  +7.50103500E-06
   +1.71559534E-04  -1.77241069E-04  -1.94540000E-08
   +4.37365500E-06  -1.58200880E-05  +3.79000000E-09
end_grad

# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]
# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...
start_hess
   +3.21181401E-01  -1.16370238E-01  +6.91689954E-01  -2.76785752E-07  +1.66101668E-07
   +8.47541808E-01  -1.62392529E-01  +1.17147771E-01  +6.98216810E-08  +7.24388623E-01
   +1.03463984E-01  -2.40682160E-01  -1.69772229E-07  +2.46912889E-02  +7.98818393E-01
   +2.53699826E-08  -8.43539954E-08  -1.28520184E-01  +4.17888595E-07  -4.10489048E-07
   +8.32163140E-01  -3.71632423E-02  +1.43103278E-03  +1.94742439E-07  -3.53956862E-01
   -1.53342678E-01  +1.91636688E-06  +6.72334415E-01  +3.38425355E-02  +1.61753981E-02
   -2.34695090E-07  -3.90167485E-02  -1.64052735E-01  +6.04771707E-07  -1.36875547E-01
   +6.28481430E-01  +8.74983038E-08  -4.92727666E-08  +4.07968756E-03  +1.42262311E-06
   +1.11683487E-06  -7.26816417E-02  -2.12476202E-06  -2.00162596E-06  +5.39782170E-02
   +1.53884139E-01  +7.01988624E-02  +6.71213348E-08  -3.11477785E-02  -1.68760869E-02
   +7.14000256E-08  +8.66067395E-04  +8.41227584E-03  -4.93729164E-08  -1.88486035E-01
   +7.73964646E-02  -2.34155354E-02  +3.33852633E-08  +1.26944076E-02  +4.73307324E-03
   +1.77781398E-09  +2.03961170E-03  -1.76112225E-03  -2.51853095E-09  -1.00869448E-01
   +3.34372891E-02  +1.27449826E-07  +5.10152297E-08  -5.31111623E-02  +1.58000837E-08
   -8.66429753E-09  -3.28220779E-04  +9.14849908E-09  -2.12596766E-09  +4.04235861E-04
   -1.77114577E-07  -3.85426227E-08  +7.02068617E-02  -5.02770815E-02  -1.54431750E-02
   +3.75380848E-02  +5.09532117E-03  +1.76886976E-02  -2.68772495E-02  +2.41424231E-03
   -1.51997063E-03  +7.59367121E-04  -1.79955269E-02  -1.96536224E-02  +2.04156669E-02
   +5.31748229E-02  -7.19028497E-03  -1.96960813E-01  +1.82919299E-01  -6.86055031E-03
   -1.43890316E-02  +2.75149922E-02  -1.12446989E-03  +1.45743856E-03  -9.34123475E-04
   -3.53318352E-03  -5.46149515E-03  +6.82049571E-03  +9.84632469E-03  +1.96454443E-01
   +2.65939748E-02  +1.76785291E-01  -3.45771465E-01  +3.31324543E-04  -5.59259198E-03
   +6.83395876E-03  +2.40931171E-04  -1.04597134E-04  +5.31985501E-05  -6.62635846E-03
   -4.38317951E-03  +2.30853555E-03  -2.82334734E-02  -1.90273202E-01  +3.51104284E-01
   -5.02770100E-02  -1.54430285E-02  -3.75380345E-02  +5.09539570E-03  +1.76886606E-02
   +2.68772656E-02  +2.41418549E-03  -1.51995175E-03  -7.59389054E-04  -1.79955782E-02
   -1.96535955E-02  -2.04157280E-02  +2.10615831E-03  +3.22976512E-03  +4.91238789E-03
   +5.31747544E-02  -7.19018828E-03  -1.96960115E-01  -1.82919216E-01  -6.86059251E-03
   -1.43890158E-02  -2.75149661E-02  -1.12442847E-03  +1.45740130E-03  +9.34119467E-04
   -3.53321923E-03  -5.46149362E-03  -6.82052168E-03  +3.22977485E-03  +1.89430862E-02
   +2.44901465E-02  +9.84622165E-03  +1.96453747E-01  -2.65939287E-02  -1.76785264E-01
   -3.45772519E-01  -3.31366552E-04  +5.59258503E-03  +6.83390437E-03  -2.40934680E-04
   +1.04609870E-04  +5.32063678E-05  +6.62635214E-03  +4.38315607E-03  +2.30852484E-03
   -4.91238111E-03  -2.44900663E-02  -1.66681281E-02  +2.82334530E-02  +1.90273089E-01
   +3.51105365E-01  +5.48421527E-03  +1.65440970E-02  +2.62981624E-02  -7.90590085E-02
   -2.78325480E-02  -4.56442319E-02  -2.40032582E-02  -3.33353273E-02  -3.31551971E-02
   +9.33894172E-04  -3.40416651E-04  +2.03950350E-04  +6.46610191E-04  -1.14711509E-03
   -3.91569712E-04  -3.68395076E-03  +4.73644893E-03  -1.55109772E-03  +9.11116762E-02
   -4.09868733E-03  -1.57418173E-02  -2.91855285E-02  -2.99172537E-02  -1.90747900E-01
   -1.60789322E-01  -1.09296864E-02  -8.11796397E-03  -1.20162372E-02  -7.77422180E-04
   +5.43982203E-04  +2.09924579E-04  -1.47603285E-03  +2.13863861E-03  +1.60697781E-04
   +4.11577084E-03  -3.38967026E-03  +1.30741562E-03  +3.70805163E-02  +1.95927391E-01
   -2.10772367E-04  +4.96445075E-03  +5.72993623E-03  -4.10545982E-02  -1.56397953E-01
   -3.20593285E-01  -2.77816680E-03  +6.12966452E-05  +8.91964288E-03  +9.95215384E-06
   +4.26573595E-05  +8.13310546E-05  +2.18283806E-04  +2.21809596E-05  +9.61361100E-04
   -1.33298051E-03  +1.06556581E-03  +8.09439187E-04  +5.01470257E-02  +1.74014342E-01
   +3.24787317E-01  +5.48403496E-03  +1.65439469E-02  -2.62982437E-02  -7.90581428E-02
   -2.78306939E-02  +4.56419207E-02  -2.40028946E-02  -3.33351354E-02  +3.31557463E-02
   +9.33905347E-04  -3.40430321E-04  -2.03936056E-04  -3.68394900E-03  +4.73642411E-03
   +1.55110876E-03  +6.46633415E-04  -1.14713427E-03  +3.91575143E-04  +4.50360151E-03
   +4.79201303E-03  -5.11896831E-03  +9.11105852E-02  -4.09846348E-03  -1.57417630E-02
   +2.91856667E-02  -2.99153949E-02  -1.90746099E-01  +1.60789036E-01  -1.09295089E-02
   -8.11790619E-03  +1.20164670E-02  -7.77381628E-04  +5.43970956E-04  -2.09925117E-04
   +4.11579886E-03  -3.38966700E-03  -1.30740852E-03  -1.47605683E-03  +2.13866165E-03
   -1.60677943E-04  +4.79178242E-03  +1.77253275E-02  -2.47636145E-02  +3.70784598E-02
   +1.95925475E-01  +2.10830573E-04  -4.96453560E-03  +5.73005681E-03  +4.10523289E-02
   +1.56397748E-01  -3.20596289E-01  +2.77855989E-03  -6.08670109E-05  +8.91918235E-03
   -9.96523351E-06  -4.26590348E-05  +8.13245909E-05  +1.33300079E-03  -1.06560464E-03
   +8.09428761E-04  -2.18278404E-04  -2.21941754E-05  +9.61379168E-04  +5.11866974E-03
   +2.47633312E-02  -2.13504759E-02  -5.01448882E-02  -1.74014145E-01  +3.24790643E-01
   -8.51531188E-03  -5.02516185E-03  +6.69815780E-08  -3.73245139E-02  +5.97108281E-02
   -1.11669355E-07  -2.38651237E-01  +2.05408392E-01  +1.21185555E-07  -2.16359759E-03
   +6.62200767E-04  +3.51296424E-09  +5.35349181E-04  -4.22686111E-04  +1.06538440E-04
   +5.35338000E-04  -4.22670173E-04  -1.06541094E-04  +3.95085875E-03  +1.00282553E-03
   +9.88166627E-05  +3.95088083E-03  +1.00281889E-03  -9.88525316E-05  +2.77182113E-01
   -2.72243300E-03  -3.41765617E-04  +2.11436005E-08  -4.49534815E-02  +9.96571230E-03
   +1.59688998E-07  +3.11253888E-01  -4.64612539E-01  +7.02302515E-07  -5.59596353E-04
   +2.54337620E-04  -1.26166360E-09  -2.25399132E-04  +1.90486269E-04  -2.10465225E-05
   -2.25395843E-04  +1.90482326E-04  +2.10417876E-05  -3.59003453E-04  +1.53840197E-03
   +9.33550569E-04  -3.59015425E-04  +1.53839730E-03  -9.33553675E-04  -2.62058998E-01
   +4.51127735E-01  +7.13852988E-08  +3.11909089E-08  +1.02419780E-03  +3.80461872E-07
   -3.25124308E-07  -2.98255063E-03  -3.89762636E-07  +1.20866166E-06  -3.73488080E-03
   +1.17846162E-08  -6.06957620E-09  +1.64247599E-04  +5.52359572E-05  -2.97592291E-05
   +4.17731389E-06  -5.52451521E-05  +2.97601713E-05  +4.17802914E-06  -9.14473178E-04
   +1.50253157E-03  +6.37352040E-04  +9.14447559E-04  -1.50255429E-03  +6.37365689E-04
   -3.71139020E-08  -8.86371036E-07  +4.26008021E-03  -1.77395532E-01  -6.95896200E-02
   -9.08901533E-08  +8.35655680E-03  +2.64074226E-03  -2.42392870E-08  -2.45280180E-04
   -2.05329549E-03  +1.60731054E-08  +1.01133417E-01  +4.80480136E-02  +6.79360423E-08
   +7.98742252E-03  +2.46918951E-03  +1.51364201E-03  +7.98744292E-03  +2.46920104E-03
   -1.51363649E-03  +1.13099927E-04  +2.08952997E-04  +2.16078849E-05  +1.13083663E-04
   +2.08941988E-04  -2.16054847E-05  +5.02314890E-04  +2.06800142E-04  -2.41532637E-09
   +5.14638771E-02  -7.30313798E-02  -1.80181744E-02  -3.75080360E-08  +2.98292814E-03
   +1.47084464E-03  -1.13656658E-08  -3.98763167E-04  -8.93715768E-04  +8.60204088E-09
   +4.83130024E-02  -3.41423962E-03  +2.60782161E-08  +3.43868900E-03  +1.01544573E-03
   +2.47867658E-04  +3.43869531E-03  +1.01544745E-03  -2.47865966E-04  -1.38807397E-04
   +1.27773772E-04  +6.02251982E-05  -1.38807866E-04  +1.27766790E-04  -6.02210785E-05
   +1.48835942E-04  +1.43861412E-04  -8.62200620E-10  +1.53923280E-02  +1.84257468E-02
   -1.38423955E-07  -5.65493901E-08  +9.07891985E-03  +5.10897935E-09  +6.75983501E-09
   -1.19954102E-04  -2.84894682E-10  -2.35823369E-09  -1.41117657E-05  +9.62476792E-08
   +3.74846800E-08  -2.21183133E-02  -2.95079025E-04  -4.85469309E-04  +3.63148000E-04
   +2.95092056E-04  +4.85472685E-04  +3.63146780E-04  -1.10264011E-04  +3.36337157E-05
   +1.77116128E-05  +1.10262297E-04  -3.36318521E-05  +1.77120901E-05  +4.55124567E-10
   +1.60823331E-10  -4.36203556E-06  +2.55931996E-08  +9.26436329E-09  +1.24163876E-02
end_hess

