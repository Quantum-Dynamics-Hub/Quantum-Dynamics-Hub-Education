

     ******************************************
     **    PROGRAM:              MCPC        **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


 original author: Daniel Robertson, FSU
 later revisions: Ron Shepard, ANL;
                  Michal Dallos, University Vienna



 This Version of Program mcpc is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



   ******  File header section  ******

 Headers form the restart file:
    Hermit Integral Program : SIFS version  compute-0-37      13:33:52.636 23-Jun-21
     title                                                                          
     title                                                                          
     title                                                                          


   ******  DRT info section  ******

 Informations for the DRT no.  1
 Header form the DRT file: 
     title                                                                          
 Molecular symmetry group:   sym1
 Total number of electrons:   42
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        55

   ***  Informations from the DRT number:   1

 
 Symmetry orbital summary:
 Symm.blocks:         1     2     3     4     5     6     7     8
 Symm.labels:         ag    b3u   b2u   b1g   b1u   b2g   b3g   au 

 List of doubly occupied orbitals:
  1 ag   2 ag   3 ag   4 ag   5 ag   6 ag   1 b3u  2 b3u  3 b3u  4 b3u  5 b3u  1 b2u
  2 b2u  3 b2u  4 b2u  1 b1g  2 b1g  3 b1g

 List of active orbitals:
  1 b1u  2 b1u  3 b1u  1 b2g  1 b3g  1 au 

 Informations for the DRT no.  2
 Header form the DRT file: 
     title                                                                          
 Molecular symmetry group:    b2u
 Total number of electrons:   42
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        40

   ***  Informations from the DRT number:   2

 
 Symmetry orbital summary:
 Symm.blocks:         1     2     3     4     5     6     7     8
 Symm.labels:         ag    b3u   b2u   b1g   b1u   b2g   b3g   au 

 List of doubly occupied orbitals:
  1 ag   2 ag   3 ag   4 ag   5 ag   6 ag   1 b3u  2 b3u  3 b3u  4 b3u  5 b3u  1 b2u
  2 b2u  3 b2u  4 b2u  1 b1g  2 b1g  3 b1g

 List of active orbitals:
  1 b1u  2 b1u  3 b1u  1 b2g  1 b3g  1 au 

 Informations for the DRT no.  3
 Header form the DRT file: 
     title                                                                          
 Molecular symmetry group:    b3u
 Total number of electrons:   42
 Spin multiplicity:            3
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        48

   ***  Informations from the DRT number:   3

 
 Symmetry orbital summary:
 Symm.blocks:         1     2     3     4     5     6     7     8
 Symm.labels:         ag    b3u   b2u   b1g   b1u   b2g   b3g   au 

 List of doubly occupied orbitals:
  1 ag   2 ag   3 ag   4 ag   5 ag   6 ag   1 b3u  2 b3u  3 b3u  4 b3u  5 b3u  1 b2u
  2 b2u  3 b2u  4 b2u  1 b1g  2 b1g  3 b1g

 List of active orbitals:
  1 b1u  2 b1u  3 b1u  1 b2g  1 b3g  1 au 


   ******  MCSCF convergence information:  ******

 MCSCF convergence criteria were satisfied.

 mcscf energy=  -230.5962673335    nuclear repulsion=   201.5958396276
 demc=             0.0000000000    wnorm=                 0.0000002903
 knorm=            0.0000000233    apxde=                 0.0000000000


 MCSCF calculation performmed for   3 symmetries.

 State averaging:
 No,  ssym, navst, wavst
  1    ag     1   0.2000
  2    b2     2   0.2000 0.2000
  3    b3     2   0.2000 0.2000

 Input the DRT No of interest: [  1]:
In the DRT No.: 1 there are  1 states.

 Which one to take? [  1]:
 The CSFs for the state No  1 of the symmetry  ag  will be printed
 according to the following print options :

 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: csfs will be printed based on coefficient magnitudes.

 input the coefficient threshold (end with 0.) [ 0.0000]:
 List of active orbitals:
  1 b1u  2 b1u  3 b1u  1 b2g  1 b3g  1 au 

   csf       coeff       coeff**2    step(*)
  -----  ------------  ------------  ------------
     15  0.9673664624  0.9357978726  300330
     16 -0.1241336931  0.0154091738  300303
      3 -0.1178522709  0.0138891578  330030
      9 -0.1121216856  0.0125712724  310122
      8  0.0766938077  0.0058819401  310212
     14 -0.0685795512  0.0047031548  301122
     17 -0.0495763723  0.0024578167  300033
      6 -0.0437523384  0.0019142671  312030
      5 -0.0363662704  0.0013225056  312300
      2 -0.0357776348  0.0012800392  330300
     55 -0.0332487640  0.0011054803  000333
     44 -0.0299390190  0.0008963449  030330
      4  0.0250322653  0.0006266143  330003
     10 -0.0204115585  0.0004166317  303300
     13  0.0142588121  0.0002033137  301212
     49 -0.0138474094  0.0001917507  012330
      7  0.0131234040  0.0001722237  312003
     36  0.0120922679  0.0001462229  102330
     30  0.0120472788  0.0001451369  120033
     22  0.0118895945  0.0001413625  130122
     29 -0.0114157199  0.0001303187  120303
     11 -0.0107399918  0.0001153474  303030
     21  0.0102920330  0.0001059259  130212
     52 -0.0096252103  0.0000926447  003330
     32  0.0077859798  0.0000606215  112122
     46  0.0064525475  0.0000416354  030033
     37 -0.0060277805  0.0000363341  102303
     45  0.0052400252  0.0000274579  030303
     12  0.0047157240  0.0000222381  303003
     27  0.0044916075  0.0000201745  121122
     35  0.0040028915  0.0000160231  103122
     50  0.0037488735  0.0000140541  012303
     33 -0.0034455896  0.0000118721  111222
     38  0.0024847134  0.0000061738  102033
     53  0.0023899247  0.0000057117  003303
     18  0.0020114361  0.0000040459  132300
     51  0.0020098765  0.0000040396  012033
     43  0.0016763578  0.0000028102  031122
     19  0.0014472542  0.0000020945  132030
     24  0.0013175483  0.0000017359  123030
      1  0.0011742253  0.0000013788  333000

 input the coefficient threshold (end with 0.) [ 0.0000]:
 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: