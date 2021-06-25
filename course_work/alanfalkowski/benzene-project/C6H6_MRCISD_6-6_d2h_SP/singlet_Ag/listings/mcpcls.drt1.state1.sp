

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
    Hermit Integral Program : SIFS version  compute-0-12      14:50:12.965 24-Jun-21
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
 demc=            -0.0000000000    wnorm=                 0.0000004468
 knorm=            0.0000000129    apxde=                 0.0000000000


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
     15 -0.9673664210  0.9357977925  300330
     16  0.1241337564  0.0154091895  300303
      3  0.1178486564  0.0138883058  330030
      9  0.1121177228  0.0125703838  310122
      8 -0.0766930814  0.0058818287  310212
     14  0.0685861581  0.0047040611  301122
     17  0.0495764706  0.0024578264  300033
      6  0.0437612979  0.0019150512  312030
      5  0.0363677591  0.0013226139  312300
      2  0.0357745075  0.0012798154  330300
     55  0.0332487607  0.0011054801  000333
     44  0.0299377210  0.0008962671  030330
      4 -0.0250311410  0.0006265580  330003
     10  0.0204148761  0.0004167672  303300
     13 -0.0142634300  0.0002034454  301212
     49  0.0138491836  0.0001917999  012330
      7 -0.0131251819  0.0001722704  312003
     36 -0.0120919379  0.0001462150  102330
     30 -0.0120471054  0.0001451327  120033
     22 -0.0118886981  0.0001413411  130122
     29  0.0114153428  0.0001303101  120303
     11  0.0107435903  0.0001154247  303030
     21 -0.0102920385  0.0001059261  130212
     52  0.0096265424  0.0000926703  003330
     32 -0.0077866741  0.0000606323  112122
     46 -0.0064523238  0.0000416325  030033
     37  0.0060284511  0.0000363422  102303
     45 -0.0052396403  0.0000274538  030303
     12 -0.0047167539  0.0000222478  303003
     27 -0.0044920999  0.0000201790  121122
     35 -0.0040037134  0.0000160297  103122
     50 -0.0037492845  0.0000140571  012303
     33  0.0034453973  0.0000118708  111222
     38 -0.0024854261  0.0000061773  102033
     53 -0.0023903085  0.0000057136  003303
     18 -0.0020113902  0.0000040457  132300
     51 -0.0020103045  0.0000040413  012033
     43 -0.0016761909  0.0000028096  031122
     19 -0.0014470659  0.0000020940  132030
     24 -0.0013176628  0.0000017362  123030
      1 -0.0011741713  0.0000013787  333000

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