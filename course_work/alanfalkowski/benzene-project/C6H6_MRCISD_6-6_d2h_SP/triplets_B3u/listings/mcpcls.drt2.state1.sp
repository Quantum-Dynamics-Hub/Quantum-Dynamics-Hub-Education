

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
    Hermit Integral Program : SIFS version  compute-0-12      14:54:08.910 24-Jun-21
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
In the DRT No.: 2 there are  2 states.

 Which one to take? [  1]:
 The CSFs for the state No  1 of the symmetry  b3  will be printed
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
     11 -0.6630346004  0.4396148813  300132
      6  0.6519019918  0.4249762069  310320
      9  0.1969446109  0.0387871798  301320
     24  0.1685950570  0.0284242932  110232
     29 -0.1412537665  0.0199526265  100323
     15 -0.1233481956  0.0152147774  130320
     20 -0.0990861295  0.0098180611  120132
     38 -0.0681340906  0.0046422543  010323
     33  0.0646963135  0.0041856130  030132
     22 -0.0594490078  0.0035341845  112320
     28  0.0492804929  0.0024285670  101232
     18 -0.0435235988  0.0018943037  121320
      2  0.0385872282  0.0014889742  330102
     27 -0.0320939027  0.0010300186  102132
     36  0.0310670756  0.0009651632  012132
     10 -0.0269478161  0.0007261848  301023
      7 -0.0268790566  0.0007224837  310023
     40 -0.0258162210  0.0006664773  001323
     25 -0.0179026694  0.0003205056  103320
      5 -0.0152922736  0.0002338536  311202
     39  0.0111661050  0.0001246819  003132
      4  0.0093715889  0.0000878267  312102
     34  0.0061723765  0.0000380982  013320
      1  0.0061689697  0.0000380562  331020
      8 -0.0047723842  0.0000227757  303102
     19 -0.0037408862  0.0000139942  121023
     31  0.0034343682  0.0000117949  031320
     14  0.0034186125  0.0000116869  131202
     21  0.0032245225  0.0000103975  113202
     37 -0.0031030619  0.0000096290  011232
     13  0.0010241469  0.0000010489  132102
     17 -0.0010105141  0.0000010211  123102

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