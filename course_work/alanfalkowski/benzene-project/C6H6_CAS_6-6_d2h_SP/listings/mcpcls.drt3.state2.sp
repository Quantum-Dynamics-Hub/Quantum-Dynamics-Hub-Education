

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
In the DRT No.: 3 there are  2 states.

 Which one to take? [  1]:
 The CSFs for the state No  2 of the symmetry  b2  will be printed
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
      7  0.6834315613  0.4670786990  310130
     12  0.6670454568  0.4449496414  300311
     29  0.1894173491  0.0358789321  110321
     10  0.1545029381  0.0238711579  301130
     17 -0.0950933668  0.0090427484  130130
     34  0.0856308183  0.0073326370  101321
      8 -0.0487142779  0.0023730809  310103
      2 -0.0473630948  0.0022432627  330011
     46 -0.0449676225  0.0020220871  010133
     40 -0.0364154939  0.0013260882  030311
     43 -0.0261591061  0.0006842988  012311
     25 -0.0227739466  0.0005186526  112130
      5  0.0222014684  0.0004929052  311021
     35 -0.0211644019  0.0004479319  101312
     18  0.0171722771  0.0002948871  130103
     20 -0.0168171083  0.0002828151  121130
     47 -0.0142629837  0.0002034327  003311
      6 -0.0137326246  0.0001885850  311012
     48 -0.0129794765  0.0001684668  001133
     33  0.0113505696  0.0001288354  102311
      1  0.0100536672  0.0001010762  331100
      3  0.0091341531  0.0000834328  313100
     36 -0.0064194873  0.0000412098  100133
     30 -0.0062884688  0.0000395448  110312
     27  0.0061036941  0.0000372551  111230
     22 -0.0059569901  0.0000354857  120311
     41  0.0057954170  0.0000335869  013130
     21  0.0048346158  0.0000233735  121103
     26  0.0041803557  0.0000174754  112103
     45  0.0032416178  0.0000105081  011312
     44 -0.0024919973  0.0000062101  011321
     14  0.0024475926  0.0000059907  132011
      9  0.0024364233  0.0000059362  303011
     19  0.0024303309  0.0000059065  123011
     11  0.0024001234  0.0000057606  301103
     15  0.0017722288  0.0000031408  131021
     16 -0.0016521505  0.0000027296  131012
     23 -0.0015910308  0.0000025314  113021
     31  0.0015543592  0.0000024160  103130
     42 -0.0014873707  0.0000022123  013103
     39 -0.0013500872  0.0000018227  031103
     38  0.0011975949  0.0000014342  031130

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