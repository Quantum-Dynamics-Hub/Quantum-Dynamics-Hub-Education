

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
      7 -0.6834223955  0.4670661707  310130
     12 -0.6670457374  0.4449500158  300311
     29 -0.1894116205  0.0358767620  110321
     10 -0.1545440652  0.0238838681  301130
     17  0.0950904392  0.0090421916  130130
     34 -0.0856415073  0.0073344678  101321
      8  0.0487143076  0.0023730838  310103
      2  0.0473626972  0.0022432251  330011
     46  0.0449668494  0.0020220175  010133
     40  0.0364130009  0.0013259066  030311
     43  0.0261612009  0.0006844084  012311
     25  0.0227808466  0.0005189670  112130
      5 -0.0222016478  0.0004929132  311021
     35  0.0211642059  0.0004479236  101312
     18 -0.0171716592  0.0002948659  130103
     20  0.0168211365  0.0002829506  121130
     47  0.0142651999  0.0002034959  003311
      6  0.0137323820  0.0001885783  311012
     48  0.0129822513  0.0001685388  001133
     33 -0.0113496530  0.0001288146  102311
      1 -0.0100528734  0.0001010603  331100
      3 -0.0091347376  0.0000834434  313100
     36  0.0064197897  0.0000412137  100133
     30  0.0062869497  0.0000395257  110312
     27 -0.0061039755  0.0000372585  111230
     22  0.0059571411  0.0000354875  120311
     41 -0.0057956813  0.0000335899  013130
     21 -0.0048352636  0.0000233798  121103
     26 -0.0041811348  0.0000174819  112103
     45 -0.0032417018  0.0000105086  011312
     44  0.0024918042  0.0000062091  011321
     14 -0.0024470819  0.0000059882  132011
      9 -0.0024362991  0.0000059356  303011
     19 -0.0024303434  0.0000059066  123011
     11 -0.0023972333  0.0000057467  301103
     15 -0.0017719703  0.0000031399  131021
     16  0.0016519756  0.0000027290  131012
     23  0.0015910203  0.0000025313  113021
     31 -0.0015527635  0.0000024111  103130
     42  0.0014874784  0.0000022126  013103
     39  0.0013499459  0.0000018224  031103
     38 -0.0011972838  0.0000014335  031130

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