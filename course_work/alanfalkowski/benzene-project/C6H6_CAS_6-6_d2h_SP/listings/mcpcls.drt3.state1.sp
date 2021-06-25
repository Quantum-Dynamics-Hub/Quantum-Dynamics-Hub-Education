

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
 The CSFs for the state No  1 of the symmetry  b2  will be printed
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
     12 -0.7038076622  0.4953452254  300311
      7  0.6373830282  0.4062571247  310130
     10  0.2236129299  0.0500027424  301130
      8 -0.1069220514  0.0114323251  310103
      2  0.1049574763  0.0110160718  330011
     30  0.0664478118  0.0044153117  110312
     17 -0.0652646750  0.0042594778  130130
     36 -0.0641349651  0.0041132937  100133
      4  0.0512165141  0.0026231313  312011
     11 -0.0443898364  0.0019704576  301103
     22 -0.0421096388  0.0017732217  120311
     35  0.0386154907  0.0014911561  101312
     34 -0.0360169458  0.0012972204  101321
     46 -0.0298380679  0.0008903103  010133
     40  0.0281162705  0.0007905247  030311
     33 -0.0246175423  0.0006060234  102311
     29 -0.0199466505  0.0003978689  110321
      9  0.0153708657  0.0002362635  303011
     43  0.0152933087  0.0002338853  012311
     27  0.0144106315  0.0002076663  111230
     25 -0.0124307584  0.0001545238  112130
     47  0.0097607669  0.0000952726  003311
     48 -0.0080130567  0.0000642091  001133
     20 -0.0071644545  0.0000513294  121130
      3  0.0059701555  0.0000356428  313100
      6 -0.0059675924  0.0000356122  311012
     15 -0.0058050025  0.0000336981  131021
      1  0.0050312682  0.0000253137  331100
     26 -0.0045046473  0.0000202918  112103
     41  0.0041713317  0.0000174000  013130
      5  0.0035934206  0.0000129127  311021
     32 -0.0035532416  0.0000126255  103103
     16  0.0034835506  0.0000121351  131012
     28 -0.0034804417  0.0000121135  111203
     14 -0.0033763919  0.0000114000  132011
     31  0.0030970343  0.0000095916  103130
     18  0.0027974509  0.0000078257  130103
     45  0.0026574361  0.0000070620  011312
     38 -0.0023239475  0.0000054007  031130
     23 -0.0022164216  0.0000049125  113021
     21 -0.0017583280  0.0000030917  121103
     44 -0.0015383626  0.0000023666  011321
     42 -0.0013075946  0.0000017098  013103
     19 -0.0012873620  0.0000016573  123011
     37 -0.0012192720  0.0000014866  033011
     24  0.0010048185  0.0000010097  113012

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