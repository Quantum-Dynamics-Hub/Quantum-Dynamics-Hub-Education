

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
    Hermit Integral Program : SIFS version  compute-0-55      14:52:49.849 24-Jun-21
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
     12 -0.7038076548  0.4953452150  300311
      7  0.6373696592  0.4062400825  310130
     10  0.2236512574  0.0500198849  301130
      8 -0.1069193877  0.0114317555  310103
      2  0.1049531889  0.0110151719  330011
     30  0.0664453854  0.0044149892  110312
     17 -0.0652629725  0.0042592556  130130
     36 -0.0641348555  0.0041132797  100133
      4  0.0512243418  0.0026239332  312011
     11 -0.0443963671  0.0019710374  301103
     22 -0.0421084083  0.0017731180  120311
     35  0.0386193556  0.0014914546  101312
     34 -0.0360179919  0.0012972957  101321
     46 -0.0298375841  0.0008902814  010133
     40  0.0281147311  0.0007904381  030311
     33 -0.0246196518  0.0006061273  102311
     29 -0.0199445885  0.0003977866  110321
      9  0.0153751630  0.0002363956  303011
     43  0.0152947520  0.0002339294  012311
     27  0.0144105575  0.0002076642  111230
     25 -0.0124359123  0.0001546519  112130
     47  0.0097620134  0.0000952969  003311
     48 -0.0080147719  0.0000642366  001133
     20 -0.0071669917  0.0000513658  121130
      3  0.0059704090  0.0000356458  313100
      6 -0.0059674759  0.0000356108  311012
     15 -0.0058049337  0.0000336973  131021
      1  0.0050307157  0.0000253081  331100
     26 -0.0045042471  0.0000202882  112103
     41  0.0041713346  0.0000174000  013130
      5  0.0035933979  0.0000129125  311021
     32 -0.0035537898  0.0000126294  103103
     16  0.0034834906  0.0000121347  131012
     28 -0.0034803609  0.0000121129  111203
     14 -0.0033761811  0.0000113986  132011
     31  0.0030957680  0.0000095838  103130
     18  0.0027978456  0.0000078279  130103
     45  0.0026574245  0.0000070619  011312
     38 -0.0023241891  0.0000054019  031130
     23 -0.0022167653  0.0000049140  113021
     21 -0.0017582551  0.0000030915  121103
     44 -0.0015382687  0.0000023663  011321
     42 -0.0013074843  0.0000017095  013103
     19 -0.0012875080  0.0000016577  123011
     37 -0.0012192341  0.0000014865  033011
     24  0.0010049557  0.0000010099  113012

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