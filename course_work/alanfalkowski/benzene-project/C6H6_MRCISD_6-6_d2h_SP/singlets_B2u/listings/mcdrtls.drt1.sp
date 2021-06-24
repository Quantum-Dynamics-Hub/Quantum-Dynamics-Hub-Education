
 program "mcdrt 4.1 a3"

 distinct row table specification and csf
 selection for mcscf wavefunction optimization.

 programmed by: ron shepard

 version date: 17-oct-91


 This Version of Program mcdrt is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCDRT       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 expanded keystroke file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets/WO
 
 input the spin multiplicity [  0]: spin multiplicity:    1    singlet 
 input the total number of electrons [  0]: nelt:     42
 input the number of irreps (1-8) [  0]: nsym:      8
 enter symmetry labels:(y,[n]) enter 8 labels (a4):
 enter symmetry label, default=   1
 enter symmetry label, default=   2
 enter symmetry label, default=   3
 enter symmetry label, default=   4
 enter symmetry label, default=   5
 enter symmetry label, default=   6
 enter symmetry label, default=   7
 enter symmetry label, default=   8
 input the molecular spatial symmetry (irrep 1:nsym) [  0]: spatial symmetry is irrep number:      1
 
 input the list of doubly-occupied orbitals (sym(i),rmo(i),i=1,ndot):
 number of doubly-occupied orbitals:     18
 number of inactive electrons:     36
 number of active electrons:      6
 level(*)        1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
                16  17  18
 symd(*)         1   1   1   1   1   1   2   2   2   2   2   3   3   3   3
                 4   4   4
 slabel(*)     ag  ag  ag  ag  ag  ag  b3u b3u b3u b3u b3u b2u b2u b2u b2u
               b1g b1g b1g
 doub(*)         1   2   3   4   5   6   1   2   3   4   5   1   2   3   4
                 1   2   3
 
 input the active orbitals (sym(i),rmo(i),i=1,nact):
 nact:      6
 level(*)        1   2   3   4   5   6
 syml(*)         5   5   5   6   7   8
 slabel(*)     b1u b1u b1u b2g b3g au 
 modrt(*)        1   2   3   1   1   1
 input the minimum cumulative occupation for each active level:
  b1u b1u b1u b2g b3g au 
    1   2   3   1   1   1
 input the maximum cumulative occupation for each active level:
  b1u b1u b1u b2g b3g au 
    1   2   3   1   1   1
 slabel(*)     b1u b1u b1u b2g b3g au 
 modrt(*)        1   2   3   1   1   1
 occmin(*)       0   0   0   0   0   6
 occmax(*)       6   6   6   6   6   6
 input the minimum b value for each active level:
  b1u b1u b1u b2g b3g au 
    1   2   3   1   1   1
 input the maximum b value for each active level:
  b1u b1u b1u b2g b3g au 
    1   2   3   1   1   1
 slabel(*)     b1u b1u b1u b2g b3g au 
 modrt(*)        1   2   3   1   1   1
 bmin(*)         0   0   0   0   0   0
 bmax(*)         6   6   6   6   6   6
 input the step masks for each active level:
 modrt:smask=
   1:1111   2:1111   3:1111   1:1111   1:1111   1:1111
 input the number of vertices to be deleted [  0]: number of vertices to be removed (a priori):      0
 number of rows in the drt:     30
 are any arcs to be manually removed?(y,[n])
 nwalk=      55
 input the range of drt levels to print (l1,l2):
 levprt(*)       0   6

 level  0 through level  6 of the drt:

 row lev a b syml lab rmo  l0  l1  l2  l3 isym xbar   y0    y1    y2    xp     z

  30   6 3 0   8 au    1    0   0   0   0   1     1     0     0     0    55     0
                                            2     0     0     0     0    40     0
                                            3     0     0     0     0    40     0
                                            4     0     0     0     0    40     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  27   5 3 0   7 b3g   1   30   0   0   0   1     1     0     0     0    19     0
                                            2     0     0     0     0    11     0
                                            3     0     0     0     0    11     0
                                            4     0     0     0     0     9     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  28   5 2 1   7 b3g   1    0   0  30   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0    22     0
                                            6     0     0     0     0    18     0
                                            7     0     0     0     0    18     0
                                            8     1     1     1     0    17    19

  29   5 2 0   7 b3g   1    0   0   0  30   1     1     1     1     1    19    36
                                            2     0     0     0     0    11     0
                                            3     0     0     0     0    11     0
                                            4     0     0     0     0     9     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  21   4 3 0   6 b2g   1   27   0   0   0   1     1     0     0     0     7     0
                                            2     0     0     0     0     3     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  22   4 2 1   6 b2g   1   28   0  27   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0    11     0
                                            6     0     0     0     0     9     0
                                            7     1     1     1     0     0     7
                                            8     1     0     0     0     0    19

  23   4 2 0   6 b2g   1   29  28   0  27   1     2     1     1     1    12     7
                                            2     1     1     0     0     8    19
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  24   4 1 2   6 b2g   1    0   0  28   0   1     0     0     0     0     6     0
                                            2     1     1     1     0     9    27
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  25   4 1 1   6 b2g   1    0   0  29  28   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0    11     0
                                            6     0     0     0     0     9     0
                                            7     1     1     1     0     0    48
                                            8     1     1     1     1     0    36

  26   4 1 0   6 b2g   1    0   0   0  29   1     1     1     1     1     7    48
                                            2     0     0     0     0     3     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  11   3 3 0   5 b1u   3   21   0   0   0   1     1     0     0     0     1     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  12   3 2 1   5 b1u   3   22   0  21   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     3     0
                                            6     1     1     1     0     0     1
                                            7     1     0     0     0     0     7
                                            8     1     0     0     0     0    19

  13   3 2 0   5 b1u   3   23  22   0  21   1     3     1     1     1     6     1
                                            2     1     0     0     0     0    19
                                            3     1     1     0     0     0    19
                                            4     1     1     0     0     0     7
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  14   3 1 2   5 b1u   3   24   0  22   0   1     0     0     0     0     3     0
                                            2     1     0     0     0     0    27
                                            3     1     1     1     0     0    19
                                            4     1     1     1     0     0     7
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  15   3 1 1   5 b1u   3   25  24  23  22   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     2     2     1     0     8    19
                                            6     2     2     2     0     0    13
                                            7     2     1     1     1     0     7
                                            8     2     1     1     1     0    19

  16   3 1 0   5 b1u   3   26  25   0  23   1     3     2     2     2     6    13
                                            2     1     1     1     1     0    27
                                            3     1     1     0     0     0    36
                                            4     1     1     0     0     0    48
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  17   3 0 3   5 b1u   3    0   0  24   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     1     1     1     0     1    35
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  18   3 0 2   5 b1u   3    0   0  25  24   1     0     0     0     0     3     0
                                            2     1     1     1     1     0    36
                                            3     1     1     1     0     0    36
                                            4     1     1     1     0     0    48
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  19   3 0 1   5 b1u   3    0   0  26  25   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     3     0
                                            6     1     1     1     0     0    54
                                            7     1     1     1     1     0    48
                                            8     1     1     1     1     0    36

  20   3 0 0   5 b1u   3    0   0   0  26   1     1     1     1     1     1    54
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

   5   2 2 0   5 b1u   2   13  12   0  11   1     4     1     1     1     1     0
                                            2     2     1     0     0     0     1
                                            3     2     1     0     0     0     7
                                            4     2     1     0     0     0    19
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   6   2 1 1   5 b1u   2   15  14  13  12   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     5     3     3     0     2     2
                                            6     5     3     2     1     0     1
                                            7     5     3     2     1     0     7
                                            8     5     3     2     1     0    19

   7   2 1 0   5 b1u   2   16  15   0  13   1     8     5     3     3     3     4
                                            2     4     3     1     1     0    19
                                            3     4     3     1     1     0    19
                                            4     4     3     1     1     0     7
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   8   2 0 2   5 b1u   2   18  17  15  14   1     3     3     2     0     1    24
                                            2     4     3     3     1     0    27
                                            3     4     3     3     1     0    19
                                            4     4     3     3     1     0     7
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   9   2 0 1   5 b1u   2   19  18  16  15   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     5     5     5     2     2    25
                                            6     5     4     3     2     0    13
                                            7     5     4     3     2     0     7
                                            8     5     4     3     2     0    19

  10   2 0 0   5 b1u   2   20  19   0  16   1     4     3     3     3     1    18
                                            2     2     2     1     1     0    27
                                            3     2     2     1     1     0    36
                                            4     2     2     1     1     0    48
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

   2   1 1 0   5 b1u   1    7   6   0   5   1    17     9     4     4     1     0
                                            2    11     7     2     2     0     1
                                            3    11     7     2     2     0     7
                                            4    11     7     2     2     0    19
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   3   1 0 1   5 b1u   1    9   8   7   6   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5    21    16    13     5     1     3
                                            6    18    13     9     5     0     1
                                            7    18    13     9     5     0     7
                                            8    18    13     9     5     0    19

   4   1 0 0   5 b1u   1   10   9   0   7   1    17    13     8     8     1     6
                                            2    11     9     4     4     0    19
                                            3    11     9     4     4     0    19
                                            4    11     9     4     4     0     7
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

   1   0 0 0   0       0    4   3   0   2   1    55    38    17    17     1     0
                                            2    40    29    11    11     0     1
                                            3    40    29    11    11     0     7
                                            4    40    29    11    11     0    19
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

 initial csf selection step:
 total number of walks in the drt, nwalk=      55
 keep all of these walks?(y,[n]) individual walks will be generated from the drt.
 apply orbital-group occupation restrictions?(y,[n]) apply reference occupation restrictions?(y,[n]) manually select individual walks?(y,[n])
 step-vector based csf selection complete.
       55 csfs selected from      55 total walks.

 beginning step-vector based csf selection.
 enter [step_vector/disposition] pairs:

 enter the active orbital step vector, (-1/ to end):

 step-vector based csf selection complete.
       55 csfs selected from      55 total walks.

 beginning numerical walk selection:
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input walk number (0 to end) [  0]:
 final csf selection complete.
       55 csfs selected from      55 total walks.
  drt construction and csf selection complete.
 
 input a title card, default=mdrt2_title
  title                                                                         
  
 input a drt file name, default=mcdrtfl
 drt and indexing arrays written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets/WO
 
 write the drt file?([y],n) include step(*) vectors?([y],n) drt file is being written...


   List of selected configurations (step vectors)


   CSF#     1    3 3 3 0 0 0
   CSF#     2    3 3 0 3 0 0
   CSF#     3    3 3 0 0 3 0
   CSF#     4    3 3 0 0 0 3
   CSF#     5    3 1 2 3 0 0
   CSF#     6    3 1 2 0 3 0
   CSF#     7    3 1 2 0 0 3
   CSF#     8    3 1 0 2 1 2
   CSF#     9    3 1 0 1 2 2
   CSF#    10    3 0 3 3 0 0
   CSF#    11    3 0 3 0 3 0
   CSF#    12    3 0 3 0 0 3
   CSF#    13    3 0 1 2 1 2
   CSF#    14    3 0 1 1 2 2
   CSF#    15    3 0 0 3 3 0
   CSF#    16    3 0 0 3 0 3
   CSF#    17    3 0 0 0 3 3
   CSF#    18    1 3 2 3 0 0
   CSF#    19    1 3 2 0 3 0
   CSF#    20    1 3 2 0 0 3
   CSF#    21    1 3 0 2 1 2
   CSF#    22    1 3 0 1 2 2
   CSF#    23    1 2 3 3 0 0
   CSF#    24    1 2 3 0 3 0
   CSF#    25    1 2 3 0 0 3
   CSF#    26    1 2 1 2 1 2
   CSF#    27    1 2 1 1 2 2
   CSF#    28    1 2 0 3 3 0
   CSF#    29    1 2 0 3 0 3
   CSF#    30    1 2 0 0 3 3
   CSF#    31    1 1 2 2 1 2
   CSF#    32    1 1 2 1 2 2
   CSF#    33    1 1 1 2 2 2
   CSF#    34    1 0 3 2 1 2
   CSF#    35    1 0 3 1 2 2
   CSF#    36    1 0 2 3 3 0
   CSF#    37    1 0 2 3 0 3
   CSF#    38    1 0 2 0 3 3
   CSF#    39    0 3 3 3 0 0
   CSF#    40    0 3 3 0 3 0
   CSF#    41    0 3 3 0 0 3
   CSF#    42    0 3 1 2 1 2
   CSF#    43    0 3 1 1 2 2
   CSF#    44    0 3 0 3 3 0
   CSF#    45    0 3 0 3 0 3
   CSF#    46    0 3 0 0 3 3
   CSF#    47    0 1 3 2 1 2
   CSF#    48    0 1 3 1 2 2
   CSF#    49    0 1 2 3 3 0
   CSF#    50    0 1 2 3 0 3
   CSF#    51    0 1 2 0 3 3
   CSF#    52    0 0 3 3 3 0
   CSF#    53    0 0 3 3 0 3
   CSF#    54    0 0 3 0 3 3
   CSF#    55    0 0 0 3 3 3
