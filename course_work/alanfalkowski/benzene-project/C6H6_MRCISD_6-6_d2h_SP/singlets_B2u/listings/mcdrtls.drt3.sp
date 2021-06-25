
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
 
 input the spin multiplicity [  0]: spin multiplicity:    3    triplet 
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
 input the molecular spatial symmetry (irrep 1:nsym) [  0]: spatial symmetry is irrep number:      2
 
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
 number of rows in the drt:     32
 are any arcs to be manually removed?(y,[n])
 nwalk=      48
 input the range of drt levels to print (l1,l2):
 levprt(*)       0   6

 level  0 through level  6 of the drt:

 row lev a b syml lab rmo  l0  l1  l2  l3 isym xbar   y0    y1    y2    xp     z

  32   6 2 2   8 au    1    0   0   0   0   1     1     0     0     0    48     0
                                            2     0     0     0     0    45     0
                                            3     0     0     0     0    48     0
                                            4     0     0     0     0    48     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  28   5 2 2   7 b3g   1   32   0   0   0   1     1     0     0     0    12     0
                                            2     0     0     0     0     9     0
                                            3     0     0     0     0    12     0
                                            4     0     0     0     0    12     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  29   5 2 1   7 b3g   1    0  32   0   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0    18     0
                                            6     0     0     0     0    22     0
                                            7     0     0     0     0    17     0
                                            8     1     1     0     0    18    12

  30   5 1 3   7 b3g   1    0   0  32   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     6     0
                                            6     0     0     0     0     2     0
                                            7     0     0     0     0    10     0
                                            8     1     1     1     0     6    30

  31   5 1 2   7 b3g   1    0   0   0  32   1     1     1     1     1    12    36
                                            2     0     0     0     0     9     0
                                            3     0     0     0     0    12     0
                                            4     0     0     0     0    12     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  19   4 2 2   6 b2g   1   28   0   0   0   1     1     0     0     0     3     0
                                            2     0     0     0     0     3     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  20   4 2 1   6 b2g   1   29  28   0   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     9     0
                                            6     0     0     0     0    11     0
                                            7     1     1     0     0     0     3
                                            8     1     0     0     0     0    12

  21   4 2 0   6 b2g   1    0  29   0   0   1     0     0     0     0     8     0
                                            2     1     1     0     0    12    12
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  22   4 1 3   6 b2g   1   30   0  28   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     3     0
                                            6     0     0     0     0     1     0
                                            7     1     1     1     0     0     3
                                            8     1     0     0     0     0    30

  23   4 1 2   6 b2g   1   31  30  29  28   1     2     1     1     1     9     3
                                            2     2     2     1     0     6    24
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  24   4 1 1   6 b2g   1    0  31   0  29   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     9     0
                                            6     0     0     0     0    11     0
                                            7     1     1     0     0     0    45
                                            8     1     1     1     1     0    30

  25   4 0 4   6 b2g   1    0   0  30   0   1     0     0     0     0     1     0
                                            2     1     1     1     0     0    36
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  26   4 0 3   6 b2g   1    0   0  31  30   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     3     0
                                            6     0     0     0     0     1     0
                                            7     1     1     1     0     0    45
                                            8     1     1     1     1     0    36

  27   4 0 2   6 b2g   1    0   0   0  31   1     1     1     1     1     3    45
                                            2     0     0     0     0     3     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

  11   3 2 1   5 b1u   3   20  19   0   0   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     1     1     0     0     3     0
                                            7     1     0     0     0     0     3
                                            8     1     0     0     0     0    12

  12   3 2 0   5 b1u   3   21  20   0   0   1     0     0     0     0     0     0
                                            2     1     0     0     0     6    12
                                            3     1     1     0     0     0    12
                                            4     1     1     0     0     0     3
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  13   3 1 2   5 b1u   3   23  22  20  19   1     3     1     1     1     0     3
                                            2     2     0     0     0     3    24
                                            3     2     2     1     0     0    12
                                            4     2     2     1     0     0     3
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  14   3 1 1   5 b1u   3   24  23  21  20   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     3     3     1     0     0    18
                                            6     2     2     0     0     8     3
                                            7     2     1     1     1     0     3
                                            8     2     1     1     1     0    12

  15   3 1 0   5 b1u   3    0  24   0  21   1     0     0     0     0     0     0
                                            2     1     1     1     1     6    18
                                            3     1     1     0     0     0    30
                                            4     1     1     0     0     0    45
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  16   3 0 3   5 b1u   3   26  25  23  22   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     3     3     2     0     0    27
                                            6     2     2     2     0     1    11
                                            7     2     1     1     1     0     3
                                            8     2     1     1     1     0    30

  17   3 0 2   5 b1u   3   27  26  24  23   1     3     2     2     2     0    12
                                            2     2     2     2     2     3    27
                                            3     2     2     1     0     0    30
                                            4     2     2     1     0     0    45
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

  18   3 0 1   5 b1u   3    0  27   0  24   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     0     0     0     0     0     0
                                            6     1     1     0     0     3    45
                                            7     1     1     1     1     0    45
                                            8     1     1     1     1     0    30
 ........................................

   5   2 2 0   5 b1u   2   12  11   0   0   1     0     0     0     0     0     0
                                            2     2     1     0     0     1     0
                                            3     2     1     0     0     0     3
                                            4     2     1     0     0     0    12
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   6   2 1 1   5 b1u   2   14  13  12  11   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     6     3     0     0     0     3
                                            6     6     4     2     1     2     1
                                            7     6     4     2     1     0     3
                                            8     6     4     2     1     0    12

   7   2 1 0   5 b1u   2   15  14   0  12   1     3     3     0     0     0    18
                                            2     4     3     1     1     3    15
                                            3     4     3     1     1     0    12
                                            4     4     3     1     1     0     3
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   8   2 0 2   5 b1u   2   17  16  14  13   1    12     9     6     3     0     3
                                            2     8     6     4     2     1    26
                                            3     8     6     4     2     0    12
                                            4     8     6     4     2     0     3
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   9   2 0 1   5 b1u   2   18  17  15  14   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5     6     6     3     3     0    18
                                            6     6     5     3     2     2     9
                                            7     6     5     3     2     0     3
                                            8     6     5     3     2     0    12

  10   2 0 0   5 b1u   2    0  18   0  15   1     0     0     0     0     0     0
                                            2     2     2     1     1     1    23
                                            3     2     2     1     1     0    30
                                            4     2     2     1     1     0    45
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

   2   1 1 0   5 b1u   1    7   6   0   5   1     9     6     0     0     0     3
                                            2    12     8     2     2     1     0
                                            3    12     8     2     2     0     3
                                            4    12     8     2     2     0    12
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0

   3   1 0 1   5 b1u   1    9   8   7   6   1     0     0     0     0     0     0
                                            2     0     0     0     0     0     0
                                            3     0     0     0     0     0     0
                                            4     0     0     0     0     0     0
                                            5    27    21     9     6     0     3
                                            6    24    18    10     6     1     2
                                            7    24    18    10     6     0     3
                                            8    24    18    10     6     0    12

   4   1 0 0   5 b1u   1   10   9   0   7   1     9     9     3     3     0    18
                                            2    12    10     4     4     1    17
                                            3    12    10     4     4     0    12
                                            4    12    10     4     4     0     3
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

   1   0 0 0   0       0    4   3   0   2   1    45    36     9     9     0     3
                                            2    48    36    12    12     1     0
                                            3    48    36    12    12     0     3
                                            4    48    36    12    12     0    12
                                            5     0     0     0     0     0     0
                                            6     0     0     0     0     0     0
                                            7     0     0     0     0     0     0
                                            8     0     0     0     0     0     0
 ........................................

 initial csf selection step:
 total number of walks in the drt, nwalk=      48
 keep all of these walks?(y,[n]) individual walks will be generated from the drt.
 apply orbital-group occupation restrictions?(y,[n]) apply reference occupation restrictions?(y,[n]) manually select individual walks?(y,[n])
 step-vector based csf selection complete.
       48 csfs selected from      48 total walks.

 beginning step-vector based csf selection.
 enter [step_vector/disposition] pairs:

 enter the active orbital step vector, (-1/ to end):

 step-vector based csf selection complete.
       48 csfs selected from      48 total walks.

 beginning numerical walk selection:
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input walk number (0 to end) [  0]:
 final csf selection complete.
       48 csfs selected from      48 total walks.
  drt construction and csf selection complete.
 
 input a title card, default=mdrt2_title
  title                                                                         
  
 input a drt file name, default=mcdrtfl
 drt and indexing arrays written to file:
 /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/MRCI_6-6/singlets/WO
 
 write the drt file?([y],n) include step(*) vectors?([y],n) drt file is being written...


   List of selected configurations (step vectors)


   CSF#     1    3 3 1 1 0 0
   CSF#     2    3 3 0 0 1 1
   CSF#     3    3 1 3 1 0 0
   CSF#     4    3 1 2 0 1 1
   CSF#     5    3 1 1 0 2 1
   CSF#     6    3 1 1 0 1 2
   CSF#     7    3 1 0 1 3 0
   CSF#     8    3 1 0 1 0 3
   CSF#     9    3 0 3 0 1 1
   CSF#    10    3 0 1 1 3 0
   CSF#    11    3 0 1 1 0 3
   CSF#    12    3 0 0 3 1 1
   CSF#    13    1 3 3 1 0 0
   CSF#    14    1 3 2 0 1 1
   CSF#    15    1 3 1 0 2 1
   CSF#    16    1 3 1 0 1 2
   CSF#    17    1 3 0 1 3 0
   CSF#    18    1 3 0 1 0 3
   CSF#    19    1 2 3 0 1 1
   CSF#    20    1 2 1 1 3 0
   CSF#    21    1 2 1 1 0 3
   CSF#    22    1 2 0 3 1 1
   CSF#    23    1 1 3 0 2 1
   CSF#    24    1 1 3 0 1 2
   CSF#    25    1 1 2 1 3 0
   CSF#    26    1 1 2 1 0 3
   CSF#    27    1 1 1 2 3 0
   CSF#    28    1 1 1 2 0 3
   CSF#    29    1 1 0 3 2 1
   CSF#    30    1 1 0 3 1 2
   CSF#    31    1 0 3 1 3 0
   CSF#    32    1 0 3 1 0 3
   CSF#    33    1 0 2 3 1 1
   CSF#    34    1 0 1 3 2 1
   CSF#    35    1 0 1 3 1 2
   CSF#    36    1 0 0 1 3 3
   CSF#    37    0 3 3 0 1 1
   CSF#    38    0 3 1 1 3 0
   CSF#    39    0 3 1 1 0 3
   CSF#    40    0 3 0 3 1 1
   CSF#    41    0 1 3 1 3 0
   CSF#    42    0 1 3 1 0 3
   CSF#    43    0 1 2 3 1 1
   CSF#    44    0 1 1 3 2 1
   CSF#    45    0 1 1 3 1 2
   CSF#    46    0 1 0 1 3 3
   CSF#    47    0 0 3 3 1 1
   CSF#    48    0 0 1 1 3 3
