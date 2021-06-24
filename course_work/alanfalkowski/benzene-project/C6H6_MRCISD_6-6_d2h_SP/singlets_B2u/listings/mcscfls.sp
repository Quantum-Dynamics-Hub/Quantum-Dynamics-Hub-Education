

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
      2621440000 of real*8 words (20000.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-4,
   tol(2)=1.e-4,
   tol(1)=1.e-8,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2,
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  5,1,40,5,2,40,5,3,40,6,1,40,7,1,40,8,1,40
   NAVST(1) = 1,
   WAVST(1,1)=1 ,
   NAVST(2) = 2,
   WAVST(2,1)=1 ,
   WAVST(2,2)=1 ,
   NAVST(3) = 2,
   WAVST(3,1)=1 ,
   WAVST(3,2)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h
 _3/

 Integral file header information:
 Hermit Integral Program : SIFS version  compute-0-55      14:52:49.849 24-Jun-21

 Core type energy values:
 energy( 1)=  2.015958396276E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =  201.595839628


   ******  Basis set information:  ******

 Number of irreps:                  8
 Total number of basis functions: 192

 irrep no.              1    2    3    4    5    6    7    8
 irrep label           Ag   B3u  B2u  B1g  B1u  B2g  B3g  Au 
 no. of bas.fcions.    39   39   30   30   16   16   11   11


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-04. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-04. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  3 DRTs.

 DRT  first state   no.of aver.states   weights
  1   ground state          1             0.200
  2   ground state          2             0.200 0.200
  3   ground state          2             0.200 0.200

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        2
    2        3
    3        3

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       5       1(139)    40
       5       2(140)    40
       5       3(141)    40
       6       1(155)    40
       7       1(171)    40
       8       1(182)    40

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0  0  0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    ag 
 Total number of electrons:   42
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        55

 Informations for the DRT no.  2

 DRT file header:
  title                                                                          
 Molecular symmetry group:    b2u
 Total number of electrons:   42
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        40

 Informations for the DRT no.  3

 DRT file header:
  title                                                                          
 Molecular symmetry group:    b3u
 Total number of electrons:   42
 Spin multiplicity:            3
 Number of active orbitals:    6
 Number of active electrons:   6
 Total number of CSFs:        48
 

 faar:   0 active-active rotations allowed out of:   3 possible.


 Number of active-double rotations:         0
 Number of active-active rotations:         0
 Number of double-virtual rotations:      553
 Number of active-virtual rotations:       74
 lenbfsdef=                131071  lenbfs=                  4244
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #          75
 class  2 (pq|ri):         #           0
 class  3 (pq|ia):         #       11405
 class  4 (pi|qa):         #       19492
 class  5 (pq|ra):         #        1522
 class  6 (pq|ij)/(pi|qj): #        2670
 class  7 (pq|ab):         #       44774
 class  8 (pa|qb):         #       88036
 class  9 p(bp,ai)         #       40922
 class 10p(ai,jp):        #           0
 class 11p(ai,bj):        #      161507

 Size of orbital-Hessian matrix B:                   205633
 Size of the orbital-state Hessian matrix C:         144837
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:         350470


 Source of the initial MO coeficients:

 Input MO coefficient file: /home3/molecula1/maplima/agf18/teste_COLUMBUS/benzene_d2h_3/
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     8, naopsy(1) =    39, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=2621346358

 inoutp: segmentation information:
 in-core transformation space,   avcinc =2621176303
 address segment size,           sizesg =2621024817
 number of in-core blocks,       nincbk =       106
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:    1795435 transformed 1/r12    array elements were written in     329 records.

 !timer: 2-e transformation              cpu_time=     2.752 walltime=     4.213

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

 trial vectors are generated internally.

 trial vector  1 is unit matrix column    15
 ciiter=   9 noldhv= 16 noldv= 16

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649519906     -432.3607916182        0.0000003143        0.0000010000
    2      -230.4272648343     -432.0231044619        0.0085709101        0.0100000000

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     6
 ciiter=  11 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534777357     -432.1493173633        0.0000002403        0.0000010000
    2      -230.4831753639     -432.0790149915        0.0000006480        0.0000010000
    3      -230.4170940945     -432.0129337222        0.0032358447        0.0100000000

 trial vectors are generated internally.

 trial vector  1 is unit matrix column    12
 ciiter=  11 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6036188788     -432.1994585064        0.0000001992        0.0000010000
    2      -230.5761126984     -432.1719523261        0.0000005122        0.0000010000
    3      -230.4791175712     -432.0749571989        0.0031485946        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.559036070285231E-008
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 4.0402E-07 rpnorm= 0.0000E+00 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.522861  -22.505885   -2.317430   -1.670076   -1.432812   -1.011848

 qvv(*) eigenvalues. symmetry block  1
     0.067001    0.112598    0.288104    0.301011    0.326554    0.377383    0.619197    0.635579    0.694843    0.870084
     0.996525    1.144862    1.178280    1.209563    1.371298    1.536205    1.759818    1.908638    1.981923    3.052423
     3.193766    3.262403    3.269203    3.541445    3.638982    3.765176    3.907539    4.381477    4.836105    5.503383
     5.525127    6.176692    7.007109

 fdd(*) eigenvalues. symmetry block  2
   -22.520750  -22.505884   -2.046481   -1.308942   -1.193050

 qvv(*) eigenvalues. symmetry block  2
     0.082453    0.150504    0.293284    0.332512    0.364538    0.595941    0.726096    0.739838    0.806926    0.845580
     1.121834    1.314156    1.341956    1.367003    1.536958    1.679076    1.966080    2.089646    2.167642    2.667882
     3.029694    3.156437    3.241435    3.868885    3.956894    4.149360    4.481802    4.719341    5.511980    5.865500
     5.870022    6.173419    7.637916    9.484363

 fdd(*) eigenvalues. symmetry block  3
   -22.522771   -2.050999   -1.258015   -1.199136

 qvv(*) eigenvalues. symmetry block  3
     0.082419    0.286451    0.292792    0.332268    0.595535    0.737980    0.762842    0.845670    1.122797    1.313348
     1.365645    1.509211    1.678901    2.169972    2.666674    2.977373    3.026598    3.238114    3.539134    3.865394
     4.146064    4.484900    4.521074    5.507870    5.866750    6.174076

 fdd(*) eigenvalues. symmetry block  4
   -22.520585   -1.669892   -1.012000

 qvv(*) eigenvalues. symmetry block  4
     0.112570    0.300875    0.377164    0.530075    0.619347    0.694701    0.869046    1.144068    1.173253    1.370581
     1.535018    1.672647    1.758152    1.908728    2.357099    3.053587    3.262277    3.539554    3.848896    3.909515
     4.377549    4.685627    4.835603    5.501795    6.176298    6.679985    7.003760
 i,qaaresolved                     1 -0.973321501659255     
 i,qaaresolved                     2  0.167311248936681     
 i,qaaresolved                     3  0.231741625532320     

 qvv(*) eigenvalues. symmetry block  5
     0.350647    0.660786    0.990369    1.188563    1.344783    1.654868    1.827603    2.324503    3.073205    3.372898
     3.756244    4.150294    4.691469
 i,qaaresolved                     1 -0.618105880023085     

 qvv(*) eigenvalues. symmetry block  6
     0.280623    0.379022    0.643792    0.844812    1.122058    1.498212    1.817466    1.893259    1.930705    2.934526
     3.401156    3.937476    4.204541    4.360737    5.274183
 i,qaaresolved                     1 -0.578801574869280     

 qvv(*) eigenvalues. symmetry block  7
     0.279337    0.843893    1.102939    1.496471    1.816586    1.928548    2.808957    3.397071    4.201253    4.363728
 i,qaaresolved                     1  0.170729933031147     

 qvv(*) eigenvalues. symmetry block  8
     0.354057    0.990374    1.348026    1.826196    2.138480    2.323007    3.376890    3.752421    4.687279    5.172127

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     2.855 walltime=     4.820

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=   -230.5962673335 demc= 2.3060E+02 wnorm= 2.8472E-07 knorm= 2.1783E-08 apxde= 2.7534E-15    *not conv.*     

               starting mcscf iteration...   2
 !timer:                                 cpu_time=     2.863 walltime=     4.854

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     8, naopsy(1) =    39, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=2621346358

 inoutp: segmentation information:
 in-core transformation space,   avcinc =2621176303
 address segment size,           sizesg =2621024817
 number of in-core blocks,       nincbk =       106
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:    1795435 transformed 1/r12    array elements were written in     329 records.

 !timer: 2-e transformation              cpu_time=     0.864 walltime=     0.871

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649519907     -432.3607916183        0.0000005491        0.0000010000
    2      -230.4272733289     -432.0231129565        0.0084327170        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   4 noldhv=  8 noldv=  8

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534777361     -432.1493173638        0.0000004188        0.0000010000
    2      -230.4831753627     -432.0790149903        0.0000002951        0.0000010000
    3      -230.4170958670     -432.0129354946        0.0030581671        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  7 noldv=  7

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6036188789     -432.1994585065        0.0000003612        0.0000010000
    2      -230.5761126991     -432.1719523267        0.0000009540        0.0000010000
    3      -230.4791184417     -432.0749580693        0.0029505552        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.585465755756817E-008
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 2.7534E-07 rpnorm= 0.0000E+00 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.522861  -22.505885   -2.317430   -1.670076   -1.432812   -1.011848

 qvv(*) eigenvalues. symmetry block  1
     0.067001    0.112598    0.288104    0.301011    0.326554    0.377383    0.619197    0.635579    0.694843    0.870084
     0.996525    1.144862    1.178280    1.209563    1.371298    1.536205    1.759819    1.908638    1.981923    3.052423
     3.193766    3.262403    3.269203    3.541445    3.638982    3.765176    3.907539    4.381477    4.836105    5.503383
     5.525127    6.176692    7.007109

 fdd(*) eigenvalues. symmetry block  2
   -22.520750  -22.505884   -2.046481   -1.308942   -1.193050

 qvv(*) eigenvalues. symmetry block  2
     0.082453    0.150504    0.293284    0.332512    0.364538    0.595941    0.726096    0.739838    0.806926    0.845580
     1.121834    1.314156    1.341956    1.367003    1.536958    1.679076    1.966080    2.089646    2.167642    2.667882
     3.029694    3.156437    3.241435    3.868885    3.956894    4.149360    4.481802    4.719341    5.511980    5.865500
     5.870022    6.173419    7.637916    9.484363

 fdd(*) eigenvalues. symmetry block  3
   -22.522771   -2.050999   -1.258015   -1.199136

 qvv(*) eigenvalues. symmetry block  3
     0.082419    0.286451    0.292792    0.332268    0.595535    0.737980    0.762842    0.845670    1.122797    1.313348
     1.365645    1.509211    1.678901    2.169972    2.666674    2.977373    3.026598    3.238114    3.539134    3.865394
     4.146064    4.484900    4.521074    5.507870    5.866750    6.174076

 fdd(*) eigenvalues. symmetry block  4
   -22.520585   -1.669892   -1.012000

 qvv(*) eigenvalues. symmetry block  4
     0.112570    0.300875    0.377164    0.530075    0.619347    0.694701    0.869046    1.144068    1.173253    1.370581
     1.535018    1.672647    1.758152    1.908728    2.357099    3.053587    3.262277    3.539554    3.848896    3.909515
     4.377549    4.685627    4.835603    5.501795    6.176298    6.679985    7.003760
 i,qaaresolved                     1 -0.973321525286263     
 i,qaaresolved                     2  0.167311263926306     
 i,qaaresolved                     3  0.231741634384711     

 qvv(*) eigenvalues. symmetry block  5
     0.350647    0.660786    0.990369    1.188563    1.344783    1.654868    1.827603    2.324503    3.073205    3.372898
     3.756244    4.150294    4.691469
 i,qaaresolved                     1 -0.618105885180668     

 qvv(*) eigenvalues. symmetry block  6
     0.280623    0.379022    0.643792    0.844812    1.122058    1.498212    1.817466    1.893259    1.930705    2.934526
     3.401156    3.937476    4.204541    4.360737    5.274183
 i,qaaresolved                     1 -0.578801588791900     

 qvv(*) eigenvalues. symmetry block  7
     0.279337    0.843893    1.102939    1.496471    1.816586    1.928548    2.808957    3.397071    4.201253    4.363728
 i,qaaresolved                     1  0.170729912516374     

 qvv(*) eigenvalues. symmetry block  8
     0.354057    0.990374    1.348026    1.826196    2.138480    2.323007    3.376890    3.752421    4.687279    5.172127

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     1.011 walltime=     1.059

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    2 emc=   -230.5962673335 demc=-2.2737E-13 wnorm= 4.4684E-07 knorm= 1.2857E-08 apxde= 2.7796E-15    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.200 total energy=     -230.764951991, rel. (eV)=   0.000000
   DRT #2 state # 1 wt 0.200 total energy=     -230.553477736, rel. (eV)=   5.754510
   DRT #2 state # 2 wt 0.200 total energy=     -230.483175363, rel. (eV)=   7.667535
   DRT #3 state # 1 wt 0.200 total energy=     -230.603618879, rel. (eV)=   4.390099
   DRT #3 state # 2 wt 0.200 total energy=     -230.576112699, rel. (eV)=   5.138581
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  Ag 
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36        MO   37        MO   38        MO   39
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  2    -  B3u
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36        MO   37        MO   38        MO   39
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  3    -  B2u
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  4    -  B1g
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  5    -  B1u
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     1.95104848     0.37766429     0.16635409     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  6    -  B2g
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     1.66016558     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  7    -  B3g
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     1.48891534     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11
  occ(*)=     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  8    -  Au 
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.35585222     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11
  occ(*)=     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        723 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        Ag  partial gross atomic populations
   ao class       1Ag        2Ag        3Ag        4Ag        5Ag        6Ag 
    C1_ s       0.024115   1.965539   0.548289   0.570135  -0.031396  -0.029880
    C1_ p       0.003468  -0.000704   0.067236   0.017042   0.323798   0.397214
    C1_ d       0.000470   0.000179   0.004670  -0.062235   0.037144  -0.079263
    C2_ s       1.965552   0.024164   1.148780   0.279751  -0.058322  -0.013963
    C2_ p       0.000801   0.007425   0.143610   0.540706   0.669714   0.894395
    C2_ d      -0.000546   0.001050   0.009672  -0.018646   0.075114  -0.021358
    H1_ s       0.001179   0.000814   0.024209   0.432803   0.323613   0.606223
    H1_ p       0.000351   0.000587  -0.000387   0.028498   0.005348  -0.030677
    H2_ s       0.003058   0.000734   0.053826   0.221114   0.644759   0.292900
    H2_ p       0.001552   0.000211   0.000096  -0.009167   0.010228  -0.015591
 
   ao class       7Ag        8Ag        9Ag       10Ag       11Ag       12Ag 
 
   ao class      13Ag       14Ag       15Ag       16Ag       17Ag       18Ag 
 
   ao class      19Ag       20Ag       21Ag       22Ag       23Ag       24Ag 
 
   ao class      25Ag       26Ag       27Ag       28Ag       29Ag       30Ag 
 
   ao class      31Ag       32Ag       33Ag       34Ag       35Ag       36Ag 
 
   ao class      37Ag       38Ag       39Ag 

                        B3u partial gross atomic populations
   ao class       1B3u       2B3u       3B3u       4B3u       5B3u       6B3u
    C1_ s       0.022924   1.952003   0.907970   0.105368  -0.119996   0.000000
    C1_ p      -0.000455   0.001641  -0.013552   0.204694   0.739348   0.000000
    C1_ d      -0.000102   0.000126   0.032980  -0.188062  -0.012470   0.000000
    C2_ s       1.958708   0.025923   0.464657   0.158806  -0.020885   0.000000
    C2_ p       0.013886   0.014358   0.277090   0.477455   0.482196   0.000000
    C2_ d       0.003923   0.004386  -0.053223  -0.386015  -0.014182   0.000000
    H1_ s       0.000617   0.000533   0.233173   0.517384   0.663779   0.000000
    H1_ p       0.000315   0.000354   0.029126  -0.000702  -0.005710   0.000000
    H2_ s       0.000067   0.000410   0.122706   1.114484   0.284368   0.000000
    H2_ p       0.000116   0.000266  -0.000926  -0.003412   0.003553   0.000000
 
   ao class       7B3u       8B3u       9B3u      10B3u      11B3u      12B3u
 
   ao class      13B3u      14B3u      15B3u      16B3u      17B3u      18B3u
 
   ao class      19B3u      20B3u      21B3u      22B3u      23B3u      24B3u
 
   ao class      25B3u      26B3u      27B3u      28B3u      29B3u      30B3u
 
   ao class      31B3u      32B3u      33B3u      34B3u      35B3u      36B3u
 
   ao class      37B3u      38B3u      39B3u

                        B2u partial gross atomic populations
   ao class       1B2u       2B2u       3B2u       4B2u       5B2u       6B2u
    C1_ p       0.016768   0.170465   0.653078   0.101007   0.000000   0.000000
    C1_ d       0.005055  -0.050968   0.015134  -0.007312   0.000000   0.000000
    C2_ s       1.965412   1.405352  -0.000037  -0.139995   0.000000   0.000000
    C2_ p       0.007796   0.067194   1.342944   1.124341   0.000000   0.000000
    C2_ d       0.001844   0.024820   0.030207  -0.018103   0.000000   0.000000
    H1_ p       0.000103  -0.011426  -0.014038   0.003649   0.000000   0.000000
    H2_ s       0.001851   0.355724   0.000546   0.942209   0.000000   0.000000
    H2_ p       0.001171   0.038838  -0.027833  -0.005799   0.000000   0.000000
 
   ao class       7B2u       8B2u       9B2u      10B2u      11B2u      12B2u
 
   ao class      13B2u      14B2u      15B2u      16B2u      17B2u      18B2u
 
   ao class      19B2u      20B2u      21B2u      22B2u      23B2u      24B2u
 
   ao class      25B2u      26B2u      27B2u      28B2u      29B2u      30B2u

                        B1g partial gross atomic populations
   ao class       1B1g       2B1g       3B1g       4B1g       5B1g       6B1g
    C1_ p       0.003407   0.331464   0.471993   0.000000   0.000000   0.000000
    C1_ d       0.000348   0.006990   0.012509   0.000000   0.000000   0.000000
    C2_ s       1.989749   0.864965  -0.043125   0.000000   0.000000   0.000000
    C2_ p       0.004512   0.200204   0.834988   0.000000   0.000000   0.000000
    C2_ d       0.001538  -0.090846  -0.111223   0.000000   0.000000   0.000000
    H1_ p       0.000044  -0.015568  -0.000572   0.000000   0.000000   0.000000
    H2_ s       0.000202   0.667375   0.881494   0.000000   0.000000   0.000000
    H2_ p       0.000200   0.035415  -0.046065   0.000000   0.000000   0.000000
 
   ao class       7B1g       8B1g       9B1g      10B1g      11B1g      12B1g
 
   ao class      13B1g      14B1g      15B1g      16B1g      17B1g      18B1g
 
   ao class      19B1g      20B1g      21B1g      22B1g      23B1g      24B1g
 
   ao class      25B1g      26B1g      27B1g      28B1g      29B1g      30B1g

                        B1u partial gross atomic populations
   ao class       1B1u       2B1u       3B1u       4B1u       5B1u       6B1u
    C1_ p       0.600857   0.200041  -0.007658   0.000000   0.000000   0.000000
    C1_ d       0.009560  -0.002809  -0.000126   0.000000   0.000000   0.000000
    C2_ p       1.309067   0.152987   0.175216   0.000000   0.000000   0.000000
    C2_ d       0.020124   0.018168   0.001101   0.000000   0.000000   0.000000
    H1_ p       0.003409   0.006192  -0.000316   0.000000   0.000000   0.000000
    H2_ p       0.008032   0.003085  -0.001863   0.000000   0.000000   0.000000
 
   ao class       7B1u       8B1u       9B1u      10B1u      11B1u      12B1u
 
   ao class      13B1u      14B1u      15B1u      16B1u

                        B2g partial gross atomic populations
   ao class       1B2g       2B2g       3B2g       4B2g       5B2g       6B2g
    C1_ p       1.075348   0.000000   0.000000   0.000000   0.000000   0.000000
    C1_ d       0.001658   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ p       0.530870   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ d       0.030542   0.000000   0.000000   0.000000   0.000000   0.000000
    H1_ p       0.014374   0.000000   0.000000   0.000000   0.000000   0.000000
    H2_ p       0.007374   0.000000   0.000000   0.000000   0.000000   0.000000
 
   ao class       7B2g       8B2g       9B2g      10B2g      11B2g      12B2g
 
   ao class      13B2g      14B2g      15B2g      16B2g

                        B3g partial gross atomic populations
   ao class       1B3g       2B3g       3B3g       4B3g       5B3g       6B3g
    C1_ d       0.016013   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ p       1.443030   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ d       0.010890   0.000000   0.000000   0.000000   0.000000   0.000000
    H2_ p       0.018982   0.000000   0.000000   0.000000   0.000000   0.000000
 
   ao class       7B3g       8B3g       9B3g      10B3g      11B3g

                        Au  partial gross atomic populations
   ao class       1Au        2Au        3Au        4Au        5Au        6Au 
    C1_ d       0.014179   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ p       0.329477   0.000000   0.000000   0.000000   0.000000   0.000000
    C2_ d       0.001931   0.000000   0.000000   0.000000   0.000000   0.000000
    H2_ p       0.010265   0.000000   0.000000   0.000000   0.000000   0.000000
 
   ao class       7Au        8Au        9Au       10Au       11Au 


                        gross atomic populations
     ao           C1_        C2_        H1_        H2_
      s         5.915071  11.975492   2.804329   5.587828
      p         5.356499  11.044263   0.012953   0.028728
      d        -0.246332  -0.478831   0.000000   0.000000
    total      11.025238  22.540924   2.817282   5.616556
 

 Total number of electrons:   42.00000000

 !timer: mcscf                           cpu_time=     3.907 walltime=     6.020
