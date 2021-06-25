

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
 Hermit Integral Program : SIFS version  compute-0-37      13:33:52.636 23-Jun-21

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
 trmain:    1795381 transformed 1/r12    array elements were written in     329 records.

 !timer: 2-e transformation              cpu_time=     2.675 walltime=     7.233

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463
 !timer: mosort                          cpu_time=     0.111 walltime=     0.866

 trial vectors are generated internally.

 trial vector  1 is unit matrix column    15
 ciiter=  10 noldhv= 16 noldv= 16

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7419868575     -432.3378264851        0.0000001454        0.0000010000
    2      -230.3789480091     -431.9747876367        0.0031369079        0.0100000000

 trial vectors are generated internally.

 trial vector  1 is unit matrix column    11
 ciiter=  10 noldhv= 18 noldv= 18

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5120792937     -432.1079189213        0.0000002908        0.0000010000
    2      -230.4499134843     -432.0457531119        0.0000002937        0.0000010000
    3      -230.4245388315     -432.0203784591        0.0015704341        0.0100000000

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     7
 ciiter=  10 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5407629343     -432.1366025619        0.0000001777        0.0000010000
    2      -230.5291308025     -432.1249704301        0.0000002582        0.0000010000
    3      -230.4517270070     -432.0475666346        0.0016104045        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.334053493826040E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.90291224 pnorm= 0.0000E+00 rznorm= 5.0015E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.668341  -22.636837   -2.428096   -1.758712   -1.510182   -1.094676

 qvv(*) eigenvalues. symmetry block  1
     0.060072    0.110085    0.262151    0.292606    0.305715    0.371842    0.598657    0.602368    0.663094    0.843714
     0.949587    1.084396    1.149484    1.160202    1.338887    1.507275    1.725175    1.847274    1.900521    2.999340
     3.130446    3.201758    3.209394    3.481116    3.562962    3.695276    3.842513    4.308726    4.767401    5.444814
     5.462322    6.091243    6.927346

 fdd(*) eigenvalues. symmetry block  2
   -22.666212  -22.636858   -2.145078   -1.371083   -1.273086

 qvv(*) eigenvalues. symmetry block  2
     0.078783    0.149388    0.284263    0.323307    0.357693    0.575393    0.700908    0.710839    0.788186    0.818335
     1.067461    1.287119    1.288386    1.336934    1.512203    1.629098    1.892337    2.048541    2.118338    2.613433
     2.954731    3.107655    3.182421    3.794917    3.908801    4.087682    4.417923    4.642462    5.444490    5.806049
     5.810675    6.094208    7.553972    9.411569

 fdd(*) eigenvalues. symmetry block  3
   -22.668311   -2.153282   -1.364508   -1.281978

 qvv(*) eigenvalues. symmetry block  3
     0.078574    0.278392    0.283440    0.322713    0.574664    0.698255    0.739898    0.818133    1.068563    1.285642
     1.335161    1.471477    1.628713    2.121931    2.611147    2.908706    2.949988    3.177200    3.462017    3.789574
     4.082658    4.422961    4.462486    5.437814    5.803747    6.095242

 fdd(*) eigenvalues. symmetry block  4
   -22.666120   -1.757413   -1.095468

 qvv(*) eigenvalues. symmetry block  4
     0.110066    0.292639    0.371456    0.528775    0.602800    0.662957    0.842010    1.082923    1.157448    1.338118
     1.505402    1.606984    1.722840    1.847493    2.303796    3.001325    3.202225    3.477935    3.795850    3.845710
     4.302825    4.638573    4.766834    5.442813    6.090686    6.605506    6.922356
 i,qaaresolved                     1  -1.06848932878413     
 i,qaaresolved                     2  9.709813873458277E-002
 i,qaaresolved                     3  0.199110989767591     

 qvv(*) eigenvalues. symmetry block  5
     0.353506    0.626198    0.976806    1.130902    1.309000    1.588382    1.787168    2.259062    3.003462    3.287813
     3.696399    4.094924    4.625613
 i,qaaresolved                     1 -0.699656292993638     

 qvv(*) eigenvalues. symmetry block  6
     0.265930    0.367904    0.591580    0.822824    1.097263    1.455734    1.768160    1.850942    1.872629    2.876921
     3.332024    3.879109    4.148793    4.275943    5.202394
 i,qaaresolved                     1 -0.665346769326246     

 qvv(*) eigenvalues. symmetry block  7
     0.264439    0.821148    1.055057    1.452425    1.765209    1.869278    2.723612    3.325408    4.144877    4.279561
 i,qaaresolved                     1  9.508977695413461E-002

 qvv(*) eigenvalues. symmetry block  8
     0.351665    0.975194    1.310637    1.784830    2.114601    2.255095    3.293955    3.690382    4.619891    5.084972

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     2.835 walltime=     8.270

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=   -230.5547746744 demc= 2.3055E+02 wnorm= 2.6672E-01 knorm= 4.2982E-01 apxde= 2.6021E-02    *not conv.*     

               starting mcscf iteration...   2
 !timer:                                 cpu_time=     2.843 walltime=     8.290

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

 !timer: 2-e transformation              cpu_time=     0.783 walltime=     0.813

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   7 noldhv= 10 noldv= 10

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7583224558     -432.3541620834        0.0000046819        0.0000100000
    2      -230.4274639577     -432.0233035853        0.0090750893        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   7 noldhv= 17 noldv= 17

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5507431076     -432.1465827352        0.0000044420        0.0000100000
    2      -230.4759339791     -432.0717736068        0.0000022419        0.0000100000
    3      -230.4299395468     -432.0257791744        0.0045179176        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   7 noldhv= 16 noldv= 16

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5938564008     -432.1896960284        0.0000088298        0.0000100000
    2      -230.5713795279     -432.1672191555        0.0000031188        0.0000100000
    3      -230.4749912814     -432.0708309090        0.0044391594        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  8.263357589383820E-003
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.97275512 pnorm= 0.0000E+00 rznorm= 7.1009E-06 rpnorm= 0.0000E+00 noldr=  7 nnewr=  7 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.559499  -22.539316   -2.345128   -1.692875   -1.453455   -1.033629

 qvv(*) eigenvalues. symmetry block  1
     0.065266    0.111898    0.282276    0.298593    0.321291    0.375886    0.615010    0.626359    0.686731    0.863150
     0.985075    1.129321    1.170574    1.196816    1.363263    1.528549    1.750633    1.892515    1.961272    3.038476
     3.177955    3.246896    3.253848    3.526063    3.619610    3.746812    3.890410    4.362980    4.817987    5.487708
     5.508227    6.154779    6.986939

 fdd(*) eigenvalues. symmetry block  2
   -22.557378  -22.539323   -2.071392   -1.325953   -1.214055

 qvv(*) eigenvalues. symmetry block  2
     0.081591    0.150113    0.291114    0.330414    0.362591    0.591174    0.722213    0.730236    0.802004    0.838444
     1.107664    1.307796    1.328505    1.358840    1.530513    1.665973    1.948063    2.078231    2.154370    2.653799
     3.011043    3.143295    3.226610    3.850296    3.944188    4.133127    4.464803    4.699566    5.494414    5.849528
     5.854236    6.152790    7.616496    9.465890

 fdd(*) eigenvalues. symmetry block  3
   -22.559430   -2.076964   -1.285276   -1.221254

 qvv(*) eigenvalues. symmetry block  3
     0.081199    0.284228    0.290301    0.329708    0.590413    0.727824    0.756595    0.838282    1.108735    1.306627
     1.356878    1.499101    1.665430    2.157388    2.652113    2.959880    3.007325    3.222373    3.519430    3.846009
     4.129162    4.468597    4.505191    5.489453    5.850031    6.153678

 fdd(*) eigenvalues. symmetry block  4
   -22.557239   -1.692572   -1.034159

 qvv(*) eigenvalues. symmetry block  4
     0.111800    0.298728    0.375362    0.529563    0.615343    0.686824    0.861581    1.128591    1.168803    1.362427
     1.526777    1.655713    1.748520    1.892888    2.342718    3.039999    3.246766    3.523566    3.835093    3.892907
     4.358377    4.672893    4.817361    5.485932    6.154299    6.660364    6.982948
 i,qaaresolved                     1 -0.996869535467501     
 i,qaaresolved                     2  0.132789899745535     
 i,qaaresolved                     3  0.215039090031067     

 qvv(*) eigenvalues. symmetry block  5
     0.369585    0.651976    0.989322    1.174862    1.337843    1.640357    1.818299    2.309100    3.055538    3.351519
     3.741163    4.135429    4.674360
 i,qaaresolved                     1 -0.638952726374150     

 qvv(*) eigenvalues. symmetry block  6
     0.277457    0.376118    0.631139    0.839304    1.115898    1.487801    1.805854    1.882765    1.915874    2.920158
     3.383951    3.922363    4.189964    4.338829    5.255453
 i,qaaresolved                     1 -0.599519368421591     

 qvv(*) eigenvalues. symmetry block  7
     0.276231    0.837661    1.090351    1.485251    1.804551    1.913231    2.787651    3.378890    4.185970    4.342159
 i,qaaresolved                     1  0.132181732009708     

 qvv(*) eigenvalues. symmetry block  8
     0.370168    0.987939    1.340535    1.816060    2.132170    2.306786    3.356125    3.736367    4.669341    5.149757

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.891 walltime=     1.270

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=   -230.5900470942 demc= 3.5272E-02 wnorm= 6.6107E-02 knorm= 2.3184E-01 apxde= 5.3673E-03    *not conv.*     

               starting mcscf iteration...   3
 !timer:                                 cpu_time=     3.734 walltime=     9.560

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

 !timer: 2-e transformation              cpu_time=     0.758 walltime=     0.787

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   8 noldhv= 12 noldv= 12

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7644468486     -432.3602864762        0.0000014669        0.0000023872
    2      -230.4285559828     -432.0243956104        0.0041333482        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   8 noldhv= 19 noldv= 19

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5540675926     -432.1499072202        0.0000006804        0.0000023872
    2      -230.4823725385     -432.0782121662        0.0000008098        0.0000023872
    3      -230.4203814871     -432.0162211147        0.0060124625        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   8 noldhv= 19 noldv= 19

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6033590925     -432.1991987201        0.0000004070        0.0000023872
    2      -230.5763727161     -432.1722123437        0.0000014942        0.0000023872
    3      -230.4798008488     -432.0756404765        0.0060497528        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.611314851475567E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99887966 pnorm= 0.0000E+00 rznorm= 7.8647E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.529190  -22.511744   -2.322017   -1.673708   -1.436025   -1.015279

 qvv(*) eigenvalues. symmetry block  1
     0.066755    0.112542    0.287272    0.300644    0.325788    0.377278    0.618634    0.634190    0.693622    0.869202
     0.994762    1.142314    1.177228    1.207574    1.370170    1.535214    1.758552    1.906001    1.978568    3.050280
     3.191343    3.260013    3.266853    3.539150    3.635886    3.762305    3.904842    4.378586    4.833297    5.501032
     5.522562    6.173153    7.003908

 fdd(*) eigenvalues. symmetry block  2
   -22.527078  -22.511743   -2.050522   -1.311486   -1.196297

 qvv(*) eigenvalues. symmetry block  2
     0.082438    0.150459    0.293066    0.332354    0.364300    0.595351    0.725574    0.738481    0.806254    0.844639
     1.119584    1.313345    1.339839    1.365919    1.536077    1.677153    1.963222    2.087821    2.165529    2.665785
     3.026771    3.154433    3.239235    3.866000    3.955042    4.146912    4.479130    4.716197    5.509266    5.863049
     5.867693    6.170106    7.634461    9.481432

 fdd(*) eigenvalues. symmetry block  3
   -22.529103   -2.055262   -1.262514   -1.202702

 qvv(*) eigenvalues. symmetry block  3
     0.082186    0.286185    0.292410    0.331840    0.594757    0.736396    0.761988    0.844597    1.120592    1.312386
     1.364240    1.507731    1.676781    2.168131    2.664400    2.974627    3.023485    3.235579    3.536005    3.862250
     4.143427    4.482425    4.518690    5.504936    5.864156    6.170851

 fdd(*) eigenvalues. symmetry block  4
   -22.526917   -1.673540   -1.015581

 qvv(*) eigenvalues. symmetry block  4
     0.112467    0.300670    0.376879    0.530019    0.618874    0.693620    0.867907    1.141624    1.172706    1.369389
     1.533745    1.669974    1.756692    1.906253    2.354839    3.051574    3.259820    3.537042    3.846788    3.906986
     4.374478    4.683795    4.832737    5.499389    6.172733    6.676883    7.000378
 i,qaaresolved                     1 -0.977071132039442     
 i,qaaresolved                     2  0.163171237752338     
 i,qaaresolved                     3  0.224719086468139     

 qvv(*) eigenvalues. symmetry block  5
     0.354834    0.659561    0.990725    1.186571    1.344047    1.652787    1.826442    2.322190    3.070448    3.369422
     3.754045    4.148083    4.688927
 i,qaaresolved                     1 -0.621902395015413     

 qvv(*) eigenvalues. symmetry block  6
     0.280188    0.378692    0.641643    0.844230    1.121136    1.496756    1.815775    1.891644    1.928394    2.932352
     3.398523    3.935166    4.202480    4.357192    5.271283
 i,qaaresolved                     1 -0.581803337430895     

 qvv(*) eigenvalues. symmetry block  7
     0.278976    0.842936    1.101086    1.494689    1.814816    1.926120    2.805485    3.394142    4.198905    4.360278
 i,qaaresolved                     1  0.164393271060986     

 qvv(*) eigenvalues. symmetry block  8
     0.356679    0.990021    1.347153    1.824587    2.137676    2.320496    3.373564    3.749890    4.684445    5.168543

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.865 walltime=     1.270

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=   -230.5961237577 demc= 6.0767E-03 wnorm= 1.2891E-02 knorm= 4.7323E-02 apxde= 1.3022E-04    *not conv.*     

               starting mcscf iteration...   4
 !timer:                                 cpu_time=     4.599 walltime=    10.830

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

 !timer: 2-e transformation              cpu_time=     0.758 walltime=     0.787

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv= 13 noldv= 13

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649130193     -432.3607526469        0.0000005438        0.0000010000
    2      -230.4271234594     -432.0229630870        0.0097714442        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534871924     -432.1493268200        0.0000003213        0.0000010000
    2      -230.4832573807     -432.0790970084        0.0000002523        0.0000010000
    3      -230.4176224086     -432.0134620363        0.0040814861        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6035885027     -432.1994281303        0.0000004916        0.0000010000
    2      -230.5760723294     -432.1719119570        0.0000003367        0.0000010000
    3      -230.4796797191     -432.0755193467        0.0062887426        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.228366838919699E-004
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995407 pnorm= 0.0000E+00 rznorm= 4.4367E-07 rpnorm= 0.0000E+00 noldr=  7 nnewr=  7 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.523243  -22.506285   -2.317735   -1.670336   -1.433059   -1.012092

 qvv(*) eigenvalues. symmetry block  1
     0.066961    0.112586    0.288046    0.300942    0.326496    0.377371    0.619150    0.635471    0.694723    0.870011
     0.996405    1.144663    1.178168    1.209405    1.371211    1.536133    1.759717    1.908424    1.981687    3.052236
     3.193580    3.262222    3.269038    3.541282    3.638758    3.764957    3.907317    4.381274    4.835893    5.503190
     5.524923    6.176443    7.006888

 fdd(*) eigenvalues. symmetry block  2
   -22.521132  -22.506283   -2.046759   -1.309159   -1.193278

 qvv(*) eigenvalues. symmetry block  2
     0.082462    0.150488    0.293275    0.332527    0.364496    0.595925    0.726042    0.739739    0.806867    0.845490
     1.121662    1.314098    1.341792    1.366935    1.536872    1.678934    1.965874    2.089512    2.167442    2.667726
     3.029498    3.156266    3.241292    3.868697    3.956736    4.149174    4.481579    4.719111    5.511784    5.865285
     5.869851    6.173169    7.637672    9.484154

 fdd(*) eigenvalues. symmetry block  3
   -22.523152   -2.051293   -1.258323   -1.199416

 qvv(*) eigenvalues. symmetry block  3
     0.082362    0.286408    0.292738    0.332203    0.595470    0.737825    0.762743    0.845542    1.122636    1.313252
     1.365488    1.509073    1.678706    2.169835    2.666475    2.977172    3.026364    3.237896    3.538903    3.865153
     4.145845    4.484710    4.520874    5.507640    5.866536    6.173844

 fdd(*) eigenvalues. symmetry block  4
   -22.520966   -1.670159   -1.012276

 qvv(*) eigenvalues. symmetry block  4
     0.112544    0.300855    0.377099    0.530050    0.619327    0.694622    0.868905    1.143900    1.173169    1.370478
     1.534871    1.672440    1.758005    1.908556    2.356914    3.053432    3.262072    3.539346    3.848723    3.909325
     4.377312    4.685459    4.835378    5.501594    6.176046    6.679750    7.003512
 i,qaaresolved                     1 -0.973570092735283     
 i,qaaresolved                     2  0.167743223737473     
 i,qaaresolved                     3  0.229428043074581     

 qvv(*) eigenvalues. symmetry block  5
     0.351265    0.660700    0.990547    1.188529    1.344846    1.654888    1.827591    2.324395    3.073000    3.372651
     3.756099    4.150104    4.691286
 i,qaaresolved                     1 -0.618565294650827     

 qvv(*) eigenvalues. symmetry block  6
     0.280568    0.379010    0.643661    0.844785    1.121967    1.498113    1.817336    1.893125    1.930554    2.934371
     3.400983    3.937289    4.204382    4.360478    5.273959
 i,qaaresolved                     1 -0.578766727749654     

 qvv(*) eigenvalues. symmetry block  7
     0.279318    0.843767    1.102807    1.496298    1.816448    1.928383    2.808714    3.396842    4.201030    4.363487
 i,qaaresolved                     1  0.170746695121317     

 qvv(*) eigenvalues. symmetry block  8
     0.354039    0.990242    1.347944    1.826011    2.138384    2.322836    3.376667    3.752197    4.687035    5.171869

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.869 walltime=     1.251

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=   -230.5962636849 demc= 1.3993E-04 wnorm= 9.8269E-04 knorm= 9.5845E-03 apxde= 2.2987E-06    *not conv.*     

               starting mcscf iteration...   5
 !timer:                                 cpu_time=     5.468 walltime=    12.081

 orbital-state coupling will be calculated this iteration.

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

 !timer: 2-e transformation              cpu_time=     0.834 walltime=     0.864

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv= 14 noldv= 14

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649413488     -432.3607809764        0.0000006326        0.0000010000
    2      -230.4272958786     -432.0231355062        0.0088755514        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534660098     -432.1493056375        0.0000004568        0.0000010000
    2      -230.4832368995     -432.0790765271        0.0000003144        0.0000010000
    3      -230.4172443949     -432.0130840226        0.0042272308        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6036079935     -432.1994476211        0.0000004541        0.0000010000
    2      -230.5760816749     -432.1719213025        0.0000003571        0.0000010000
    3      -230.4793320276     -432.0751716553        0.0064672615        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.757143374829148E-005
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99998355 pnorm= 2.4029E-02 rznorm= 5.2076E-07 rpnorm= 7.4646E-07 noldr=  9 nnewr=  9 nolds=  8 nnews=  8
 

 fdd(*) eigenvalues. symmetry block  1
   -22.522852  -22.505953   -2.317456   -1.670107   -1.432850   -1.011878

 qvv(*) eigenvalues. symmetry block  1
     0.066988    0.112593    0.288095    0.300987    0.326544    0.377378    0.619187    0.635561    0.694815    0.870067
     0.996511    1.144832    1.178250    1.209537    1.371282    1.536191    1.759797    1.908603    1.981896    3.052387
     3.193736    3.262377    3.269180    3.541420    3.638953    3.765143    3.907503    4.381448    4.836073    5.503347
     5.525092    6.176663    7.007080

 fdd(*) eigenvalues. symmetry block  2
   -22.520741  -22.505951   -2.046517   -1.308983   -1.193084

 qvv(*) eigenvalues. symmetry block  2
     0.082456    0.150497    0.293283    0.332521    0.364523    0.595945    0.726083    0.739821    0.806913    0.845558
     1.121809    1.314145    1.341931    1.366993    1.536937    1.679054    1.966055    2.089628    2.167603    2.667857
     3.029669    3.156405    3.241416    3.868863    3.956863    4.149328    4.481766    4.719310    5.511949    5.865459
     5.869992    6.173385    7.637887    9.484337

 fdd(*) eigenvalues. symmetry block  3
   -22.522761   -2.051022   -1.258045   -1.199173

 qvv(*) eigenvalues. symmetry block  3
     0.082398    0.286436    0.292775    0.332249    0.595523    0.737947    0.762814    0.845635    1.122774    1.313326
     1.365606    1.509181    1.678862    2.169949    2.666637    2.977345    3.026568    3.238077    3.539103    3.865362
     4.146028    4.484867    4.521036    5.507839    5.866715    6.174046

 fdd(*) eigenvalues. symmetry block  4
   -22.520575   -1.669926   -1.012039

 qvv(*) eigenvalues. symmetry block  4
     0.112560    0.300869    0.377140    0.530064    0.619347    0.694688    0.869008    1.144051    1.173224    1.370561
     1.534980    1.672616    1.758119    1.908708    2.357068    3.053561    3.262241    3.539519    3.848865    3.909484
     4.377517    4.685590    4.835568    5.501759    6.176270    6.679952    7.003732
 i,qaaresolved                     1 -0.973343928024979     
 i,qaaresolved                     2  0.167586212969606     
 i,qaaresolved                     3  0.230843841522769     

 qvv(*) eigenvalues. symmetry block  5
     0.350852    0.660773    0.990441    1.188594    1.344826    1.654927    1.827622    2.324502    3.073175    3.372872
     3.756223    4.150256    4.691439
 i,qaaresolved                     1 -0.618237640724614     

 qvv(*) eigenvalues. symmetry block  6
     0.280606    0.379023    0.643785    0.844808    1.122037    1.498196    1.817442    1.893235    1.930685    2.934504
     3.401133    3.937443    4.204510    4.360706    5.274147
 i,qaaresolved                     1 -0.578699480203627     

 qvv(*) eigenvalues. symmetry block  7
     0.279337    0.843857    1.102920    1.496436    1.816563    1.928532    2.808930    3.397040    4.201207    4.363697
 i,qaaresolved                     1  0.170892760345402     

 qvv(*) eigenvalues. symmetry block  8
     0.354004    0.990317    1.347998    1.826147    2.138451    2.322987    3.376863    3.752381    4.687237    5.172094

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.958 walltime=     1.367

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=   -230.5962667853 demc= 3.1004E-06 wnorm= 4.6057E-04 knorm= 5.7337E-03 apxde= 5.4477E-07    *not conv.*     

               starting mcscf iteration...   6
 !timer:                                 cpu_time=     6.426 walltime=    13.448

 orbital-state coupling will be calculated this iteration.

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

 !timer: 2-e transformation              cpu_time=     0.759 walltime=     0.788

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv= 15 noldv= 15

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649518397     -432.3607914673        0.0000002366        0.0000010000
    2      -230.4272693585     -432.0231089861        0.0090228192        0.0100000000

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  10 noldhv=  5 noldv=  5

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534776916     -432.1493173192        0.0000002496        0.0000010000
    2      -230.4831758381     -432.0790154657        0.0000004532        0.0000010000
    3      -230.4170554843     -432.0128951120        0.0057215131        0.0100000000

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6036187784     -432.1994584060        0.0000003975        0.0000010000
    2      -230.5761125195     -432.1719521471        0.0000004093        0.0000010000
    3      -230.4790901638     -432.0749297914        0.0057904752        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.936847670926039E-007
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-1.00000000 pnorm= 1.5700E-04 rznorm= 2.6075E-07 rpnorm= 9.2649E-07 noldr=  5 nnewr=  5 nolds=  2 nnews=  2
 

 fdd(*) eigenvalues. symmetry block  1
   -22.522862  -22.505886   -2.317430   -1.670076   -1.432812   -1.011848

 qvv(*) eigenvalues. symmetry block  1
     0.067001    0.112598    0.288104    0.301010    0.326554    0.377383    0.619197    0.635578    0.694842    0.870084
     0.996525    1.144861    1.178280    1.209563    1.371298    1.536205    1.759818    1.908638    1.981923    3.052422
     3.193766    3.262403    3.269203    3.541445    3.638982    3.765176    3.907538    4.381477    4.836104    5.503383
     5.525127    6.176692    7.007109

 fdd(*) eigenvalues. symmetry block  2
   -22.520750  -22.505884   -2.046482   -1.308942   -1.193050

 qvv(*) eigenvalues. symmetry block  2
     0.082453    0.150504    0.293284    0.332512    0.364538    0.595941    0.726096    0.739837    0.806925    0.845579
     1.121833    1.314155    1.341956    1.367003    1.536958    1.679076    1.966080    2.089646    2.167642    2.667881
     3.029693    3.156437    3.241435    3.868885    3.956894    4.149360    4.481801    4.719340    5.511980    5.865500
     5.870021    6.173419    7.637916    9.484363

 fdd(*) eigenvalues. symmetry block  3
   -22.522771   -2.050999   -1.258015   -1.199137

 qvv(*) eigenvalues. symmetry block  3
     0.082419    0.286451    0.292792    0.332268    0.595535    0.737980    0.762842    0.845669    1.122796    1.313348
     1.365645    1.509211    1.678900    2.169972    2.666674    2.977372    3.026597    3.238114    3.539133    3.865394
     4.146063    4.484899    4.521073    5.507870    5.866750    6.174076

 fdd(*) eigenvalues. symmetry block  4
   -22.520585   -1.669893   -1.012000

 qvv(*) eigenvalues. symmetry block  4
     0.112570    0.300875    0.377164    0.530075    0.619347    0.694701    0.869045    1.144068    1.173252    1.370581
     1.535018    1.672647    1.758151    1.908727    2.357099    3.053587    3.262277    3.539554    3.848896    3.909514
     4.377548    4.685627    4.835602    5.501795    6.176298    6.679985    7.003760
 i,qaaresolved                     1 -0.973321840127040     
 i,qaaresolved                     2  0.167313231377829     
 i,qaaresolved                     3  0.231733618726633     

 qvv(*) eigenvalues. symmetry block  5
     0.350649    0.660785    0.990370    1.188563    1.344783    1.654868    1.827603    2.324503    3.073204    3.372898
     3.756244    4.150294    4.691469
 i,qaaresolved                     1 -0.618107034859375     

 qvv(*) eigenvalues. symmetry block  6
     0.280623    0.379022    0.643792    0.844812    1.122058    1.498212    1.817466    1.893259    1.930705    2.934525
     3.401156    3.937475    4.204541    4.360736    5.274183
 i,qaaresolved                     1 -0.578800960686420     

 qvv(*) eigenvalues. symmetry block  7
     0.279337    0.843892    1.102939    1.496471    1.816585    1.928547    2.808957    3.397071    4.201252    4.363728
 i,qaaresolved                     1  0.170730902264471     

 qvv(*) eigenvalues. symmetry block  8
     0.354057    0.990374    1.348026    1.826195    2.138480    2.323006    3.376889    3.752421    4.687279    5.172126

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.876 walltime=     1.270

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=   -230.5962673334 demc= 5.4815E-07 wnorm= 3.1495E-06 knorm= 5.0348E-05 apxde= 4.2304E-11    *not conv.*     

               starting mcscf iteration...   7
 !timer:                                 cpu_time=     7.302 walltime=    14.719

 orbital-state coupling will be calculated this iteration.

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

 !timer: 2-e transformation              cpu_time=     0.769 walltime=     0.799

 mosort: allocated sort2 space, avc2is=  2621208211 available sort2 space, avcisx=  2621208463

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv= 15 noldv= 15

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.7649519906     -432.3607916182        0.0000002160        0.0000010000
    2      -230.4272524771     -432.0230921047        0.0096206759        0.0100000000

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  10 noldhv=  5 noldv=  5

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.5534777357     -432.1493173633        0.0000001939        0.0000010000
    2      -230.4831753639     -432.0790149915        0.0000007067        0.0000010000
    3      -230.4170526631     -432.0128922907        0.0053803284        0.0100000000

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -230.6036188788     -432.1994585064        0.0000004153        0.0000010000
    2      -230.5761126984     -432.1719523261        0.0000004117        0.0000010000
    3      -230.4790879443     -432.0749275719        0.0057447888        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.628665117667171E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 4.1979E-07 rpnorm= 1.2241E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

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
 i,qaaresolved                     1 -0.973321493940772     
 i,qaaresolved                     2  0.167311238809862     
 i,qaaresolved                     3  0.231741623592518     

 qvv(*) eigenvalues. symmetry block  5
     0.350647    0.660786    0.990369    1.188563    1.344783    1.654868    1.827603    2.324503    3.073205    3.372898
     3.756244    4.150294    4.691469
 i,qaaresolved                     1 -0.618105893904590     

 qvv(*) eigenvalues. symmetry block  6
     0.280623    0.379022    0.643792    0.844812    1.122058    1.498212    1.817466    1.893259    1.930705    2.934526
     3.401156    3.937476    4.204541    4.360737    5.274183
 i,qaaresolved                     1 -0.578801587704251     

 qvv(*) eigenvalues. symmetry block  7
     0.279337    0.843893    1.102939    1.496471    1.816586    1.928548    2.808957    3.397071    4.201253    4.363728
 i,qaaresolved                     1  0.170729918484340     

 qvv(*) eigenvalues. symmetry block  8
     0.354057    0.990374    1.348026    1.826196    2.138480    2.323007    3.376890    3.752421    4.687279    5.172127

 restrt: restart information saved on the restart file (unit= 13).
 !timer: mcscf iteration                 cpu_time=     0.934 walltime=     1.359

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    7 emc=   -230.5962673335 demc= 4.3912E-11 wnorm= 2.9029E-07 knorm= 2.3343E-08 apxde= 3.9171E-15    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.200 total energy=     -230.764951991, rel. (eV)=   0.000000
   DRT #2 state # 1 wt 0.200 total energy=     -230.553477736, rel. (eV)=   5.754510
   DRT #2 state # 2 wt 0.200 total energy=     -230.483175364, rel. (eV)=   7.667535
   DRT #3 state # 1 wt 0.200 total energy=     -230.603618879, rel. (eV)=   4.390099
   DRT #3 state # 2 wt 0.200 total energy=     -230.576112698, rel. (eV)=   5.138581
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
  occ(*)=     1.95104834     0.37766436     0.16635412     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  6    -  B2g
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     1.66016563     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  7    -  B3g
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     1.48891533     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11
  occ(*)=     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  8    -  Au 
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.35585221     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO    9        MO   10        MO   11
  occ(*)=     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        723 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        Ag  partial gross atomic populations
   ao class       1Ag        2Ag        3Ag        4Ag        5Ag        6Ag 
    C1_ s       0.024115   1.965539   0.548289   0.570135  -0.031396  -0.029880
    C1_ p       0.003468  -0.000704   0.067236   0.017042   0.323798   0.397213
    C1_ d       0.000470   0.000179   0.004670  -0.062235   0.037144  -0.079263
    C2_ s       1.965551   0.024164   1.148780   0.279751  -0.058322  -0.013963
    C2_ p       0.000801   0.007425   0.143610   0.540706   0.669714   0.894395
    C2_ d      -0.000546   0.001050   0.009672  -0.018646   0.075114  -0.021358
    H1_ s       0.001179   0.000814   0.024209   0.432804   0.323613   0.606223
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
    C1_ p      -0.000455   0.001641  -0.013552   0.204695   0.739347   0.000000
    C1_ d      -0.000102   0.000126   0.032980  -0.188062  -0.012470   0.000000
    C2_ s       1.958708   0.025923   0.464657   0.158806  -0.020885   0.000000
    C2_ p       0.013886   0.014358   0.277090   0.477455   0.482196   0.000000
    C2_ d       0.003923   0.004386  -0.053223  -0.386015  -0.014182   0.000000
    H1_ s       0.000617   0.000533   0.233173   0.517384   0.663779   0.000000
    H1_ p       0.000315   0.000354   0.029126  -0.000702  -0.005710   0.000000
    H2_ s       0.000067   0.000410   0.122706   1.114484   0.284369   0.000000
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
    C2_ d       0.001844   0.024820   0.030207  -0.018102   0.000000   0.000000
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
    H2_ s       0.000202   0.667375   0.881495   0.000000   0.000000   0.000000
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

 !timer: mcscf                           cpu_time=     8.284 walltime=    16.332
