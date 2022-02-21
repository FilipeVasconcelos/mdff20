#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "io.h"
#include "cell.h"
#include "pbc.h"
#include "config.h"
#include "radial_distribution.h"
#include "tools.h"

int read_rdf (char * controlfn)
{
    int c;
    char buffer[MAX_LEN+1];
    FILE * fp;

    fp = fopen (controlfn, "r");
    if (NULL == fp ) {
        pError("opening control file (reading rdf) ");
        return (-1);
    }

    while (EOF != fscanf(fp, "%s\n", buffer)) {
        if (strcmp(buffer,"resg") == 0 ) {
            c=fscanf(fp,"%lf",&resg);
        }
        if (strcmp(buffer,"nconf") == 0 ) {
            c=fscanf(fp,"%d",&nconf);
        }
   }
   fclose(fp);
   return 0;
}
void default_rdf(){
    resg  = 0.1;
    nconf = 0;
    cutgr = 0.5 * dmin3(simuCell.w[0],simuCell.w[1],simuCell.w[2]);
    nbins=(int) (cutgr/resg)+1;
}
void check_rdf(){
}
void info_rdf(){
    if (ionode) {
        SEPARATOR;
        printf("radial distribution function info \n");
        LSEPARATOR;
        printf("resolution of g(r) function resg        = %.5f \n",resg);
        printf("number of config. in TRAJFF             = %d \n", nconf );
        printf("number of points in g(r)                = %d \n", nbins );
        printf("save radial distribution function to file RDFFF\n");
        printf("read configuration from TRAJFF");
        putchar('\n');
    }
}
void alloc_rdf(){
    int igr,it1,it2;
    gr=malloc(nbins*sizeof *gr);
    for (igr=0;igr<nbins;igr++){
        //printf("%d %d\n",igr,ntype);
        for(it1=0;it1<ntype;it1++){
            for(it2=0;it2<ntype;it2++){
            gr[igr][it1][it2]=0.0;
            }
        }
    }
}

void free_rdf(){
    free(gr);
}

int rdf_readtraj(){

    int ic,igr;
    FILE *fp;
    fp = fopen ("TRAJFF","r");
    for(ic=0;ic<nconf;ic++) {
        printf("read conf %d\n",ic);
        read_next_traj(fp);
        sample_config(1);
        get_gr();
    }
    //closing TRAJFF 
    if (fclose(fp))     {
       io_node printf("error closing file.");
       return -1;
    }
    /* debug
    for (igr=0;igr<nbins;igr++){
        printf("%d\n",gr[igr][0][0]);
    }
    */
    write_rdf();
    printf("write RDFFF");
    return 0;
}

int write_rdf(){

    int igr,mp;
    int npairs;
    int it1,it2;
    double rr,vol;
    double (*grr)[NPAIRSMAX];
    double average_volume;

    npairs=ntype*(ntype+1)/2;

    grr=malloc(nbins*sizeof *grr);
    FILE * fp;
    fp = fopen ("RDFFF", "w");
    if (NULL == fp) {
        pError("opening file :");
        printf(" %s\n","RDFFF");
        return -1;
    }

    for (igr=0;igr<nbins;igr++)
    {
        rr  = ((double) igr+0.5)*resg;
        vol = 4 * PI * ( resg * rr * rr + ( pow(resg,3) ) / 12 );
        //vol = vol / average_volume;  //! only  make sense for trajectory in NPT otherwise the average volume is the current volume ) 
        vol = vol / simuCell.Omega;  //! only  make sense for trajectory in NPT otherwise the average volume is the current volume ) 
        // all - all pairs 
        grr[igr][mp] = ((double) gr[igr][0][0]) / ( nconf * vol * nion * nion );
        //grr ( 0 , igr ) = REAL ( gr ( 0 , 0 , igr ) ) / ( ngr * vol * natm * natm )
        // type pairs
        mp = 0;
        for(it1=0;it1<ntype;it1++){
            for(it2=0;it2<ntype;it2++){
//#ifdef debug
//       io_node WRITE ( stdout , '(a,3i5)' ) 'debug ( pair ) : ', mp , it1 , it2
//#endif        
            if ( mp < 0 || mp > npairs ) {
                pError("ERROR out of bound of gr in write_rdf");
            }
            //nr[mp] = it1;          
            grr[igr][mp] = ((double) gr[igr][it1][it2]) / (double) ( nconf * vol * nionit[it1] * nionit[it2] ); 
            //grr ( mp , igr ) = REAL ( gr ( it1 , it2 , igr ) ) / REAL ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
            mp++;
            }
        }
        for (mp=0;mp<npairs;mp++) fprintf(fp,EE2,rr,grr[igr][mp]);fprintf(fp,"\n");
/*
      if ( ionode ) then
#ifdef GFORTRAN
        WRITE ( kunit_GRTFF ,'(8e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(8e20.10)') rr , ( REAL ( grr ( mp , igr ) , kind = dp ) * 4.0_dp * pi * rr * rr * REAL ( natmi(nr(mp)) * vol , kind = dp ) , mp = 0 , npairs )
#else
        WRITE ( kunit_GRTFF ,'(<npairs+2>e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(<npairs+2>e20.10)') rr , ( REAL ( grr ( mp , igr ) , kind = dp ) * 4.0_dp * pi * rr * rr * REAL ( natmi(nr(mp)) * vol , kind = dp) , mp = 0 , npairs )
#endif
*/
    }
    free(grr);
    return 0;
}

void init_rdf(char * controlfn ){
    printf("here 0\n");
    /* default values */
    default_rdf();
    printf("here 2\n");
    /* read parameters */
    read_rdf(controlfn);
    printf("here 3\n");
    /* check parameters */
    //check_rdf();
    /* alloc memory */
    alloc_rdf();
    printf("here 4\n");
    /* print information */
    info_rdf();
    rdf_readtraj();
    printf("here 5\n");
}
void get_gr(){
  
    int  igr ; 
    int  ia , ja , ierr , ita , jta ;
    double cut2 , rijsq , rr ;
    double rxi , ryi , rzi ;
    double rxij , ryij , rzij;
    /***************************************
            cut-off half box
    ****************************************/
    cut2 = cutgr * cutgr; 
    /***************************************
            cartesian to direct
     ***************************************/
    kardir ( nion , rx , ry , rz , simuCell.B ) ;

    for(ia=atomDec.iaStart;ia<atomDec.iaEnd;ia++) {
        ita = typia[ia];
        rxi = rx[ia];
        ryi = ry[ia];
        rzi = rz[ia];
        for (ja=0;ja<nion;ja++) {
            if ( ja != ia ) {  
                jta = typia[ja];
                rxij = rxi - rx[ja];
                ryij = ryi - ry[ja];
                rzij = rzi - rz[ja];
                //printf("%f %f %f\n",rx[ja],ry[ja],rz[ja]);
                pbc(&rxij,&ryij,&rzij);
                rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                if ( rijsq < cut2 ) {
                    rr = sqrt(rijsq);
                    igr = (int) (rr / resg ) ;
                    if ( igr < 0 && igr > nbins-1 ) {
                        printf("ERROR out of bound of gr in gr_main");
                        exit(-1);   
                    } 
                    // all pairs
                    //printf("%d %d %d %d %d\n",igr,ita,jta,gr[igr][0][0],gr[igr][ita][jta]);
                    gr[igr][0][0]++;
                    gr[igr][ita][jta]++;
//!          gr ( 0    , 0   , igr ) = gr ( 0   ,   0 , igr ) + 1 
//!          gr ( ita  , jta , igr ) = gr ( ita , jta , igr ) + 1
                }
            }
        }
    }
    /***************************************
            direct to cartesian
     ***************************************/
    dirkar ( nion , rx , ry , rz , simuCell.A ) ;
}

/*MODULE radial_distrib 

  USE io,               ONLY :  ionode, stdout, stdin, stderr
  USE constants,        ONLY :  dp
  USE mpimdff

  implicit none

  integer :: nbins            !< (internal) number of bins in g(r) distribution
  integer :: npairs            !< (internal) number of bins in g(r) distribution
  integer :: nconf            !< number of configurations used in g(r) calculation
  real(kind=dp) :: cutgr      !< radial cut-off 
  real(kind=dp) :: resg       !< resolution in g(r) distribution 
  !> g(r) function ( bin x ntype x ntype )
  integer, dimension(:,:,:), allocatable :: gr  

CONTAINS

! *********************** SUBROUTINE grcalc_init *******************************
!
!> \brief
!! initialize radial distribution calculation parameters
!
! ******************************************************************************
SUBROUTINE gr_init

  USE config,                   ONLY :  simu_cell
  USE control,                  ONLY :  calc

  implicit none

  ! local
  integer            :: ioerr 
!28/05/13  integer            :: npangr, i
  character(len=132) :: filename

  namelist /grtag/   nconf , &
                     cutgr , &
                     resg  

  if ( calc .ne. 'gr' ) return

  CALL gr_default_tag
  
  ! =================
  !  read grtag tags
  ! =================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , grtag , iostat=ioerr )
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : grtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : grtag wrong tag'
   STOP
  endif
  CLOSE ( stdin )

  ! ==========================================
  ! define a new resolution to be 2^N points
  ! ==========================================
 nbins=int(cutgr/resg)+1
! print*,cutgr,resg,nbins
!28/05/13 ! i = 1
!28/05/13 ! do while ( 2**i .lt. nbins )
!28/05/13 !    i = i + 1
!28/05/13 ! enddo
!28/05/13 ! npangr = i
!28/05/13 ! nbins = 2** npangr
!28/05/13 ! resg = cutgr / REAL ( nbins  )

  CALL gr_print_info(stdout)

  return 
 
END SUBROUTINE gr_init

! *********************** SUBROUTINE gr_alloc **********************************
!
!> \brief
!! Allocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_alloc

  USE control,                  ONLY :  calc
  USE config,                   ONLY :  ntype

  implicit none

  if ( calc .ne. 'gr' .and. calc .ne. 'rmc' ) return
  allocate(  gr ( 0 : nbins - 1 , 0 : ntype , 0 : ntype ) ) 
  gr = 0      

  return 
 
END SUBROUTINE gr_alloc


! *********************** SUBROUTINE gr_dealloc ********************************
!
!> \brief
!! Deallocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_dealloc

  USE control,                  ONLY :  calc

  implicit none

  if ( calc .ne. 'gr' ) return
  
  deallocate( gr )

  return 
 
END SUBROUTINE gr_dealloc


! *********************** SUBROUTINE gr_default_tag ****************************
!
!> \brief
!! set default values to gr tag
!
! ******************************************************************************
SUBROUTINE gr_default_tag

  USE config,           ONLY : simu_cell

  implicit none

  ! ===============
  !  default value
  ! ===============
  resg = 0.1_dp
  nconf = 0
  cutgr=0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC)

  return

END SUBROUTINE gr_default_tag


! *********************** SUBROUTINE gr_print_info *****************************
!
!> \brief
!! print infog on g(r) calculation
!
! ******************************************************************************
SUBROUTINE gr_print_info(kunit)

  USe control,                  ONLY :  calc

  implicit none
 
  ! local
  integer :: kunit

   if ( ionode ) then
                  WRITE ( kunit ,'(a,f10.5,a)')         'resolution of g(r) function resg     = ',resg,' new value to have 2^N points in g(r)'
                  WRITE ( kunit ,'(a,i5)')              'number of points in g(r)             = ',nbins
                  WRITE ( kunit ,'(a)')                 'save radial_distribution in file     :   GRTFF' 
      if ( calc .eq. 'gr' )     then 
                  WRITE ( kunit ,'(a)')                 'read configuration from file         :   TRAJFF'
                  blankline(kunit)
                  WRITE ( kunit ,'(a,i5)')              'number of config. in TRAJFF          = ',nconf        
                  blankline(kunit)
      endif
   endif 
  return

END SUBROUTINE gr_print_info

! *********************** SUBROUTINE grcalc ************************************
!
!> \brief
!! main driver of radial distribution function calculation
!! this subroutine read the trajectory, Allocate, call the  
!
! ******************************************************************************
SUBROUTINE grcalc

  USE control,                  ONLY :  itraj_format , trajff_data
  USE config,                   ONLY :  system , natm , ntype , rx , ry , rz , atype , &
                                        rhoN , config_alloc , simu_cell , atypei , itype, natmi, &
                                        coord_format_allowed , atom_dec , read_traj , read_traj_header
  USE io,                       ONLY :  kunit_TRAJFF , kunit_GRTFF , kunit_NRTFF
  USE constants,                ONLY :  pi 
  USE cell,                     ONLY :  lattice , dirkar , periodicbc, kardir
  USE time,                     ONLY :  grtimetot_comm

  implicit none

  ! local 
  integer                                              :: ic , ngr , igr 
  integer                                              :: it1 , it2 , mp , ierr 
  real(kind=dp),     dimension ( : , : ) , allocatable :: grr 
  integer,           dimension ( : )     , allocatable :: nr 
  character(len=15), dimension ( : )     , allocatable :: cint
  real(kind=dp)                                        :: rr , vol
  real(kind=dp)                                        :: ttt1 , ttt2      
  real(kind=dp)                                        :: average_volume
  real(kind=dp)                                        :: rhoN_av 


  rhoN_av = 0.0_dp
  average_volume = 0.0_dp      

  OPEN ( kunit_GRTFF , file = 'GRTFF' )
  OPEN ( kunit_NRTFF , file = 'NRTFF' )

  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')
  CALL read_traj_header ( kunit_TRAJFF , itraj_format )
  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')

  CALL lattice ( simu_cell ) 
  rhoN = REAL ( natm , kind = dp )  / simu_cell%omega 

#ifdef debug
  write(*,*) simu_cell%A
  write(*,*) simu_cell%omega
  write(*,*) simu_cell%ANORM
#endif

  CALL gr_init

  CALL print_general_info( stdout )

#ifdef debug
  write(*,*)
  write(*,*) simu_cell%A
  write(*,*) simu_cell%omega
  write(*,*) simu_cell%ANORM
#endif
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms' )
  CALL gr_alloc

  npairs =  ntype * ( ntype + 1 ) / 2
  allocate ( grr ( 0 : npairs , 0 : nbins-1 ) , nr ( 0 : npairs ) , cint ( 0 : npairs  ))
  grr  = 0.0_dp
  nr   = 0
  cint = ''

#ifdef debug
  if ( ionode ) then 
    WRITE ( stdout , '(a,2i6)' ) 'debug : atom decomposition istart, iend ', atom_dec%istart , atom_dec%iend
    WRITE ( stdout , '(a,2i6)' ) 'debug : number of type npairs ', npairs
  endif
#endif

  CALL typeinfo_init

  ngr = 0
  do ic = 1, nconf

    CALL read_traj ( kunit_TRAJFF , itraj_format , trajff_data ) 

    CALL lattice ( simu_cell )
    rhoN_av = rhoN_av + ( REAL ( natm ,kind=dp )  / simu_cell%omega )
    average_volume = average_volume + simu_cell%omega    
    io_node WRITE ( stdout , '(a,i6,a,i6,a,f12.3)' ) 'config : [ ',ic,' / ',nconf,' ]   vol : ',simu_cell%omega
#ifdef debug
    CALL distance_tab
!    print*,simu_cell%omega,average_volume/ REAL(ic,kind=dp)
#endif
    CALL kardir     ( natm , rx , ry , rz , simu_cell%B ) 
    CALL periodicbc ( natm , rx , ry , rz  )
    CALL dirkar     ( natm , rx , ry , rz , simu_cell%A ) 

#ifdef debug
    CALL distance_tab
#endif

    ngr=ngr+1 
    ! ==========================
    !  calc radial_distribution 
    ! ==========================  
    call gr_main 

  enddo !nconf 
  rhoN_av = rhoN_av / REAL(nconf,kind=dp)      
  average_volume = average_volume / REAL(nconf,kind=dp)
  if ( ionode .and. average_volume .ne. simu_cell%omega ) write(stdout,'(a,e16.8)') 'average volume : ',average_volume

  ! ===========================================
  !        merge results  
  ! ===========================================
#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif

  CALL MPI_ALL_REDUCE_INTEGER ( gr(:,0,0), nbins )
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      CALL MPI_ALL_REDUCE_INTEGER ( gr(:,it1,it2), nbins )
    enddo  
  enddo
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  grtimetot_comm = grtimetot_comm + ( ttt2 - ttt1 )
#endif


#ifdef debug
  do igr=0, nbins
    io_node WRITE (stdout , '(a,5i6)') 'debug ( total ) : ',igr,gr(igr,1,1)
  enddo
#endif

  ! ======================================================= 
  !  write output files GRTFF , NRTFF
  ! ======================================================= 
  cint ( 0) = atypei ( 0 )//' - '//atypei ( 0 )
  mp = 1
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      cint(mp) = atypei(it1)//' - '//atypei(it2)
      mp = mp + 1 
    enddo
  enddo

#ifdef GFORTRAN
  io_node WRITE ( kunit_GRTFF , '(8a)' )          '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
  io_node WRITE ( kunit_NRTFF , '(8a)' )          '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
#else
  io_node WRITE ( kunit_GRTFF , '(<npairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
  io_node WRITE ( kunit_NRTFF , '(<npairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
#endif

  nr ( 0 ) = 0 

  do igr = 0 , nbins-1
    rr  = ( REAL ( igr ,kind=dp )+0.5d0)*resg
    vol = 4.d0 * pi * ( resg * rr* rr + ( resg**3 ) / 12.d0 )
    vol = vol / average_volume  ! only  make sense for trajectory in NPT otherwise the average volume is the current volume ) 
    ! all - all pairs 
    grr ( 0 , igr ) = REAL ( gr ( igr , 0 , 0 ) , kind = dp ) / ( ngr * vol * natm * natm )
!    grr ( 0 , igr ) = REAL ( gr ( 0 , 0 , igr ) ) / ( ngr * vol * natm * natm )
    ! type pairs
    mp = 1
    do it1 = 1 , ntype
      do it2 = it1 , ntype
#ifdef debug
       io_node WRITE ( stdout , '(a,3i5)' ) 'debug ( pair ) : ', mp , it1 , it2
#endif        
        if ( mp .lt. 0 .and. mp .gt. npairs ) then
          WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
          STOP
        endif 
        nr ( mp ) = it1          
        grr ( mp , igr ) = REAL ( gr ( igr , it1 , it2 ) , kind = dp) / REAL ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) , kind = dp )  
!        grr ( mp , igr ) = REAL ( gr ( it1 , it2 , igr ) ) / REAL ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
        mp = mp + 1
      enddo
    enddo
      if ( ionode ) then
#ifdef GFORTRAN
        WRITE ( kunit_GRTFF ,'(8e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(8e20.10)') rr , ( REAL ( grr ( mp , igr ) , kind = dp ) * 4.0_dp * pi * rr * rr * REAL ( natmi(nr(mp)) * vol , kind = dp ) , mp = 0 , npairs )
#else
        WRITE ( kunit_GRTFF ,'(<npairs+2>e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(<npairs+2>e20.10)') rr , ( REAL ( grr ( mp , igr ) , kind = dp ) * 4.0_dp * pi * rr * rr * REAL ( natmi(nr(mp)) * vol , kind = dp) , mp = 0 , npairs )
#endif
      endif
  enddo


  CLOSE ( kunit_NRTFF )
  CLOSE ( kunit_GRTFF )

  CLOSE( kunit_TRAJFF )

  !CALL static_struc_fac ( grr , nbins , npairs ) 

  deallocate ( grr ,nr , cint)
  CALL gr_dealloc

  return

END SUBROUTINE grcalc

! *********************** SUBROUTINE gr_main ***********************************
!
!> \brief
!! based on Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE gr_main 

  USE control,                  ONLY :  myrank
  USE config,                   ONLY :  natm , natmi , rx , ry , rz , atype , simu_cell , ntype , itype, atom_dec
  USE time,                     ONLY :  grtimetot
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer :: ia , ja , ierr , ita , jta 
  integer :: igr 
  real(kind=dp) :: cut2 , rijsq , rr 
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: ttt1 , ttt2      

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
  ! =======================
  !  cut-off half box
  ! =======================
  cut2 = cutgr * cutgr 
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend

    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ita = itype( ia )
    do ja = 1, natm
      if ( ja .ne. ia ) then  
        jta = itype( ja )
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        if ( rijsq.lt.cut2 ) then
          rr = SQRT ( rijsq )
          igr = INT ( rr / resg ) 
          if ( igr .lt. 0 .and. igr .gt. nbins-1 ) then
            WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
            STOP
          endif 
          ! all pairs
          gr ( igr , 0    , 0   ) = gr ( igr , 0   ,   0 ) + 1 
          gr ( igr , ita  , jta ) = gr ( igr , ita , jta  ) + 1
!          gr ( 0    , 0   , igr ) = gr ( 0   ,   0 , igr ) + 1 
!          gr ( ita  , jta , igr ) = gr ( ita , jta , igr ) + 1
        endif
      endif
    enddo
  enddo
  
#ifdef debug2
  do igr=0, nbins-1
    WRITE (stdout , '(a,5i12)') 'debug: ',myrank,igr,gr(igr,0,0)
  enddo
#endif 

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  grtimetot = grtimetot + ( ttt2 - ttt1 )
#endif

  return
 
END SUBROUTINE gr_main

! *********************** SUBROUTINE static_struc_fac **************************
!
!
! ******************************************************************************
!SUBROUTINE static_struc_fac ( gr , nbins , npairs )

!  USE io,                       ONLY :  ionode , kunit_STRFACFF , stdout 
!  USE config,                   ONLY :  rhoN
!  USE constants,                ONLY :  pi , tpi , imag

!  implicit none
!  ! global
!  integer :: nbins , npairs
!  real(kind=dp) :: gr ( 0 : nbins-1 , 0 : npairs ) , rr
!
!  ! local
!  integer :: i , j , is , NN
!  real(kind=dp) :: q , qj , ri , rip
!  complex(kind=dp) :: Sk
!  real(kind=dp) , dimension (:) , allocatable  :: stat_str
!  real(kind=dp) , dimension (:) , allocatable :: in 
!  real(kind=dp) , dimension (:,:) , allocatable :: Uji 
!  complex(kind=dp)   , dimension (:) , allocatable :: out
!  real(kind=dp) :: res , shift
!  real(kind=dp) :: x , k
!
!  io_node WRITE ( stdout , '(a)' ) 'in static_struc_fac'
!
!  allocate ( in ( nbins ) , out ( nbins /2 + 1 ) )
!
!  in  = ( 0.0,0.0)
!  out = ( 0.0,0.0)
!  ! ========
!  !   FFT
!  ! ========
!  do i=0,nbins-1
!    in ( i + 1 ) = gr ( i , 0 ) 
!  enddo
!
!!  CALL fft_1D_complex ( in , out , nbins )
!  CALL fft_1D_real(in,out,nbins)
!
!  do i= 1 , nbins/2+1
!    q = ( dble ( i )  + 0.5_dp ) / REAL ( nbins ) / resg
!    Sk = 1.0_dp + rhoN * out( i + 1 )  
!    io_node WRITE ( 20000 , '(3e16.8)' )  q , Sk  
!  enddo
!
!  deallocate ( in , out )
!
!! other version
!  ! Uji (eq 12) J Phys Cndens Matter V17 2005 )
!!  allocate ( Uji ( nbins , nbins ) ) 
!
!  do i = 1 , nbins
!!    ri  = ( dble ( i )  + 0.5_dp )  * res
!    rip = ( dble ( i + 1 )  + 0.5_dp )  * res
!!    do j = 1 , nbins
!      qj = tpi * REAL ( j ) + 0.5_dp / REAL ( nbins / 2 + 1 ) / resg
!      Uji ( j , i ) = SR ( ri , qj ) - SR ( rip , qj ) 
! !     Uji ( j , i ) = Uji ( j , i ) / qj  
!    enddo
! ! enddo
!  Uji = 2.0_dp * tpi * Uji
!!
!  do i= 1 , nbins/2+1
!    q = tpi * ( dble ( i )  + 0.5_dp ) / REAL ( nbins / 2 + 1) / resg
!!    do j = 1 , nbins
!      Sk =  Sk + Uji ( j , i ) * ( gr ( j , 0 ) -1.0_dp )  
!    enddo
!    io_node WRITE ( 30000 , '(3e16.8)' )  q , Sk
!  enddo
!
!
!  deallocate ( Uji )

! test purpose because I'm dumb
! I was enable to get k in "real life" 
! I did a simple discret case  ( N = 4 ) 
! to find the relation between k and q 

!  NN = 4
!  allocate ( in ( NN ) , out  ( NN ) )
!
!  in = ( 0.0,0.0)
!  is = 3 
!  in(is) = ( 1.0,0.0)
!  print*,'in in '
!  CALL fft_1D_complex ( in , out , NN )
!
!  res = 2.0_dp
!  shift= (is-1) * res
!  write( stdout , '(8a)' ) '            x       Re in(i)        Im in(i)            k          Re out(i)       Im out(i)        Re phase        Im phase'
!  do i=1,NN
!    x = REAL(i-1)*res
!    k = REAL(i-1) / REAL ( NN ) 
!    q = REAL(i-1) / REAL ( NN ) / res
!    write( stdout , '(8e16.8)' ) x , in(i) , k , out(i) , exp ( -imag * tpi * k * (is-1) )
!    write( stdout , '(8e16.8)' ) x , in(i) , q , out(i) , exp ( -imag * tpi * q * shift )
!  enddo
!
!  deallocate ( in ,out )

!  return
!CONTAINS

!real(kind=dp) function SR(r,qj)
!  implicit none
!  real(kind=dp) :: r , qj  
!  SR = SIN ( qj * r ) / qj / qj 
!  SR = SR - r * COS ( qj * r ) / qj
!end function 

!END SUBROUTINE static_struc_fac


SUBROUTINE read_GRTFF( grr , filename)

  USE constants,         only : dp
  implicit none

  integer :: mp , bin
  real(kind=dp) , dimension ( 0:npairs , 0: nbins-1 ) :: grr
  real(kind=dp) :: rr
  character(*) :: filename

  OPEN(UNIT=1000,FILE=filename)
    read(1000,*)      
    do bin=0,nbins-1
      read(1000,*) rr, ( grr ( mp , bin ), mp = 0 , npairs )
    enddo 
  CLOSE(1000)

END SUBROUTINE

END MODULE radial_distrib 
! ===== fmV =====
*/
