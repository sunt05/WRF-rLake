! MODULE module_sf_lake
!
! The lake scheme was retrieved from the Community Land Model version 4.5
! (Oleson et al. 2013) with some modifications by Gu et al. (2013). It is a
! one-dimensional mass and energy balance scheme with 20-25 model layers,
! including up to 5 snow layers on the lake ice, 10 water layers, and 10 soil
! layers on the lake bottom. The lake scheme is used with actual lake points and
! lake depth derived from the WPS, and it also can be used with user defined
! lake points and lake depth in WRF (lake_min_elev and lakedepth_default).
! The lake scheme is independent of a land surface scheme and therefore
! can be used with any land surface scheme embedded in WRF. The lake scheme
! developments and evaluations were included in Subin et al. (2012) and Gu et al. (2013)
!
!   Subin et al. 2012: Improved lake model for climate simulations, J. Adv. Model.
!   Earth Syst., 4, M02001. DOI:10.1029/2011MS000072;
!   Gu et al. 2013: Calibration and validation of lake surface temperature simulations
!   with the coupled WRF-Lake model. Climatic Change, 1-13, 10.1007/s10584-013-0978-y.
!
! USE module_wrf_error
! USE module_model_constants, ONLY : rcp
!
! IMPLICIT NONE
!
! SUBROUTINE Lake( t_phy        ,p8w            ,dz8w         ,qvcurr          ,&  !i
!      u_phy        ,v_phy          , glw         ,emiss           ,&
!      rainbl       ,dtbl           ,swdown       ,albedo          ,&
!      xlat_urb2d   ,z_lake3d       ,dz_lake3d    ,lakedepth2d     ,&
!      watsat3d     ,csol3d         ,tkmg3d       ,tkdry3d         ,&
!      tksatu3d     ,ivgtyp         ,ht           ,xland           ,&
!      iswater, xice, xice_threshold, lake_min_elev                ,&
!      ids          ,ide            ,jds          ,jde             ,&
!      kds          ,kde            ,ims          ,ime             ,&
!      jms          ,jme            ,kms          ,kme             ,&
!      its          ,ite            ,jts          ,jte             ,&
!      kts          ,kte                                           ,&
!      h2osno2d     ,snowdp2d       ,snl2d        ,z3d             ,&  !h
!      dz3d         ,zi3d           ,h2osoi_vol3d ,h2osoi_liq3d    ,&
!      h2osoi_ice3d ,t_grnd2d       ,t_soisno3d   ,t_lake3d        ,&
!      savedtke12d  ,lake_icefrac3d                                ,&
!      lakemask                                          ,&
!      hfx          ,lh             ,grdflx       ,tsk             ,&  !o
!      qfx          ,t2             ,th2          ,q2 )
!
!   !==============================================================================
!   ! This subroutine was first edited by Hongping Gu and Jiming Jin for coupling
!   ! 07/20/2010
!   !==============================================================================
!   IMPLICIT NONE
!
!   REAL    , PARAMETER :: r_d          = 287.
!   REAL    , PARAMETER :: cp           = 7.*r_d/2.
!   REAL    , PARAMETER :: rcp          = r_d/cp
!
!   INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)
!   INTEGER,PARAMETER  ::     lbp = 1
!   INTEGER,PARAMETER  ::     ubp = 1
!   integer,parameter  ::     lbc = 1                        ! column-index bounds
!   integer,parameter  ::     ubc = 1
!   INTEGER,PARAMETER  ::     num_shlakec       = 1
!   INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
!   INTEGER,PARAMETER  ::     num_shlakep       = 1
!   INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
!   INTEGER,PARAMETER  ::     pcolumn(1)        = 1
!   INTEGER,PARAMETER  ::     pgridcell(1)      = 1
!   INTEGER,PARAMETER  ::     cgridcell(1)      = 1
!   INTEGER,PARAMETER  ::     clandunit(1)      = 1
!
!   INTEGER,PARAMETER  ::     begg = 1
!   INTEGER,PARAMETER  ::     endg = 1
!   INTEGER,PARAMETER  ::     begl = 1
!   INTEGER,PARAMETER  ::     endl = 1
!   INTEGER,PARAMETER  ::     begc = 1
!   INTEGER,PARAMETER  ::     endc = 1
!   INTEGER,PARAMETER  ::     begp = 1
!   INTEGER,PARAMETER  ::     endp = 1
!
!   INTEGER,PARAMETER  ::     column    =1
!   LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.
!
!
!
!   REAL(r8), PARAMETER :: vkc    = 0.4_r8
!   REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
!   REAL(r8), PARAMETER :: grav   = 9.80616_r8
!   REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
!   REAL(r8), PARAMETER :: tfrz   = 273.16_r8
!   REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
!   REAL(r8), PARAMETER :: denice = 0.917e3_r8
!   REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
!   REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
!   REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
!   REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
!   REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
!   REAL(r8), PARAMETER :: rair   = 287.0423_r8
!   REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
!   REAL(r8), PARAMETER :: tcrit  = 2.5
!   REAL(r8), PARAMETER :: tkwat  = 0.6
!   REAL(r8), PARAMETER :: tkice  = 2.290
!   REAL(r8), PARAMETER :: tkairc = 0.023
!   REAL(r8), PARAMETER :: bdsno = 250.
!
!   REAL(r8),  PARAMETER :: spval = 1.e36
!
!   REAL(r8), PARAMETER  ::     depth_c = 300.
!
!   REAL(r8), PARAMETER :: wimp   = 0.05
!   REAL(r8), PARAMETER :: ssi    = 0.033
!   REAL(r8), PARAMETER :: cnfac  = 0.5
!
!   INTEGER,PARAMETER :: istsoil = 1
!   REAL(r8) :: dtime = 60
!
!   REAL(r8) :: zlak(1:nlevlake)     !lake z  (layers)
!   REAL(r8) :: dzlak(1:nlevlake)    !lake dz (thickness)
!   REAL(r8) :: zsoi(1:nlevsoil)     !soil z  (layers)
!   REAL(r8) :: dzsoi(1:nlevsoil)    !soil dz (thickness)
!   REAL(r8) :: zisoi(0:nlevsoil)    !soil zi (interfaces)
!
!
!   REAL(r8) :: sand(19)
!   REAL(r8) :: clay(19)
!   INTEGER :: i
!
!   DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
!        10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
!
!   DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
!        33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./
!
!   REAL(r8) :: watsat(1,nlevsoil)      ! volumetric soil water at saturation (porosity)
!   REAL(r8) :: tksatu(1,nlevsoil)      ! thermal conductivity, saturated soil [W/m-K]
!   REAL(r8) :: tkmg(1,nlevsoil)        ! thermal conductivity, soil minerals  [W/m-K]
!   REAL(r8) :: tkdry(1,nlevsoil)       ! thermal conductivity, dry soil (W/m/Kelvin)
!   REAL(r8) :: csol(1,nlevsoil)        ! heat capacity, soil solids (J/m**3/Kelvin)
!
!   INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &
!        ims,ime, jms,jme, kms,kme,  &
!        its,ite, jts,jte, kts,kte
!   INTEGER , INTENT (IN) :: iswater
!   REAL(r8),     INTENT(IN)  :: xice_threshold
!   REAL(r8), DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   XICE
!   REAL(r8), DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   LAKEMASK
!
!   REAL(r8),           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: t_phy
!   REAL(r8),           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: p8w
!   REAL(r8),           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: dz8w
!   REAL(r8),           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: qvcurr
!   REAL(r8), OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: U_PHY
!   REAL(r8), OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)  :: V_PHY
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: glw
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: emiss
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: rainbl
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: swdown
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(INOUT)  :: albedo
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: XLAND
!   REAL(r8), OPTIONAL, DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: XLAT_URB2D
!   INTEGER,        DIMENSION( ims:ime, jms:jme )         ,INTENT(INOUT)  :: IVGTYP
!   REAL   :: dtbl
!
!   REAL(r8),           DIMENSION( ims:ime,1:nlevlake,jms:jme ),INTENT(IN)  :: z_lake3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevlake,jms:jme ),INTENT(IN)  :: dz_lake3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)  :: watsat3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)  :: csol3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)  :: tkmg3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)  :: tkdry3d
!   REAL(r8),           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)  :: tksatu3d
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: lakedepth2d
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)  :: ht
!   REAL                                                  ,INTENT(IN)  :: lake_min_elev
!
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: HFX
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: LH
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: GRDFLX
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: TSK
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: QFX
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: T2
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: TH2
!   REAL(r8),           DIMENSION( ims:ime, jms:jme )         ,INTENT(OUT) :: Q2
!
!   REAL(r8),           DIMENSION(ims:ime,jms:jme )                ,INTENT(inout)  :: savedtke12d
!   REAL(r8),           DIMENSION(ims:ime,jms:jme )                ,INTENT(inout)  :: snowdp2d,       &
!        h2osno2d,       &
!        snl2d,          &
!        t_grnd2d
!
!   REAL(r8),    DIMENSION( ims:ime,1:nlevlake, jms:jme )           ,INTENT(inout)  :: t_lake3d,       &
!        lake_icefrac3d
!   REAL(r8),    DIMENSION( ims:ime,-nlevsnow+1:nlevsoil, jms:jme )  ,INTENT(inout)  :: t_soisno3d,     &
!        h2osoi_ice3d,   &
!        h2osoi_liq3d,   &
!        h2osoi_vol3d,   &
!        z3d,            &
!        dz3d
!   REAL(r8),    DIMENSION( ims:ime,-nlevsnow+0:nlevsoil, jms:jme )  ,INTENT(inout)  :: zi3d
!
!
!
!   REAL(r8)     :: SFCTMP,PBOT,PSFC,ZLVL,Q2K,EMISSI,LWDN,PRCP,SOLDN,SOLNET
!   INTEGER  :: C,i,j,k
!
!   REAL(r8)  :: forc_t(1)
!   REAL(r8)  :: forc_pbot(1)
!   REAL(r8)  :: forc_psrf(1)
!   REAL(r8)  :: forc_hgt(1)
!   REAL(r8)  :: forc_hgt_q(1)
!   REAL(r8)  :: forc_hgt_t(1)
!   REAL(r8)  :: forc_hgt_u(1)
!   REAL(r8)  :: forc_q(1)
!   REAL(r8)  :: forc_u(1)
!   REAL(r8)  :: forc_v(1)
!   REAL(r8)  :: forc_lwrad(1)
!   REAL(r8)  :: prec(1)
!   REAL(r8)  :: sabg(1)
!   REAL(r8)  :: lat(1)
!   REAL(r8)  :: z_lake(1,nlevlake)
!   REAL(r8)  :: dz_lake(1,nlevlake)
!
!   REAL(r8)  :: lakedepth(1)
!   LOGICAL   :: do_capsnow(1)
!
!   REAL(r8)  :: h2osoi_vol(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: t_grnd(1)
!   REAL(r8)  :: h2osno(1)
!   REAL(r8)  :: snowdp(1)
!   REAL(r8)  :: z(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: dz(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: t_soisno(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: t_lake(1,nlevlake)
!   INTEGER   :: snl(1)
!   REAL(r8)  :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: savedtke1(1)
!   REAL(r8)  :: zi(1,-nlevsnow+1:nlevsoil)
!   REAL(r8)  :: lake_icefrac(1,nlevlake)
!
!   REAL(r8)  :: eflx_gnet(1)
!   REAL(r8)  :: eflx_lwrad_net(1)
!   REAL(r8)  :: eflx_sh_tot(1)
!   REAL(r8)  :: eflx_lh_tot(1)
!   REAL(r8)  :: t_ref2m(1)
!   REAL(r8)  :: q_ref2m(1)
!   REAL(r8)  :: taux(1)
!   REAL(r8)  :: tauy(1)
!   REAL(r8)  :: ram1(1)
!
!   REAL(r8)  :: z0mg(1)
!
!   dtbl = 60
!   dtime = dtbl
!
!   DO J = jts,jte
!      DO I = its,ite
!
!         SFCTMP  = t_phy(i,1,j)
!         PBOT    = p8w(i,2,j)
!         PSFC    = P8w(i,1,j)
!         ZLVL    = 0.5 * dz8w(i,1,j)
!         Q2K     = qvcurr(i,1,j)/(1.0 + qvcurr(i,1,j))
!         EMISSI  = EMISS(I,J)
!         LWDN    = GLW(I,J)*EMISSI
!         PRCP    = RAINBL(i,j)/dtbl
!         SOLDN   = SWDOWN(I,J)                        ! SOLDN is total incoming solar
!         SOLNET  = SOLDN*(1.-ALBEDO(I,J))             ! use mid-day albedo to determine net downward solar
!         ! (no solar zenith angle correction)
!         !        IF (XLAND(I,J).GT.1.5) THEN
!
!         !  if ( xice(i,j).gt.xice_threshold) then
!         !   ivgtyp(i,j) = iswater
!         !   xland(i,j) = 2.
!         !   lake_icefrac3d(i,1,j) = xice(i,j)
!         !   endif
!
!         !#if (EM_CORE==1)
!         IF (lakemask(i,j).EQ.1) THEN
!            !#else
!            !        if (ivgtyp(i,j)==iswater.and.ht(i,j)>= lake_min_elev ) THEN
!            !#endif
!
!            DO c = 1,column
!
!               forc_t(c)          = SFCTMP           ! [K]
!               forc_pbot(c)       = PBOT
!               forc_psrf(c)       = PSFC
!               forc_hgt(c)        = ZLVL             ! [m]
!               forc_hgt_q(c)      = ZLVL             ! [m]
!               forc_hgt_t(c)      = ZLVL             ! [m]
!               forc_hgt_u(c)      = ZLVL             ! [m]
!               forc_q(c)          = Q2K              ! [kg/kg]
!               forc_u(c)          = U_PHY(I,1,J)
!               forc_v(c)          = V_PHY(I,1,J)
!               ! forc_rho(c)        = SFCPRS / (287.04 * SFCTMP * (1.0+ 0.61 * Q2K)) ![kg/m/m/m]
!               forc_lwrad(c)      = LWDN             ! [W/m/m]
!               prec(c)            = PRCP             ! [mm/s]
!               sabg(c)            = SOLNET
!               lat(c)             = XLAT_URB2D(I,J)*pie/180  ! [radian]
!               do_capsnow(c)      = .FALSE.
!
!               lakedepth(c)           = lakedepth2d(i,j)
!               savedtke1(c)           = savedtke12d(i,j)
!               snowdp(c)              = snowdp2d(i,j)
!               h2osno(c)              = h2osno2d(i,j)
!               snl(c)                 = snl2d(i,j)
!               t_grnd(c)              = t_grnd2d(i,j)
!               DO k = 1,nlevlake
!                  t_lake(c,k)        = t_lake3d(i,k,j)
!                  lake_icefrac(c,k)  = lake_icefrac3d(i,k,j)
!                  z_lake(c,k)        = z_lake3d(i,k,j)
!                  dz_lake(c,k)       = dz_lake3d(i,k,j)
!               ENDDO
!               DO k = -nlevsnow+1,nlevsoil
!                  t_soisno(c,k)      = t_soisno3d(i,k,j)
!                  h2osoi_ice(c,k)    = h2osoi_ice3d(i,k,j)
!                  h2osoi_liq(c,k)    = h2osoi_liq3d(i,k,j)
!                  h2osoi_vol(c,k)    = h2osoi_vol3d(i,k,j)
!                  z(c,k)             = z3d(i,k,j)
!                  dz(c,k)            = dz3d(i,k,j)
!               ENDDO
!               DO k = -nlevsnow+0,nlevsoil
!                  zi(c,k)            = zi3d(i,k,j)
!               ENDDO
!               DO k = 1,nlevsoil
!                  watsat(c,k)        = watsat3d(i,k,j)
!                  csol(c,k)          = csol3d(i,k,j)
!                  tkmg(c,k)          = tkmg3d(i,k,j)
!                  tkdry(c,k)         = tkdry3d(i,k,j)
!                  tksatu(c,k)        = tksatu3d(i,k,j)
!               ENDDO
!
!            ENDDO
!
!            CALL LakeMain(watsat,tksatu,tkmg,tkdry,csol, &
!              forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q,  & !I
!                 forc_hgt_t,forc_hgt_u,forc_q, forc_u,         &
!                 forc_v,forc_lwrad,prec, sabg,lat,             &
!                 z_lake,dz_lake,lakedepth,do_capsnow,          &
!                 h2osno,snowdp,snl,z,dz,zi,                    & !H
!                 h2osoi_vol,h2osoi_liq,h2osoi_ice,             &
!                 t_grnd,t_soisno,t_lake,                       &
!                 savedtke1,lake_icefrac,                       &
!                 eflx_lwrad_net,eflx_gnet,                     & !O
!                 eflx_sh_tot,eflx_lh_tot,                      &
!                 t_ref2m,q_ref2m,                              &
!                 taux,tauy,ram1,z0mg)
!
!
!            DO c = 1,column
!               HFX(I,J)          = eflx_sh_tot(c)            ![W/m/m]
!               LH(I,J)           = eflx_lh_tot(c)            !W/m/m]
!               GRDFLX(I,J)       = eflx_gnet(c)              !W/m/m]
!               TSK(I,J)          = t_grnd(c)                 ![K]
!               T2(I,J)           = t_ref2m(c)
!               TH2(I,J)          = T2(I,J)*(1.E5/PSFC)**RCP
!               Q2(I,J)           = q_ref2m(c)
!               albedo(i,j)       = ( 0.6 * lake_icefrac(c,1) ) + ( (1.0-lake_icefrac(c,1)) * 0.08)
!
!               IF( tsk(i,j) >= tfrz ) THEN
!                  qfx(i,j)      = eflx_lh_tot(c)/hvap
!               ELSE
!                  qfx(i,j)      = eflx_lh_tot(c)/hsub       ! heat flux (W/m^2)=>mass flux(kg/(sm^2))
!               ENDIF
!            ENDDO
!
!            ! Renew Lake State Varialbes:(14)
!            DO c = 1,column
!
!               savedtke12d(i,j)         = savedtke1(c)
!               snowdp2d(i,j)            = snowdp(c)
!               h2osno2d(i,j)            = h2osno(c)
!               snl2d(i,j)               = snl(c)
!               t_grnd2d(i,j)            = t_grnd(c)
!               DO k = = 1,nlevlake
!                  t_lake3d(i,k,j)       = t_lake(c,k)
!                  lake_icefrac3d(i,k,j) = lake_icefrac(c,k)
!               ENDDO
!               DO k = -nlevsnow+1,nlevsoil
!                  z3d(i,k,j)            = z(c,k)
!                  dz3d(i,k,j)           = dz(c,k)
!                  t_soisno3d(i,k,j)     = t_soisno(c,k)
!                  h2osoi_liq3d(i,k,j)   = h2osoi_liq(c,k)
!                  h2osoi_ice3d(i,k,j)   = h2osoi_ice(c,k)
!                  h2osoi_vol3d(i,k,j)   = h2osoi_vol(c,k)
!               ENDDO
!               DO k = -nlevsnow+0,nlevsoil
!                  zi3d(i,k,j)           = zi(c,k)
!               ENDDO
!
!            ENDDO
!
!         ENDIF
!         !        ENDIF    ! if xland = 2
!      ENDDO
!   ENDDO
!
! END SUBROUTINE Lake

SUBROUTINE LakeMain(watsat,tksatu,tkmg,tkdry,csol,& !I by wfs
     forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q,& !I
     forc_hgt_t,forc_hgt_u,forc_q, forc_u,         &
     forc_v,forc_lwrad,prec, sabg,lat,             &
     z_lake,dz_lake,lakedepth,do_capsnow,          &
     h2osno,snowdp,snl,z,dz,zi,                    & !H
     h2osoi_vol,h2osoi_liq,h2osoi_ice,             &
     t_grnd,t_soisno,t_lake,                       &
     savedtke1,lake_icefrac,                       &
     kme,eflx_lwrad_net,eflx_gnet,          & !O by wfs
     eflx_sh_tot,eflx_lh_tot,                      &
     t_ref2m,q_ref2m,                              &
     taux,tauy,ram1,z0mg,z0hg,z0qg)                   !by wfs
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  ! REAL(r8),INTENT(in) :: times(1,nlevlake)  !times for kme by wfs
  ! REAL(r8),INTENT(in) :: infv(1,nlevlake)   !volume of inflow by wfs： m/s
  ! REAL(r8),INTENT(in) :: inft(1,nlevlake)   !temperature of inflow by wfs： K
  REAL(r8),INTENT(in) :: forc_t(1)
  REAL(r8),INTENT(in) :: forc_pbot(1)
  REAL(r8),INTENT(in) :: forc_psrf(1)
  REAL(r8),INTENT(in) :: forc_hgt(1)
  REAL(r8),INTENT(in) :: forc_hgt_q(1)
  REAL(r8),INTENT(in) :: forc_hgt_t(1)
  REAL(r8),INTENT(in) :: forc_hgt_u(1)
  REAL(r8),INTENT(in) :: forc_q(1)
  REAL(r8),INTENT(in) :: forc_u(1)
  REAL(r8),INTENT(in) :: forc_v(1)
  REAL(r8),INTENT(in) :: forc_lwrad(1)
  REAL(r8),INTENT(in) :: prec(1)
  REAL(r8),INTENT(in) :: sabg(1)
  REAL(r8),INTENT(in) :: lat(1)
  REAL(r8),INTENT(in) :: z_lake(1,nlevlake)
  REAL(r8),INTENT(in) :: dz_lake(1,nlevlake)

  REAL(r8), INTENT(in) :: lakedepth(1)

  LOGICAL , INTENT(in) :: do_capsnow(1)



  REAL(r8),INTENT(inout) :: h2osoi_vol(1,-nlevsnow+1:nlevsoil) ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
  REAL(r8),INTENT(inout) :: t_grnd(1)          ! ground temperature (Kelvin)
  REAL(r8),INTENT(inout) :: h2osno(1)          ! snow water (mm H2O)
  REAL(r8),INTENT(inout) :: snowdp(1)          ! snow height (m)
  REAL(r8),INTENT(inout) :: z(1,-nlevsnow+1:nlevsoil)             ! layer depth for snow & soil (m)
  REAL(r8),INTENT(inout) :: dz(1,-nlevsnow+1:nlevsoil)            ! layer thickness for soil or snow (m)
  REAL(r8),INTENT(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)      ! soil (or snow) temperature (Kelvin)
  REAL(r8),INTENT(inout) :: t_lake(1,nlevlake)                   ! lake temperature (Kelvin)
  INTEGER ,INTENT(inout) :: snl(1)                              ! number of snow layers
  REAL(r8),INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)    ! liquid water (kg/m2)
  REAL(r8),INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)    ! ice lens (kg/m2)
  REAL(r8),INTENT(inout) :: savedtke1(1)       ! top level eddy conductivity from previous timestep (W/m.K)
  REAL(r8),INTENT(inout) :: zi(1,-nlevsnow+0:nlevsoil)            ! interface level below a "z" level (m)
  REAL(r8),INTENT(inout) :: lake_icefrac(1,nlevlake)  ! mass fraction of lake layer that is frozen

  ! REAL(r8),INTENT(out) :: tkix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)! interface conductivity by wfs
  ! REAL(r8),INTENT(out) :: extraE(1,nlevlake)   !extraE from inflow by wfs： W
  ! REAL(r8),INTENT(out) :: phi(1,nlevlake) !heat source from solar by wfs： W
  REAL(r8),INTENT(out) :: kme(lbc:ubc,nlevlake)!kme as output by wfs
  ! REAL(r8),INTENT(out) :: fin(lbc:ubc)  ! fin as output by wfs
  ! REAL(r8),INTENT(out) :: ocvts(lbc:ubc) ! ocvts as output by wfs
  ! REAL(r8),INTENT(out) :: ncvts(lbc:ubc) ! ncvts as output by wfs
  REAL(r8),INTENT(out) :: eflx_gnet(1)
  REAL(r8),INTENT(out) :: eflx_lwrad_net(1)
  REAL(r8),INTENT(out) :: eflx_sh_tot(1)
  REAL(r8),INTENT(out) :: eflx_lh_tot(1)
  REAL(r8),INTENT(out) :: t_ref2m(1)
  REAL(r8),INTENT(out) :: q_ref2m(1)
  REAL(r8),INTENT(out) :: taux(1)
  REAL(r8),INTENT(out) :: tauy(1)
  REAL(r8),INTENT(out) :: ram1(1)


  REAL(r8),INTENT(out) :: z0mg(1)
  REAL(r8),INTENT(out) :: z0hg(lbp:ubp) ! by wfs
  REAL(r8),INTENT(out) :: z0qg(lbp:ubp) ! by wfs



  ! REAL(r8) :: begwb(1)=0.
  REAL(r8) :: t_veg(1)
  REAL(r8) :: eflx_soil_grnd(1)
  REAL(r8) :: eflx_lh_grnd(1)
  REAL(r8) :: eflx_sh_grnd(1)
  REAL(r8) :: eflx_lwrad_out(1)
  REAL(r8) :: qflx_evap_tot(1)
  REAL(r8) :: qflx_evap_soi(1)
  REAL(r8) :: qflx_prec_grnd(1)
  REAL(r8) :: forc_snow(1)
  REAL(r8) :: forc_rain(1)
  REAL(r8) :: ws(1)
  REAL(r8) :: ks(1)
  REAL(r8) :: qflx_snomelt(1)
  INTEGER  :: imelt(1,-nlevsnow+1:nlevsoil)
  REAL(r8) :: endwb(1)
  REAL(r8) :: snowage(1)
  REAL(r8) :: snowice(1)
  REAL(r8) :: snowliq(1)
  REAL(r8) :: t_snow(1)
  REAL(r8) :: qflx_drain(1)
  REAL(r8) :: qflx_surf(1)
  REAL(r8) :: qflx_infl(1)
  REAL(r8) :: qflx_qrgwl(1)
  REAL(r8) :: qcharge(1)
  REAL(r8) :: qflx_snowcap(1)
  REAL(r8) :: qflx_snowcap_col(1)
  REAL(r8) :: qflx_snow_grnd_pft(1)
  REAL(r8) :: qflx_snow_grnd_col(1)
  REAL(r8) :: qflx_rain_grnd(1)
  REAL(r8) :: frac_iceold(1,-nlevsnow+1:nlevsoil)
  REAL(r8) :: qflx_evap_tot_col(1)
  REAL(r8) :: soilalpha(1)
  REAL(r8) :: zwt(1)
  REAL(r8) :: fcov(1)
  REAL(r8) :: rootr_column(1,1:nlevsoil)
  REAL(r8) :: qflx_evap_grnd(1)
  REAL(r8) :: qflx_sub_snow(1)
  REAL(r8) :: qflx_dew_snow(1)
  REAL(r8) :: qflx_dew_grnd(1)
  REAL(r8) :: qflx_rain_grnd_col(1)



  !print*,'good4'
  !    print *,tfrz
  !    lat  = lat*pie/180  ! [radian]

  IF (prec(1)> 0.) THEN
     IF ( forc_t(1) > (tfrz + tcrit)) THEN
        forc_rain(1) = prec(1)
        forc_snow(1) = 0.
        !   flfall(1) = 1.
     ELSE
        forc_rain(1) = 0.
        forc_snow(1) = prec(1)

        !  if ( forc_t(1) <= tfrz) then
        !      flfall(1) = 0.
        !  else if ( forc_t(1) <= tfrz+2.) then
        !      flfall(1) = -54.632 + 0.2 *  forc_t(1)
        !  else
        !      flfall(1) = 0.4
     ENDIF
  ELSE
     forc_rain(1) = 0.
     forc_snow(1) = 0.
     !  flfall(1) = 1.
  ENDIF




  CALL ShalLakeFluxes(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q,   &  !i
       forc_hgt_t,forc_hgt_u,forc_q,                   &
       forc_u,forc_v,forc_lwrad,forc_snow,             &
       forc_rain,t_grnd,h2osno,snowdp,sabg,lat,        &
       dz,dz_lake,t_soisno,t_lake,snl,h2osoi_liq,      &
       h2osoi_ice,savedtke1,lakedepth,                 &  !by wfs
       qflx_prec_grnd,qflx_evap_soi,qflx_evap_tot,     &  !o
       eflx_sh_grnd,eflx_lwrad_out,eflx_lwrad_net,     &
       eflx_soil_grnd,eflx_sh_tot,eflx_lh_tot,         &
       eflx_lh_grnd,t_veg,t_ref2m,q_ref2m,taux,tauy,   &
       ram1,ws,ks,eflx_gnet,z0mg,z0hg,z0qg)             ! by wfs
  ! PRINT*,''
  ! PRINT*,'after ShalLakeFluxes:'
  ! PRINT *,'eflx_gnet ',eflx_gnet
  ! PRINT *,'eflx_sh_tot',eflx_sh_tot
  ! PRINT*,''
    ! print *,'t_grnd3=',t_grnd
  CALL ShalLakeTemperature(watsat,tksatu,tkmg,tkdry,csol,&
       t_grnd,h2osno,sabg,dz,dz_lake,z,zi,& !i
       z_lake,ws,ks,snl,eflx_gnet,lakedepth,       &
       lake_icefrac,snowdp,                        & !i&o
       kme,eflx_sh_grnd,eflx_sh_tot,eflx_soil_grnd,    & !o
       t_lake,t_soisno,h2osoi_liq,                 &
       h2osoi_ice,savedtke1,                       &
       frac_iceold,qflx_snomelt,imelt)
  ! PRINT*,''
  ! PRINT*,'after ShalLakeTemperature:'
  ! PRINT *,'eflx_gnet ',eflx_gnet
  ! PRINT *,'eflx_sh_tot',eflx_sh_tot
  ! PRINT*,''
    !  print *,'t_grnd4=',t_grnd
  ! begwb removed: this varialbe was used without initialization, TS, 20170816
  !  begwb,qflx_evap_tot,forc_t,do_capsnow,            &
  CALL ShalLakeHydrology(dz_lake,forc_rain,forc_snow,                          & !i
       qflx_evap_tot,forc_t,do_capsnow,            &
       t_grnd,qflx_evap_soi,                             &
       qflx_snomelt,imelt,frac_iceold,                   & !i add by guhp
       z,dz,zi,snl,h2osno,snowdp,lake_icefrac,t_lake,      & !i&o
       endwb,snowage,snowice,snowliq,t_snow,             & !o
       t_soisno,h2osoi_ice,h2osoi_liq,h2osoi_vol,        &
       qflx_drain,qflx_surf,qflx_infl,qflx_qrgwl,        &
       qcharge,qflx_prec_grnd,qflx_snowcap,              &
       qflx_snowcap_col,qflx_snow_grnd_pft,              &
       qflx_snow_grnd_col,qflx_rain_grnd,                &
       qflx_evap_tot_col,soilalpha,zwt,fcov,             &
       rootr_column,qflx_evap_grnd,qflx_sub_snow,        &
       qflx_dew_snow,qflx_dew_grnd,qflx_rain_grnd_col)
  !print*,'good8'
      !  print *,'t_grnd5=',t_grnd


  !==================================================================================
  ! !DESCRIPTION:
  ! Calculation of Shallow Lake Hydrology. Full hydrology of snow layers is
  ! done. However, there is no infiltration, and the water budget is balanced with

END SUBROUTINE LakeMain


SUBROUTINE ShalLakeFluxes(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q,           &  !i
     forc_hgt_t,forc_hgt_u,forc_q,                   &
     forc_u,forc_v,forc_lwrad,forc_snow,             &
     forc_rain,t_grnd,h2osno,snowdp,sabg,lat,        &
     dz,dz_lake,t_soisno,t_lake,snl,h2osoi_liq,      &
     h2osoi_ice,savedtke1,lakedepth,                 &  !by wfs
     qflx_prec_grnd,qflx_evap_soi,qflx_evap_tot,     &  !o
     eflx_sh_grnd,eflx_lwrad_out,eflx_lwrad_net,     &
     eflx_soil_grnd,eflx_sh_tot,eflx_lh_tot,         &
     eflx_lh_grnd,t_veg,t_ref2m,q_ref2m,taux,tauy,   &
     ram1,ws,ks,eflx_gnet,z0mg,z0hg,z0qg)               ! by wfs
  !==============================================================================
  ! DESCRIPTION:
  ! Calculates lake temperatures and surface fluxes for shallow lakes.
  !
  ! Shallow lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
  ! and soil layers with active temperature and gas diffusion below.
  !
  ! WARNING: This subroutine assumes lake columns have one and only one pft.
  !
  ! REVISION HISTORY:
  ! Created by Zack Subin, 2009
  ! Reedited by Hongping Gu, 2010
  !==============================================================================

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  ! INTEGER :: i
  !
  ! DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
  !      10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
  !
  ! DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
  !      33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)


  REAL(r8),INTENT(in) :: forc_t(1)
  REAL(r8),INTENT(in) :: forc_pbot(1)
  REAL(r8),INTENT(in) :: forc_psrf(1)
  REAL(r8),INTENT(in) :: forc_hgt(1)
  REAL(r8),INTENT(in) :: forc_hgt_q(1)
  REAL(r8),INTENT(in) :: forc_hgt_t(1)
  REAL(r8),INTENT(in) :: forc_hgt_u(1)
  REAL(r8),INTENT(in) :: forc_q(1)
  REAL(r8),INTENT(in) :: forc_u(1)
  REAL(r8),INTENT(in) :: forc_v(1)
  REAL(r8),INTENT(in) :: forc_lwrad(1)
  ! real(r8),intent(in) :: forc_rho(1)
  REAL(r8),INTENT(in) :: forc_snow(1)
  REAL(r8),INTENT(in) :: forc_rain(1)
  REAL(r8),INTENT(in) :: h2osno(1)
  REAL(r8),INTENT(in) :: snowdp(1)
  REAL(r8),INTENT(in) :: sabg(1)
  REAL(r8),INTENT(in) :: lat(1)
  REAL(r8),INTENT(in) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8),INTENT(in) :: dz_lake(1,nlevlake)
  REAL(r8),INTENT(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8),INTENT(in) :: t_lake(1,nlevlake)
  INTEGER ,INTENT(in) :: snl(1)
  REAL(r8),INTENT(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
  REAL(r8),INTENT(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8),INTENT(in) :: savedtke1(1)
  REAL(r8),INTENT(in) :: lakedepth(1)


  REAL(r8),INTENT(inout) :: t_grnd(1)

  REAL(r8),INTENT(out):: qflx_prec_grnd(1)
  REAL(r8),INTENT(out):: qflx_evap_soi(1)
  REAL(r8),INTENT(out):: qflx_evap_tot(1)
  REAL(r8),INTENT(out):: eflx_sh_grnd(1)
  REAL(r8),INTENT(out):: eflx_lwrad_out(1)
  REAL(r8),INTENT(out):: eflx_lwrad_net(1)
  REAL(r8),INTENT(out):: eflx_soil_grnd(1)
  REAL(r8),INTENT(out):: eflx_sh_tot(1)
  REAL(r8),INTENT(out):: eflx_lh_tot(1)
  REAL(r8),INTENT(out):: eflx_lh_grnd(1)
  REAL(r8),INTENT(out):: t_veg(1)
  REAL(r8),INTENT(out):: t_ref2m(1)
  REAL(r8),INTENT(out):: q_ref2m(1)
  REAL(r8),INTENT(out):: taux(1)
  REAL(r8),INTENT(out):: tauy(1)
  REAL(r8),INTENT(out):: ram1(1)
  REAL(r8),INTENT(out):: ws(1)
  REAL(r8),INTENT(out):: ks(1)

  REAL(r8),INTENT(out):: eflx_gnet(1)

  REAL(r8),INTENT(out):: z0mg(1)
  REAL(r8),INTENT(out):: z0hg(lbp:ubp) ! by wfs
  REAL(r8),INTENT(out):: z0qg(lbp:ubp) ! by wfs



  INTEGER , PARAMETER :: islak  = 2
  INTEGER , PARAMETER :: niters = 4
  REAL(r8), PARAMETER :: beta1  = 1._r8
  REAL(r8), PARAMETER :: emg    = 0.97_r8
  REAL(r8), PARAMETER :: zii    = 1000._r8
  REAL(r8), PARAMETER :: tdmax  = 277._r8
  REAL(r8) :: forc_th(1)
  REAL(r8) :: forc_vp(1)
  REAL(r8) :: forc_rho(1)
  INTEGER  :: i,fc,fp,g,c,p
  INTEGER  :: fncopy
  INTEGER  :: fnold
  INTEGER  :: fpcopy(num_shlakep)
  INTEGER  :: iter
  INTEGER  :: nmozsgn(lbp:ubp)
  INTEGER  :: jtop(lbc:ubc)
  REAL(r8) :: ax
  REAL(r8) :: bx
  REAL(r8) :: degdT
  REAL(r8) :: dqh(lbp:ubp)
  REAL(r8) :: dth(lbp:ubp)
  REAL(r8) :: dthv
  REAL(r8) :: dzsur(lbc:ubc)
  REAL(r8) :: eg
  REAL(r8) :: htvp(lbc:ubc)
  REAL(r8) :: obu(lbp:ubp)
  REAL(r8) :: obuold(lbp:ubp)
  REAL(r8) :: qsatg(lbc:ubc)
  REAL(r8) :: qsatgdT(lbc:ubc)
  REAL(r8) :: qstar
  REAL(r8) :: ram(lbp:ubp)
  REAL(r8) :: rah(lbp:ubp)
  REAL(r8) :: raw(lbp:ubp)
  REAL(r8) :: stftg3(lbp:ubp)
  REAL(r8) :: temp1(lbp:ubp)
  REAL(r8) :: temp12m(lbp:ubp)
  REAL(r8) :: temp2(lbp:ubp)
  REAL(r8) :: temp22m(lbp:ubp)
  REAL(r8) :: tgbef(lbc:ubc)
  REAL(r8) :: thm(lbc:ubc)
  REAL(r8) :: thv(lbc:ubc)
  REAL(r8) :: thvstar
  REAL(r8) :: tksur
  REAL(r8) :: tsur
  REAL(r8) :: tstar
  REAL(r8) :: um(lbp:ubp)
  REAL(r8) :: ur(lbp:ubp)
  REAL(r8) :: ustar(lbp:ubp)
  REAL(r8) :: wc
  REAL(r8) :: zeta
  REAL(r8) :: zldis(lbp:ubp)
  REAL(r8) :: displa(lbp:ubp)
  ! REAL(r8) :: z0hg(lbp:ubp) !by wfs
  ! REAL(r8) :: z0qg(lbp:ubp) !by wfs
  REAL(r8) :: beta(2)
  REAL(r8) :: u2m
  REAL(r8) :: u10(1)
  REAL(r8) :: fv(1)

  REAL(r8) :: fm(lbp:ubp)
  REAL(r8) :: bw
  REAL(r8) :: t_grnd_temp
  REAL(r8) :: betaprime(lbc:ubc)
  REAL(r8) :: miu !kinematic viscosity of air by wfs
  REAL(r8) :: R0  !near-surface atmospheric roughness Reynolds number by wfs
  REAL(r8) :: cc  !Charnock coefficient by wfs
  REAL(r8) :: fl  !fetch-limitation for Charnock coefficient by wfs
  REAL(r8) :: dl  !depth-limitation for Charnock coefficient by wfs

  CHARACTER*256 :: message


  DATA beta/0.4_r8, 0.4_r8/

  ! (deep lake, shallow lake)
  ! This is the energy absorbed at the lake surface if no snow.
  !    data za  /0.6_r8, 0.5_r8/
  !    data eta /0.1_r8, 0.5_r8/
  !-----------------------------------------------------------------------


  !    dtime = get_step_size()

  ! Begin calculations

  !dir$ concurrent
  !cdir nodep
  forc_th(1)  = forc_t(1) * (forc_psrf(1)/ forc_pbot(1))**(rair/cpair)
  forc_vp(1)  = forc_q(1) * forc_pbot(1)/ (0.622 + 0.378 * forc_q(1))
  forc_rho(1) = (forc_pbot(1) - 0.378 * forc_vp(1)) / (rair * forc_t(1))

  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)
     g = cgridcell(c)

     ! Surface temperature and fluxes

     ! Find top layer
     !       if (snl(c) > 0 .or. snl(c) < -5) then
     !         WRITE(message,*)  'snl is not defined in ShalLakeFluxesMod'
     !         CALL wrf_message(message)
     !         CALL wrf_error_fatal("snl: out of range value")
     !       end if
     !       if (snl(c) /= 0) then
     !           write(6,*)'snl is not equal to zero in ShalLakeFluxesMod'
     !           call endrun()
     !       end if
     jtop(c) = snl(c) + 1


     IF (snl(c) < 0) THEN
        betaprime(c) = 1._r8  !Assume all solar rad. absorbed at the surface of the top snow layer.
        dzsur(c) = dz(c,jtop(c))/2._r8
     ELSE
        betaprime(c) = beta(islak)
        dzsur(c) = dz_lake(c,1)/2._r8
     END IF
     ! Originally this was 1*dz, but shouldn't it be 1/2?

     ! Saturated vapor pressure, specific humidity and their derivatives
     ! at lake surface

     CALL QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

     ! Potential, virtual potential temperature, and wind speed at the
     ! reference height

     thm(c) = forc_t(g) + 0.0098_r8*forc_hgt_t(g)   ! intermediate variable
     thv(c) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))     ! virtual potential T
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     g = pgridcell(p)

     nmozsgn(p) = 0
     obuold(p) = 0._r8
     displa(p) = 0._r8

     ! Roughness lengths


     ! changed by Hongping Gu
     !   if (t_grnd(c) >= tfrz) then   ! for unfrozen lake
     !      z0mg(p) = 0.01_r8
     !   else                          ! for frozen lake
     !   ! Is this okay even if it is snow covered?  What is the roughness over
     !   non-veg. snow?
     !      z0mg(p) = 0.04_r8
     !   end if

     IF (t_grnd(c) >= tfrz) THEN   ! for unfrozen lake
        z0mg(p) = 0.001_r8        !original 0.01, Gu's version is 0.001
     ELSE IF(snl(c) == 0 ) THEN                         ! for frozen lake
        ! Is this okay even if it is snow covered?  What is the roughness over
        ! non-veg. snow?
        z0mg(p) = 0.001_r8          !original 0.04, now for frozen lake without snow, Gu's version was 0.05
     ELSE                          ! for frozen lake with snow
        z0mg(p) = 0.0024_r8       ! Gu's version was 0.024
     END IF


     ! z0hg(p) = 0.0003_r8
     ! z0qg(p) = 0.0003_r8!by wfs

     ! Latent heat

     !#define PERGRO
     !#if (defined PERGRO)
     htvp(c) = hvap
     !#else
     !       if (t_grnd(c) > tfrz) then
     !          htvp(c) = hvap
     !       else
     !          htvp(c) = hsub
     !       end if
     !#endif
     ! Zack Subin, 3/26/09: Shouldn't this be the ground temperature rather than the air temperature above?
     ! I'll change it for now.

     ! Initialize stability variables

     ur(p)    = MAX(1.0_r8,SQRT(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
     dth(p)   = thm(c)-t_grnd(c)
     dqh(p)   = forc_q(g)-qsatg(c)
     dthv     = dth(p)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(p)
     zldis(p) = forc_hgt_u(g) - 0._r8

     ! Initialize Monin-Obukhov length and wind speed

     CALL MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg(p), um(p), obu(p))

  END DO

  iter = 1
  fncopy = num_shlakep
  fpcopy(1:num_shlakep) = filter_shlakep(1:num_shlakep)

  ! Begin stability iteration

  ITERATION : DO WHILE (iter <= niters .AND. fncopy > 0)

     ! Determine friction velocity, and potential temperature and humidity
     ! profiles of the surface boundary layer

    !following Subin's work by wfs
    IF (iter == 1)  ustar(p)=0.06_r8 !initial ustar
    miu = 1.51e-5_r8*(t_grnd(c)/293.15_r8)**1.5_r8*1.013e+5_r8/forc_psrf(1)
    R0 = z0mg(p)*ustar(p)/miu

    cc = 0.01_r8 + (0.11_r8 - 0.01_r8) * exp(-min(fl,dl))
    fl = (25._r8*lakedepth(p)*grav/ustar(p)**2._r8)**(1._r8/3._r8)/100._r8
    dl = 1._r8*sqrt(lakedepth(p)*grav)/ur(p)

    IF (t_grnd(c) < tfrz) THEN
      IF (snl(c) == 0) THEN
        z0mg(p) = 0.001_r8   !for frozen lake without snow
      ELSE
        z0mg(p) = 0.0024_r8  !for frozen lake with snow
      END IF
      !z0hg and z0qg for frozen lake
      z0hg(p) = z0mg(p)*exp(-0.13_r8*R0**0.45_r8)
      z0qg(p) = z0mg(p)*exp(-0.13_r8*R0**0.45_r8)
    ELSE IF (iter == 1) THEN
      z0mg(p) = 0.001_r8     !for unfrozen lake, first iteration
      z0hg(p) = z0mg(p)*exp(-vkc/0.71_r8*(4_r8*sqrt(R0)-3.2_r8))
      z0qg(p) = z0mg(p)*exp(-vkc/0.66_r8*(4_r8*sqrt(R0)-4.2_r8))
    ELSE                    !for unfrozen lake, 2~4 iteration
      z0mg(p) = max(0.1_r8*miu/ustar(p),cc*ustar(p)/grav)
      z0hg(p) = z0mg(p)*exp(-vkc/0.71_r8*(4_r8*sqrt(R0)-3.2_r8))
      z0qg(p) = z0mg(p)*exp(-vkc/0.66_r8*(4_r8*sqrt(R0)-4.2_r8))
    END IF

    IF (z0mg(p) < 0.00001_r8)  z0mg(p) = 0.00001_r8
    IF (z0hg(p) < 0.00001_r8)  z0hg(p) = 0.00001_r8
    IF (z0qg(p) < 0.00001_r8)  z0qg(p) = 0.00001_r8

    !end by wfs

     CALL FrictionVelocity(pgridcell,forc_hgt,forc_hgt_u,          & !i
          forc_hgt_t,forc_hgt_q,                  & !i
          lbp, ubp, fncopy, fpcopy,               & !i
          displa, z0mg, z0hg, z0qg,               & !i
          obu, iter, ur, um,                      & !i
          ustar,temp1, temp2, temp12m, temp22m,   & !o
          u10,fv,                                 & !o
          fm)  !i&o

     !dir$ concurrent
     !cdir nodep
     DO fp = 1, fncopy
        p = fpcopy(fp)
        c = pcolumn(p)
        g = pgridcell(p)

        tgbef(c) = t_grnd(c)
        IF (t_grnd(c) > tfrz .AND. t_lake(c,1) > tfrz .AND. snl(c) == 0) THEN
           tksur = savedtke1(c)
           ! Set this to the eddy conductivity from the last
           ! timestep, as the molecular conductivity will be orders of magnitude too small.
           ! Will have to deal with first timestep.
           tsur = t_lake(c,1)
        ELSE IF (snl(c) == 0) THEN  !frozen but no snow layers
           tksur = tkice
           tsur = t_lake(c,1)
        ELSE
           !Need to calculate thermal conductivity of the top snow layer
           bw = (h2osoi_ice(c,jtop(c))+h2osoi_liq(c,jtop(c)))/dz(c,jtop(c))
           tksur = tkairc + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkairc)
           tsur = t_soisno(c,jtop(c))
        END IF

        ! Determine aerodynamic resistances

        ram(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
        rah(p)  = 1._r8/(temp1(p)*ustar(p))
        raw(p)  = 1._r8/(temp2(p)*ustar(p))
        ram1(p) = ram(p)   !pass value to global variable

        ! Get derivative of fluxes with respect to ground temperature

        stftg3(p) = emg*sb*tgbef(c)*tgbef(c)*tgbef(c)

        ! Changed surface temperature from t_lake(c,1) to tsur.
        ! Also adjusted so that if there are snow layers present, all radiation is absorbed in the top layer.
        ax  = betaprime(c)*sabg(p) + emg*forc_lwrad(g) + 3._r8*stftg3(p)*tgbef(c) &
             + forc_rho(g)*cpair/rah(p)*thm(c) &
             - htvp(c)*forc_rho(g)/raw(p)*(qsatg(c)-qsatgdT(c)*tgbef(c) - forc_q(g)) &
             + tksur*tsur/dzsur(c)
        !Changed sabg(p) and to betaprime(c)*sabg(p).
        bx  = 4._r8*stftg3(p) + forc_rho(g)*cpair/rah(p) &
             + htvp(c)*forc_rho(g)/raw(p)*qsatgdT(c) + tksur/dzsur(c)

        t_grnd(c) = ax/bx
        ! print *,'t_grnd1=',t_grnd(c)

        ! Update htvp
        !#define PERGRO
        !#ifndef PERGRO
        !       if (t_grnd(c) > tfrz) then
        !          htvp(c) = hvap
        !       else
        !          htvp(c) = hsub
        !       end if
        !#endif

        ! Surface fluxes of momentum, sensible and latent heat
        ! using ground temperatures from previous time step

        eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
        qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c))-forc_q(g))/raw(p)

        ! Re-calculate saturated vapor pressure, specific humidity and their
        ! derivatives at lake surface

        CALL QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

        dth(p)=thm(c)-t_grnd(c)
        dqh(p)=forc_q(g)-qsatg(c)

        tstar = temp1(p)*dth(p)
        qstar = temp2(p)*dqh(p)

        thvstar=tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
        zeta=zldis(p)*vkc * grav*thvstar/(ustar(p)**2*thv(c))

        IF (zeta >= 0._r8) THEN     !stable
           zeta = MIN(2._r8,MAX(zeta,0.01_r8))
           um(p) = MAX(ur(p),0.1_r8)
        ELSE                     !unstable
           zeta = MAX(-100._r8,MIN(zeta,-0.01_r8))
           wc = beta1*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
           um(p) = SQRT(ur(p)*ur(p)+wc*wc)
        END IF
        obu(p) = zldis(p)/zeta

        IF (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1

        obuold(p) = obu(p)

     END DO   ! end of filtered pft loop

     iter = iter + 1
     IF (iter <= niters ) THEN
        ! Rebuild copy of pft filter for next pass through the ITERATION loop

        fnold = fncopy
        fncopy = 0
        DO fp = 1, fnold
           p = fpcopy(fp)
           IF (nmozsgn(p) < 3) THEN
              fncopy = fncopy + 1
              fpcopy(fncopy) = p
           END IF
        END DO   ! end of filtered pft loop
     END IF

  END DO ITERATION   ! end of stability iteration

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     g = pgridcell(p)

     ! If there is snow on the ground and t_grnd > tfrz: reset t_grnd = tfrz.
     ! Re-evaluate ground fluxes.
     ! h2osno > 0.5 prevents spurious fluxes.
     ! note that qsatg and qsatgdT should be f(tgbef) (PET: not sure what this
     ! comment means)
     ! Zack Subin, 3/27/09: Since they are now a function of whatever t_grnd was before cooling
     !    to freezing temperature, then this value should be used in the derivative correction term.
     ! Should this happen if the lake temperature is below freezing, too? I'll assume that for now.
     ! Also, allow convection if ground temp is colder than lake but warmer than 4C, or warmer than
     !    lake which is warmer than freezing but less than 4C.
     !#ifndef SHLAKETEST
     IF ( (h2osno(c) > 0.5_r8 .OR. t_lake(c,1) <= tfrz) .AND. t_grnd(c) > tfrz) THEN
        !#else
        !       if ( t_lake(c,1) <= tfrz .and. t_grnd(c) > tfrz) then
        !#endif
        t_grnd_temp = t_grnd(c)
        t_grnd(c) = tfrz
        eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
        qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
     ELSE IF ( (t_lake(c,1) > t_grnd(c) .AND. t_grnd(c) > tdmax) .OR. &
          (t_lake(c,1) < t_grnd(c) .AND. t_lake(c,1) > tfrz .AND. t_grnd(c) < tdmax) ) THEN
        ! Convective mixing will occur at surface
        t_grnd_temp = t_grnd(c)
        t_grnd(c) = t_lake(c,1)
        eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
        qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
     END IF
    !  print *,'t_grnd2=',t_grnd(c)

     ! Update htvp
     !#define PERGRO
     !#ifndef PERGRO
     !       if (t_grnd(c) > tfrz) then
     !          htvp(c) = hvap
     !       else
     !          htvp(c) = hsub
     !       end if
     !#endif

     ! Net longwave from ground to atmosphere

     !       eflx_lwrad_out(p) = (1._r8-emg)*forc_lwrad(g) + stftg3(p)*(-3._r8*tgbef(c)+4._r8*t_grnd(c))
     ! What is tgbef doing in this equation? Can't it be exact now? --Zack Subin, 4/14/09
     eflx_lwrad_out(p) = (1._r8-emg)*forc_lwrad(g) + emg*sb*t_grnd(c)**4

     ! Ground heat flux

     eflx_soil_grnd(p) = sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) - &
          eflx_sh_grnd(p) - htvp(c)*qflx_evap_soi(p)
     !Why is this sabg(p) and not beta*sabg(p)??
     !I've kept this as the incorrect sabg so that the energy balance check will be correct.
     !This is the effective energy flux into the ground including the lake solar absorption
     !below the surface.  The variable eflx_gnet will be used to pass the actual heat flux
     !from the ground interface into the lake.

     taux(p) = -forc_rho(g)*forc_u(g)/ram(p)
     tauy(p) = -forc_rho(g)*forc_v(g)/ram(p)

     eflx_sh_tot(p)   = eflx_sh_grnd(p)
     qflx_evap_tot(p) = qflx_evap_soi(p)
     eflx_lh_tot(p)   = htvp(c)*qflx_evap_soi(p)
     eflx_lh_grnd(p)  = htvp(c)*qflx_evap_soi(p)
     !#define LAKEDEBUG
     !#if (defined LAKEDEBUG)
     !       write(message,*) 'c, sensible heat = ', c, eflx_sh_tot(p), 'latent heat = ', eflx_lh_tot(p) &
     !              , 'ground temp = ', t_grnd(c), 'h2osno = ', h2osno(c)
     !       CALL wrf_message(message)
     !       if (abs(eflx_sh_tot(p)) > 1500 .or. abs(eflx_lh_tot(p)) > 1500) then
     !           write(message,*)'WARNING: SH, LH = ', eflx_sh_tot(p), eflx_lh_tot(p)
     !           CALL wrf_message(message)
     !       end if
     !       if (abs(eflx_sh_tot(p)) > 10000 .or. abs(eflx_lh_tot(p)) > 10000 &
     !             .or. abs(t_grnd(c)-288)>200 ) CALL wrf_error_fatal ( 't_grnd is out of range' )
     !#endif
     ! 2 m height air temperature
     t_ref2m(p) = thm(c) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

     ! 2 m height specific humidity
     q_ref2m(p) = forc_q(g) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

     ! Energy residual used for melting snow
     ! Effectively moved to ShalLakeTemp

     ! Prepare for lake layer temperature calculations below
     ! fin(c) = betaprime * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
     !          eflx_sh_tot(p) + eflx_lh_tot(p))
     ! NOW this is just the net ground heat flux calculated below.

     eflx_gnet(p) = betaprime(c) * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
          eflx_sh_tot(p) + eflx_lh_tot(p))
     !            print *,'eflx_gnet=',eflx_gnet(p),eflx_gnet
     ! This is the actual heat flux from the ground interface into the lake, not including
     ! the light that penetrates the surface.

     !       u2m = max(1.0_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))
     ! u2 often goes below 1 m/s; it seems like the only reason for this minimum is to
     ! keep it from being zero in the ks equation below; 0.1 m/s is a better limit for
     ! stable conditions --ZS
     u2m = MAX(0.1_r8,ustar(p)/vkc*LOG(2._r8/z0mg(p)))

     ws(c) = 1.2e-03_r8 * u2m
     ks(c) = 6.6_r8*SQRT(ABS(SIN(lat(g))))*(u2m**(-1.84_r8))

  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of surface flux relevant code in original BiogeophysicsLakeMod until history loop.

  ! The following are needed for global average on history tape.

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     g = pgridcell(p)
     !       t_veg(p) = forc_t(g)
     !This is an odd choice, since elsewhere t_veg = t_grnd for bare ground.
     !Zack Subin, 4/09
     t_veg(p) = t_grnd(c)
     eflx_lwrad_net(p)  = eflx_lwrad_out(p) - forc_lwrad(g)
     qflx_prec_grnd(p) = forc_rain(g) + forc_snow(g)
  END DO

END SUBROUTINE ShalLakeFluxes

SUBROUTINE ShalLakeTemperature(watsat,tksatu,tkmg,tkdry,csol,&
        t_grnd,h2osno,sabg,dz,dz_lake,z,zi, & !i by wfs
        z_lake,ws,ks,snl,eflx_gnet,lakedepth,       &
        lake_icefrac,snowdp,                        & !i&o
        kme,eflx_sh_grnd,eflx_sh_tot,eflx_soil_grnd,    & !o by wfs
        t_lake,t_soisno,h2osoi_liq,                 &
        h2osoi_ice,savedtke1,                       &
        frac_iceold,qflx_snomelt,imelt)
  !=======================================================================================================
  ! !DESCRIPTION:
  ! Calculates temperatures in the 20-25 layer column of (possible) snow,
  ! lake water, and soil beneath lake.
  ! Snow and soil temperatures are determined as in SoilTemperature, except
  ! for appropriate boundary conditions at the top of the snow (the flux is fixed
  ! to be the ground heat flux calculated in ShalLakeFluxes), the bottom of the snow
  ! (adjacent to top lake layer), and the top of the soil (adjacent to the bottom
  ! lake layer). Also, the soil is assumed to be always fully saturated (ShalLakeHydrology
  ! will have to insure this). The whole column is solved simultaneously as one tridiagonal matrix.
  ! Lake temperatures are determined from the Hostetler model as before, except now:
  !    i) Lake water layers can freeze by any fraction and release latent heat; thermal
  !       and mechanical properties are adjusted for ice fraction.
  !   ii) Convective mixing (though not eddy diffusion) still occurs for frozen lakes.
  !  iii) No sunlight is absorbed in the lake if there are snow layers.
  !   iv) Light is allowed to reach the top soil layer (where it is assumed to be completely absorbed).
  !    v) Lakes have variable depth, set ultimately in surface data set but now in initShalLakeMod.
  !
  ! Eddy + molecular diffusion:
  ! d ts    d            d ts     1 ds
  ! ---- = -- [(km + ke) ----] + -- --
  !  dt    dz             dz     cw dz
  !
  ! where: ts = temperature (kelvin)
  !         t = time (s)
  !         z = depth (m)
  !        km = molecular diffusion coefficient (m**2/s)
  !        ke = eddy diffusion coefficient (m**2/s)
  !        cw = heat capacity (j/m**3/kelvin)
  !         s = heat source term (w/m**2)
  !
  !   Shallow lakes are allowed to have variable depth, set in _____.
  !
  !   For shallow lakes:    ke > 0 if unfrozen,
  !       and convective mixing occurs WHETHER OR NOT frozen. (See e.g. Martynov...)
  !
  ! Use the Crank-Nicholson method to set up tridiagonal system of equations to
  ! solve for ts at time n+1, where the temperature equation for layer i is
  ! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
  !
  ! The solution conserves energy as:
  !
  ! [For lake layers]
  ! cw*([ts(      1)] n+1 - [ts(      1)] n)*dz(      1)/dt + ... +
  ! cw*([ts(10)] n+1 - [ts(10)] n)*dz(10)/dt = fin
  ! But now there is phase change, so cv is not constant and there is
  ! latent heat.
  !
  ! where:
  ! [ts] n   = old temperature (kelvin)
  ! [ts] n+1 = new temperature (kelvin)
  ! fin      = heat flux into lake (w/m**2)
  !          = betaprime*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot
  !          (This is now the same as the ground heat flux.)
  !            + phi(1) + ... + phi(10) + phi(top soil level)
  ! betaprime = beta(islak) for no snow layers, and 1 for snow layers.
  ! This assumes all radiation is absorbed in the top snow layer and will need
  ! to be changed for CLM 4.
  !
  ! WARNING: This subroutine assumes lake columns have one and only one pft.
  !
  ! Outline:
  ! 1!) Initialization
  ! 2!) Lake density
  ! 3!) Diffusivity
  ! 4!) Heat source term from solar radiation penetrating lake
  ! 5!) Set thermal props and find initial energy content
  ! 6!) Set up vectors for tridiagonal matrix solution
  ! 7!) Solve tridiagonal and back-substitute
  ! 8!) (Optional) Do first energy check using temperature change at constant heat capacity.
  ! 9!) Phase change
  ! 9.5!) (Optional) Do second energy check using temperature change and latent heat, considering changed heat capacity.
  !                  Also do soil water balance check.
  !10!) Convective mixing
  !11!) Do final energy check to detect small numerical errors (especially from convection)
  !     and dump small imbalance into sensible heat, or pass large errors to BalanceCheckMod for abort.
  !
  ! REVISION HISTORY:
  ! Created by Zack Subin, 2009.
  ! Reedited by Hongping Gu, 2010.
  !=========================================================================================================


  !    use TridiagonalMod     , only : Tridiagonal

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  ! INTEGER :: i
  !
  ! DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
  !      10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
  !
  ! DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
  !      33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  ! REAL(r8), INTENT(in) :: times(1,nlevlake)  !times for kme by wfs
  ! REAL(r8), INTENT(in) :: infv(1,nlevlake)   !volume of inflow by wfs： m/s
  ! REAL(r8), INTENT(in) :: inft(1,nlevlake)   !temperature of inflow by wfs： K
  REAL(r8), INTENT(in) :: t_grnd(1)
  REAL(r8), INTENT(inout) :: h2osno(1)
  REAL(r8), INTENT(in) :: sabg(1)
  REAL(r8), INTENT(in) :: dz(1,-nlevsnow + 1:nlevsoil)
  REAL(r8), INTENT(in) :: dz_lake(1,nlevlake)
  REAL(r8), INTENT(in) :: z(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: zi(1,-nlevsnow+0:nlevsoil)

  REAL(r8), INTENT(in) :: z_lake(1,nlevlake)
  REAL(r8), INTENT(in) :: ws(1)
  REAL(r8), INTENT(in) :: ks(1)

  INTEGER , INTENT(in) :: snl(1)
  REAL(r8), INTENT(inout) :: eflx_gnet(1)
  REAL(r8), INTENT(in) :: lakedepth(1)

  REAL(r8), INTENT(inout) :: snowdp(1)

  ! REAL(r8), INTENT(out) :: tkix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil) !interface conductivity by wfs
  ! REAL(r8), INTENT(out) :: extraE(1,nlevlake)   !extraE from inflow by wfs： W
  ! REAL(r8), INTENT(out) :: phi(1,nlevlake) !heat source from solar by wfs： W
  REAL(r8), INTENT(out) :: kme(lbc:ubc,nlevlake)!kme as output by wfs
  ! REAL(r8), INTENT(out) :: fin(lbc:ubc)  ! fin as output by wfs
  ! REAL(r8), INTENT(out) :: ocvts(lbc:ubc) ! ocvts as output by wfs
  ! REAL(r8), INTENT(out) :: ncvts(lbc:ubc) ! ncvts as output by wfs
  REAL(r8), INTENT(out) :: eflx_sh_grnd(1)
  REAL(r8), INTENT(out) :: eflx_sh_tot(1)
  REAL(r8), INTENT(out) :: eflx_soil_grnd(1)

  REAL(r8), INTENT(inout) :: t_lake(1,nlevlake)
  REAL(r8), INTENT(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: lake_icefrac(1,nlevlake)
  REAL(r8), INTENT(out) :: savedtke1(1)
  REAL(r8), INTENT(out) :: frac_iceold(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: qflx_snomelt(1)
  INTEGER, INTENT(out)  :: imelt(1,-nlevsnow+1:nlevsoil)



  INTEGER , PARAMETER  :: islak = 2
  REAL(r8), PARAMETER  :: p0 = 1._r8
  INTEGER  :: i,j,fc,fp,g,c,p
  REAL(r8) :: beta(2)
  REAL(r8) :: za(2)
  REAL(r8) :: eta(2)
  REAL(r8) :: cwat
  REAL(r8) :: cice_eff

  REAL(r8) :: cfus
  ! REAL(r8) :: times(1,nlevlake)  !times for kme by wfs
  REAL(r8) :: km
  REAL(r8) :: tkice_eff
  REAL(r8) :: a(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: b(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: c1(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: r(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: rhow(lbc:ubc,nlevlake)
  REAL(r8) :: phi(lbc:ubc,nlevlake)
  ! REAL(r8) :: kme(lbc:ubc,nlevlake) !by wfs
  REAL(r8) :: rsfin
  REAL(r8) :: rsfout
  REAL(r8) :: phi_soil(lbc:ubc)
  REAL(r8) :: ri
  REAL(r8) :: fin(lbc:ubc)
  REAL(r8) :: ocvts(lbc:ubc)
  REAL(r8) :: ncvts(lbc:ubc)
  REAL(r8) :: ke
  REAL(r8) :: zin
  REAL(r8) :: zout
  REAL(r8) :: drhodz
  REAL(r8) :: n2
  REAL(r8) :: num
  REAL(r8) :: den
  REAL(r8) :: tav_froz(lbc:ubc)
  REAL(r8) :: tav_unfr(lbc:ubc)
  REAL(r8) :: nav(lbc:ubc)
  REAL(r8) :: phidum
  REAL(r8) :: iceav(lbc:ubc)
  REAL(r8) :: qav(lbc:ubc)
  INTEGER  :: jtop(lbc:ubc)
  REAL(r8) :: cv (lbc:ubc,-nlevsnow+1:nlevsoil)
  REAL(r8) :: tk (lbc:ubc,-nlevsnow+1:nlevsoil)

  REAL(r8) :: cv_lake (lbc:ubc,1:nlevlake)
  REAL(r8) :: tk_lake (lbc:ubc,1:nlevlake)
  REAL(r8) :: cvx (lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: tkix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: tx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil) !
  REAL(r8) :: tktopsoillay(lbc:ubc)
  REAL(r8) :: fnx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: phix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: zx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: dzm
  REAL(r8) :: dzp
  INTEGER  :: jprime
  REAL(r8) :: factx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)
  REAL(r8) :: t_lake_bef(lbc:ubc,1:nlevlake)
  REAL(r8) :: t_soisno_bef(lbc:ubc,-nlevsnow+1:nlevsoil) !
  REAL(r8) :: lhabs(lbc:ubc)
  REAL(r8) :: esum1(lbc:ubc)
  REAL(r8) :: esum2(lbc:ubc)
  REAL(r8) :: zsum(lbc:ubc)
  REAL(r8) :: wsum(lbc:ubc)
  REAL(r8) :: wsum_end(lbc:ubc)
  REAL(r8) :: errsoi(1)
  REAL(r8) :: eflx_snomelt(1)
  REAL(r8) :: n22 ! inter for enhanced diffusivity by wfs
  REAL(r8) :: ded ! enhanced diffusivity by wfs

  CHARACTER*256 :: message



  DATA beta/0.4_r8, 0.4_r8/
  DATA za  /0.6_r8, 0.6_r8/
  ! DATA za  /0.1_r8, 0.1_r8/ ! by wfs

  !   For now, keep beta and za for shallow lake the same as deep lake, until better data is found.
  !   It looks like eta is key and that larger values give better results for shallow lakes.  Use
  !   empirical expression from Hakanson (below). This is still a very unconstrained parameter
  !   that deserves more attention.
  !   Some radiation will be allowed to reach the soil.
  !-----------------------------------------------------------------------



  ! 1!) Initialization
  ! Determine step size

  !    dtime = get_step_size()

  ! Initialize constants
  cwat = cpliq*denh2o ! water heat capacity per unit volume
  cice_eff = cpice*denh2o !use water density because layer depth is not adjusted
  !for freezing
  cfus = hfus*denh2o  ! latent heat per unit volume
  tkice_eff = tkice * denice/denh2o !effective conductivity since layer depth is constant
  km = tkwat/cwat     ! a constant (molecular diffusivity)

  ! Begin calculations

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)

     ! Initialize Ebal quantities computed below

     ocvts(c) = 0._r8
     ncvts(c) = 0._r8
     esum1(c) = 0._r8
     esum2(c) = 0._r8

  END DO

  ! Initialize set of previous time-step variables as in DriverInit,
  ! which is currently not called over lakes. This has to be done
  ! here because phase change will occur in this routine.
  ! Ice fraction of snow at previous time step

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        IF (j >= snl(c) + 1) THEN
           frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
        END IF
     END DO
  END DO

  ! Sum soil water.
  DO j = 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        IF (j == 1) wsum(c) = 0._r8
        wsum(c) = wsum(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)


     ! Prepare for lake layer temperature calculations below

     ! fin(c) = betaprime * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
     !     eflx_sh_tot(p) + eflx_lh_tot(p))
     ! fin(c) now passed from ShalLakeFluxes as eflx_gnet
    !  PRINT*, 'eflx_gnet',eflx_gnet
     fin(c) = eflx_gnet(p)
    !  PRINT*, 'eflx_gnet after use 1',eflx_gnet
    !    PRINT*, 't_lake after use 1',t_lake

  END DO

  ! 2!) Lake density

  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * &
             1000._r8*( 1.0_r8 - 1.9549e-05_r8*(ABS(t_lake(c,j)-277._r8))**1.68_r8 ) &
             + lake_icefrac(c,j)*denice
        ! Allow for ice fraction; assume constant ice density.
        ! Is this the right weighted average?
        ! Using this average will make sure that surface ice is treated properly during
        ! convective mixing.
     END DO
  END DO

  ! 3!) Diffusivity and implied thermal "conductivity" = diffusivity * cwat
  DO j = 1, nlevlake-1
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        drhodz = (rhow(c,j+1)-rhow(c,j)) / (z_lake(c,j+1)-z_lake(c,j))
        n2 = grav / rhow(c,j) * drhodz
        ! Fixed sign error here: our z goes up going down into the lake, so no negative
        ! sign is needed to make this positive unlike in Hostetler. --ZS
        num = 40._r8 * n2 * (vkc*z_lake(c,j))**2
        den = MAX( (ws(c)**2) * EXP(-2._r8*ks(c)*z_lake(c,j)), 1.e-10_r8 )
        ri = ( -1._r8 + SQRT( MAX(1._r8+num/den, 0._r8) ) ) / 20._r8
        n22 = n2
        IF (n22 < 0.000075_r8) n22 = 0.000075_r8
        ded = 1.04e-8_r8*n22**(-0.43_r8)
        IF (t_grnd(c) > tfrz .AND. t_lake(c,1) > tfrz .AND. snl(c) == 0) THEN
           ! ke = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
           IF( t_lake(c,1) > 277.15_r8 ) THEN
              IF (lakedepth(c) > 15.0 ) THEN
                 ! ke = 1.e+2_r8*vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
                ke = vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri) !by wfs
              ELSE
                ke = vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
              ENDIF
           ELSE
              IF (lakedepth(c) > 15.0 ) THEN
                 IF (lakedepth(c) > 150.0 ) THEN
                    ke = 1.e+5_r8*vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
                 ELSE
                    ke =1.e+4_r8*vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
                 END IF
              ELSE
                 ke = vkc*ws(c)*z_lake(c,j)/p0 * EXP(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
              ENDIF
           END IF

           IF (z_lake(c,j) > 20.0 ) THEN
             IF (z_lake(c,j) > 40.0 ) THEN
               kme(c,j) = 30*ded + 100*ke
             ELSE
               kme(c,j) = 100*ded + 100*ke
             END IF
           ELSE
             kme(c,j) = 100*ded + ke
           ENDIF
           ! end ke times by wfs

           IF ( j==1 ) kme(c,j) = km + ke
           ! end ke times by wfs

           tk_lake(c,j) = kme(c,j)*cwat
           ! If there is some ice in this layer (this should rarely happen because the surface
           ! is unfrozen and it will be unstable), still use the cwat to get out the tk b/c the eddy
           ! diffusivity equation assumes water.
        ELSE
           kme(c,j) = km
           tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                + tkwat*lake_icefrac(c,j) )
           ! Assume the resistances add as for the calculation of conductivities at layer interfaces.
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)

     j = nlevlake
     kme(c,nlevlake) = kme(c,nlevlake-1)

     IF (t_grnd(c) > tfrz .AND. t_lake(c,1) > tfrz .AND. snl(c) == 0) THEN
        tk_lake(c,j) = tk_lake(c,j-1)
     ELSE
        tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
             + tkwat*lake_icefrac(c,j) )
     END IF

     ! Use in surface flux calculation for next timestep.
     savedtke1(c) = kme(c,1)*cwat ! Will only be used if unfrozen
     ! set number of column levels for use by Tridiagonal below
     jtop(c) = snl(c) + 1
  END DO

  ! 4!) Heat source term: unfrozen lakes only
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fp = 1, num_shlakep
        p = filter_shlakep(fp)
        c = pcolumn(p)

        ! Set eta(:), the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
        ! (regression of Secchi Depth with lake depth for small glacial basin lakes), and the
        ! Poole & Atkins expression for extinction coeffient of 1.7 / Secchi Depth (m).
        !#ifndef ETALAKE
        ! eta(:) = 1.1925_r8*lakedepth(c)**(-0.424) !by wfs
        eta(:) = 1
        !#else
        !          eta(:) = ETALAKE
        !#endif

        zin  = z_lake(c,j) - 0.5_r8*dz_lake(c,j)
        zout = z_lake(c,j) + 0.5_r8*dz_lake(c,j)
        rsfin  = EXP( -eta(islak)*MAX(  zin-za(islak),0._r8 ) )
        rsfout = EXP( -eta(islak)*MAX( zout-za(islak),0._r8 ) )
        !extra energy from inflow /W by wfs
        ! extraE(c,j) = infv(c,j)*dtime*cwat*(inft(c,j)-t_lake(c,j))
        ! PRINT*, j,'t_lake=',t_lake(c,j),'extraE=',extraE(c,j)/dz_lake(c,j)
        ! Let rsfout for bottom layer go into soil.
        ! This looks like it should be robust even for pathological cases,
        ! like lakes thinner than za.
        ! inflow energy added as another heat source by wfs

        IF (t_grnd(c) > tfrz .AND. t_lake(c,1) > tfrz .AND. snl(c) == 0) THEN
           phidum = (rsfin-rsfout) * sabg(p) * (1._r8-beta(islak))
          !  PRINT*, 'phidum=',phidum/dz_lake(c,j)
           IF (j == nlevlake) THEN
              phi_soil(c) = rsfout * sabg(p) * (1._r8-beta(islak))
           END IF
        ELSE IF (j == 1 .AND. snl(c) == 0) THEN !if frozen but no snow layers
           phidum = sabg(p) * (1._r8-beta(islak))
        ELSE !radiation absorbed at surface
           phidum = 0._r8
           IF (j == nlevlake) phi_soil(c) = 0._r8
        END IF
        phi(c,j) = phidum

     END DO
  END DO

  ! 5!) Set thermal properties and check initial energy content.

  ! For lake
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
     END DO
  END DO
  ! PRINT*, 'eflx_gnet after use 2',eflx_gnet
  ! For snow / soil
  CALL SoilThermProp_Lake (watsat,tksatu,tkmg,tkdry,csol, &
       snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
       tk, cv, tktopsoillay)
  ! PRINT*, 'eflx_gnet after use 3',eflx_gnet
  ! Sum cv*t_lake for energy check
  ! Include latent heat term, and correction for changing heat capacity with phase change.

  ! This will need to be over all soil / lake / snow layers. Lake is below.
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        !          ocvts(c) = ocvts(c) + cv_lake(c,j)*t_lake(c,j) &
        ocvts(c) = ocvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
             + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) !&
        !                   + (cwat-cice_eff)*lake_icefrac(c)*tfrz*dz_lake(c,j) !enthalpy reconciliation term
        t_lake_bef(c,j) = t_lake(c,j)
     END DO
  END DO

  ! Now do for soil / snow layers
  DO j = -nlevsnow + 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        IF (j >= jtop(c)) THEN
           !             ocvts(c) = ocvts(c) + cv(c,j)*t_soisno(c,j) &
           ocvts(c) = ocvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                + hfus*h2osoi_liq(c,j) !&
           !                      + (cpliq-cpice)*h2osoi_ice(c,j)*tfrz !enthalpy reconciliation term
           IF (j == 1 .AND. h2osno(c) > 0._r8 .AND. j == jtop(c)) THEN
              ocvts(c) = ocvts(c) - h2osno(c)*hfus
           END IF
           t_soisno_bef(c,j) = t_soisno(c,j)
           !             if(abs(t_soisno(c,j)-288) > 150)   then
           !                WRITE( message,* ) 'WARNING: Extreme t_soisno at c, level',c, j
           !                CALL wrf_error_fatal ( message )
           !             endif
        END IF
     END DO
  END DO

!!!!!!!!!!!!!!!!!!!
  ! 6!) Set up vector r and vectors a, b, c1 that define tridiagonal matrix

  ! Heat capacity and resistance of snow without snow layers (<1cm) is ignored during diffusion,
  ! but its capacity to absorb latent heat may be used during phase change.

  ! Set up interface depths, zx, heat capacities, cvx, solar source terms, phix, and temperatures, tx.
  DO j = -nlevsnow+1, nlevlake+nlevsoil
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)

        jprime = j - nlevlake

        IF (j >= jtop(c)) THEN
           IF (j < 1) THEN !snow layer
              zx(c,j) = z(c,j)
              cvx(c,j) = cv(c,j)
              phix(c,j) = 0._r8
              tx(c,j) = t_soisno(c,j)

           ELSE IF (j <= nlevlake) THEN !lake layer
              zx(c,j) = z_lake(c,j)
              cvx(c,j) = cv_lake(c,j)
              phix(c,j) = phi(c,j)
              tx(c,j) = t_lake(c,j)

           ELSE !soil layer
              zx(c,j) = zx(c,nlevlake) + dz_lake(c,nlevlake)/2._r8 + z(c,jprime)
              cvx(c,j) = cv(c,jprime)
              IF (j == nlevlake + 1) THEN !top soil layer
                 phix(c,j) = phi_soil(c)
              ELSE !middle or bottom soil layer
                 phix(c,j) = 0._r8
              END IF
              tx(c,j) = t_soisno(c,jprime)

           END IF
        END IF

     END DO
  END DO

  ! Determine interface thermal conductivities, tkix

  DO j = -nlevsnow+1, nlevlake+nlevsoil
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)

        jprime = j - nlevlake

        IF (j >= jtop(c)) THEN
           IF (j < 0) THEN !non-bottom snow layer
              tkix(c,j) = tk(c,j)
           ELSE IF (j == 0) THEN !bottom snow layer
              dzp = zx(c,j+1) - zx(c,j)
              tkix(c,j) = tk_lake(c,1)*tk(c,j)*dzp / &
                   (tk(c,j)*z_lake(c,1) + tk_lake(c,1)*(-z(c,j)) )
              ! tk(c,0) is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
           ELSE IF (j < nlevlake) THEN !non-bottom lake layer
              tkix(c,j) = ( tk_lake(c,j)*tk_lake(c,j+1) * (dz_lake(c,j+1)+dz_lake(c,j)) ) &
                   / ( tk_lake(c,j)*dz_lake(c,j+1) + tk_lake(c,j+1)*dz_lake(c,j) )
           ELSE IF (j == nlevlake) THEN !bottom lake layer
              dzp = zx(c,j+1) - zx(c,j)
              tkix(c,j) = (tktopsoillay(c)*tk_lake(c,j)*dzp / &
                   (tktopsoillay(c)*dz_lake(c,j)/2._r8 + tk_lake(c,j)*z(c,1) ) )
              ! tktopsoillay is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
           ELSE !soil layer
              tkix(c,j) = tk(c,jprime)
           END IF
        END IF

     END DO
  END DO


  ! Determine heat diffusion through the layer interface and factor used in computing
  ! tridiagonal matrix and set up vector r and vectors a, b, c1 that define tridiagonal
  ! matrix and solve system

  DO j = -nlevsnow+1, nlevlake+nlevsoil
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        IF (j >= jtop(c)) THEN
           IF (j < nlevlake+nlevsoil) THEN !top or interior layer
              factx(c,j) = dtime/cvx(c,j)
              fnx(c,j) = tkix(c,j)*(tx(c,j+1)-tx(c,j))/(zx(c,j+1)-zx(c,j))
           ELSE !bottom soil layer
              factx(c,j) = dtime/cvx(c,j)
              fnx(c,j) = 0._r8 !not used
           END IF
        END IF
     ENDDO
  END DO

  !    print *,10,10
  DO j =  -nlevsnow+1,nlevlake+nlevsoil
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        IF (j >= jtop(c)) THEN
           IF (j == jtop(c)) THEN !top layer
              dzp    = zx(c,j+1)-zx(c,j)
              a(c,j) = 0._r8
              b(c,j) = 1+(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
              c1(c,j) =  -(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
              r(c,j) = tx(c,j) + factx(c,j)*( fin(c) + phix(c,j) + cnfac*fnx(c,j) )
           ELSE IF (j < nlevlake+nlevsoil) THEN !middle layer
              dzm    = (zx(c,j)-zx(c,j-1))
              dzp    = (zx(c,j+1)-zx(c,j))
              a(c,j) =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j-1)/dzm
              b(c,j) = 1._r8+ (1._r8-cnfac)*factx(c,j)*(tkix(c,j)/dzp + tkix(c,j-1)/dzm)
              c1(c,j) =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j)/dzp
              r(c,j) = tx(c,j) + cnfac*factx(c,j)*( fnx(c,j) - fnx(c,j-1) ) + factx(c,j)*phix(c,j)
           ELSE  !bottom soil layer
              dzm     = (zx(c,j)-zx(c,j-1))
              a(c,j) =   - (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
              b(c,j) = 1._r8+ (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
              c1(c,j) = 0._r8
              r(c,j) = tx(c,j) - cnfac*factx(c,j)*fnx(c,j-1)
           END IF
        END IF
     ENDDO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    print *,'t1',tx(1,1)
  ! PRINT*, 'eflx_gnet after use 3',eflx_gnet
  ! PRINT*, 't_lake after use 3',t_lake
  ! PRINT*, 'tx after use 3',tx
  ! 7!) Solve for tdsolution

  CALL Tridiagonal(lbc, ubc, -nlevsnow + 1, nlevlake + nlevsoil, jtop, num_shlakec, filter_shlakec, &
       a, b, c1, r, tx)  !wfs
  !   print *,a(1,1), b(1,1), c1(1,1), r(1,1), tx(1,1)
  !    print *,'t2', tx(1,1)
  ! PRINT*, 'eflx_gnet after use 3.5',eflx_gnet
  ! PRINT*, 't_lake after use 3.5',t_lake
  ! PRINT*, 'tx after use 3.5',tx
  ! Set t_soisno and t_lake
  DO j = -nlevsnow+1, nlevlake + nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        jprime = j - nlevlake

        ! Don't do anything with invalid snow layers.
        IF (j >= jtop(c)) THEN
           IF (j < 1) THEN !snow layer
              t_soisno(c,j) = tx(c,j)
           ELSE IF (j <= nlevlake) THEN !lake layer
              t_lake(c,j)   = tx(c,j)
           ELSE !soil layer
              t_soisno(c,jprime) = tx(c,j)
           END IF
        END IF
     END DO
  END DO

!!!!!!!!!!!!!!!!!!!!!!!
  ! PRINT*, 'eflx_gnet after use 4',eflx_gnet
  ! PRINT*, 't_lake after use 4',t_lake

  ! 8!) Sum energy content and total energy into lake for energy check. Any errors will be from the
  !     Tridiagonal solution.
  !#define LAKEDEBUG
  !#if (defined LAKEDEBUG)
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        esum1(c) = esum1(c) + (t_lake(c,j)-t_lake_bef(c,j))*cv_lake(c,j)
        esum2(c) = esum2(c) + (t_lake(c,j)-tfrz)*cv_lake(c,j)
     END DO
  END DO

  DO j = -nlevsnow+1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        IF (j >= jtop(c)) THEN
           esum1(c) = esum1(c) + (t_soisno(c,j)-t_soisno_bef(c,j))*cv(c,j)
           esum2(c) = esum2(c) + (t_soisno(c,j)-tfrz)*cv(c,j)
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     ! Again assuming only one pft per column
     !          esum1(c) = esum1(c) + lhabs(c)
     errsoi(c) = esum1(c)/dtime - eflx_soil_grnd(p)
     ! eflx_soil_grnd includes all the solar radiation absorbed in the lake,
     ! unlike eflx_gnet
     !          if(abs(errsoi(c)) > 1.e-5_r8) then
     !             WRITE( message,* )'Primary soil energy conservation error in shlake &
     !                                column during Tridiagonal Solution,', 'error (W/m^2):', c, errsoi(c)
     !             CALL wrf_error_fatal ( message )
     !          end if
  END DO
  ! This has to be done before convective mixing because the heat capacities for each layer
  ! will get scrambled.

  !#endif

!!!!!!!!!!!!!!!!!!!!!!!
  ! PRINT*, 'eflx_gnet after use 5',eflx_gnet

  ! 9!) Phase change
  CALL PhaseChange_Lake (snl,h2osno,dz,dz_lake,                            & !i
       t_soisno,h2osoi_liq,h2osoi_ice,               & !i&o
       lake_icefrac,t_lake, snowdp,                  & !i&o
       qflx_snomelt,eflx_snomelt,imelt,              & !o
       cv, cv_lake,                                  & !i&o
       lhabs)                                          !o
  !print *,'t3',t_lake(1,1)
!!!!!!!!!!!!!!!!!!!!!!!
  ! PRINT*, 'eflx_gnet after use 6',eflx_gnet

  ! 9.5!) Second energy check and water check.  Now check energy balance before and after phase
  !       change, considering the possibility of changed heat capacity during phase change, by
  !       using initial heat capacity in the first step, final heat capacity in the second step,
  !       and differences from tfrz only to avoid enthalpy correction for (cpliq-cpice)*melt*tfrz.
  !       Also check soil water sum.
  !#define LAKEDEBUG
  !#if (defined LAKEDEBUG)
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        esum2(c) = esum2(c) - (t_lake(c,j)-tfrz)*cv_lake(c,j)
     END DO
  END DO

  DO j = -nlevsnow+1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        IF (j >= jtop(c)) THEN
           esum2(c) = esum2(c) - (t_soisno(c,j)-tfrz)*cv(c,j)
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     ! Again assuming only one pft per column
     esum2(c) = esum2(c) - lhabs(c)
     errsoi(c) = esum2(c)/dtime
     !          if(abs(errsoi(c)) > 1.e-5_r8) then
     !             write(message,*)'Primary soil energy conservation error in shlake column during Phase Change, error (W/m^2):', &
     !                       c, errsoi(c)
     !             CALL wrf_error_fatal ( message )
     !          end if
  END DO

  ! Check soil water
  ! Sum soil water.
  DO j = 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        IF (j == 1) wsum_end(c) = 0._r8
        wsum_end(c) = wsum_end(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
        !          if (j == nlevsoil) then
        !!             if (abs(wsum(c)-wsum_end(c))>1.e-7_r8) then
        !!                write(message,*)'Soil water balance error during phase change in ShalLakeTemperature.', &
        !!                          'column, error (kg/m^2):', c, wsum_end(c)-wsum(c)
        !!                CALL wrf_error_fatal ( message )
        !!             end if
        !          end if
     END DO
  END DO

  !#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 10!) Convective mixing: make sure fracice*dz is conserved, heat content c*dz*T is conserved, and
  ! all ice ends up at the top. Done over all lakes even if frozen.
  ! Either an unstable density profile or ice in a layer below an incompletely frozen layer will trigger.

  !Recalculate density
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * &
             1000._r8*( 1.0_r8 - 1.9549e-05_r8*(ABS(t_lake(c,j)-277._r8))**1.68_r8 ) &
             + lake_icefrac(c,j)*denice
     END DO
  END DO

  ! DO j = 1, nlevlake-1
  !    !dir$ concurrent
  !    !cdir nodep
  !    DO fc = 1, num_shlakec
  !       c = filter_shlakec(fc)
  !       IF (rhow(c,j) > rhow(c,j+1)+0.0001_r8 .OR. &
  !             (lake_icefrac(c,j) < 1._r8 .AND. lake_icefrac(c,j+1) > 0._r8) ) THEN
  !             print*, 'Mix',j,rhow(c,j),rhow(c,j+1)
  !       END IF
  !    END DO
  ! END DO

  DO j = 1, nlevlake-1
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        qav(c) = 0._r8
        nav(c) = 0._r8
        iceav(c) = 0._r8
     END DO

     ! IF (rhow(c,j) > rhow(c,j+1) .OR. &
     !      (lake_icefrac(c,j) < 1._r8 .AND. lake_icefrac(c,j+1) > 0._r8) ) THEN
     !    print *, 'Mix',j
     ! END IF

     DO i = 1, j+1
        !dir$ concurrent
        !cdir nodep
        DO fc = 1, num_shlakec
           c = filter_shlakec(fc)
           IF (rhow(c,j) > rhow(c,j+1)+0.0001_r8 .OR. &
                (lake_icefrac(c,j) < 1._r8 .AND. lake_icefrac(c,j+1) > 0._r8) ) THEN
              !#define LAKEDEBUG
              !#if (defined LAKEDEBUG)
              !                if (i==1)  then
              !                  write(message,*), 'Convective Mixing in column ', c, '.'
              !                  CALL wrf_message(message)
              !                endif
              !#endif
              qav(c) = qav(c) + dz_lake(c,i)*(t_lake(c,i)-tfrz) * &
                   ((1._r8 - lake_icefrac(c,i))*cwat + lake_icefrac(c,i)*cice_eff)
              !                tav(c) = tav(c) + t_lake(c,i)*dz_lake(c,i)
              iceav(c) = iceav(c) + lake_icefrac(c,i)*dz_lake(c,i)
              nav(c) = nav(c) + dz_lake(c,i)
           END IF
        END DO
     END DO

     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        IF (rhow(c,j) > rhow(c,j+1)+0.0001_r8 .OR. &
             (lake_icefrac(c,j) < 1._r8 .AND. lake_icefrac(c,j+1) > 0._r8) ) THEN
           qav(c) = qav(c)/nav(c)
           iceav(c) = iceav(c)/nav(c)
           !If the average temperature is above freezing, put the extra energy into the water.
           !If it is below freezing, take it away from the ice.
           IF (qav(c) > 0._r8) THEN
              tav_froz(c) = 0._r8 !Celsius
              tav_unfr(c) = qav(c) / ((1._r8 - iceav(c))*cwat)
           ELSE IF (qav(c) < 0._r8) THEN
              tav_froz(c) = qav(c) / (iceav(c)*cice_eff)
              tav_unfr(c) = 0._r8 !Celsius
           ELSE
              tav_froz(c) = 0._r8
              tav_unfr(c) = 0._r8
           END IF
        END IF
     END DO

     DO i = 1, j+1
        !dir$ concurrent
        !cdir nodep
        DO fc = 1, num_shlakec
           c = filter_shlakec(fc)
           IF (nav(c) > 0._r8) THEN
              !             if(0==1) then

              !Put all the ice at the top.!
              !If the average temperature is above freezing, put the extra energy into the water.
              !If it is below freezing, take it away from the ice.
              !For the layer with both ice & water, be careful to use the average temperature
              !that preserves the correct total heat content given what the heat capacity of that
              !layer will actually be.
              IF (i == 1) zsum(c) = 0._r8

              IF ((zsum(c)+dz_lake(c,i))/nav(c) <= iceav(c)) THEN
                 lake_icefrac(c,i) = 1._r8
                 t_lake(c,i) = tav_froz(c) + tfrz
              ELSE IF (zsum(c)/nav(c) < iceav(c)) THEN
                 lake_icefrac(c,i) = (iceav(c)*nav(c) - zsum(c)) / dz_lake(c,i)
                 ! Find average value that preserves correct heat content.
                 t_lake(c,i) = ( lake_icefrac(c,i)*tav_froz(c)*cice_eff &
                      + (1._r8 - lake_icefrac(c,i))*tav_unfr(c)*cwat ) &
                      / ( lake_icefrac(c,i)*cice_eff + (1-lake_icefrac(c,i))*cwat ) + tfrz
              ELSE
                 lake_icefrac(c,i) = 0._r8
                 t_lake(c,i) = tav_unfr(c) + tfrz
              END IF
              zsum(c) = zsum(c) + dz_lake(c,i)

              rhow(c,i) = (1._r8 - lake_icefrac(c,i)) * &
                   1000._r8*( 1.0_r8 - 1.9549e-05_r8*(ABS(t_lake(c,i)-277._r8))**1.68_r8 ) &
                   + lake_icefrac(c,i)*denice
           END IF
        END DO
     END DO
  END DO
! PRINT*, 'eflx_gnet after use 6',eflx_gnet
! PRINT*, 't_lake after use 6',t_lake
!!!!!!!!!!!!!!!!!!!!!!!
  ! 11!) Re-evaluate thermal properties and sum energy content.
  ! For lake
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
        !#define LAKEDEBUG
        !#if (defined LAKEDEBUG)
        !          write(message,*)'Lake Ice Fraction, c, level:', c, j, lake_icefrac(c,j)
        !          CALL wrf_message(message)
        !#endif
     END DO
  END DO




  ! For snow / soil
  !  call SoilThermProp_Lake(lbc, ubc, num_shlakec, filter_shlakec, tk, cv, tktopsoillay)
  CALL SoilThermProp_Lake (watsat,tksatu,tkmg,tkdry,csol, &
        snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
        tk, cv, tktopsoillay)


  ! PRINT*, 'eflx_gnet after use 7',eflx_gnet
  ! PRINT*, 'ncvts after use 7',ncvts
  ! PRINT*, 't_lake after use 7',t_lake
  ! PRINT*, 'cv',cv
  ! Do as above to sum energy content
  DO j = 1, nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        ! PRINT*, 'components of ncvts after use 7.1.0',ncvts(c) , cv_lake(c,j),t_lake(c,j),tfrz,&
        !  cfus,dz_lake(c,j),(1._r8-lake_icefrac(c,j))

        !          ncvts(c) = ncvts(c) + cv_lake(c,j)*t_lake(c,j) &
        ncvts(c) = ncvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
             + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) !&
        !                   + (cwat-cice_eff)*lake_icefrac(c)*tfrz*dz_lake(c,j) !enthalpy reconciliation term
        fin(c) = fin(c) + phi(c,j)

        ! PRINT*, 'ncvts after use 7.1.0',ncvts
     END DO
  END DO
  !print *,tfrz
! PRINT*, 'ncvts after use 7.1',ncvts
  DO j = -nlevsnow + 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        ! PRINT*, 'components of ncvts after use 7.2.0',ncvts(c), cv(c,j),t_soisno(c,j),tfrz,&
        !  hfus,h2osoi_liq(c,j), h2osno(c)

        IF (j >= jtop(c)) THEN
           !             ncvts(c) = ncvts(c) + cv(c,j)*t_soisno(c,j) &
           ncvts(c) = ncvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                + hfus*h2osoi_liq(c,j) !&
           !                      + (cpliq-cpice)*h2osoi_ice(c,j)*tfrz !enthalpy reconciliation term
           IF (j == 1 .AND. h2osno(c) > 0._r8 .AND. j == jtop(c)) THEN
              ncvts(c) = ncvts(c) - h2osno(c)*hfus
           END IF
        END IF
        IF (j == 1) fin(c) = fin(c) + phi_soil(c)
     END DO
  END DO
! PRINT*, 'ncvts after use 7.2',ncvts

  ! Check energy conservation.
  !  print*, 'before checking:'
  !  print*, 'errsoi(c)',errsoi(c)
  !  print*, 'eflx_sh_tot',eflx_sh_tot(p)
  !  print*, 'eflx_gnet',eflx_gnet(p)
  ! PRINT*, 'eflx_gnet after use 8',eflx_gnet

  DO fp = 1, num_shlakep

     p = filter_shlakep(fp)
     c = pcolumn(p)
     errsoi(c) = (ncvts(c)-ocvts(c)) / dtime - fin(c)   !CUOWU
     !  #ifndef LAKEDEBUG
     !  !       if (abs(errsoi(c)) < 0.10_r8) then ! else send to Balance Check and abort
     !        if (abs(errsoi(c)) < 10._r8) then ! else send to Balance Check and abort
     !  #else
     !        if (abs(errsoi(c)) < 1._r8) then
     !  #endif
    !  PRINT*, 'sub of errsoi',ncvts(c),ocvts(c) , dtime , fin(c)
    !  PRINT*, 'errsoi(c)',errsoi(c)
    !  PRINT*, 'eflx_gnet after use 8.1',eflx_gnet
     eflx_sh_tot(p) = eflx_sh_tot(p) - errsoi(c)
     eflx_sh_grnd(p) = eflx_sh_grnd(p) - errsoi(c)
     eflx_soil_grnd(p) = eflx_soil_grnd(p) + errsoi(c)
     eflx_gnet(p) = eflx_gnet(p) + errsoi(c)   !CUOWU

     !          if (abs(errsoi(c)) > 1.e-3_r8) then
     !          if (abs(errsoi(c)) > 1.e-1_r8) then
     !             write(message,*)'errsoi incorporated into sensible heat in ShalLakeTemperature: c, (W/m^2):', c, errsoi(c)
     !             CALL wrf_message(message)
     !          end if
     errsoi(c) = 0._r8
     !  print*, 'after checking:'
     !  PRINT*, 'errsoi(c)',errsoi(c)
     !  print*, 'eflx_sh_tot',eflx_sh_tot(p)
     !  print*, 'eflx_gnet',eflx_gnet(p)
     !  print*, ''
     !#if (defined LAKEDEBUG)
     !       else
     !          write(message,*)'Soil Energy Balance Error at column, ', c, 'G, fintotal, column E tendency = ', &
     !             eflx_gnet(p), fin(c), (ncvts(c)-ocvts(c)) / dtime
     !          CALL wrf_message(message)
     !#endif
     !       end if
    !  PRINT*, 'eflx_gnet after use 8.2',eflx_gnet
  END DO
  ! This loop assumes only one point per column.
  ! PRINT*, 'eflx_gnet after use 9',eflx_gnet

END SUBROUTINE ShalLakeTemperature

! CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! ROUTINE: PhaseChange_Lake
!
! !INTERFACE:
SUBROUTINE PhaseChange_Lake (snl,h2osno,dz,dz_lake,                        & !i
     t_soisno,h2osoi_liq,h2osoi_ice,               & !i&o
     lake_icefrac,t_lake, snowdp,                  & !i&o
     qflx_snomelt,eflx_snomelt,imelt,              & !o
     cv, cv_lake,                                  & !i&o
     lhabs)                                          !o
  !=============================================================================================
  ! !DESCRIPTION:
  ! Calculation of the phase change within snow, soil, & lake layers:
  ! (1) Check the conditions for which the phase change may take place,
  !     i.e., the layer temperature is great than the freezing point
  !     and the ice mass is not equal to zero (i.e. melting),
  !     or the layer temperature is less than the freezing point
  !     and the liquid water mass is greater than the allowable supercooled
  !    (i.e. freezing).
  ! (2) Assess the amount of phase change from the energy excess (or deficit)
  !     after setting the layer temperature to freezing point, depending on
  !     how much water or ice is available.
  ! (3) Re-adjust the ice and liquid mass, and the layer temperature: either to
  !     the freezing point if enough water or ice is available to fully compensate,
  !     or to a remaining temperature.
  ! The specific heats are assumed constant. Potential cycling errors resulting from
  ! this assumption will be trapped at the end of ShalLakeTemperature.
  ! !CALLED FROM:
  ! subroutine ShalLakeTemperature in this module
  !
  ! !REVISION HISTORY:
  ! 04/2009 Zack Subin: Initial code
  !==============================================================================================
  ! !USES:
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER , INTENT(in) :: snl(1)
  REAL(r8), INTENT(inout) :: h2osno(1)
  REAL(r8), INTENT(in) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: dz_lake(1,nlevlake)


  REAL(r8), INTENT(inout) :: snowdp(1)
  REAL(r8), INTENT(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: lake_icefrac(1,nlevlake)
  REAL(r8), INTENT(inout) :: t_lake(1,nlevlake)

  REAL(r8), INTENT(out) :: qflx_snomelt(1)
  REAL(r8), INTENT(out) :: eflx_snomelt(1)
  INTEGER, INTENT(out)  :: imelt(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: cv(lbc:ubc,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: cv_lake (lbc:ubc,1:nlevlake)
  REAL(r8), INTENT(out):: lhabs(1:1)



  INTEGER  :: j,c,g
  INTEGER  :: fc
  REAL(r8) :: heatavail
  REAL(r8) :: heatrem
  REAL(r8) :: melt
  REAL(r8), PARAMETER :: smallnumber = 1.e-7_r8
  LOGICAL  :: dophasechangeflag
  !-----------------------------------------------------------------------

  !    dtime = get_step_size()

  ! Initialization

  !dir$ concurrent
  !cdir nodep
  DO fc = 1,num_shlakec
     c = filter_shlakec(fc)

     qflx_snomelt(c) = 0._r8
     eflx_snomelt(c) = 0._r8
     lhabs(c)        = 0._r8
  END DO

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)

        IF (j >= snl(c) + 1) imelt(c,j) = 0
     END DO
  END DO

  ! Check for case of snow without snow layers and top lake layer temp above freezing.

  !dir$ concurrent
  !cdir nodep
  DO fc = 1,num_shlakec
     c = filter_shlakec(fc)

     IF (snl(c) == 0 .AND. h2osno(c) > 0._r8 .AND. t_lake(c,1) > tfrz) THEN
        heatavail = (t_lake(c,1) - tfrz) * cv_lake(c,1)
        melt = MIN(h2osno(c), heatavail/hfus)
        heatrem = MAX(heatavail - melt*hfus, 0._r8)
        !catch small negative value to keep t at tfrz
        t_lake(c,1) = tfrz + heatrem/(cv_lake(c,1))
        snowdp(c) = snowdp(c)*(1._r8 - melt/h2osno(c))
        h2osno(c) = h2osno(c) - melt
        lhabs(c) = lhabs(c) + melt*hfus
        qflx_snomelt(c) = qflx_snomelt(c) + melt
        ! Prevent tiny residuals
        IF (h2osno(c) < smallnumber) h2osno(c) = 0._r8
        IF (snowdp(c) < smallnumber) snowdp(c) = 0._r8
     END IF
  END DO

  ! Lake phase change

  DO j = 1,nlevlake
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)

        dophasechangeflag = .FALSE.
        IF (t_lake(c,j) > tfrz .AND. lake_icefrac(c,j) > 0._r8) THEN ! melting
           dophasechangeflag = .TRUE.
           heatavail = (t_lake(c,j) - tfrz) * cv_lake(c,j)
           melt = MIN(lake_icefrac(c,j)*denh2o*dz_lake(c,j), heatavail/hfus)
           !denh2o is used because layer thickness is not adjusted for freezing
           heatrem = MAX(heatavail - melt*hfus, 0._r8)
           !catch small negative value to keep t at tfrz
        ELSE IF (t_lake(c,j) < tfrz .AND. lake_icefrac(c,j) < 1._r8) THEN !freezing
           dophasechangeflag = .TRUE.
           heatavail = (t_lake(c,j) - tfrz) * cv_lake(c,j)
           melt = MAX(-(1._r8-lake_icefrac(c,j))*denh2o*dz_lake(c,j), heatavail/hfus)
           !denh2o is used because layer thickness is not adjusted for freezing
           heatrem = MIN(heatavail - melt*hfus, 0._r8)
           !catch small positive value to keep t at tfrz
        END IF
        ! Update temperature and ice fraction.
        IF (dophasechangeflag) THEN
           lake_icefrac(c,j) = lake_icefrac(c,j) - melt/(denh2o*dz_lake(c,j))
           lhabs(c) = lhabs(c) + melt*hfus
           ! Update heat capacity
           cv_lake(c,j) = cv_lake(c,j) + melt*(cpliq-cpice)
           t_lake(c,j) = tfrz + heatrem/cv_lake(c,j)
           ! Prevent tiny residuals
           IF (lake_icefrac(c,j) > 1._r8 - smallnumber) lake_icefrac(c,j) = 1._r8
           IF (lake_icefrac(c,j) < smallnumber)         lake_icefrac(c,j) = 0._r8
        END IF
     END DO
  END DO

  ! Snow & soil phase change

  DO j = -nlevsnow+1,nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        dophasechangeflag = .FALSE.

        IF (j >= snl(c) + 1) THEN

           IF (t_soisno(c,j) > tfrz .AND. h2osoi_ice(c,j) > 0._r8) THEN ! melting
              dophasechangeflag = .TRUE.
              heatavail = (t_soisno(c,j) - tfrz) * cv(c,j)
              melt = MIN(h2osoi_ice(c,j), heatavail/hfus)
              heatrem = MAX(heatavail - melt*hfus, 0._r8)
              !catch small negative value to keep t at tfrz
              IF (j <= 0) THEN !snow
                 imelt(c,j) = 1
                 qflx_snomelt(c) = qflx_snomelt(c) + melt
              END IF
           ELSE IF (t_soisno(c,j) < tfrz .AND. h2osoi_liq(c,j) > 0._r8) THEN !freezing
              dophasechangeflag = .TRUE.
              heatavail = (t_soisno(c,j) - tfrz) * cv(c,j)
              melt = MAX(-h2osoi_liq(c,j), heatavail/hfus)
              heatrem = MIN(heatavail - melt*hfus, 0._r8)
              !catch small positive value to keep t at tfrz
              IF (j <= 0) THEN !snow
                 imelt(c,j) = 2
                 qflx_snomelt(c) = qflx_snomelt(c) + melt
                 ! Does this works for both signs of melt in SnowHydrology? I think
                 ! qflx_snomelt(c) is just output.
              END IF
           END IF

           ! Update temperature and soil components.
           IF (dophasechangeflag) THEN
              h2osoi_ice(c,j) = h2osoi_ice(c,j) - melt
              h2osoi_liq(c,j) = h2osoi_liq(c,j) + melt
              lhabs(c) = lhabs(c) + melt*hfus
              ! Update heat capacity
              cv(c,j) = cv(c,j) + melt*(cpliq-cpice)
              t_soisno(c,j) = tfrz + heatrem/cv(c,j)
              ! Prevent tiny residuals
              IF (h2osoi_ice(c,j) < smallnumber) h2osoi_ice(c,j) = 0._r8
              IF (h2osoi_liq(c,j) < smallnumber) h2osoi_liq(c,j) = 0._r8
           END IF

        END IF
     END DO
  END DO

  ! Update eflx_snomelt(c)
  !dir$ concurrent
  !cdir nodep
  DO fc = 1,num_shlakec
     c = filter_shlakec(fc)
     eflx_snomelt(c) = qflx_snomelt(c)*hfus
  END DO
!!!

END SUBROUTINE PhaseChange_Lake

SUBROUTINE lakeini(IVGTYP,         ISLTYP,          HT,              SNOW, & !i
     lake_min_elev,     restart,        lakedepth_default, lake_depth,     &
     lakedepth2d,    savedtke12d,     snowdp2d,        h2osno2d,       & !o
     snl2d,          t_grnd2d,        t_lake3d,        lake_icefrac3d, &
     z_lake3d,       dz_lake3d,       t_soisno3d,      h2osoi_ice3d,   &
     h2osoi_liq3d,   h2osoi_vol3d,    z3d,             dz3d,           &
     zi3d,           watsat3d,        csol3d,          tkmg3d,         &
     iswater,        xice,            xice_threshold,  xland,   tsk,   &
     lakemask,       lakeflag,                                         &
     lake_depth_flag, use_lakedepth,                                   &
     tkdry3d,        tksatu3d,        lake,            its, ite, jts, jte, &
     ims,ime, jms,jme)

  !==============================================================================
  ! This subroutine was first edited by Hongping Gu for coupling
  ! 07/20/2010
  !==============================================================================

  !  USE module_wrf_error
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  ! INTEGER :: i
  !
  ! DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
  !      10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
  !
  ! DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
  !      33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER , INTENT (IN) :: iswater
  REAL(r8),     INTENT(IN)  :: xice_threshold
  REAL(r8), DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   XICE
  REAL(r8), DIMENSION( ims:ime , jms:jme ), INTENT(IN)::      TSK
  REAL(r8), DIMENSION( ims:ime, jms:jme )  ,INTENT(INOUT)  :: XLAND

  REAL(r8), DIMENSION( ims:ime , jms:jme ) ::   LAKEMASK
  INTEGER , INTENT (IN) :: lakeflag

  INTEGER , INTENT (INOUT) :: lake_depth_flag
  INTEGER , INTENT (IN) ::   use_lakedepth

  LOGICAL , INTENT(IN)      ::     restart
  INTEGER,  INTENT(IN   )   ::     ims,ime, jms,jme
  INTEGER,  INTENT(IN   )   ::     its,ite, jts,jte
  INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)       :: IVGTYP,       &
       ISLTYP
  REAL(r8),    DIMENSION( ims:ime, jms:jme ), INTENT(IN)       :: HT
  REAL(r8),    DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)    :: SNOW
  REAL(r8),    INTENT(in)         :: lakedepth_default,lake_min_elev

  REAL(r8),    DIMENSION(ims:ime,jms:jme ),INTENT(out) :: lakedepth2d,    &
       savedtke12d
  REAL(r8),    DIMENSION(ims:ime,jms:jme ),INTENT(out)    :: snowdp2d,       &
       h2osno2d,       &
       snl2d,          &
       t_grnd2d

  REAL(r8),    DIMENSION( ims:ime,1:nlevlake, jms:jme ),INTENT(out)   :: t_lake3d,       &
       lake_icefrac3d, &
       z_lake3d,       &
       dz_lake3d
  REAL(r8),    DIMENSION( ims:ime,-nlevsnow+1:nlevsoil, jms:jme ),INTENT(out)   :: t_soisno3d,     &
       h2osoi_ice3d,   &
       h2osoi_liq3d,   &
       h2osoi_vol3d,   &
       z3d,            &
       dz3d
  REAL(r8),    DIMENSION( ims:ime,1:nlevsoil, jms:jme ),INTENT(out)            :: watsat3d,       &
       csol3d,         &
       tkmg3d,         &
       tkdry3d,        &
       tksatu3d
  REAL(r8),    DIMENSION( ims:ime,-nlevsnow+0:nlevsoil, jms:jme ),INTENT(out)   :: zi3d

  LOGICAL, DIMENSION( ims:ime, jms:jme ),INTENT(out)                      :: lake
  REAL(r8), OPTIONAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)    ::  lake_depth

  REAL(r8),    DIMENSION( ims:ime,1:nlevsoil, jms:jme )   :: bsw3d,    &
       bsw23d,   &
       psisat3d, &
       vwcsat3d, &
       watdry3d, &
       watopt3d, &
       hksat3d,  &
       sucsat3d, &
       clay3d,   &
       sand3d
  INTEGER  :: n,i,j,k,ib,lev,bottom
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: bd2d
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: tkm2d
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: xksat2d
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: depthratio2d
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: clay2d
  REAL(r8),DIMENSION(ims:ime,jms:jme )    :: sand2d

  REAL(r8)                 :: scalez  = 0.025_r8
  LOGICAL,PARAMETER        :: arbinit = .TRUE.
  REAL(r8),PARAMETER           :: defval  = -999.0
  INTEGER                  :: isl
  INTEGER                  :: numb_lak
  CHARACTER*256 :: message

  IF ( RESTART ) RETURN

  DO j = jts,jte
     DO i = its,ite
        snowdp2d(i,j)         = snow(i,j)*0.005               ! SNOW in kg/m^2 and snowdp in m
        h2osno2d(i,j)         = snow(i,j) ! mm
     ENDDO
  ENDDO

  ! initialize all the grid with default value
  DO j = jts,jte
     DO i = its,ite

        lakedepth2d(i,j)             = defval
        snl2d(i,j)                   = defval
        DO k = -nlevsnow+1,nlevsoil
           h2osoi_liq3d(i,k,j)      = defval
           h2osoi_ice3d(i,k,j)      = defval
           t_soisno3d(i,k,j)        = defval
           z3d(i,k,j)               = defval
           dz3d(i,k,j)              = defval
        ENDDO
        DO k = 1,nlevlake
           t_lake3d(i,k,j)          = defval
           lake_icefrac3d(i,k,j)    = defval
           z_lake3d(i,k,j)          = defval
           dz_lake3d(i,k,j)         = defval
        ENDDO

     ENDDO
  ENDDO

  ! judge whether the grid is lake grid
  numb_lak = 0
  DO i=its,ite
     DO j=jts,jte
        !#if (EM_CORE==1)
        IF (lakeflag.EQ.0) THEN
           IF(ht(i,j)>=lake_min_elev) THEN
              IF ( xice(i,j).GT.xice_threshold) THEN   !mchen
                 ivgtyp(i,j) = iswater
                 xland(i,j) = 2.
                 lake_icefrac3d(i,1,j) = xice(i,j)
                 xice(i,j)=0.0
              ENDIF
           ENDIF

           IF(ivgtyp(i,j)==iswater.AND.ht(i,j)>=lake_min_elev) THEN
              lake(i,j)  = .TRUE.
              lakemask(i,j) = 1
              numb_lak   = numb_lak + 1
           ELSE
              lake(i,j)  = .FALSE.
              lakemask(i,j) = 0
           END IF
        ELSE
           IF(lakemask(i,j).EQ.1) THEN
              lake(i,j)  = .TRUE.
              numb_lak   = numb_lak + 1
              IF ( xice(i,j).GT.xice_threshold) THEN   !mchen
                 ivgtyp(i,j) = iswater
                 xland(i,j) = 2.
                 lake_icefrac3d(i,1,j) = xice(i,j)
                 xice(i,j)=0.0
              ENDIF
           ELSE
              lake(i,j)  = .FALSE.
           ENDIF
        ENDIF   ! end if lakeflag=0

        !!           lake = .true.   !wfs

        !#else
        !            if(ht(i,j)>=lake_min_elev) then
        !              if ( xice(i,j).gt.xice_threshold) then   !mchen
        !                   ivgtyp(i,j) = iswater
        !                   xland(i,j) = 2.
        !                   lake_icefrac3d(i,1,j) = xice(i,j)
        !                   xice(i,j)=0.0
        !               endif
        !            endif
        !            if(ivgtyp(i,j)==iswater.and.ht(i,j)>=lake_min_elev) then
        !                lake(i,j)  = .true.
        !                numb_lak   = numb_lak + 1
        !            else
        !                lake(i,j)  = .false.
        !            end if
        !
        !#endif
     END DO
  END DO
  !    write(message,*) "the total number of lake grid is :", numb_lak
  !    CALL wrf_message(message)
  !    CALL LakeDebug(msg)
  ! initialize lake grid

  DO j = jts,jte
     DO i = its,ite

        IF ( lake(i,j) ) THEN

           !	t_soisno3d(i,:,j)      = tsk(i,j)
           !        t_lake3d(i,:,j)        = tsk(i,j)
           !        t_grnd2d(i,j)          = tsk(i,j)

           z3d(i,:,j)             = 0.0
           dz3d(i,:,j)            = 0.0
           zi3d(i,:,j)            = 0.0
           h2osoi_liq3d(i,:,j)    = 0.0
           h2osoi_ice3d(i,:,j)    = 0.0
           lake_icefrac3d(i,:,j)  = 0.0
           h2osoi_vol3d(i,:,j)    = 0.0
           snl2d(i,j)             = 0.0
           !          if ( use_lakedepth.eq.1 .and.lake_depth_flag.eq.0 ) then !mchen
           !          call wrf_error_fatal ( 'STOP: You need lake-depth information. Rerun WPS or set use_lakedepth = 0')
           !          end if
           IF ( use_lakedepth.EQ.0 .AND.lake_depth_flag.EQ.1 ) THEN !mchen
              lake_depth_flag = 0
           END IF
           IF ( lake_depth_flag.EQ.1 ) THEN

              IF (lake_depth(i,j) > 0.0) THEN
                 lakedepth2d(i,j)   = lake_depth(i,j)
              ELSE
                 IF ( lakedepth_default  > 0.0 ) THEN
                    lakedepth2d(i,j)   = lakedepth_default
                 ELSE
                    lakedepth2d(i,j)   = spval
                 ENDIF
              ENDIF

           ELSE
              IF ( lakedepth_default  > 0.0 ) THEN
                 lakedepth2d(i,j)   = lakedepth_default
              ELSE
                 lakedepth2d(i,j)   = spval
              ENDIF
           ENDIF
        ENDIF

     ENDDO
  ENDDO


  !#ifndef EXTRALAKELAYERS
  !  dzlak(1) = 0.1_r8
  !  dzlak(2) = 1._r8
  !  dzlak(3) = 2._r8
  !  dzlak(4) = 3._r8
  !  dzlak(5) = 4._r8
  !  dzlak(6) = 5._r8
  !  dzlak(7) = 7._r8
  !  dzlak(8) = 7._r8
  !  dzlak(9) = 10.45_r8
  !  dzlak(10)= 10.45_r8
  !
  !  zlak(1) =  0.05_r8
  !  zlak(2) =  0.6_r8
  !  zlak(3) =  2.1_r8
  !  zlak(4) =  4.6_r8
  !  zlak(5) =  8.1_r8
  !  zlak(6) = 12.6_r8
  !  zlak(7) = 18.6_r8
  !  zlak(8) = 25.6_r8
  !  zlak(9) = 34.325_r8
  !  zlak(10)= 44.775_r8
  dzlak(1) = 0.1_r8
  dzlak(2) = 0.1_r8
  dzlak(3) = 0.1_r8
  dzlak(4) = 0.1_r8
  dzlak(5) = 0.1_r8
  dzlak(6) = 0.1_r8
  dzlak(7) = 0.1_r8
  dzlak(8) = 0.1_r8
  dzlak(9) = 0.1_r8
  dzlak(10)= 0.1_r8

  zlak(1) =  0.05_r8
  zlak(2) =  0.15_r8
  zlak(3) =  0.25_r8
  zlak(4) =  0.35_r8
  zlak(5) =  0.45_r8
  zlak(6) = 0.55_r8
  zlak(7) = 0.65_r8
  zlak(8) = 0.75_r8
  zlak(9) = 0.85_r8
  zlak(10)= 0.95_r8
  !#else
  !  dzlak(1) =0.1_r8
  !  dzlak(2) =0.25_r8
  !  dzlak(3) =0.25_r8
  !  dzlak(4) =0.25_r8
  !  dzlak(5) =0.25_r8
  !  dzlak(6) =0.5_r8
  !  dzlak(7) =0.5_r8
  !  dzlak(8) =0.5_r8
  !  dzlak(9) =0.5_r8
  !  dzlak(10) =0.75_r8
  !  dzlak(11) =0.75_r8
  !  dzlak(12) =0.75_r8
  !  dzlak(13) =0.75_r8
  !  dzlak(14) =2_r8
  !  dzlak(15) =2_r8
  !  dzlak(16) =2.5_r8
  !  dzlak(17) =2.5_r8
  !  dzlak(18) =3.5_r8
  !  dzlak(19) =3.5_r8
  !  dzlak(20) =3.5_r8
  !  dzlak(21) =3.5_r8
  !  dzlak(22) =5.225_r8
  !  dzlak(23) =5.225_r8
  !  dzlak(24) =5.225_r8
  !  dzlak(25) =5.225_r8
  !
  !  zlak(1) = dzlak(1)/2._r8
  !  do k = 2,10
  !     zlak(k) = zlak(k-1) + (dzlak(k-1)+dzlak(k))/2._r8
  !  end do
  !#endif

  ! "0" refers to soil surface and "10" refers to the bottom of model soil

  DO j = 1, nlevsoil
     zsoi(j) = scalez*(EXP(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
  ENDDO

  dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
  DO j = 2,nlevsoil-1
     dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
  ENDDO
  dzsoi(nlevsoil) = zsoi(nlevsoil)-zsoi(nlevsoil-1)

  zisoi(0) = 0._r8
  DO j = 1, nlevsoil-1
     zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
  ENDDO
  zisoi(nlevsoil) = zsoi(nlevsoil) + 0.5_r8*dzsoi(nlevsoil)


!!!!!!!!!!!!!!!!!!begin to initialize lake variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO j = jts,jte
     DO i = its,ite

        IF ( lake(i,j) ) THEN

           ! Soil hydraulic and thermal properties
           isl = ISLTYP(i,j)
           IF (isl == 14 ) isl = isl + 1
           DO k = 1,nlevsoil
              sand3d(i,k,j)  = sand(isl)
              clay3d(i,k,j)  = clay(isl)
           ENDDO

           DO k = 1,nlevsoil
              clay2d(i,j) = clay3d(i,k,j)
              sand2d(i,j) = sand3d(i,k,j)
              watsat3d(i,k,j) = 0.489_r8 - 0.00126_r8*sand2d(i,j)
              bd2d(i,j)    = (1._r8-watsat3d(i,k,j))*2.7e3_r8
              xksat2d(i,j) = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand2d(i,j)) ) ! mm/s
              tkm2d(i,j) = (8.80_r8*sand2d(i,j)+2.92_r8*clay2d(i,j))/(sand2d(i,j)+clay2d(i,j))          ! W/(m K)

              bsw3d(i,k,j) = 2.91_r8 + 0.159_r8*clay2d(i,j)
              bsw23d(i,k,j) = -(3.10_r8 + 0.157_r8*clay2d(i,j) - 0.003_r8*sand2d(i,j))
              psisat3d(i,k,j) = -(EXP((1.54_r8 - 0.0095_r8*sand2d(i,j) + 0.0063_r8*(100.0_r8-sand2d(i,j)  &
                   -clay2d(i,j)))*LOG(10.0_r8))*9.8e-5_r8)
              vwcsat3d(i,k,j) = (50.5_r8 - 0.142_r8*sand2d(i,j) - 0.037_r8*clay2d(i,j))/100.0_r8
              hksat3d(i,k,j) = xksat2d(i,j)
              sucsat3d(i,k,j) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand2d(i,j)) )
              tkmg3d(i,k,j) = tkm2d(i,j) ** (1._r8- watsat3d(i,k,j))
              tksatu3d(i,k,j) = tkmg3d(i,k,j)*0.57_r8**watsat3d(i,k,j)
              tkdry3d(i,k,j) = (0.135_r8*bd2d(i,j) + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd2d(i,j))
              csol3d(i,k,j) = (2.128_r8*sand2d(i,j)+2.385_r8*clay2d(i,j)) / (sand2d(i,j)+clay2d(i,j))*1.e6_r8  ! J/(m3 K)
              watdry3d(i,k,j) = watsat3d(i,k,j) * (316230._r8/sucsat3d(i,k,j)) ** (-1._r8/bsw3d(i,k,j))
              watopt3d(i,k,j) = watsat3d(i,k,j) * (158490._r8/sucsat3d(i,k,j)) ** (-1._r8/bsw3d(i,k,j))
           END DO
           IF (lakedepth2d(i,j) == spval) THEN
              lakedepth2d(i,j) = zlak(nlevlake) + 0.5_r8*dzlak(nlevlake)
              z_lake3d(i,1:nlevlake,j) = zlak(1:nlevlake)
              dz_lake3d(i,1:nlevlake,j) = dzlak(1:nlevlake)
           ELSE
              depthratio2d(i,j) = lakedepth2d(i,j) / (zlak(nlevlake) + 0.5_r8*dzlak(nlevlake))
              z_lake3d(i,1,j) = zlak(1)
              dz_lake3d(i,1,j) = dzlak(1)
              dz_lake3d(i,2:nlevlake,j) = dzlak(2:nlevlake)*depthratio2d(i,j)
              z_lake3d(i,2:nlevlake,j) = zlak(2:nlevlake)*depthratio2d(i,j) + dz_lake3d(i,1,j)*(1._r8 - depthratio2d(i,j))
           END IF
           ! initial t_lake3d here
           t_soisno3d(i,1,j)      = tsk(i,j)
           t_lake3d(i,1,j)        = tsk(i,j)
           t_grnd2d(i,j)          = 277.0
           DO k = 2, nlevlake
              IF(z_lake3d(i,k,j).LE.depth_c) THEN
                 t_soisno3d(i,k,j)=tsk(i,j)+(277.0-tsk(i,j))/depth_c*z_lake3d(i,k,j)
                 t_lake3d(i,k,j)=tsk(i,j)+(277.0-tsk(i,j))/depth_c*z_lake3d(i,k,j)
              ELSE
                 t_soisno3d(i,k,j)      = 277.0
                 t_lake3d(i,k,j)        = 277.0
              END IF
           ENDDO
           !end initial t_lake3d here
           z3d(i,1:nlevsoil,j) = zsoi(1:nlevsoil)
           zi3d(i,0:nlevsoil,j) = zisoi(0:nlevsoil)
           dz3d(i,1:nlevsoil,j) = dzsoi(1:nlevsoil)
           savedtke12d(i,j) = tkwat ! Initialize for first timestep.


           IF (snowdp2d(i,j) < 0.01_r8) THEN
              snl2d(i,j) = 0
              dz3d(i,-nlevsnow+1:0,j) = 0._r8
              z3d (i,-nlevsnow+1:0,j) = 0._r8
              zi3d(i,-nlevsnow+0:0,j) = 0._r8
           ELSE
              IF ((snowdp2d(i,j) >= 0.01_r8) .AND. (snowdp2d(i,j) <= 0.03_r8)) THEN
                 snl2d(i,j) = -1
                 dz3d(i,0,j)  = snowdp2d(i,j)
              ELSE IF ((snowdp2d(i,j) > 0.03_r8) .AND. (snowdp2d(i,j) <= 0.04_r8)) THEN
                 snl2d(i,j) = -2
                 dz3d(i,-1,j) = snowdp2d(i,j)/2._r8
                 dz3d(i, 0,j) = dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.04_r8) .AND. (snowdp2d(i,j) <= 0.07_r8)) THEN
                 snl2d(i,j) = -2
                 dz3d(i,-1,j) = 0.02_r8
                 dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.07_r8) .AND. (snowdp2d(i,j) <= 0.12_r8)) THEN
                 snl2d(i,j) = -3
                 dz3d(i,-2,j) = 0.02_r8
                 dz3d(i,-1,j) = (snowdp2d(i,j) - 0.02_r8)/2._r8
                 dz3d(i, 0,j) = dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.12_r8) .AND. (snowdp2d(i,j) <= 0.18_r8)) THEN
                 snl2d(i,j) = -3
                 dz3d(i,-2,j) = 0.02_r8
                 dz3d(i,-1,j) = 0.05_r8
                 dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-2,j) - dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.18_r8) .AND. (snowdp2d(i,j) <= 0.29_r8)) THEN
                 snl2d(i,j) = -4
                 dz3d(i,-3,j) = 0.02_r8
                 dz3d(i,-2,j) = 0.05_r8
                 dz3d(i,-1,j) = (snowdp2d(i,j) - dz3d(i,-3,j) - dz3d(i,-2,j))/2._r8
                 dz3d(i, 0,j) = dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.29_r8) .AND. (snowdp2d(i,j) <= 0.41_r8)) THEN
                 snl2d(i,j) = -4
                 dz3d(i,-3,j) = 0.02_r8
                 dz3d(i,-2,j) = 0.05_r8
                 dz3d(i,-1,j) = 0.11_r8
                 dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-3,j) - dz3d(i,-2,j) - dz3d(i,-1,j)
              ELSE IF ((snowdp2d(i,j) > 0.41_r8) .AND. (snowdp2d(i,j) <= 0.64_r8)) THEN
                 snl2d(i,j) = -5
                 dz3d(i,-4,j) = 0.02_r8
                 dz3d(i,-3,j) = 0.05_r8
                 dz3d(i,-2,j) = 0.11_r8
                 dz3d(i,-1,j) = (snowdp2d(i,j) - dz3d(i,-4,j) - dz3d(i,-3,j) - dz3d(i,-2,j))/2._r8
                 dz3d(i, 0,j) = dz3d(i,-1,j)
              ELSE IF (snowdp2d(i,j) > 0.64_r8) THEN
                 snl2d(i,j) = -5
                 dz3d(i,-4,j) = 0.02_r8
                 dz3d(i,-3,j) = 0.05_r8
                 dz3d(i,-2,j) = 0.11_r8
                 dz3d(i,-1,j) = 0.23_r8
                 dz3d(i, 0,j)=snowdp2d(i,j)-dz3d(i,-4,j)-dz3d(i,-3,j)-dz3d(i,-2,j)-dz3d(i,-1,j)
              ENDIF
           END IF

           DO k = 0, snl2d(i,j)+1, -1
              z3d(i,k,j)    = zi3d(i,k,j) - 0.5_r8*dz3d(i,k,j)
              zi3d(i,k-1,j) = zi3d(i,k,j) - dz3d(i,k,j)
           END DO

           ! 3:subroutine makearbinit

           IF (snl2d(i,j) < 0) THEN
              DO k = snl2d(i,j)+1, 0
                 ! Be careful because there may be new snow layers with bad temperatures like 0 even if
                 ! coming from init. con. file.
                 IF(arbinit .OR. t_soisno3d(i,k,j) > 300 .OR. t_soisno3d(i,k,j) < 200) t_soisno3d(i,k,j) = 250._r8
              ENDDO
           END IF

           DO k = 1, nlevsoil
              IF(arbinit .OR. t_soisno3d(i,k,j) > 1000 .OR. t_soisno3d(i,k,j) < 0) t_soisno3d(i,k,j) = t_lake3d(i,nlevlake,j)
           END DO

           DO k = 1, nlevlake
              IF(arbinit .OR. lake_icefrac3d(i,k,j) > 1._r8 .OR. lake_icefrac3d(i,k,j) < 0._r8) THEN
                 IF(t_lake3d(i,k,j) >= tfrz) THEN
                    lake_icefrac3d(i,k,j) = 0._r8
                 ELSE
                    lake_icefrac3d(i,k,j) = 1._r8
                 END IF
              END IF
           END DO

           DO k = 1,nlevsoil
              IF (arbinit .OR. h2osoi_vol3d(i,k,j) > 10._r8 .OR. h2osoi_vol3d(i,k,j) < 0._r8) h2osoi_vol3d(i,k,j) = 1.0_r8
              h2osoi_vol3d(i,k,j) = MIN(h2osoi_vol3d(i,k,j),watsat3d(i,k,j))

              ! soil layers
              IF (t_soisno3d(i,k,j) <= tfrz) THEN
                 h2osoi_ice3d(i,k,j)  = dz3d(i,k,j)*denice*h2osoi_vol3d(i,k,j)
                 h2osoi_liq3d(i,k,j) = 0._r8
              ELSE
                 h2osoi_ice3d(i,k,j) = 0._r8
                 h2osoi_liq3d(i,k,j) = dz3d(i,k,j)*denh2o*h2osoi_vol3d(i,k,j)
              ENDIF
           ENDDO

           DO k = -nlevsnow+1, 0
              IF (k > snl2d(i,j)) THEN
                 h2osoi_ice3d(i,k,j) = dz3d(i,k,j)*bdsno
                 h2osoi_liq3d(i,k,j) = 0._r8
              END IF
           END DO

        END IF   !lake(i,j)
     ENDDO
  ENDDO
  ! print *,'good9'
END SUBROUTINE lakeini

SUBROUTINE SoilThermProp_Lake (watsat,tksatu,tkmg,tkdry,csol, &
      snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
      tk, cv, tktopsoillay)


  !
  ! !DESCRIPTION:
  ! Calculation of thermal conductivities and heat capacities of
  ! snow/soil layers
  ! (1) The volumetric heat capacity is calculated as a linear combination
  !     in terms of the volumetric fraction of the constituent phases.
  !
  ! (2) The thermal conductivity of soil is computed from the algorithm of
  !     Johansen (as reported by Farouki 1981), and of snow is from the
  !     formulation used in SNTHERM (Jordan 1991).
  ! The thermal conductivities at the interfaces between two neighboring
  ! layers (j, j+1) are derived from an assumption that the flux across
  ! the interface is equal to that from the node j to the interface and the
  ! flux from the interface to the node j+1.
  !
  ! For lakes, the proper soil layers (not snow) should always be saturated.
  !
  ! !USES:

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER , INTENT(in) :: snl(1)
  REAL(r8), INTENT(in) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: zi(1,-nlevsnow+0:nlevsoil)
  REAL(r8), INTENT(in) :: z(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)

  REAL(r8), INTENT(out) :: cv(lbc:ubc,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: tk(lbc:ubc,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: tktopsoillay(1:1)

  INTEGER  :: l,c,j
  INTEGER  :: fc
  REAL(r8) :: bw
  REAL(r8) :: dksat
  REAL(r8) :: dke
  REAL(r8) :: fl
  REAL(r8) :: satw
  REAL(r8) :: thk(lbc:ubc,-nlevsnow+1:nlevsoil)
  CHARACTER*256 :: message

  ! Thermal conductivity of soil from Farouki (1981)
  ! PRINT*, 'components of cv(c,j)', watsat,tksatu,tkmg,tkdry,csol

  DO j = -nlevsnow+1,nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        ! Only examine levels from 1->nlevsoil
        IF (j >= 1) THEN
           !             l = clandunit(c)
           !             if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
           ! This could be altered later for allowing this to be over glaciers.

           ! Soil should be saturated.
           !#define LAKEDEBUG
           !#if (defined LAKEDEBUG)
           satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
           !                satw = min(1._r8, satw)
           !                if (satw < 0.999_r8) then
           !                   write(message,*)'WARNING: soil layer unsaturated in SoilThermProp_Lake, satw, j = ', satw, j
           !                   CALL wrf_error_fatal ( message )
           !                end if
           ! Could use denice because if it starts out frozen, the volume of water will go below sat.,
           ! since we're not yet doing excess ice.
           ! But take care of this in HydrologyLake.
           !#endif
           satw = 1._r8
           fl = h2osoi_liq(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j))
           IF (t_soisno(c,j) >= tfrz) THEN       ! Unfrozen soil
              dke = MAX(0._r8, LOG10(satw) + 1.0_r8)
              dksat = tksatu(c,j)
           ELSE                               ! Frozen soil
              dke = satw
              dksat = tkmg(c,j)*0.249_r8**(fl*watsat(c,j))*2.29_r8**watsat(c,j)
           ENDIF
           thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
           !             else
           !                thk(c,j) = tkwat
           !                if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
           !             endif
        ENDIF

        ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
        ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
        IF (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) THEN
           bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dz(c,j)
           thk(c,j) = tkairc + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkairc)
        END IF

     END DO
  END DO
  ! PRINT*, 'thk',thk
  ! Thermal conductivity at the layer interface

  ! Have to correct for the fact that bottom snow layer and top soil layer border lake.
  ! For the first case, the snow layer conductivity for the middle of the layer will be returned.
  ! Because the interfaces are below the soil layers, the conductivity for the top soil layer
  ! will have to be returned separately.
  DO j = -nlevsnow+1,nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        IF (j >= snl(c)+1 .AND. j <= nlevsoil-1 .AND. j /= 0) THEN
           tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
        ELSE IF (j == 0) THEN
           tk(c,j) = thk(c,j)
        ELSE IF (j == nlevsoil) THEN
           tk(c,j) = 0._r8
        END IF
        ! For top soil layer.
        IF (j == 1) tktopsoillay(c) = thk(c,j)
     END DO
  END DO

  ! Soil heat capacity, from de Vires (1963)
  DO j = 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        !          l = clandunit(c)
        !          if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
        cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) +   &
             (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
        !          else
        !             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
        !          endif
        !          if (j == 1) then
        !             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
        !                cv(c,j) = cv(c,j) + cpice*h2osno(c)
        !             end if
        !          end if
        ! Won't worry about heat capacity for thin snow on lake with no snow layers.
     ENDDO
  END DO

  ! Snow heat capacity

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,num_shlakec
        c = filter_shlakec(fc)
        IF (snl(c)+1 < 1 .AND. j >= snl(c)+1) THEN
           cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
        END IF
     END DO
  END DO

  !print *,'thk=',thk
  !print *,'tk=',tk

END SUBROUTINE SoilThermProp_Lake

SUBROUTINE ShalLakeHydrology(dz_lake,forc_rain,forc_snow,                      & !i  qflx_evap_tot,forc_t,do_capsnow,            &
     qflx_evap_tot,forc_t,do_capsnow,            &
     t_grnd,qflx_evap_soi,                             &
     qflx_snomelt,imelt,frac_iceold,                   & !i add by guhp
     z,dz,zi,snl,h2osno,snowdp,lake_icefrac,t_lake,      & !i&o
     endwb,snowage,snowice,snowliq,t_snow,             & !o
     t_soisno,h2osoi_ice,h2osoi_liq,h2osoi_vol,        &
     qflx_drain,qflx_surf,qflx_infl,qflx_qrgwl,        &
     qcharge,qflx_prec_grnd,qflx_snowcap,              &
     qflx_snowcap_col,qflx_snow_grnd_pft,              &
     qflx_snow_grnd_col,qflx_rain_grnd,                &
     qflx_evap_tot_col,soilalpha,zwt,fcov,             &
     rootr_column,qflx_evap_grnd,qflx_sub_snow,        &
     qflx_dew_snow,qflx_dew_grnd,qflx_rain_grnd_col)

  !==================================================================================
  ! !DESCRIPTION:
  ! Calculation of Shallow Lake Hydrology. Full hydrology of snow layers is
  ! done. However, there is no infiltration, and the water budget is balanced with
  ! qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained at
  ! volumetric saturation if ice melting frees up pore space. Likewise, if the water
  ! portion alone at some point exceeds pore capacity, it is reduced. This is consistent
  ! with the possibility of initializing the soil layer with excess ice. The only
  ! real error with that is that the thermal conductivity will ignore the excess ice
  ! (and accompanying thickness change).
  !
  ! If snow layers are present over an unfrozen lake, and the top layer of the lake
  ! is capable of absorbing the latent heat without going below freezing,
  ! the snow-water is runoff and the latent heat is subtracted from the lake.
  !
  ! WARNING: This subroutine assumes lake columns have one and only one pft.
  !
  ! Sequence is:
  !  ShalLakeHydrology:
  !    Do needed tasks from Hydrology1, Biogeophysics2, & top of Hydrology2.
  !    -> SnowWater:             change of snow mass and snow water onto soil
  !    -> SnowCompaction:        compaction of snow layers
  !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
  !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
  !    Add water to soil if melting has left it with open pore space.
  !    Cleanup and do water balance.
  !    If snow layers are found above a lake with unfrozen top layer, whose top
  !    layer has enough heat to melt all the snow ice without freezing, do so
  !    and eliminate the snow layers.
  !
  ! !REVISION HISTORY:
  ! Created by Zack Subin, 2009
  !
  !============================================================================================

  ! USES:
  !
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in) :: dz_lake(1,nlevlake)
  REAL(r8), INTENT(in) :: forc_rain(1)
  REAL(r8), INTENT(in) :: forc_snow(1)
  REAL(r8), INTENT(in) :: qflx_evap_tot(1)
  REAL(r8), INTENT(in) :: forc_t(1)
  LOGICAL , INTENT(in) :: do_capsnow(1)
  REAL(r8), INTENT(in) :: t_grnd(1)
  REAL(r8), INTENT(in) :: qflx_evap_soi(1)
  REAL(r8), INTENT(in) :: qflx_snomelt(1)
  INTEGER,  INTENT(in) :: imelt(1,-nlevsnow+1:nlevsoil)

  ! REAL(r8), INTENT(inout) :: begwb(1)
  REAL(r8) :: begwb(1) =0. ! change to a LOCAL varialbe


  REAL(r8), INTENT(inout) :: z(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: zi(1,-nlevsnow+0:nlevsoil)
  INTEGER , INTENT(inout) :: snl(1)
  REAL(r8), INTENT(inout) :: h2osno(1)
  REAL(r8), INTENT(inout) :: snowdp(1)
  REAL(r8), INTENT(inout) :: lake_icefrac(1,nlevlake)
  REAL(r8), INTENT(inout) :: t_lake(1,nlevlake)

  REAL(r8), INTENT(inout) :: frac_iceold(1,-nlevsnow+1:nlevsoil)


  REAL(r8), INTENT(out) :: endwb(1)
  REAL(r8), INTENT(out) :: snowage(1)
  REAL(r8), INTENT(out) :: snowice(1)
  REAL(r8), INTENT(out) :: snowliq(1)
  REAL(r8), INTENT(out) :: t_snow(1)
  REAL(r8), INTENT(out) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: h2osoi_vol(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(out) :: qflx_drain(1)
  REAL(r8), INTENT(out) :: qflx_surf(1)
  REAL(r8), INTENT(out) :: qflx_infl(1)
  REAL(r8), INTENT(out) :: qflx_qrgwl(1)
  REAL(r8), INTENT(out) :: qcharge(1)
  REAL(r8), INTENT(out) :: qflx_prec_grnd(1)
  REAL(r8), INTENT(out) :: qflx_snowcap(1)
  REAL(r8), INTENT(out) :: qflx_snowcap_col(1)
  REAL(r8), INTENT(out) :: qflx_snow_grnd_pft(1)
  REAL(r8), INTENT(out) :: qflx_snow_grnd_col(1)
  REAL(r8), INTENT(out) :: qflx_rain_grnd(1)
  REAL(r8), INTENT(out) :: qflx_evap_tot_col(1)
  REAL(r8) ,INTENT(out) :: soilalpha(1)
  REAL(r8), INTENT(out) :: zwt(1)
  REAL(r8), INTENT(out) :: fcov(1)
  REAL(r8), INTENT(out) :: rootr_column(1,1:nlevsoil)
  REAL(r8), INTENT(out) :: qflx_evap_grnd(1)
  REAL(r8), INTENT(out) :: qflx_sub_snow(1)
  REAL(r8), INTENT(out) :: qflx_dew_snow(1)
  REAL(r8), INTENT(out) :: qflx_dew_grnd(1)
  REAL(r8), INTENT(out) :: qflx_rain_grnd_col(1)
  REAL(r8) :: watdry
  REAL(r8) :: rwat(lbc:ubc)
  REAL(r8) :: swat(lbc:ubc)
  REAL(r8) :: rz(lbc:ubc)
  REAL(r8) :: tsw
  REAL(r8) :: stsw

  INTEGER  :: p,fp,g,l,c,j,fc,jtop
  INTEGER  :: num_shlakesnowc
  INTEGER  :: filter_shlakesnowc(ubc-lbc+1)
  INTEGER  :: num_shlakenosnowc
  INTEGER  :: filter_shlakenosnowc(ubc-lbc+1)

  INTEGER  :: newnode
  REAL(r8) :: dz_snowf
  REAL(r8) :: bifall
  REAL(r8) :: fracsnow(lbp:ubp)
  REAL(r8) :: fracrain(lbp:ubp)
  REAL(r8) :: qflx_prec_grnd_snow(lbp:ubp)
  REAL(r8) :: qflx_prec_grnd_rain(lbp:ubp)
  REAL(r8) :: qflx_evap_soi_lim
  REAL(r8) :: h2osno_temp
  REAL(r8), PARAMETER :: snow_bd = 250._r8
  REAL(r8) :: sumsnowice(lbc:ubc)
  LOGICAL  :: unfrozen(lbc:ubc)
  REAL(r8) :: heatrem
  REAL(r8) :: heatsum(lbc:ubc)
  REAL(r8) :: qflx_top_soil(1)
  CHARACTER*256 :: message

  REAL(r8) :: snow_water(lbc:ubc)

  ! Determine step size

  !    dtime = get_step_size()

  ! Add soil water to water balance.
  DO j = 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
     END DO
  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Do precipitation onto ground, etc., from Hydrology1.

  !dir$ concurrent
  !cdir nodep
  DO fp = 1, num_shlakep
     p = filter_shlakep(fp)
     g = pgridcell(p)
     !       l = plandunit(p)
     c = pcolumn(p)

     ! Precipitation onto ground (kg/(m2 s))
     !       ! PET, 1/18/2005: Added new terms for mass balance correction
     !       ! due to dynamic pft weight shifting (column-level h2ocan_loss)
     !       ! Because the fractionation between rain and snow is indeterminate if
     !       ! rain + snow = 0, I am adding this very small flux only to the rain
     !       ! components.
     ! Not relevant unless PFTs are added to lake later.
     !       if (frac_veg_nosno(p) == 0) then
     qflx_prec_grnd_snow(p) = forc_snow(g)
     qflx_prec_grnd_rain(p) = forc_rain(g) !+ h2ocan_loss(c)
     !       else
     !          qflx_prec_grnd_snow(p) = qflx_through_snow(p) + (qflx_candrip(p) * fracsnow(p))
     !          qflx_prec_grnd_rain(p) = qflx_through_rain(p) + (qflx_candrip(p) * fracrain(p)) + h2ocan_loss(c)
     !       end if
     qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

     IF (do_capsnow(c)) THEN
        qflx_snowcap(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)
        qflx_snow_grnd_pft(p) = 0._r8
        qflx_rain_grnd(p) = 0._r8
     ELSE
        qflx_snowcap(p) = 0._r8
        !#define OFFLINE
        !#if (defined OFFLINE)
        !          qflx_snow_grnd_pft(p) = qflx_prec_grnd(p)*(1._r8-flfall(g)) ! ice onto ground (mm/s)
        !          qflx_rain_grnd(p)     = qflx_prec_grnd(p)*flfall(g)      ! liquid water onto ground (mm/s)
        !#else
        !          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
        !          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
        !#endif
     END IF
     ! Assuming one PFT; needed for below
     qflx_snow_grnd_col(c) = qflx_snow_grnd_pft(p)
     qflx_rain_grnd_col(c) = qflx_rain_grnd(p)

  END DO ! (end pft loop)

  ! Determine snow height and snow water

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)
     !       l = clandunit(c)
     g = cgridcell(c)

     ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
     ! U.S.Department of Agriculture Forest Service, Project F,
     ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

     IF (do_capsnow(c)) THEN
        dz_snowf = 0._r8
     ELSE
        IF (forc_t(g) > tfrz + 2._r8) THEN
           bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
        ELSE IF (forc_t(g) > tfrz - 15._r8) THEN
           bifall=50._r8 + 1.7_r8*(forc_t(g) - tfrz + 15._r8)**1.5_r8
        ELSE
           bifall=50._r8
        END IF
        dz_snowf = qflx_snow_grnd_col(c)/bifall
        snowdp(c) = snowdp(c) + dz_snowf*dtime
        h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)
     END IF

     !       if (itype(l)==istwet .and. t_grnd(c)>tfrz) then
     !          h2osno(c)=0._r8
     !          snowdp(c)=0._r8
     !          snowage(c)=0._r8
     !       end if
     ! Take care of this later in function.

     ! When the snow accumulation exceeds 10 mm, initialize snow layer
     ! Currently, the water temperature for the precipitation is simply set
     ! as the surface air temperature

     newnode = 0    ! flag for when snow node will be initialized
     IF (snl(c) == 0 .AND. qflx_snow_grnd_col(c) > 0.0_r8 .AND. snowdp(c) >= 0.01_r8) THEN
        newnode = 1
        snl(c) = -1
        dz(c,0) = snowdp(c)                       ! meter
        z(c,0) = -0.5_r8*dz(c,0)
        zi(c,-1) = -dz(c,0)
        snowage(c) = 0._r8                        ! snow age
        t_soisno(c,0) = MIN(tfrz, forc_t(g))      ! K
        h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
        h2osoi_liq(c,0) = 0._r8                   ! kg/m2
        frac_iceold(c,0) = 1._r8
     END IF

     ! The change of ice partial density of surface node due to precipitation.
     ! Only ice part of snowfall is added here, the liquid part will be added
     ! later.

     IF (snl(c) < 0 .AND. newnode == 0) THEN
        h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
        dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
     END IF

  END DO

  ! Calculate sublimation and dew, adapted from HydrologyLake and Biogeophysics2.

  !dir$ concurrent
  !cdir nodep
  DO fp = 1,num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     jtop = snl(c)+1

     ! Use column variables here
     qflx_evap_grnd(c) = 0._r8
     qflx_sub_snow(c) = 0._r8
     qflx_dew_snow(c) = 0._r8
     qflx_dew_grnd(c) = 0._r8

     IF (jtop <= 0) THEN ! snow layers
        j = jtop
        ! Assign ground evaporation to sublimation from soil ice or to dew
        ! on snow or ground

        IF (qflx_evap_soi(p) >= 0._r8) THEN
           ! for evaporation partitioning between liquid evap and ice sublimation,
           ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
           ! Since we're not limiting evap over lakes, but still can't remove more from top
           ! snow layer than there is there, create temp. limited evap_soi.
           qflx_evap_soi_lim = MIN(qflx_evap_soi(p), (h2osoi_liq(c,j)+h2osoi_ice(c,j))/dtime)
           IF ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) THEN
              qflx_evap_grnd(c) = MAX(qflx_evap_soi_lim*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
           ELSE
              qflx_evap_grnd(c) = 0._r8
           END IF
           qflx_sub_snow(c) = qflx_evap_soi_lim - qflx_evap_grnd(c)
        ELSE
           IF (t_grnd(c) < tfrz) THEN
              qflx_dew_snow(c) = ABS(qflx_evap_soi(p))
           ELSE
              qflx_dew_grnd(c) = ABS(qflx_evap_soi(p))
           END IF
        END IF
        ! Update the pft-level qflx_snowcap
        ! This was moved in from Hydrology2 to keep all pft-level
        ! calculations out of Hydrology2
        IF (do_capsnow(c)) qflx_snowcap(p) = qflx_snowcap(p) + qflx_dew_snow(c) + qflx_dew_grnd(c)

     ELSE ! No snow layers: do as in HydrologyLake but with actual clmtype variables
        IF (qflx_evap_soi(p) >= 0._r8) THEN
           ! Sublimation: do not allow for more sublimation than there is snow
           ! after melt.  Remaining surface evaporation used for infiltration.
           qflx_sub_snow(c) = MIN(qflx_evap_soi(p), h2osno(c)/dtime)
           qflx_evap_grnd(c) = qflx_evap_soi(p) - qflx_sub_snow(c)
        ELSE
           IF (t_grnd(c) < tfrz-0.1_r8) THEN
              qflx_dew_snow(c) = ABS(qflx_evap_soi(p))
           ELSE
              qflx_dew_grnd(c) = ABS(qflx_evap_soi(p))
           END IF
        END IF

        ! Update snow pack for dew & sub.
        h2osno_temp = h2osno(c)
        IF (do_capsnow(c)) THEN
           h2osno(c) = h2osno(c) - qflx_sub_snow(c)*dtime
           qflx_snowcap(p) = qflx_snowcap(p) + qflx_dew_snow(c) + qflx_dew_grnd(c)
        ELSE
           h2osno(c) = h2osno(c) + (-qflx_sub_snow(c)+qflx_dew_snow(c))*dtime
        END IF
        IF (h2osno_temp > 0._r8) THEN
           snowdp(c) = snowdp(c) * h2osno(c) / h2osno_temp
        ELSE
           snowdp(c) = h2osno(c)/snow_bd !Assume a constant snow bulk density = 250.
        END IF
        !
        !#define PERGRO
        !#if (defined PERGRO)
        IF (ABS(h2osno(c)) < 1.e-10_r8) h2osno(c) = 0._r8
        !#else
        !          h2osno(c) = max(h2osno(c), 0._r8)
        !#endif

     END IF

     qflx_snowcap_col(c) = qflx_snowcap(p)

  END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Determine initial snow/no-snow filters (will be modified possibly by
  ! routines CombineSnowLayers and DivideSnowLayers below

  CALL BuildSnowFilter(lbc, ubc, num_shlakec, filter_shlakec,snl,       &            !i
       num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc) !o

  ! Determine the change of snow mass and the snow water onto soil

  CALL SnowWater(lbc, ubc, num_shlakesnowc, filter_shlakesnowc,         & !i
       num_shlakenosnowc, filter_shlakenosnowc,               & !i
       snl,do_capsnow,qflx_snomelt,qflx_rain_grnd,            & !i
       qflx_sub_snow,qflx_evap_grnd,                          & !i
       qflx_dew_snow,qflx_dew_grnd,dz,                        & !i
       h2osoi_ice,h2osoi_liq,                                 & !i&o
       qflx_top_soil)                                           !o


  ! Determine soil hydrology
  ! Here this consists only of making sure that soil is saturated even as it melts and 10%
  ! of pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
  ! saturation value, then remove water.

  DO j = 1,nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        IF (h2osoi_vol(c,j) < watsat(c,j)) THEN
           h2osoi_liq(c,j) = (watsat(c,j)*dz(c,j) - h2osoi_ice(c,j)/denice)*denh2o
           ! h2osoi_vol will be updated below, and this water addition will come from qflx_qrgwl
        ELSE IF (h2osoi_liq(c,j) > watsat(c,j)*denh2o*dz(c,j)) THEN
           h2osoi_liq(c,j) = watsat(c,j)*denh2o*dz(c,j)
        END IF

     END DO
  END DO
!!!!!!!!!!

  !    if (.not. is_perpetual()) then
  IF (1==1) THEN

     ! Natural compaction and metamorphosis.

     CALL SnowCompaction(lbc, ubc, num_shlakesnowc, filter_shlakesnowc,   &!i
          snl,imelt,frac_iceold,t_soisno,                  &!i
          h2osoi_ice,h2osoi_liq,                           &!i
          dz)                                               !&o

     ! Combine thin snow elements

     CALL CombineSnowLayers(lbc, ubc,                            & !i
          num_shlakesnowc, filter_shlakesnowc, & !i&o
          snl,h2osno,snowdp,dz,zi,             & !i&o
          t_soisno,h2osoi_ice,h2osoi_liq,      & !i&o
          z)  !o


     ! Divide thick snow elements

     CALL DivideSnowLayers(lbc, ubc,                             & !i
          num_shlakesnowc, filter_shlakesnowc,  & !i&o
          snl,dz,zi,t_soisno,                   & !i&o
          h2osoi_ice,h2osoi_liq,                & !i&o
          z)  !o


  ELSE

     DO fc = 1, num_shlakesnowc
        c = filter_shlakesnowc(fc)
        h2osno(c) = 0._r8
     END DO
     DO j = -nlevsnow+1,0
        DO fc = 1, num_shlakesnowc
           c = filter_shlakesnowc(fc)
           IF (j >= snl(c)+1) THEN
              h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
           END IF
        END DO
     END DO

  END IF

  ! Check for snow layers above lake with unfrozen top layer.  Mechanically,
  ! the snow will fall into the lake and melt or turn to ice.  If the top layer has
  ! sufficient heat to melt the snow without freezing, then that will be done.
  ! Otherwise, the top layer will undergo freezing, but only if the top layer will
  ! not freeze completely.  Otherwise, let the snow layers persist and melt by diffusion.
  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)

     IF (t_lake(c,1) > tfrz .AND. lake_icefrac(c,1) == 0._r8 .AND. snl(c) < 0) THEN
        unfrozen(c) = .TRUE.
     ELSE
        unfrozen(c) = .FALSE.
     END IF
  END DO

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        IF (unfrozen(c)) THEN
           IF (j == -nlevsnow+1) THEN
              sumsnowice(c) = 0._r8
              heatsum(c) = 0._r8
           END IF
           IF (j >= snl(c)+1) THEN
              sumsnowice(c) = sumsnowice(c) + h2osoi_ice(c,j)
              heatsum(c) = heatsum(c) + h2osoi_ice(c,j)*cpice*(tfrz - t_soisno(c,j)) &
                   + h2osoi_liq(c,j)*cpliq*(tfrz - t_soisno(c,j))
           END IF
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec
     c = filter_shlakec(fc)

     IF (unfrozen(c)) THEN
        heatsum(c) = heatsum(c) + sumsnowice(c)*hfus
        heatrem = (t_lake(c,1) - tfrz)*cpliq*denh2o*dz_lake(c,1) - heatsum(c)

        IF (heatrem + denh2o*dz_lake(c,1)*hfus > 0._r8) THEN
           ! Remove snow and subtract the latent heat from the top layer.
           h2osno(c) = 0._r8
           snl(c) = 0
           ! The rest of the bookkeeping for the removed snow will be done below.
           !#define LAKEDEBUG
           !#if (defined LAKEDEBUG)
           !                write(message,*)'Snow layers removed above unfrozen lake for column, snowice:', &
           !                          c, sumsnowice(c)
           !                CALL wrf_message(message)
           !#endif
           IF (heatrem > 0._r8) THEN ! simply subtract the heat from the layer
              t_lake(c,1) = t_lake(c,1) - heatrem/(cpliq*denh2o*dz_lake(c,1))
           ELSE !freeze part of the layer
              t_lake(c,1) = tfrz
              lake_icefrac(c,1) = -heatrem/(denh2o*dz_lake(c,1)*hfus)
           END IF
        END IF
     END IF
  END DO
!!!!!!!!!!!!

  ! Set snow age to zero if no snow

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakesnowc
     c = filter_shlakesnowc(fc)
     IF (snl(c) == 0) THEN
        snowage(c) = 0._r8
     END IF
  END DO

  ! Set empty snow layers to zero

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakesnowc
        c = filter_shlakesnowc(fc)
        IF (j <= snl(c) .AND. snl(c) > nlevsnow) THEN
           h2osoi_ice(c,j) = 0._r8
           h2osoi_liq(c,j) = 0._r8
           t_soisno(c,j) = 0._r8
           dz(c,j) = 0._r8
           z(c,j) = 0._r8
           zi(c,j-1) = 0._r8
        END IF
     END DO
  END DO

  ! Build new snow filter

  CALL BuildSnowFilter(lbc, ubc, num_shlakec, filter_shlakec, snl,&   !i
       num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc) !o

  ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
  ! over all snow layers for history output

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakesnowc
     c = filter_shlakesnowc(fc)
     t_snow(c)  = 0._r8
     snowice(c) = 0._r8
     snowliq(c) = 0._r8
  END DO
  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakenosnowc
     c = filter_shlakenosnowc(fc)
     t_snow(c)  = spval
     snowice(c) = spval
     snowliq(c) = spval
  END DO

  DO j = -nlevsnow+1, 0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakesnowc
        c = filter_shlakesnowc(fc)
        IF (j >= snl(c)+1) THEN
           t_snow(c)  = t_snow(c) + t_soisno(c,j)
           snowice(c) = snowice(c) + h2osoi_ice(c,j)
           snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
        END IF
     END DO
  END DO

  ! Determine ending water balance and volumetric soil water

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_shlakec

     c = filter_shlakec(fc)
     IF (snl(c) < 0) t_snow(c) = t_snow(c)/ABS(snl(c))
     endwb(c) = h2osno(c)
  END DO

  DO j = 1, nlevsoil
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)
        endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
        h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
     END DO
  END DO

  !#define LAKEDEBUG
  !#if (defined LAKEDEBUG)
  ! Check to make sure snow water adds up correctly.
  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_shlakec
        c = filter_shlakec(fc)

        jtop = snl(c)+1
        IF(j == jtop) snow_water(c) = 0._r8
        IF(j >= jtop) THEN
           snow_water(c) = snow_water(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
           !            if(j == 0 .and. abs(snow_water(c)-h2osno(c))>1.e-7_r8) then
           !               write(message,*)'h2osno does not equal sum of snow layers in ShalLakeHydrology:', &
           !                         'column, h2osno, sum of snow layers =', c, h2osno(c), snow_water(c)
           !               CALL c( message )
           !            end if
        END IF
     END DO
  END DO
  !#endif

!!!!!!!!!!!!!
  ! Do history variables and set special landunit runoff (adapted from end of HydrologyLake)
  !dir$ concurrent
  !cdir nodep
  DO fp = 1,num_shlakep
     p = filter_shlakep(fp)
     c = pcolumn(p)
     g = pgridcell(p)

     qflx_infl(c)      = 0._r8
     qflx_surf(c)      = 0._r8
     qflx_drain(c)     = 0._r8
     rootr_column(c,:) = spval
     soilalpha(c)      = spval
     zwt(c)            = spval
     fcov(c)           = spval
     qcharge(c)        = spval
     !       h2osoi_vol(c,:)   = spval

     ! Insure water balance using qflx_qrgwl
     qflx_qrgwl(c)     = forc_rain(g) + forc_snow(g) - qflx_evap_tot(p) - (endwb(c)-begwb(c))/dtime
     !#define LAKEDEBUG
     !#if (defined LAKEDEBUG)
     !    write(message,*)'c, rain, snow, evap, endwb, begwb, qflx_qrgwl:', &
     !       c, forc_rain(g), forc_snow(g), qflx_evap_tot(p), endwb(c), begwb(c), qflx_qrgwl(c)
     !    CALL wrf_message(message)
     !#endif

     ! The pft average must be done here for output to history tape
     qflx_evap_tot_col(c) = qflx_evap_tot(p)
  END DO

!!!!!!!!!!!!!
  !For now, bracket off the remaining biogeochem code.  May need to bring it back
  !to do soil carbon and methane beneath lakes.
  !#defined CN
  !#if (defined CN)
  !#ifndef SHLAKE
  !    do j = 1, nlevsoil
  !!dir$ concurrent
  !!cdir nodep
  !       do fc = 1, nlevsoil
  !          c = filter_soilc(fc)
  !
  !          if (h2osoi_liq(c,j) > 0._r8) then
  !             vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
  !
  !             ! the following limit set to catch very small values of
  !             ! fractional saturation that can crash the calculation of psi
  !
  !             fsat = max(vwc/vwcsat(c,j), 0.001_r8)
  !             psi = psisat(c,j) * (fsat)**bsw2(c,j)
  !             soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)
  !          else
  !             soilpsi(c,j) = -15.0_r8
  !          end if
  !       end do
  !    end do
  !#endif
  !#endif

  !#defined DGVM
  !#defined CN
  !#if (defined DGVM) || (defined CN)
  !#ifndef SHLAKE
  ! Available soil water up to a depth of 0.5 m.
  ! Potentially available soil water (=whc) up to a depth of 0.5 m.
  ! Water content as fraction of whc up to a depth of 0.5 m.

  !!dir$ concurrent
  !!cdir nodep
  !    do c = lbc,ubc
  !       l = clandunit(c)
  !!       if (ityplun(l) == istsoil) then
  !          rwat(c) = 0._r8
  !          swat(c) = 0._r8
  !          rz(c)   = 0._r8
  !!       end if
  !    end do
  !
  !    do j = 1, nlevsoil
  !!dir$ concurrent
  !!cdir nodep
  !       do c = lbc,ubc
  !          l = clandunit(c)
  !!          if (ityplun(l) == istsoil) then
  !             if (z(c,j)+0.5_r8*dz(c,j) <= 0.5_r8) then
  !                watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
  !                rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
  !                swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
  !                rz(c) = rz(c) + dz(c,j)
  !!             end if
  !          end if
  !       end do
  !    end do
  !
  !!dir$ concurrent
  !!cdir nodep
  !    do c = lbc,ubc
  !       l = clandunit(c)
  !!       if (ityplun(l) == istsoil) then
  !          if (rz(c) /= 0._r8) then
  !             tsw  = rwat(c)/rz(c)
  !             stsw = swat(c)/rz(c)
  !          else
  !             watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
  !             tsw = h2osoi_vol(c,1) - watdry
  !             stsw = watsat(c,1) - watdry
  !          end if
  !          wf(c) = tsw/stsw
  !!       else
  !!          wf(c) = 1.0_r8
  !!       end if
  !    end do

  !#endif
  !#endif

END SUBROUTINE ShalLakeHydrology

SUBROUTINE QSat (T, p, es, esdT, qs, qsdT)
  !
  ! !DESCRIPTION:
  ! Computes saturation mixing ratio and the change in saturation
  ! mixing ratio with respect to temperature.
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in)  :: T
  REAL(r8), INTENT(in)  :: p
  REAL(r8), INTENT(out) :: es
  REAL(r8), INTENT(out) :: esdT
  REAL(r8), INTENT(out) :: qs
  REAL(r8), INTENT(out) :: qsdT

  REAL(r8) :: T_limit
  REAL(r8) :: td,vp,vp1,vp2

  REAL(r8), PARAMETER :: a0 =  6.11213476
  REAL(r8), PARAMETER :: a1 =  0.444007856
  REAL(r8), PARAMETER :: a2 =  0.143064234e-01
  REAL(r8), PARAMETER :: a3 =  0.264461437e-03
  REAL(r8), PARAMETER :: a4 =  0.305903558e-05
  REAL(r8), PARAMETER :: a5 =  0.196237241e-07
  REAL(r8), PARAMETER :: a6 =  0.892344772e-10
  REAL(r8), PARAMETER :: a7 = -0.373208410e-12
  REAL(r8), PARAMETER :: a8 =  0.209339997e-15

  REAL(r8), PARAMETER :: b0 =  0.444017302
  REAL(r8), PARAMETER :: b1 =  0.286064092e-01
  REAL(r8), PARAMETER :: b2 =  0.794683137e-03
  REAL(r8), PARAMETER :: b3 =  0.121211669e-04
  REAL(r8), PARAMETER :: b4 =  0.103354611e-06
  REAL(r8), PARAMETER :: b5 =  0.404125005e-09
  REAL(r8), PARAMETER :: b6 = -0.788037859e-12
  REAL(r8), PARAMETER :: b7 = -0.114596802e-13
  REAL(r8), PARAMETER :: b8 =  0.381294516e-16

  REAL(r8), PARAMETER :: c0 =  6.11123516
  REAL(r8), PARAMETER :: c1 =  0.503109514
  REAL(r8), PARAMETER :: c2 =  0.188369801e-01
  REAL(r8), PARAMETER :: c3 =  0.420547422e-03
  REAL(r8), PARAMETER :: c4 =  0.614396778e-05
  REAL(r8), PARAMETER :: c5 =  0.602780717e-07
  REAL(r8), PARAMETER :: c6 =  0.387940929e-09
  REAL(r8), PARAMETER :: c7 =  0.149436277e-11
  REAL(r8), PARAMETER :: c8 =  0.262655803e-14

  REAL(r8), PARAMETER :: d0 =  0.503277922
  REAL(r8), PARAMETER :: d1 =  0.377289173e-01
  REAL(r8), PARAMETER :: d2 =  0.126801703e-02
  REAL(r8), PARAMETER :: d3 =  0.249468427e-04
  REAL(r8), PARAMETER :: d4 =  0.313703411e-06
  REAL(r8), PARAMETER :: d5 =  0.257180651e-08
  REAL(r8), PARAMETER :: d6 =  0.133268878e-10
  REAL(r8), PARAMETER :: d7 =  0.394116744e-13
  REAL(r8), PARAMETER :: d8 =  0.498070196e-16
  !-----------------------------------------------------------------------

  T_limit = T - tfrz
  IF (T_limit > 100.0) T_limit=100.0
  IF (T_limit < -75.0) T_limit=-75.0

  td       = T_limit
  IF (td >= 0.0) THEN
     es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
          + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
     esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
          + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
  ELSE
     es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
          + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
     esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
          + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
  ENDIF

  es    = es    * 100.            ! pa
  esdT  = esdT  * 100.            ! pa/K

  vp    = 1.0   / (p - 0.378*es)
  vp1   = 0.622 * vp
  vp2   = vp1   * vp

  qs    = es    * vp1             ! kg/kg
  qsdT  = esdT  * vp2 * p         ! 1 / K

END SUBROUTINE QSat

SUBROUTINE Tridiagonal (lbc, ubc, lbj, ubj, jtop, numf, filter, &
     a, b, c, r, u)
  !
  ! !DESCRIPTION:
  ! Tridiagonal matrix solution
  !
  ! !USES:
  !  use shr_kind_mod, only: r8 => shr_kind_r8
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  INTEGER , INTENT(in)    :: lbc, ubc
  INTEGER , INTENT(in)    :: lbj, ubj
  INTEGER , INTENT(in)    :: jtop(lbc:ubc)
  INTEGER , INTENT(in)    :: numf
  INTEGER , INTENT(in)    :: filter(1:numf)
  REAL(r8), INTENT(in)    :: a(lbc:ubc, lbj:ubj)
  REAL(r8), INTENT(in)    :: b(lbc:ubc, lbj:ubj)
  REAL(r8), INTENT(in)    :: c(lbc:ubc, lbj:ubj)
  REAL(r8), INTENT(in)    :: r(lbc:ubc, lbj:ubj)
  REAL(r8), INTENT(inout) :: u(lbc:ubc, lbj:ubj)


  INTEGER  :: j,ci,fc
  REAL(r8) :: gam(lbc:ubc,lbj:ubj)
  REAL(r8) :: bet(lbc:ubc)
  !-----------------------------------------------------------------------

  ! Solve the matrix
  ! print *,'a=',a
  ! print *,'b=',b
  ! print *,'c=',c

  !dir$ concurrent
  !cdir nodep
  DO fc = 1,numf
     ci = filter(fc)
     bet(ci) = b(ci,jtop(ci))
  END DO

  DO j = lbj, ubj
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,numf
        ci = filter(fc)
        IF (j >= jtop(ci)) THEN
           IF (j == jtop(ci)) THEN
              u(ci,j) = r(ci,j) / bet(ci)

           ELSE
              gam(ci,j) = c(ci,j-1) / bet(ci)
              bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
              u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)

           END IF
        END IF
     END DO
  END DO


  !Cray X1 unroll directive used here as work-around for compiler issue 2003/10/20
  !dir$ unroll 0
  DO j = ubj-1,lbj,-1
     !dir$ prefervector
     !dir$ concurrent
     !cdir nodep
     DO fc = 1,numf
        ci = filter(fc)
        IF (j >= jtop(ci)) THEN
           u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
        END IF
     END DO
  END DO
                !  print *,'u=',u
END SUBROUTINE Tridiagonal

SUBROUTINE SnowWater(lbc, ubc, num_snowc, filter_snowc,         & !i
     num_nosnowc, filter_nosnowc,               & !i
     snl,do_capsnow,qflx_snomelt,qflx_rain_grnd,            & !i
     qflx_sub_snow,qflx_evap_grnd,                          & !i
     qflx_dew_snow,qflx_dew_grnd,dz,                        & !i
     h2osoi_ice,h2osoi_liq,                                 & !i&o
     qflx_top_soil)                                           !o
  !===============================================================================
  ! !DESCRIPTION:
  ! Evaluate the change of snow mass and the snow water onto soil.
  ! Water flow within snow is computed by an explicit and non-physical
  ! based scheme, which permits a part of liquid water over the holding
  ! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
  ! percolate into the underlying layer.  Except for cases where the
  ! porosity of one of the two neighboring layers is less than 0.05, zero
  ! flow is assumed. The water flow out of the bottom of the snow pack will
  ! participate as the input of the soil water and runoff.  This subroutine
  ! uses a filter for columns containing snow which must be constructed prior
  ! to being called.
  !
  ! !REVISION HISTORY:
  ! 15 September 1999: Yongjiu Dai; Initial code
  ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
  ! 15 November 2000: Mariana Vertenstein
  ! 2/26/02, Peter Thornton: Migrated to new data structures.
  !=============================================================================
  ! !USES:
  !  use clmtype

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER, INTENT(in) :: lbc, ubc
  INTEGER, INTENT(in) :: num_snowc
  INTEGER, INTENT(in) :: filter_snowc(ubc-lbc+1)
  INTEGER, INTENT(in) :: num_nosnowc
  INTEGER, INTENT(in) :: filter_nosnowc(ubc-lbc+1)

  INTEGER , INTENT(in) :: snl(1)
  LOGICAL , INTENT(in) :: do_capsnow(1)
  REAL(r8), INTENT(in) :: qflx_snomelt(1)
  REAL(r8), INTENT(in) :: qflx_rain_grnd(1)
  REAL(r8), INTENT(in) :: qflx_sub_snow(1)
  REAL(r8), INTENT(in) :: qflx_evap_grnd(1)
  REAL(r8), INTENT(in) :: qflx_dew_snow(1)
  REAL(r8), INTENT(in) :: qflx_dew_grnd(1)
  REAL(r8), INTENT(in) :: dz(1,-nlevsnow+1:nlevsoil)

  REAL(r8), INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)


  REAL(r8), INTENT(out) :: qflx_top_soil(1)

  INTEGER  :: c, j, fc
  REAL(r8) :: qin(lbc:ubc)
  REAL(r8) :: qout(lbc:ubc)
  REAL(r8) :: wgdif
  REAL(r8) :: vol_liq(lbc:ubc,-nlevsnow+1:0)
  REAL(r8) :: vol_ice(lbc:ubc,-nlevsnow+1:0)
  REAL(r8) :: eff_porosity(lbc:ubc,-nlevsnow+1:0)
  !-----------------------------------------------------------------------
  ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
  ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

  !dir$ concurrent
  !cdir nodep
  DO fc = 1,num_snowc
     c = filter_snowc(fc)
     IF (do_capsnow(c)) THEN
        wgdif = h2osoi_ice(c,snl(c)+1) - qflx_sub_snow(c)*dtime
        h2osoi_ice(c,snl(c)+1) = wgdif
        IF (wgdif < 0.) THEN
           h2osoi_ice(c,snl(c)+1) = 0.
           h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
        END IF
        h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) - qflx_evap_grnd(c) * dtime
     ELSE
        wgdif = h2osoi_ice(c,snl(c)+1) + (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime
        h2osoi_ice(c,snl(c)+1) = wgdif
        IF (wgdif < 0.) THEN
           h2osoi_ice(c,snl(c)+1) = 0.
           h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
        END IF
        h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
             (qflx_rain_grnd(c) + qflx_dew_grnd(c) - qflx_evap_grnd(c)) * dtime
     END IF
     h2osoi_liq(c,snl(c)+1) = MAX(0._r8, h2osoi_liq(c,snl(c)+1))
  END DO

  ! Porosity and partial volume

  DO j = -nlevsnow+1, 0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN
           vol_ice(c,j) = MIN(1._r8, h2osoi_ice(c,j)/(dz(c,j)*denice))
           eff_porosity(c,j) = 1. - vol_ice(c,j)
           vol_liq(c,j) = MIN(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*denh2o))
        END IF
     END DO
  END DO

  ! Capillary forces within snow are usually two or more orders of magnitude
  ! less than those of gravity. Only gravity terms are considered.
  ! the genernal expression for water flow is "K * ss**3", however,
  ! no effective parameterization for "K".  Thus, a very simple consideration
  ! (not physically based) is introduced:
  ! when the liquid water of layer exceeds the layer's holding
  ! capacity, the excess meltwater adds to the underlying neighbor layer.

  qin(:) = 0._r8

  DO j = -nlevsnow+1, 0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN
           h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(c)
           IF (j <= -1) THEN
              ! No runoff over snow surface, just ponding on surface
              IF (eff_porosity(c,j) < wimp .OR. eff_porosity(c,j+1) < wimp) THEN
                 qout(c) = 0._r8
              ELSE
                 qout(c) = MAX(0._r8,(vol_liq(c,j)-ssi*eff_porosity(c,j))*dz(c,j))
                 qout(c) = MIN(qout(c),(1.-vol_ice(c,j+1)-vol_liq(c,j+1))*dz(c,j+1))
              END IF
           ELSE
              qout(c) = MAX(0._r8,(vol_liq(c,j) - ssi*eff_porosity(c,j))*dz(c,j))
           END IF
           qout(c) = qout(c)*1000.
           h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(c)
           qin(c) = qout(c)
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_snowc
     c = filter_snowc(fc)
     ! Qout from snow bottom
     qflx_top_soil(c) = qout(c) / dtime
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_nosnowc
     c = filter_nosnowc(fc)
     qflx_top_soil(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
  END DO

END SUBROUTINE SnowWater

SUBROUTINE SnowCompaction(lbc, ubc, num_snowc, filter_snowc,   &!i
     snl,imelt,frac_iceold,t_soisno,                  &!i
     h2osoi_ice,h2osoi_liq,                           &!i
     dz)                                               !i&o


  !================================================================================
  ! !DESCRIPTION:
  ! Determine the change in snow layer thickness due to compaction and
  ! settling.
  ! Three metamorphisms of changing snow characteristics are implemented,
  ! i.e., destructive, overburden, and melt. The treatments of the former
  ! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
  ! due to melt metamorphism is simply taken as a ratio of snow ice
  ! fraction after the melting versus before the melting.
  !
  ! CALLED FROM:
  ! subroutine Hydrology2 in module Hydrology2Mod
  !
  ! REVISION HISTORY:
  ! 15 September 1999: Yongjiu Dai; Initial code
  ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
  ! 2/28/02, Peter Thornton: Migrated to new data structures
  !==============================================================================
  ! USES:
  !  use clmtype
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER, INTENT(in) :: lbc, ubc
  INTEGER, INTENT(in) :: num_snowc
  INTEGER, INTENT(in) :: filter_snowc(ubc-lbc+1)
  INTEGER,  INTENT(in) :: snl(1)
  INTEGER,  INTENT(in) :: imelt(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: frac_iceold(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)


  REAL(r8), INTENT(inout) :: dz(1,-nlevsnow+1:nlevsoil)

  INTEGER :: j, c, fc
  REAL(r8), PARAMETER :: c2 = 23.e-3
  REAL(r8), PARAMETER :: c3 = 2.777e-6
  REAL(r8), PARAMETER :: c4 = 0.04
  REAL(r8), PARAMETER :: c5 = 2.0
  REAL(r8), PARAMETER :: dm = 100.0
  REAL(r8), PARAMETER :: eta0 = 9.e+5
  REAL(r8) :: burden(lbc:ubc)
  REAL(r8) :: ddz1
  REAL(r8) :: ddz2
  REAL(r8) :: ddz3
  REAL(r8) :: dexpf
  REAL(r8) :: fi
  REAL(r8) :: td
  REAL(r8) :: pdzdtc
  REAL(r8) :: void
  REAL(r8) :: wx
  REAL(r8) :: bi

  !-----------------------------------------------------------------------


  ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

  burden(:) = 0._r8

  DO j = -nlevsnow+1, 0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN

           wx = h2osoi_ice(c,j) + h2osoi_liq(c,j)
           void = 1. - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o) / dz(c,j)

           ! Allow compaction only for non-saturated node and higher ice lens node.
           IF (void > 0.001 .AND. h2osoi_ice(c,j) > .1) THEN
              bi = h2osoi_ice(c,j) / dz(c,j)
              fi = h2osoi_ice(c,j) / wx
              td = tfrz-t_soisno(c,j)
              dexpf = EXP(-c4*td)

              ! Settling as a result of destructive metamorphism

              ddz1 = -c3*dexpf
              IF (bi > dm) ddz1 = ddz1*EXP(-46.0e-3*(bi-dm))

              ! Liquid water term

              IF (h2osoi_liq(c,j) > 0.01*dz(c,j)) ddz1=ddz1*c5

              ! Compaction due to overburden

              ddz2 = -burden(c)*EXP(-0.08*td - c2*bi)/eta0

              ! Compaction occurring during melt

              IF (imelt(c,j) == 1) THEN
                 ddz3 = - 1./dtime * MAX(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
              ELSE
                 ddz3 = 0._r8
              END IF

              ! Time rate of fractional change in dz (units of s-1)

              pdzdtc = ddz1 + ddz2 + ddz3

              ! The change in dz due to compaction

              dz(c,j) = dz(c,j) * (1.+pdzdtc*dtime)
           END IF

           ! Pressure of overlying snow

           burden(c) = burden(c) + wx

        END IF
     END DO
  END DO

END SUBROUTINE SnowCompaction

SUBROUTINE CombineSnowLayers(lbc, ubc,                            & !i
     num_snowc, filter_snowc, & !i&o
     snl,h2osno,snowdp,dz,zi,             & !i&o
     t_soisno,h2osoi_ice,h2osoi_liq,      & !i&o
     z)  !o
  !==========================================================================
  ! !DESCRIPTION:
  ! Combine snow layers that are less than a minimum thickness or mass
  ! If the snow element thickness or mass is less than a prescribed minimum,
  ! then it is combined with a neighboring element.  The subroutine
  ! clm\_combo.f90 then executes the combination of mass and energy.
  ! !CALLED FROM:
  ! subroutine Hydrology2 in module Hydrology2Mod
  !
  ! !REVISION HISTORY:
  ! 15 September 1999: Yongjiu Dai; Initial code
  ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
  ! 2/28/02, Peter Thornton: Migrated to new data structures.
  !=========================================================================
  ! !USES:
  !  use clmtype
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  ! INTEGER :: i
  !
  ! DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
  !      10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
  !
  ! DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
  !      33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER, INTENT(in)    :: lbc, ubc
  INTEGER, INTENT(inout) :: num_snowc
  INTEGER, INTENT(inout) :: filter_snowc(ubc-lbc+1)
  INTEGER , INTENT(inout) :: snl(1)
  REAL(r8), INTENT(inout) :: h2osno(1)
  REAL(r8), INTENT(inout) :: snowdp(1)
  REAL(r8), INTENT(inout) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: zi(1,-nlevsnow+0:nlevsoil)
  REAL(r8), INTENT(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)

  REAL(r8), INTENT(out) :: z(1,-nlevsnow+1:nlevsoil)
  INTEGER :: c, fc
  INTEGER :: i,k
  INTEGER :: j,l
  INTEGER :: msn_old(lbc:ubc)
  INTEGER :: mssi(lbc:ubc)
  INTEGER :: neibor
  REAL(r8):: zwice(lbc:ubc)
  REAL(r8):: zwliq (lbc:ubc)
  REAL(r8):: dzmin(5)

  DATA dzmin /0.010, 0.015, 0.025, 0.055, 0.115/
  !-----------------------------------------------------------------------

  ! Check the mass of ice lens of snow, when the total is less than a small value,
  ! combine it with the underlying neighbor.

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_snowc
     c = filter_snowc(fc)
     msn_old(c) = snl(c)
  END DO

  ! The following loop is NOT VECTORIZED

  DO fc = 1, num_snowc
     c = filter_snowc(fc)
     !    l = clandunit(c)
     DO j = msn_old(c)+1,0
        IF (h2osoi_ice(c,j) <= .1) THEN
           !  if (ityplun(l) == istsoil) then
           !     h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
           !     h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
           !  else if (ityplun(l) /= istsoil .and. j /= 0) then
           h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
           h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
           !  end if

           ! shift all elements above this down one.
           IF (j > snl(c)+1 .AND. snl(c) < -1) THEN
              DO i = j, snl(c)+2, -1
                 t_soisno(c,i)   = t_soisno(c,i-1)
                 h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                 h2osoi_ice(c,i) = h2osoi_ice(c,i-1)
                 dz(c,i)         = dz(c,i-1)
              END DO
           END IF
           snl(c) = snl(c) + 1
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_snowc
     c = filter_snowc(fc)
     h2osno(c) = 0._r8
     snowdp(c) = 0._r8
     zwice(c)  = 0._r8
     zwliq(c)  = 0._r8
  END DO

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN
           h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
           snowdp(c) = snowdp(c) + dz(c,j)
           zwice(c)  = zwice(c) + h2osoi_ice(c,j)
           zwliq(c)  = zwliq(c) + h2osoi_liq(c,j)
        END IF
     END DO
  END DO

  ! Check the snow depth - all snow gone
  ! The liquid water assumes ponding on soil surface.

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_snowc
     c = filter_snowc(fc)
     ! l = clandunit(c)
     IF (snowdp(c) < 0.01 .AND. snowdp(c) > 0.) THEN
        snl(c) = 0
        h2osno(c) = zwice(c)
        IF (h2osno(c) <= 0.) snowdp(c) = 0._r8
        !    if (ityplun(l) == istsoil) h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(c)    !change by guhp
     END IF
  END DO

  ! Check the snow depth - snow layers combined
  ! The following loop IS NOT VECTORIZED

  DO fc = 1, num_snowc
     c = filter_snowc(fc)

     ! Two or more layers

     IF (snl(c) < -1) THEN

        msn_old(c) = snl(c)
        mssi(c) = 1

        DO i = msn_old(c)+1,0
           IF (dz(c,i) < dzmin(mssi(c))) THEN

              IF (i == snl(c)+1) THEN
                 ! If top node is removed, combine with bottom neighbor.
                 neibor = i + 1
              ELSE IF (i == 0) THEN
                 ! If the bottom neighbor is not snow, combine with the top neighbor.
                 neibor = i - 1
              ELSE
                 ! If none of the above special cases apply, combine with the thinnest neighbor
                 neibor = i + 1
                 IF ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
              END IF

              ! Node l and j are combined and stored as node j.
              IF (neibor > i) THEN
                 j = neibor
                 l = i
              ELSE
                 j = i
                 l = neibor
              END IF

              CALL Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                   t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

              ! Now shift all elements above this down one.
              IF (j-1 > snl(c)+1) THEN
                 DO k = j-1, snl(c)+2, -1
                    t_soisno(c,k) = t_soisno(c,k-1)
                    h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                    h2osoi_liq(c,k) = h2osoi_liq(c,k-1)
                    dz(c,k) = dz(c,k-1)
                 END DO
              END IF

              ! Decrease the number of snow layers
              snl(c) = snl(c) + 1
              IF (snl(c) >= -1) EXIT

           ELSE

              ! The layer thickness is greater than the prescribed minimum value
              mssi(c) = mssi(c) + 1

           END IF
        END DO

     END IF

  END DO

  ! Reset the node depth and the depth of layer interface

  DO j = 0, -nlevsnow+1, -1
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c) + 1) THEN
           z(c,j) = zi(c,j) - 0.5*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
        END IF
     END DO
  END DO

END SUBROUTINE CombineSnowLayers

SUBROUTINE DivideSnowLayers(lbc, ubc,                             & !i
     num_snowc, filter_snowc,  & !i&o
     snl,dz,zi,t_soisno,                   & !i&o
     h2osoi_ice,h2osoi_liq,                & !i&o
     z)  !o


  !============================================================================
  ! !DESCRIPTION:
  ! Subdivides snow layers if they exceed their prescribed maximum thickness.
  ! !CALLED FROM:
  ! subroutine Hydrology2 in module Hydrology2Mod
  !
  ! !REVISION HISTORY:
  ! 15 September 1999: Yongjiu Dai; Initial code
  ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
  ! 2/28/02, Peter Thornton: Migrated to new data structures.
  !============================================================================
  ! !USES:
  !   use clmtype
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER, INTENT(in)    :: lbc, ubc
  INTEGER, INTENT(inout) :: num_snowc
  INTEGER, INTENT(inout) :: filter_snowc(ubc-lbc+1)
  INTEGER , INTENT(inout) :: snl(1)
  REAL(r8), INTENT(inout) :: dz(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: zi(1,-nlevsnow+0:nlevsoil)
  REAL(r8), INTENT(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)
  REAL(r8), INTENT(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)


  REAL(r8), INTENT(out) :: z(1,-nlevsnow+1:nlevsoil)

  INTEGER  :: j, c, fc
  REAL(r8) :: drr
  INTEGER  :: msno
  REAL(r8) :: dzsno(lbc:ubc,nlevsnow)
  REAL(r8) :: swice(lbc:ubc,nlevsnow)
  REAL(r8) :: swliq(lbc:ubc,nlevsnow)
  REAL(r8) :: tsno(lbc:ubc ,nlevsnow)
  REAL(r8) :: zwice
  REAL(r8) :: zwliq
  REAL(r8) :: propor
  !-----------------------------------------------------------------------

  ! Begin calculation - note that the following column loops are only invoked
  ! for snow-covered columns

  DO j = 1,5
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j <= ABS(snl(c))) THEN
           dzsno(c,j) = dz(c,j+snl(c))
           swice(c,j) = h2osoi_ice(c,j+snl(c))
           swliq(c,j) = h2osoi_liq(c,j+snl(c))
           tsno(c,j)  = t_soisno(c,j+snl(c))
        END IF
     END DO
  END DO

  !dir$ concurrent
  !cdir nodep
  DO fc = 1, num_snowc
     c = filter_snowc(fc)

     msno = ABS(snl(c))

     IF (msno == 1) THEN
        ! Specify a new snow layer
        IF (dzsno(c,1) > 0.03) THEN
           msno = 2
           dzsno(c,1) = dzsno(c,1)/2.
           swice(c,1) = swice(c,1)/2.
           swliq(c,1) = swliq(c,1)/2.
           dzsno(c,2) = dzsno(c,1)
           swice(c,2) = swice(c,1)
           swliq(c,2) = swliq(c,1)
           tsno(c,2)  = tsno(c,1)
        END IF
     END IF

     IF (msno > 1) THEN
        IF (dzsno(c,1) > 0.02) THEN
           drr = dzsno(c,1) - 0.02
           propor = drr/dzsno(c,1)
           zwice = propor*swice(c,1)
           zwliq = propor*swliq(c,1)
           propor = 0.02/dzsno(c,1)
           swice(c,1) = propor*swice(c,1)
           swliq(c,1) = propor*swliq(c,1)
           dzsno(c,1) = 0.02

           CALL Combo (dzsno(c,2), swliq(c,2), swice(c,2), tsno(c,2), drr, &
                zwliq, zwice, tsno(c,1))

           ! Subdivide a new layer
           IF (msno <= 2 .AND. dzsno(c,2) > 0.07) THEN
              msno = 3
              dzsno(c,2) = dzsno(c,2)/2.
              swice(c,2) = swice(c,2)/2.
              swliq(c,2) = swliq(c,2)/2.
              dzsno(c,3) = dzsno(c,2)
              swice(c,3) = swice(c,2)
              swliq(c,3) = swliq(c,2)
              tsno(c,3)  = tsno(c,2)
           END IF
        END IF
     END IF

     IF (msno > 2) THEN
        IF (dzsno(c,2) > 0.05) THEN
           drr = dzsno(c,2) - 0.05
           propor = drr/dzsno(c,2)
           zwice = propor*swice(c,2)
           zwliq = propor*swliq(c,2)
           propor = 0.05/dzsno(c,2)
           swice(c,2) = propor*swice(c,2)
           swliq(c,2) = propor*swliq(c,2)
           dzsno(c,2) = 0.05

           CALL Combo (dzsno(c,3), swliq(c,3), swice(c,3), tsno(c,3), drr, &
                zwliq, zwice, tsno(c,2))

           ! Subdivided a new layer
           IF (msno <= 3 .AND. dzsno(c,3) > 0.18) THEN
              msno =  4
              dzsno(c,3) = dzsno(c,3)/2.
              swice(c,3) = swice(c,3)/2.
              swliq(c,3) = swliq(c,3)/2.
              dzsno(c,4) = dzsno(c,3)
              swice(c,4) = swice(c,3)
              swliq(c,4) = swliq(c,3)
              tsno(c,4)  = tsno(c,3)
           END IF
        END IF
     END IF

     IF (msno > 3) THEN
        IF (dzsno(c,3) > 0.11) THEN
           drr = dzsno(c,3) - 0.11
           propor = drr/dzsno(c,3)
           zwice = propor*swice(c,3)
           zwliq = propor*swliq(c,3)
           propor = 0.11/dzsno(c,3)
           swice(c,3) = propor*swice(c,3)
           swliq(c,3) = propor*swliq(c,3)
           dzsno(c,3) = 0.11

           CALL Combo (dzsno(c,4), swliq(c,4), swice(c,4), tsno(c,4), drr, &
                zwliq, zwice, tsno(c,3))

           ! Subdivided a new layer
           IF (msno <= 4 .AND. dzsno(c,4) > 0.41) THEN
              msno = 5
              dzsno(c,4) = dzsno(c,4)/2.
              swice(c,4) = swice(c,4)/2.
              swliq(c,4) = swliq(c,4)/2.
              dzsno(c,5) = dzsno(c,4)
              swice(c,5) = swice(c,4)
              swliq(c,5) = swliq(c,4)
              tsno(c,5)  = tsno(c,4)
           END IF
        END IF
     END IF

     IF (msno > 4) THEN
        IF (dzsno(c,4) > 0.23) THEN
           drr = dzsno(c,4) - 0.23
           propor = drr/dzsno(c,4)
           zwice = propor*swice(c,4)
           zwliq = propor*swliq(c,4)
           propor = 0.23/dzsno(c,4)
           swice(c,4) = propor*swice(c,4)
           swliq(c,4) = propor*swliq(c,4)
           dzsno(c,4) = 0.23

           CALL Combo (dzsno(c,5), swliq(c,5), swice(c,5), tsno(c,5), drr, &
                zwliq, zwice, tsno(c,4))
        END IF
     END IF

     snl(c) = -msno

  END DO

  DO j = -nlevsnow+1,0
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN
           dz(c,j)         = dzsno(c,j-snl(c))
           h2osoi_ice(c,j) = swice(c,j-snl(c))
           h2osoi_liq(c,j) = swliq(c,j-snl(c))
           t_soisno(c,j)   = tsno(c,j-snl(c))
        END IF
     END DO
  END DO

  DO j = 0, -nlevsnow+1, -1
     !dir$ concurrent
     !cdir nodep
     DO fc = 1, num_snowc
        c = filter_snowc(fc)
        IF (j >= snl(c)+1) THEN
           z(c,j)    = zi(c,j) - 0.5*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
        END IF
     END DO
  END DO

END SUBROUTINE DivideSnowLayers

SUBROUTINE Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
  !
  ! !DESCRIPTION:
  ! Combines two elements and returns the following combined
  ! variables: dz, t, wliq, wice.
  ! The combined temperature is based on the equation:
  ! the sum of the enthalpies of the two elements =
  ! that of the combined element.
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in)    :: dz2
  REAL(r8), INTENT(in)    :: wliq2
  REAL(r8), INTENT(in)    :: wice2
  REAL(r8), INTENT(in)    :: t2
  REAL(r8), INTENT(inout) :: dz
  REAL(r8), INTENT(inout) :: wliq
  REAL(r8), INTENT(inout) :: wice
  REAL(r8), INTENT(inout) :: t

  REAL(r8) :: dzc
  REAL(r8) :: wliqc
  REAL(r8) :: wicec
  REAL(r8) :: tc
  REAL(r8) :: h
  REAL(r8) :: h2
  REAL(r8) :: hc
  !-----------------------------------------------------------------------

  dzc = dz+dz2
  wicec = (wice+wice2)
  wliqc = (wliq+wliq2)
  h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
  h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

  hc = h + h2
  IF(hc < 0.)THEN
     tc = tfrz + hc/(cpice*wicec + cpliq*wliqc)
  ELSE IF (hc.LE.hfus*wliqc) THEN
     tc = tfrz
  ELSE
     tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)
  END IF

  dz = dzc
  wice = wicec
  wliq = wliqc
  t = tc

END SUBROUTINE Combo

SUBROUTINE BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec,snl, & !i
     num_snowc, filter_snowc, &                   !o
     num_nosnowc, filter_nosnowc)                 !o
  !
  ! !DESCRIPTION:
  ! Constructs snow filter for use in vectorized loops for snow hydrology.
  !
  ! !USES:
  !    use clmtype
  !
  ! !ARGUMENTS:

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER, INTENT(in)  :: lbc, ubc
  INTEGER, INTENT(in)  :: num_nolakec
  INTEGER, INTENT(in)  :: filter_nolakec(ubc-lbc+1)
  INTEGER, INTENT(in)  :: snl(1)
  INTEGER, INTENT(out) :: num_snowc
  INTEGER, INTENT(out) :: filter_snowc(ubc-lbc+1)
  INTEGER, INTENT(out) :: num_nosnowc
  INTEGER, INTENT(out) :: filter_nosnowc(ubc-lbc+1)
  !
  ! !CALLED FROM:
  ! subroutine Hydrology2 in Hydrology2Mod
  ! subroutine CombineSnowLayers in this module
  !
  ! !REVISION HISTORY:
  ! 2003 July 31: Forrest Hoffman
  !
  ! !LOCAL VARIABLES:
  ! local pointers to implicit in arguments
  !
  !EOP
  !
  ! !OTHER LOCAL VARIABLES:
  INTEGER  :: fc, c
  !-----------------------------------------------------------------------


  ! Build snow/no-snow filters for other subroutines

  num_snowc = 0
  num_nosnowc = 0
  DO fc = 1, num_nolakec
     c = filter_nolakec(fc)
     IF (snl(c) < 0) THEN
        num_snowc = num_snowc + 1
        filter_snowc(num_snowc) = c
     ELSE
        num_nosnowc = num_nosnowc + 1
        filter_nosnowc(num_nosnowc) = c
     END IF
  END DO

END SUBROUTINE BuildSnowFilter

SUBROUTINE FrictionVelocity(pgridcell,forc_hgt,forc_hgt_u,        & !i
     forc_hgt_t,forc_hgt_q,                  & !i
     lbp, ubp, fn, filterp,                  & !i
     displa, z0m, z0h, z0q,                  & !i
     obu, iter, ur, um,                      & !i
     ustar,temp1, temp2, temp12m, temp22m,   & !o
     u10,fv,                                 & !o
     fm)  !i&o

  !=============================================================================
  ! !DESCRIPTION:
  ! Calculation of the friction velocity, relation for potential
  ! temperature and humidity profiles of surface boundary layer.
  ! The scheme is based on the work of Zeng et al. (1998):
  ! Intercomparison of bulk aerodynamic algorithms for the computation
  ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
  ! Vol. 11, 2628-2644.
  !
  ! !REVISION HISTORY:
  ! 15 September 1999: Yongjiu Dai; Initial code
  ! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
  ! 12/19/01, Peter Thornton
  ! Added arguments to eliminate passing clm derived type into this function.
  ! Created by Mariana Vertenstein
  !============================================================================
  ! !USES:
  ! use clmtype
  !!use clm_atmlnd, only : clm_a2l
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  ! INTEGER,PARAMETER  ::     lbp = 1
  ! INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  ! INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  INTEGER , INTENT(in) :: pgridcell(1)
  REAL(r8), INTENT(in) :: forc_hgt(1)
  REAL(r8), INTENT(in) :: forc_hgt_u(1)
  REAL(r8), INTENT(in) :: forc_hgt_t(1)
  REAL(r8), INTENT(in) :: forc_hgt_q(1)
  INTEGER , INTENT(in)  :: lbp, ubp
  INTEGER , INTENT(in)  :: fn
  INTEGER , INTENT(in)  :: filterp(fn)
  REAL(r8), INTENT(in)  :: displa(lbp:ubp)
  REAL(r8), INTENT(in)  :: z0m(lbp:ubp)
  REAL(r8), INTENT(in)  :: z0h(lbp:ubp)
  REAL(r8), INTENT(in)  :: z0q(lbp:ubp)
  REAL(r8), INTENT(in)  :: obu(lbp:ubp)
  INTEGER,  INTENT(in)  :: iter
  REAL(r8), INTENT(in)  :: ur(lbp:ubp)
  REAL(r8), INTENT(in)  :: um(lbp:ubp)


  REAL(r8), INTENT(out) :: ustar(lbp:ubp)
  REAL(r8), INTENT(out) :: temp1(lbp:ubp)
  REAL(r8), INTENT(out) :: temp12m(lbp:ubp)
  REAL(r8), INTENT(out) :: temp2(lbp:ubp)
  REAL(r8), INTENT(out) :: temp22m(lbp:ubp)
  REAL(r8), INTENT(out) :: u10(1)
  REAL(r8), INTENT(out) :: fv(1)

  REAL(r8), INTENT(inout) :: fm(lbp:ubp)

  REAL(r8), PARAMETER :: zetam = 1.574_r8
  REAL(r8), PARAMETER :: zetat = 0.465_r8
  INTEGER :: f
  INTEGER :: p
  INTEGER :: g
  REAL(r8):: zldis(lbp:ubp)
  REAL(r8):: zeta(lbp:ubp)

  REAL(r8) :: tmp1,tmp2,tmp3,tmp4
  REAL(r8) :: fmnew
  REAL(r8) :: fm10
  REAL(r8) :: zeta10
  REAL(r8),EXTERNAL ::StabilityFunc1,StabilityFunc2
  !#endif
  !------------------------------------------------------------------------------


  ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

  !#define PERGRO
  !#if (!defined PERGRO)

  !dir$ concurrent
  !cdir nodep
  DO f = 1, fn
     p = filterp(f)
     g = pgridcell(p)

     ! Wind profile

     zldis(p) = forc_hgt_u(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetam) THEN
        ustar(p) = vkc*um(p)/(LOG(-zetam*obu(p)/z0m(p))&
             - StabilityFunc1(-zetam) &
             + StabilityFunc1(z0m(p)/obu(p)) &
             + 1.14_r8*((-zeta(p))**0.333_r8-(zetam)**0.333_r8))
     ELSE IF (zeta(p) < 0._r8) THEN
        ustar(p) = vkc*um(p)/(LOG(zldis(p)/z0m(p))&
             - StabilityFunc1(zeta(p))&
             + StabilityFunc1(z0m(p)/obu(p)))
     ELSE IF (zeta(p) <=  1._r8) THEN
        ustar(p) = vkc*um(p)/(LOG(zldis(p)/z0m(p)) + 5._r8*zeta(p) -5._r8*z0m(p)/obu(p))
     ELSE
        ustar(p) = vkc*um(p)/(LOG(obu(p)/z0m(p))+5._r8-5._r8*z0m(p)/obu(p) &
             +(5._r8*LOG(zeta(p))+zeta(p)-1._r8))
     END IF

     ! Temperature profile

     zldis(p) = forc_hgt_t(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp1(p) = vkc/(LOG(-zetat*obu(p)/z0h(p))&
             - StabilityFunc2(-zetat) &
             + StabilityFunc2(z0h(p)/obu(p)) &
             + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp1(p) = vkc/(LOG(zldis(p)/z0h(p)) &
             - StabilityFunc2(zeta(p)) &
             + StabilityFunc2(z0h(p)/obu(p)))
     ELSE IF (zeta(p) <=  1._r8) THEN
        temp1(p) = vkc/(LOG(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
     ELSE
        temp1(p) = vkc/(LOG(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
             + (5._r8*LOG(zeta(p))+zeta(p)-1._r8))
     END IF

     ! Humidity profile

     IF (forc_hgt_q(g) == forc_hgt_t(g) .AND. z0q(p) == z0h(p)) THEN
        temp2(p) = temp1(p)
     ELSE
        zldis(p) = forc_hgt_q(g)-displa(p)
        zeta(p) = zldis(p)/obu(p)
        IF (zeta(p) < -zetat) THEN
           temp2(p) = vkc/(LOG(-zetat*obu(p)/z0q(p)) &
                - StabilityFunc2(-zetat) &
                + StabilityFunc2(z0q(p)/obu(p)) &
                + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
        ELSE IF (zeta(p) < 0._r8) THEN
           temp2(p) = vkc/(LOG(zldis(p)/z0q(p)) &
                - StabilityFunc2(zeta(p)) &
                + StabilityFunc2(z0q(p)/obu(p)))
        ELSE IF (zeta(p) <=  1._r8) THEN
           temp2(p) = vkc/(LOG(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
        ELSE
           temp2(p) = vkc/(LOG(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p) &
                + (5._r8*LOG(zeta(p))+zeta(p)-1._r8))
        END IF
     ENDIF

     ! Temperature profile applied at 2-m

     zldis(p) = 2.0_r8 + z0h(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp12m(p) = vkc/(LOG(-zetat*obu(p)/z0h(p))&
             - StabilityFunc2(-zetat) &
             + StabilityFunc2(z0h(p)/obu(p)) &
             + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp12m(p) = vkc/(LOG(zldis(p)/z0h(p)) &
             - StabilityFunc2(zeta(p))  &
             + StabilityFunc2(z0h(p)/obu(p)))
     ELSE IF (zeta(p) <=  1._r8) THEN
        temp12m(p) = vkc/(LOG(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
     ELSE
        temp12m(p) = vkc/(LOG(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
             + (5._r8*LOG(zeta(p))+zeta(p)-1._r8))
     END IF

     ! Humidity profile applied at 2-m

     IF (z0q(p) == z0h(p)) THEN
        temp22m(p) = temp12m(p)
     ELSE
        zldis(p) = 2.0_r8 + z0q(p)
        zeta(p) = zldis(p)/obu(p)
        IF (zeta(p) < -zetat) THEN
           temp22m(p) = vkc/(LOG(-zetat*obu(p)/z0q(p)) - &
                StabilityFunc2(-zetat) + StabilityFunc2(z0q(p)/obu(p)) &
                + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
        ELSE IF (zeta(p) < 0._r8) THEN
           temp22m(p) = vkc/(LOG(zldis(p)/z0q(p)) - &
                StabilityFunc2(zeta(p))+StabilityFunc2(z0q(p)/obu(p)))
        ELSE IF (zeta(p) <=  1._r8) THEN
           temp22m(p) = vkc/(LOG(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
        ELSE
           temp22m(p) = vkc/(LOG(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p) &
                + (5._r8*LOG(zeta(p))+zeta(p)-1._r8))
        END IF
     END IF

     !#defined DGVM
     !#defined DUST
     !#if (defined DGVM) || (defined DUST)
     ! diagnose 10-m wind for dust model (dstmbl.F)
     ! Notes from C. Zender's dst.F:
     ! According to Bon96 p. 62, the displacement height d (here displa) is
     ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
     ! Therefore d <= 0.034*z1 and may safely be neglected.
     ! Code from LSM routine SurfaceTemperature was used to obtain u10

     zldis(p) = forc_hgt_u(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (MIN(zeta(p), 1._r8) < 0._r8) THEN
        tmp1 = (1._r8 - 16._r8*MIN(zeta(p),1._r8))**0.25_r8
        tmp2 = LOG((1._r8+tmp1*tmp1)/2._r8)
        tmp3 = LOG((1._r8+tmp1)/2._r8)
        fmnew = 2._r8*tmp3 + tmp2 - 2._r8*ATAN(tmp1) + 1.5707963_r8
     ELSE
        fmnew = -5._r8*MIN(zeta(p),1._r8)
     ENDIF
     IF (iter == 1) THEN
        fm(p) = fmnew
     ELSE
        fm(p) = 0.5_r8 * (fm(p)+fmnew)
     END IF
     zeta10 = MIN(10._r8/obu(p), 1._r8)
     IF (zeta(p) == 0._r8) zeta10 = 0._r8
     IF (zeta10 < 0._r8) THEN
        tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
        tmp2 = LOG((1.0_r8 + tmp1*tmp1)/2.0_r8)
        tmp3 = LOG((1.0_r8 + tmp1)/2.0_r8)
        fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*ATAN(tmp1) + 1.5707963_r8
     ELSE                ! not stable
        fm10 = -5.0_r8 * zeta10
     END IF
     tmp4 = LOG(forc_hgt(g) / 10._r8)
     u10(p) = ur(p) - ustar(p)/vkc * (tmp4 - fm(p) + fm10)
     fv(p)  = ustar(p)
     !#endif

  END DO
  !#endif

  !#define PERGRO
  !#if (defined PERGRO)

  !===============================================================================
  ! The following only applies when PERGRO is defined
  !===============================================================================

  !dir$ concurrent
  !cdir nodep
  DO f = 1, fn
     p = filterp(f)
     g = pgridcell(p)

     zldis(p) = forc_hgt_u(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetam) THEN           ! zeta < -1
        ustar(p) = vkc * um(p) / LOG(-zetam*obu(p)/z0m(p))
     ELSE IF (zeta(p) < 0._r8) THEN         ! -1 <= zeta < 0
        ustar(p) = vkc * um(p) / LOG(zldis(p)/z0m(p))
     ELSE IF (zeta(p) <= 1._r8) THEN        !  0 <= ztea <= 1
        ustar(p)=vkc * um(p)/LOG(zldis(p)/z0m(p))
     ELSE                             !  1 < zeta, phi=5+zeta
        ustar(p)=vkc * um(p)/LOG(obu(p)/z0m(p))
     ENDIF

     zldis(p) = forc_hgt_t(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp1(p)=vkc/LOG(-zetat*obu(p)/z0h(p))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp1(p)=vkc/LOG(zldis(p)/z0h(p))
     ELSE IF (zeta(p) <= 1._r8) THEN
        temp1(p)=vkc/LOG(zldis(p)/z0h(p))
     ELSE
        temp1(p)=vkc/LOG(obu(p)/z0h(p))
     END IF

     zldis(p) = forc_hgt_q(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp2(p)=vkc/LOG(-zetat*obu(p)/z0q(p))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp2(p)=vkc/LOG(zldis(p)/z0q(p))
     ELSE IF (zeta(p) <= 1._r8) THEN
        temp2(p)=vkc/LOG(zldis(p)/z0q(p))
     ELSE
        temp2(p)=vkc/LOG(obu(p)/z0q(p))
     END IF

     zldis(p) = 2.0_r8 + z0h(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp12m(p)=vkc/LOG(-zetat*obu(p)/z0h(p))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp12m(p)=vkc/LOG(zldis(p)/z0h(p))
     ELSE IF (zeta(p) <= 1._r8) THEN
        temp12m(p)=vkc/LOG(zldis(p)/z0h(p))
     ELSE
        temp12m(p)=vkc/LOG(obu(p)/z0h(p))
     END IF

     zldis(p) = 2.0_r8 + z0q(p)
     zeta(p) = zldis(p)/obu(p)
     IF (zeta(p) < -zetat) THEN
        temp22m(p)=vkc/LOG(-zetat*obu(p)/z0q(p))
     ELSE IF (zeta(p) < 0._r8) THEN
        temp22m(p)=vkc/LOG(zldis(p)/z0q(p))
     ELSE IF (zeta(p) <= 1._r8) THEN
        temp22m(p)=vkc/LOG(zldis(p)/z0q(p))
     ELSE
        temp22m(p)=vkc/LOG(obu(p)/z0q(p))
     END IF

     !#defined DGVM
     !#defined DUST
     !#if (defined DGVM) || (defined DUST)
     ! diagnose 10-m wind for dust model (dstmbl.F)
     ! Notes from C. Zender's dst.F:
     ! According to Bon96 p. 62, the displacement height d (here displa) is
     ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
     ! Therefore d <= 0.034*z1 and may safely be neglected.
     ! Code from LSM routine SurfaceTemperature was used to obtain u10

     zldis(p) = forc_hgt_u(g)-displa(p)
     zeta(p) = zldis(p)/obu(p)
     IF (MIN(zeta(p), 1._r8) < 0._r8) THEN
        tmp1 = (1._r8 - 16._r8*MIN(zeta(p),1._r8))**0.25_r8
        tmp2 = LOG((1._r8+tmp1*tmp1)/2._r8)
        tmp3 = LOG((1._r8+tmp1)/2._r8)
        fmnew = 2._r8*tmp3 + tmp2 - 2._r8*ATAN(tmp1) + 1.5707963_r8
     ELSE
        fmnew = -5._r8*MIN(zeta(p),1._r8)
     ENDIF
     IF (iter == 1) THEN
        fm(p) = fmnew
     ELSE
        fm(p) = 0.5_r8 * (fm(p)+fmnew)
     END IF
     zeta10 = MIN(10._r8/obu(p), 1._r8)
     IF (zeta(p) == 0._r8) zeta10 = 0._r8
     IF (zeta10 < 0._r8) THEN
        tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
        tmp2 = LOG((1.0_r8 + tmp1*tmp1)/2.0_r8)
        tmp3 = LOG((1.0_r8 + tmp1)/2.0_r8)
        fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*ATAN(tmp1) + 1.5707963_r8
     ELSE                ! not stable
        fm10 = -5.0_r8 * zeta10
     END IF
     tmp4 = LOG(forc_hgt(g) / 10._r8)
     u10(p) = ur(p) - ustar(p)/vkc * (tmp4 - fm(p) + fm10)
     fv(p)  = ustar(p)
     !#endif
  END DO

  !#endif

END SUBROUTINE FrictionVelocity


! !IROUTINE: StabilityFunc
!
! !INTERFACE:
REAL(SELECTED_REAL_KIND(12)) FUNCTION StabilityFunc1(zeta)

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in) :: zeta

  REAL(r8) :: chik, chik2
  !------------------------------------------------------------------------------

  chik2 = SQRT(1._r8-16._r8*zeta)
  chik = SQRT(chik2)
  StabilityFunc1 = 2._r8*LOG((1._r8+chik)*0.5_r8) &
                                !Changed to pie, Zack Subin, 7/9/08
       + LOG((1._r8+chik2)*0.5_r8)-2._r8*ATAN(chik)+pie*0.5_r8

END FUNCTION StabilityFunc1


!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc2
!
! !INTERFACE:
REAL(SELECTED_REAL_KIND(12)) FUNCTION StabilityFunc2(zeta)

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in) :: zeta
  REAL(r8) :: chik2
  !------------------------------------------------------------------------------

  chik2 = SQRT(1._r8-16._r8*zeta)
  StabilityFunc2 = 2._r8*LOG((1._r8+chik2)*0.5_r8)

END FUNCTION StabilityFunc2

!-----------------------------------------------------------------------
!BOP

!
! !IROUTINE: MoninObukIni
!
! !INTERFACE:
SUBROUTINE MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
  !
  ! !DESCRIPTION:
  ! Initialization of the Monin-Obukhov length.
  ! The scheme is based on the work of Zeng et al. (1998):
  ! Intercomparison of bulk aerodynamic algorithms for the computation
  ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
  ! Vol. 11, 2628-2644.
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  REAL(r8), INTENT(in)  :: ur
  REAL(r8), INTENT(in)  :: thv
  REAL(r8), INTENT(in)  :: dthv
  REAL(r8), INTENT(in)  :: zldis
  REAL(r8), INTENT(in)  :: z0m
  REAL(r8), INTENT(out) :: um
  REAL(r8), INTENT(out) :: obu

  REAL(r8) :: wc
  REAL(r8) :: rib
  REAL(r8) :: zeta
  REAL(r8) :: ustar


  ! Initial values of u* and convective velocity

  ustar=0.06_r8
  wc=0.5_r8
  IF (dthv >= 0._r8) THEN
     um=MAX(ur,0.1_r8)
  ELSE
     um=SQRT(ur*ur+wc*wc)
  ENDIF

  rib=grav*zldis*dthv/(thv*um*um)
  !#define PERGRO
  !#if (defined PERGRO)
  rib = 0._r8
  !#endif

  IF (rib >= 0._r8) THEN      ! neutral or stable
     zeta = rib*LOG(zldis/z0m)/(1._r8-5._r8*MIN(rib,0.19_r8))
     zeta = MIN(2._r8,MAX(zeta,0.01_r8 ))
  ELSE                     ! unstable
     zeta=rib*LOG(zldis/z0m)
     zeta = MAX(-100._r8,MIN(zeta,-0.01_r8 ))
  ENDIF

  obu=zldis/zeta

END SUBROUTINE MoninObukIni

SUBROUTINE LakeDebug( str )

  IMPLICIT NONE

  REAL    , PARAMETER :: r_d          = 287.
  REAL    , PARAMETER :: cp           = 7.*r_d/2.
  REAL    , PARAMETER :: rcp          = r_d/cp

  INTEGER, PARAMETER ::  r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: nlevsoil     =  10   ! number of soil layers
  INTEGER, PARAMETER :: nlevlake     =  25   ! number of lake layers
  INTEGER, PARAMETER :: nlevsnow     =   5   ! maximum number of snow layers

  INTEGER,PARAMETER  ::     lbp = 1
  INTEGER,PARAMETER  ::     ubp = 1
  INTEGER,PARAMETER  ::     lbc = 1                        ! column-index bounds
  INTEGER,PARAMETER  ::     ubc = 1

  INTEGER,PARAMETER  ::     num_shlakec       = 1
  INTEGER,PARAMETER  ::     filter_shlakec(1) = 1
  INTEGER,PARAMETER  ::     num_shlakep       = 1
  INTEGER,PARAMETER  ::     filter_shlakep(1) = 1
  INTEGER,PARAMETER  ::     pcolumn(1)        = 1
  INTEGER,PARAMETER  ::     pgridcell(1)      = 1
  INTEGER,PARAMETER  ::     cgridcell(1)      = 1
  INTEGER,PARAMETER  ::     clandunit(1)      = 1

  INTEGER,PARAMETER  ::     begg = 1
  INTEGER,PARAMETER  ::     endg = 1
  INTEGER,PARAMETER  ::     begl = 1
  INTEGER,PARAMETER  ::     endl = 1
  INTEGER,PARAMETER  ::     begc = 1
  INTEGER,PARAMETER  ::     endc = 1
  INTEGER,PARAMETER  ::     begp = 1
  INTEGER,PARAMETER  ::     endp = 1

  INTEGER,PARAMETER  ::     column    =1
  LOGICAL,PARAMETER  ::     lakpoi(1) = .TRUE.




  REAL(r8), PARAMETER :: vkc    = 0.4_r8
  REAL(r8), PARAMETER :: pie    = 3.141592653589793_r8
  REAL(r8), PARAMETER :: grav   = 9.80616_r8
  REAL(r8), PARAMETER :: sb     = 5.67e-8_r8
  REAL(r8), PARAMETER :: tfrz   = 273.16_r8
  REAL(r8), PARAMETER :: denh2o = 1.000e3_r8
  REAL(r8), PARAMETER :: denice = 0.917e3_r8
  REAL(r8), PARAMETER :: cpice  = 2.11727e3_r8
  REAL(r8), PARAMETER :: cpliq  = 4.188e3_r8
  REAL(r8), PARAMETER :: hfus   = 3.337e5_r8
  REAL(r8), PARAMETER :: hvap   = 2.501e6_r8
  REAL(r8), PARAMETER :: hsub   = 2.501e6_r8+3.337e5_r8
  REAL(r8), PARAMETER :: rair   = 287.0423_r8
  REAL(r8), PARAMETER :: cpair  = 1.00464e3_r8
  REAL(r8), PARAMETER :: tcrit  = 2.5
  REAL(r8), PARAMETER :: tkwat  = 0.6
  REAL(r8), PARAMETER :: tkice  = 2.290
  REAL(r8), PARAMETER :: tkairc = 0.023
  REAL(r8), PARAMETER :: bdsno = 250.

  REAL(r8),  PARAMETER :: spval = 1.e36

  REAL(r8), PARAMETER  ::     depth_c = 300.


  REAL(r8), PARAMETER :: wimp   = 0.05
  REAL(r8), PARAMETER :: ssi    = 0.033
  REAL(r8), PARAMETER :: cnfac  = 0.5


  INTEGER,PARAMETER :: istsoil = 1

  REAL(r8) :: dtime = 60

  REAL(r8) :: zlak(1:nlevlake)
  REAL(r8) :: dzlak(1:nlevlake)
  REAL(r8) :: zsoi(1:nlevsoil)
  REAL(r8) :: dzsoi(1:nlevsoil)
  REAL(r8) :: zisoi(0:nlevsoil)


  REAL(r8) :: sand(19)
  REAL(r8) :: clay(19)
  INTEGER :: i

  DATA(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
       10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./

  DATA(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
       33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./


  REAL(r8) :: watsat(1,nlevsoil)
  REAL(r8) :: tksatu(1,nlevsoil)
  REAL(r8) :: tkmg(1,nlevsoil)
  REAL(r8) :: tkdry(1,nlevsoil)
  REAL(r8) :: csol(1,nlevsoil)

  CHARACTER*(*), str

!!!  CALL wrf_debug( 0 , TRIM(str) )

END SUBROUTINE LakeDebug
