module globals
    
    implicit none
    public
    integer, parameter          ::  dp = selected_real_kind(p=13, r=200),                                   &
                                    Noutvarsmax = 20,                                                       &
                                    Noutmax = 1000,                                                         &
                                    Nhpmax = 20,                                                            &
                                    Ngrainmax = 10000,                                                      &
                                    Nslipsmax = 100
    
    real(kind=8), parameter     ::  tol = 1.d-12,                                                           &
                                    tolinv = 1.d12,                                                         &
                                    pi = 4.d0*atan(1.d0)
    
    integer                     ::  Nsteps, Nout, Nclus, output_steps(Noutmax), outputvars(Noutvarsmax),    &
                                    Nhv, Nhp, outINC, Ngrains, grain_interaction, hardening_law,            &
                                    sslookup(12288,5), rotate_boundary, Nslips
                                    
    real(kind=8)                ::  R1(3,3), R2(3,3), R3(3,3), symR1(3,3), 		                            &
									symR2(3,3), symR3(3,3), crss_mean, dt, relax_penalties(3), totalslip,   &
									eps4jacobian, Fglob(3,3), hparam(Nhpmax), totalslipG1, totalslipG2
									
    real(kind=8), allocatable   ::  P(:,:,:), Omega(:,:,:) 
    
    logical                     ::  stress_only, is_output_step


! mixed BC defining variables
    integer                     ::  ind_d(5), ind_sdir(5), ind_sabs(5)
    real(kind=8)                ::  W_prescribed(3,3), d_prescribed(5), sdir_prescribed(5), sabs_prescribed(5)
    
    
! output variables
    integer, allocatable        ::  out_activesIDG1(:,:,:),                                                 &
                                    out_num_activesG1(:,:),                                                 &
                                    out_activesIDG2(:,:,:),                                                 &
                                    out_num_activesG2(:,:)
    
    real(kind=8), allocatable   ::  out_average_stress(:,:)     ,                                           &
                                    out_average_slip(:)         ,                                           &
                                    out_stress_locG1(:,:,:)     ,                                           &
                                    out_stress_globG1(:,:,:)    ,                                           &
                                    out_euler_anglesG1(:,:,:)   ,                                           &
                                    out_crssG1(:,:,:)           ,                                           &
                                    out_slipratesG1(:,:,:)      ,                                           &
                                    out_statevarG1(:,:,:)       ,                                           &
                                    out_stress_locG2(:,:,:)     ,                                           &
                                    out_stress_globG2(:,:,:)    ,                                           &
                                    out_euler_anglesG2(:,:,:)   ,                                           &
                                    out_crssG2(:,:,:)           ,                                           &
                                    out_slipratesG2(:,:,:)      ,                                           &
                                    out_statevarG2(:,:,:)       ,                                           &
                                    out_relaxsliprates(:,:,:)
                                    
! orientations and hardening variables
    real(kind=8), allocatable   ::  Qg1_list(:,:,:), Qg2_list(:,:,:), Qb_list(:,:,:), Rg1(:,:,:),           &
                                    Q0g1_list(:,:,:), Q0g2_list(:,:,:), Q0b_list(:,:,:), Rg2(:,:,:),        &
                                    hvarG1_list(:,:), hvarG2_list(:,:), rhs(:)
                                    
    end module globals
    
module SCYLglobals

    use globals, only: Nslipsmax

    implicit none
    public
    integer                     ::  SCYLon
    real(kind=8)                ::  SCYL_CRSS(Nslipsmax), SCYL_CRSSg1(Nslipsmax), SCYL_CRSSg2(Nslipsmax),  &
                                    SCYLexp, SCYL_Dp(5), SCYL_S1solved(3,3), SCYL_S2solved(3,3),           &
                                    SCYL_QbQg1T(3,3), SCYL_QbQg2T(3,3)

end module SCYLglobals