module crystal_plasticity

    use globals
    use SCYLglobals
    use utils
    use minpack
    use linesearch
    
    implicit none
    
    contains
    
    subroutine init_slipsystems(cs)
    ! Subroutine defines slip systems, build Schmid matrix,its symmetric part "P" 
    ! and skew part "Omega", both global variables defined in module "globals"
    !
    ! Inputs: cs    -   specifies the set of slip systems used for the calculations as:
    !                   dictionary of possible crystalographic slip systems specified by "cs"
    !                              1              2               3               4     
    !                   1,   '{111} <110>'
    !                   12,  '{111} <110>'  '{011} <011>' 
    !                   13,  '{111} <110>'                  '{100} <011>'
    !                   14,  '{111} <110>'                                  '{112} <011>'
    !                   123, '{111} <110>'  '{011} <011>'   '{100} <011>'
    !                   124, '{111} <110>'  '{011} <011>'                   '{112} <011>'
    !                   134, '{111} <110>'                  '{100} <011>'   '{112} <011>'
    !                   1234,'{111} <110>'  '{011} <011>'   '{100} <011>'   '{112} <011>'


        real(kind=8)        :: n(36,3), b(36,3), schmid(3,3,36), schmid_sym(3,3,36), schmid_skw(3,3,36)
        integer             :: i, cs, s
        real(kind=8), parameter     ::  inv3 = 1.d0/sqrt(3.d0),  &
                                        inv2 = 1.d0/sqrt(2.d0),  &
                                        inv6 = 1.d0/sqrt(6.d0)
        
        intent(in)          :: cs
    
        ! slip systems described by slip plane normal "n" and slip direction "b"
        ! slip plane normals "n"
        ! cs = 1 (FCC octahedral slips)
        n(1,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(2,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(3,:) = (/ 1.d0,  1.d0, -1.d0/)
        n(4,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(5,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(6,:) = (/ 1.d0, -1.d0, -1.d0/)
        n(7,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(8,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(9,:) = (/ 1.d0, -1.d0,  1.d0/)
        n(10,:)= (/ 1.d0,  1.d0,  1.d0/)
        n(11,:)= (/ 1.d0,  1.d0,  1.d0/)
        n(12,:)= (/ 1.d0,  1.d0,  1.d0/)
        ! cs = 2
        n(13,:)= (/ 1.d0,  1.d0,  0.d0/)
        n(14,:)= (/ 1.d0, -1.d0,  0.d0/)
        n(15,:)= (/ 1.d0,  0.d0,  1.d0/)
        n(16,:)= (/ 1.d0,  0.d0, -1.d0/)
        n(17,:)= (/ 0.d0,  1.d0,  1.d0/)
        n(18,:)= (/ 0.d0,  1.d0, -1.d0/)
        ! cs = 3
        n(19,:)= (/ 1.d0,  0.d0,  0.d0/)
        n(20,:)= (/ 1.d0,  0.d0,  0.d0/)
        n(21,:)= (/ 0.d0,  1.d0,  0.d0/)
        n(22,:)= (/ 0.d0,  1.d0,  0.d0/)
        n(23,:)= (/ 0.d0,  0.d0,  1.d0/)
        n(24,:)= (/ 0.d0,  0.d0,  1.d0/)
        ! cs = 4
        n(25,:)= (/-1.d0,  1.d0,  2.d0/)
        n(26,:)= (/-1.d0,  1.d0, -2.d0/)
        n(27,:)= (/ 1.d0,  1.d0,  2.d0/)
        n(28,:)= (/ 1.d0,  1.d0, -2.d0/)
        n(29,:)= (/ 1.d0,  2.d0, -1.d0/)
        n(30,:)= (/ 1.d0, -2.d0, -1.d0/)
        n(31,:)= (/ 1.d0,  2.d0,  1.d0/)
        n(32,:)= (/ 1.d0, -2.d0,  1.d0/)
        n(33,:)= (/ 2.d0,  1.d0, -1.d0/)
        n(34,:)= (/-2.d0,  1.d0, -1.d0/)
        n(35,:)= (/ 2.d0,  1.d0,  1.d0/)
        n(36,:)= (/-2.d0,  1.d0,  1.d0/)
        
        ! slip directions "b"
        ! cs = 1
        b(1,:) = (/ 0.d0,  1.d0,  1.d0/)
        b(2,:) = (/ 1.d0,  0.d0,  1.d0/)
        b(3,:) = (/ 1.d0, -1.d0,  0.d0/)
        b(4,:) = (/ 0.d0,  1.d0, -1.d0/)
        b(5,:) = (/ 1.d0,  0.d0,  1.d0/)
        b(6,:) = (/ 1.d0,  1.d0,  0.d0/)
        b(7,:) = (/ 0.d0,  1.d0,  1.d0/)
        b(8,:) = (/ 1.d0,  0.d0, -1.d0/)
        b(9,:) = (/ 1.d0,  1.d0,  0.d0/)
        b(10,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(11,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(12,:)= (/ 1.d0, -1.d0,  0.d0/)
        ! cs = 2
        b(13,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(14,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(15,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(16,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(17,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(18,:)= (/ 0.d0,  1.d0,  1.d0/)
        ! cs = 3
        b(19,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(20,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(21,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(22,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(23,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(24,:)= (/ 1.d0, -1.d0,  0.d0/)
        ! cas = 4
        b(25,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(26,:)= (/ 1.d0,  1.d0,  0.d0/)
        b(27,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(28,:)= (/ 1.d0, -1.d0,  0.d0/)
        b(29,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(30,:)= (/ 1.d0,  0.d0,  1.d0/)
        b(31,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(32,:)= (/ 1.d0,  0.d0, -1.d0/)
        b(33,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(34,:)= (/ 0.d0,  1.d0,  1.d0/)
        b(35,:)= (/ 0.d0,  1.d0, -1.d0/)
        b(36,:)= (/ 0.d0,  1.d0, -1.d0/)
        
        
        n(1:12,:)  = n(1:12,:) *inv3
        n(13:18,:) = n(13:18,:)*inv2
        n(25:36,:) = n(25:36,:)*inv6
        b = b*inv2
    
        do i=1, 36
            schmid(:,:,i)     = tensordot(b(i,:), n(i,:))
            schmid_sym(:,:,i) = getsym(schmid(:,:,i))
            schmid_skw(:,:,i) = getskw(schmid(:,:,i))
        end do
        
        if (allocated(P)) deallocate(P, Omega, STAT=s)
        if (s /= 0) write(*,*) 'Deallocation of P and Omega unsuccsessful.'

        select case (cs)
        case (1)
            allocate(P(3,3,12*2), Omega(3,3,12*2), STAT=s)
            P(:,:,:12)       =  schmid_sym(:,:,:12)
            P(:,:,13:24)     = -schmid_sym(:,:,:12)
            Omega(:,:,:12)   =  schmid_skw(:,:,:12)
            Omega(:,:,13:24) = -schmid_skw(:,:,:12)
            Nslips = 24
        case (12)
            allocate(P(3,3,18*2), Omega(3,3,18*2), STAT=s)
            P(:,:,:12)       =  schmid_sym(:,:,:12)
            P(:,:,13:18)     =  schmid_sym(:,:,13:18)
            P(:,:,19:30)     = -schmid_sym(:,:,:12)
            P(:,:,31:36)     = -schmid_sym(:,:,13:18)
            Omega(:,:,:12)   =  schmid_skw(:,:,:12)
            Omega(:,:,13:18) =  schmid_skw(:,:,13:18)
            Omega(:,:,19:30) = -schmid_skw(:,:,:12)
            Omega(:,:,31:36) = -schmid_skw(:,:,13:18)
            Nslips = 36
        case (13)
            allocate(P(3,3,18*2), Omega(3,3,18*2), STAT=s)
            P(:,:,:12)       =  schmid_sym(:,:,:12)
            P(:,:,13:18)     =  schmid_sym(:,:,19:24)
            P(:,:,19:30)     = -schmid_sym(:,:,:12)
            P(:,:,31:36)     = -schmid_sym(:,:,19:24)
            Omega(:,:,:12)   =  schmid_skw(:,:,:12)
            Omega(:,:,13:18) =  schmid_skw(:,:,19:24)
            Omega(:,:,19:30) = -schmid_skw(:,:,:12)
            Omega(:,:,31:36) = -schmid_skw(:,:,19:24)
            Nslips = 36
        case (14)
            allocate(P(3,3,24*2), Omega(3,3,24*2), STAT=s)
            P(:,:,:12)       =  schmid_sym(:,:,:12)
            P(:,:,13:24)     =  schmid_sym(:,:,25:36)
            P(:,:,25:36)     = -schmid_sym(:,:,:12)
            P(:,:,37:48)     = -schmid_sym(:,:,25:36)
            Omega(:,:,:12)   =  schmid_skw(:,:,:12)
            Omega(:,:,13:24) =  schmid_skw(:,:,25:36)
            Omega(:,:,25:36) = -schmid_skw(:,:,:12)
            Omega(:,:,37:48) = -schmid_skw(:,:,25:36)
            Nslips = 48
        case (1234)
            allocate(P(3,3,36*2), Omega(3,3,36*2), STAT=s)
            P(:,:,:36)       =  schmid_sym
            P(:,:,37:)       = -schmid_sym
            Omega(:,:,:36)   =  schmid_skw
            Omega(:,:,37:72) = -schmid_skw   
            Nslips = 72
        end select
        if (s /= 0) write(*,*) 'Allocation of P and Omega unsuccessful.'
    
        ! definition of relaxation
        R1       = 0.d0
        R1(1,3)  = 1.d0
        R2       = 0.d0
        R2(2,3)  = 1.d0
		R3       = 0.d0
        R3(1,2)  = 1.d0

        symR1 = getsym(R1)
        symR2 = getsym(R2)
		symR3 = getsym(R3)
    
    end subroutine init_slipsystems
    

    subroutine init_orientation_taylor(pyQ_all)
    
    integer             ::  s(3), tmp(3)   
    real(kind=8)        ::  pyQ_all(:,:,:)
    
    if (allocated(Qg1_list)) deallocate(Qg1_list, STAT=s(1))
    if (allocated(Qg2_list)) deallocate(Qg2_list, STAT=s(2))
    
    allocate(Qg1_list, SOURCE=pyQ_all, STAT=s(3))
    
    if (any(s /= 0)) then
        write(*,*) 'Deallocation/allocation of initial orientations not successful.'
        read(*,*)
        stop
    end if
	
	tmp = shape(Qg1_list)
	Ngrains = tmp(3)
    
    end subroutine init_orientation_taylor
    
    
    subroutine init_hardening_taylor(pyinithvar_all, pyhparam, pyhardlaw)
    
    integer             ::  s(3), tmp(2), pyhardlaw
    real(kind=8)        ::  pyinithvar_all(:,:), pyhparam(:)
    
    if (allocated(hvarG1_list)) deallocate(hvarG1_list, STAT=s(1))
    if (allocated(hvarG2_list)) deallocate(hvarG2_list, STAT=s(2))
    
    allocate(hvarG1_list, SOURCE=pyinithvar_all, STAT=s(3))
    
    if (any(s /= 0)) then
        write(*,*) 'Deallocation/Allocation of initial hardening variables not successful.'
        read(*,*)
        stop
    end if
	
	hardening_law = pyhardlaw
	Nhp           = size(pyhparam)
	tmp           = shape(hvarG1_list)
	Nhv           = tmp(1)
	hparam(1:Nhp) = pyhparam
    
    end subroutine init_hardening_taylor
    
    
    subroutine alloc_outputvars_taylor()
    
        integer     ::  s(11), k
        
        s = 0
        do k = 1, size(outputvars)
            select case(outputvars(k))
				case(-1)
					exit
                case(1)
                    if (.not. allocated(out_average_stress))                & 
                    allocate(out_average_stress(5,Nout),            STAT=s(1))
                    out_average_stress = 0.d0
                case(2)
                    if (.not. allocated(out_average_slip))                  & 
                    allocate(out_average_slip(Nout),                STAT=s(2))
                    out_average_slip   = 0.d0
                case(3)
                    if (.not. allocated(out_stress_locG1))                  &
                    allocate(out_stress_locG1(5,Ngrains,Nout),       STAT=s(3))
                case(4)
                    if (.not. allocated(out_stress_globG1))                 &
                    allocate(out_stress_globG1(5,Ngrains,Nout),      STAT=s(4))
                case(5)
                    if (.not. allocated(out_euler_anglesG1))                &
                    allocate(out_euler_anglesG1(3,Ngrains,Nout),     STAT=s(5))
                case(6)
                    if (.not. allocated(out_crssG1))                        &
                    allocate(out_crssG1(Nslips,Ngrains,Nout),            STAT=s(6))
                case(7)
                    if (.not. allocated(out_slipratesG1))                   &
                    allocate(out_slipratesG1(Nslips,Ngrains,Nout),       STAT=s(7))
                case(8)
                    if (.not. allocated(out_activesIDG1))                   &
                    allocate(out_activesIDG1(Nslips,Ngrains,Nout),       STAT=s(8))
                case(9)
                    if (.not. allocated(out_num_activesG1))                 &
                    allocate(out_num_activesG1(Ngrains,Nout),        STAT=s(9))
                case(10)
                    if (.not. allocated(out_statevarG1))                    &
                    allocate(out_statevarG1(Nhv,Ngrains,Nout),       STAT=s(10))
            end select
        end do
                       
        if (any(s /= 0) == .true.) then
            write(*,*) 'Allocation of output variables was unsuccesfull.'
            read(*,*)
            stop
        end if
        
    end subroutine alloc_outputvars_taylor
    
    
    subroutine dealloc_outputvars_taylor()
    
        integer   :: s(11)

        if (allocated(out_average_stress)) deallocate(out_average_stress, STAT=s(1))
        if (allocated(out_average_slip)) deallocate(out_average_slip, STAT=s(2))
        if (allocated(out_stress_locG1)) deallocate(out_stress_locG1, STAT=s(3))
        if (allocated(out_stress_globG1)) deallocate(out_stress_globG1, STAT=s(4))
        if (allocated(out_euler_anglesG1)) deallocate(out_euler_anglesG1, STAT=s(5))
        if (allocated(out_crssG1)) deallocate(out_crssG1, STAT=s(6))
        if (allocated(out_slipratesG1)) deallocate(out_slipratesG1, STAT=s(7))
        if (allocated(out_activesIDG1)) deallocate(out_activesIDG1, STAT=s(8))
        if (allocated(out_num_activesG1)) deallocate(out_num_activesG1, STAT=s(9))
        if (allocated(out_statevarG1)) deallocate(out_statevarG1, STAT=s(10))
              
        if (any(s /= 0) == .true.) then
            write(*,*) 'Deallocation of output variables was unsuccesfull.'
            read(*,*)
            stop
        end if    
        
    end subroutine dealloc_outputvars_taylor
    
    
    subroutine init_orientation_alamel(pyQ1_all, pyQ2_all, pyQb_all)
    
    integer             ::  s(6), tmp(3)   
    real(kind=8)        ::  pyQ1_all(:,:,:), pyQ2_all(:,:,:), pyQb_all(:,:,:)
    
    if (allocated(Qg1_list)) deallocate(Qg1_list, STAT=s(1))
    if (allocated(Qg2_list)) deallocate(Qg2_list, STAT=s(2))
    if (allocated(Qb_list)) deallocate(Qb_list, STAT=s(3))
    
    allocate(Qg1_list, SOURCE=pyQ1_all, STAT=s(4))
    allocate(Qg2_list, SOURCE=pyQ2_all, STAT=s(5))
    allocate(Qb_list,  SOURCE=pyQb_all, STAT=s(6))
    
    if (any(s /= 0)) then
        write(*,*) 'Deallocation/allocation of initial orientations not successful.'
        read(*,*)
        stop
    end if
	
	tmp = shape(Qb_list)
	Nclus = tmp(3)
    
    end subroutine init_orientation_alamel
    
    
    subroutine init_hardening_alamel(pyinithvarG1_all, pyinithvarG2_all, pyhparam, pyhardlaw)
    
    integer             ::  s(4), tmp(2), pyhardlaw
    real(kind=8)        ::  pyinithvarG1_all(:,:), pyinithvarG2_all(:,:), pyhparam(:)
    
    if (allocated(hvarG1_list)) deallocate(hvarG1_list, STAT=s(1))
    if (allocated(hvarG2_list)) deallocate(hvarG2_list, STAT=s(2))
    
    allocate(hvarG1_list, SOURCE=pyinithvarG1_all, STAT=s(3))
    allocate(hvarG2_list, SOURCE=pyinithvarG2_all, STAT=s(4))
    
    if (any(s /= 0)) then
        write(*,*) 'Deallocation/Allocation of initial hardening variables not successful.'
        read(*,*)
        stop
    end if
	
	hardening_law = pyhardlaw
	Nhp           = size(pyhparam)
	tmp           = shape(hvarG1_list)
	Nhv           = tmp(1)
	hparam(1:Nhp) = pyhparam
    
    end subroutine init_hardening_alamel
    
    
    subroutine alloc_outputvars_alamel()
    
        integer     ::  s(11), k
        
        s = 0
        do k = 1, size(outputvars)
            select case(outputvars(k))
				case(-1)
					exit
                case(1)
                    if (.not. allocated(out_average_stress))                & 
                    allocate(out_average_stress(5,Nout),            STAT=s(1))
                    out_average_stress = 0.d0
                case(2)
                    if (.not. allocated(out_average_slip))                  & 
                    allocate(out_average_slip(Nout),                STAT=s(2))
                    out_average_slip   = 0.d0
                case(3)
                    if (.not. allocated(out_stress_locG1))                  &
                    allocate(out_stress_locG1(5,Nclus,Nout),                &
                             out_stress_locG2(5,Nclus,Nout),        STAT=s(3))
                case(4)
                    if (.not. allocated(out_stress_globG1))                 &
                    allocate(out_stress_globG1(5,Nclus,Nout),               &
                             out_stress_globG2(5,Nclus,Nout),       STAT=s(4))
                case(5)
                    if (.not. allocated(out_euler_anglesG1))                &
                    allocate(out_euler_anglesG1(3,Nclus,Nout),              &
                             out_euler_anglesG2(3,Nclus,Nout),      STAT=s(5))
                case(6)
                    if (.not. allocated(out_crssG1))                        &
                    allocate(out_crssG1(Nslips,Nclus,Nout),                 &
                             out_crssG2(Nslips,Nclus,Nout),         STAT=s(6))
                case(7)
                    if (.not. allocated(out_slipratesG1))                   &
                    allocate(out_slipratesG1(Nslips,Nclus,Nout),            &
                             out_slipratesG2(Nslips,Nclus,Nout),    STAT=s(7))
                case(8)
                    if (.not. allocated(out_activesIDG1))                   &
                    allocate(out_activesIDG1(Nslips,Nclus,Nout),            &
                             out_activesIDG2(Nslips,Nclus,Nout),    STAT=s(8))
                case(9)
                    if (.not. allocated(out_num_activesG1))                 &
                    allocate(out_num_activesG1(Nclus,Nout),                 &
                             out_num_activesG2(Nclus,Nout),         STAT=s(9))
                case(10)
                    if (.not. allocated(out_statevarG1))                    &
                    allocate(out_statevarG1(Nhv,Nclus,Nout),                &
                             out_statevarG2(Nhv,Nclus,Nout),        STAT=s(10))
                case(11)
                    if (.not. allocated(out_relaxsliprates))                &
                    allocate(out_relaxsliprates(2,Nclus,Nout),      STAT=s(11))
            end select
        end do
                       
        if (any(s /= 0) == .true.) then
            write(*,*) 'Allocation of output variables was unsuccesfull.'
            read(*,*)
            stop
        end if
        
    end subroutine alloc_outputvars_alamel
    
    
    subroutine dealloc_outputvars_alamel()
    
        integer   :: s(11)

        if (allocated(out_average_stress)) deallocate(out_average_stress, STAT=s(1))
        if (allocated(out_average_slip))   deallocate(out_average_slip, STAT=s(2))
        if (allocated(out_stress_locG1))   deallocate(out_stress_locG1, out_stress_locG2, STAT=s(3))
        if (allocated(out_stress_globG1))  deallocate(out_stress_globG1, out_stress_globG2, STAT=s(4))
        if (allocated(out_euler_anglesG1)) deallocate(out_euler_anglesG1, out_euler_anglesG2, STAT=s(5))
        if (allocated(out_crssG1))         deallocate(out_crssG1, out_crssG2, STAT=s(6))
        if (allocated(out_slipratesG1))    deallocate(out_slipratesG1, out_slipratesG2, STAT=s(7))
        if (allocated(out_activesIDG1))    deallocate(out_activesIDG1, out_activesIDG2, STAT=s(8))
        if (allocated(out_num_activesG1))  deallocate(out_num_activesG1, out_num_activesG2, STAT=s(9))
        if (allocated(out_relaxsliprates)) deallocate(out_relaxsliprates, STAT=s(10))
        if (allocated(out_statevarG1))     deallocate(out_statevarG1, out_statevarG2, STAT=s(11))
              
        if (any(s /= 0) == .true.) then
            write(*,*) 'Deallocation of output variables was unsuccesfull.'
            read(*,*)
            stop
        end if    
        
    end subroutine dealloc_outputvars_alamel
    
    
    function rescaleD(d, ind) result(dn)
    
    ! Function rescales those elements d[i] from 1D array d for which ind[i] = 0 so that
    ! VonMises norm of whole d array is 1. Function will not change elements d[i] for which
    ! ind[i] = 1.
    
        integer         ::  i, ind(5)
        real(kind=8)    ::  d(5), tmp, dn(5), k, vm, A, B, C, dis
        logical         ::  success = .true.
        
        vm = 3.d0/2.d0
        tmp = 0.d0
        do i=1, 5
            if (ind(i) == 1) then
                if (i > 2) then
                    vm = vm - 2.d0*d(i)**2
                else
                    vm = vm - d(i)**2
                end if
            else
                if (i > 2) then
                    tmp = tmp + 2.d0*d(i)**2
                else
                    tmp = tmp + d(i)**2
                end if
            end if
        end do
        
        if (tmp > tol) then
            ! if at least one of ind(i) is 0, something to be rescaled
            if (ind(1) == ind(2)) then
                if (ind(1) == 1) then
                    vm = vm - (-d(1)-d(2))**2
                else 
                    tmp = tmp + (-d(1)-d(2))**2
                end if
                if (vm > 0.d0) then
                    k = sqrt(vm/tmp)
                else
                    success = .false.
                end if
            else
                if (ind(1)) then
                    vm = vm - d(1)**2
                else
                    vm = vm - d(2)**2
                end if
                A = tmp
                B = 2.d0*d(1)*d(2)
                C = -vm
                dis = (B**2 - 4.d0*A*C)
                if (dis > 0.d0) then
                    k = (-B + sqrt(dis))/(2.d0*A)
                    if (k < 0.d0) success = .false.
                end if
            end if
                
            dn = d
            if (success) then
                do i=1, 5
                    if (ind(i)==0) then
                        dn(i) = dn(i)*k
                    end if
                end do  
            else
                write(*,*) 'Rescaling unsuccessfull.'
            end if    
        end if
            
    end function rescaleD
    
    function euler2rotm(ANG1d, ANG2d, ANG3d) result(Q0)
    ! ANG three Euler angles in degrees
        real(kind=8)    ::  ANG(3), Q0(3,3), ANG1, ANG2, ANG3, ANG1d, ANG2d, ANG3d

        ANG1 = deg2rad(ANG1d)
        ANG2 = deg2rad(ANG2d)
        ANG3 = deg2rad(ANG3d)

        Q0(1,1) = COS(ANG1)*COS(ANG3)-SIN(ANG1)*COS(ANG2)*SIN(ANG3)
        Q0(1,2) = SIN(ANG1)*COS(ANG3)+COS(ANG1)*COS(ANG2)*SIN(ANG3)
        Q0(1,3) = SIN(ANG2)*SIN(ANG3)
        Q0(2,1) =-COS(ANG1)*SIN(ANG3)-SIN(ANG1)*COS(ANG2)*COS(ANG3)
        Q0(2,2) =-SIN(ANG1)*SIN(ANG3)+COS(ANG1)*COS(ANG2)*COS(ANG3)
        Q0(2,3) = SIN(ANG2)*COS(ANG3)
        Q0(3,1) = SIN(ANG1)*SIN(ANG2)
        Q0(3,2) =-COS(ANG1)*SIN(ANG2)
        Q0(3,3) = COS(ANG2)

    end function euler2rotm
    
    
!C**********************************************************************
!C                         FUNCTION rotm2euler                         *
!C**********************************************************************
!C Computes the Euler angles associated with the orientation matrix    *
!C in terms of the three Euler angles: phi1 (rotation about Z1), PHI   *
!C (rotation about X2) and phi2 (rotation about Z3) - Bunge notation   *
!C**********************************************************************
      function rotm2euler(Q) result(ANG)

      real(kind=8)      ::  ANG(3), Q(3,3), ANG1, ANG2, ANG3, STH

      if (abs(Q(3,3)) < 1.D0) then
        ANG2 = acos(Q(3,3))
        STH  = sin(ANG2)
        ANG1 = atan2(Q(3,1)/STH,-Q(3,2)/STH)
        ANG3 = atan2(Q(1,3)/STH,Q(2,3)/STH)
      else
        ANG1 = atan2(Q(1,2),Q(1,1))
        ANG2 = 0.D0
        ANG3 = 0.D0
      end if
      ANG(1) = rad2deg(ANG1)
      ANG(2) = rad2deg(ANG2)
      ANG(3) = rad2deg(ANG3)

      end function rotm2euler
    
    
    function m2voigt(M) result(v)
        
        real(kind=8)    ::  M(3,3), v(5)
    
        v = (/ M(1,1), M(2,2), M(2,3), M(1,3), M(1,2) /)
    
    end function m2voigt
    
    function voigt2m(x) result(M)
        
        real(kind=8)    ::  M(3,3), x(5)
    
        M = reshape( (/ x(1), x(5), x(4), x(5), x(2), x(3), x(4), x(3), -x(1)-x(2) /), (/3,3/) )
    
    end function voigt2m

    
    function vonMises(A, m) result(VM)
    ! von Mises norm of stress or strain
        real(kind=8)   ::  A(3,3), VM, fVM, trace
        character(6)   ::  m
        
        trace = A(1,1) + A(2,2) + A(3,3)
        A = A - trace*eye(3)/3.d0
        
        if (m == 'strain') then
            fVM = sqrt(2.d0/3.d0)
        elseif (m == 'stress') then
            fVM = sqrt(3.d0/2.d0)
        end if
        
        VM = fVM*sqrt(A(1,1)**2 + A(2,2)**2 + A(3,3)**2 + 2.d0*(A(1,2)**2 + A(1,3)**2 + A(2,3)**2))
        
    end function vonMises
    
    
    subroutine vonMises_dvec(A, m, MA, VM, An) 
    ! von Mises norm of deviatoric stress or dev strain given in Voigt notation
    ! masked array MA
        real(kind=8)    ::  A(5), An(5), VM, fVM
        integer         ::  i, MA(5)
        character(6)    ::  m
        
        intent(in)      ::  A, m, MA
        intent(out)     ::  VM, An
        
        if (m == 'strain') then
            fVM = sqrt(2.d0/3.d0)
        elseif (m == 'stress') then
            fVM = sqrt(3.d0/2.d0)
        end if
        
        An = 0.d0
        forall (i=1:5, MA(i)==1) An(i) = A(i)
			
        VM = fVM*sqrt(An(1)**2 + An(2)**2 + (-An(1)-An(2))**2 + 2.d0*An(3)**2 + 2.d0*An(4)**2 + 2.d0*An(5)**2)
        
        An = A
        forall (i=1:5, MA(i)==1) An(i) = An(i)/VM
        
    end subroutine vonMises_dvec
    

!C**********************************************************************
!C                         SUBROUTINE findslips                        *
!C**********************************************************************
!C By means of simplex algorithm and linear programming finds slip     *
!C slip rates which minimizes the internal energy and fullfill the     *
!C the geometrical constrains given by strain rate tensor              *
!C**********************************************************************
    subroutine findslips(CRSS, D, sliprates, activesID, Nactive)
        
        ! input
        real(kind=8)    ::  D(3,3), CRSS(:)
        ! output
        integer         ::  Nactive, activesID(size(CRSS))
        real(kind=8)    ::  sliprates(size(CRSS))
        
        integer         ::  I, STAT, IZROV(24), IPOSV(5)
        real(kind=8)    ::  Msimplex(7,25), tmp1, tmp2
        
        intent(in)      ::  CRSS, D
        intent(out)     ::  activesID, sliprates, Nactive

        
        !***  PREPARE FOR SIMPLX METHOD INPUT ***
        !CONSTANT COLUMM
        Msimplex(1,1)= 0.d0
        Msimplex(2,1)= D(1,1)
        Msimplex(3,1)= D(2,2)
        Msimplex(4,1)= D(2,3)
        Msimplex(5,1)= D(1,3)
        Msimplex(6,1)= D(1,2)      
        !SECOND TO SIX ROW, dE=sum(dr)    !CURRENT COLUM 
        do I=1,24
            !2-13 columns
            Msimplex(1,I+1)= -CRSS(I)
            Msimplex(2,I+1)= -P(1,1,I)
            Msimplex(3,I+1)= -P(2,2,I)
            Msimplex(4,I+1)= -P(2,3,I)
            Msimplex(5,I+1)= -P(1,3,I)
            Msimplex(6,I+1)= -P(1,2,I)
        end do
        !IF D(I,J)<0, MULTIPLY BY -1.0
        do I=2,6
            if (Msimplex(I,1) < 0.D0) then
                Msimplex(I,:) = -Msimplex(I,:)
            end if
        end do
        !CALL THE SIMPLEX METHOD SUBROUTINE
        call simplx(Msimplex,5,24,7,25,0,0,5,STAT,IZROV,IPOSV)

        if (STAT /= 0) then
            WRITE(*,*) 'ERROR: No solution found by simplex!' 
            read(*,*)
            stop
        end if

        !***  RETRIVE THE OPTIMUM SLIP SYSTEM ID AND ITS SHEAR VALUE  ***  
        Nactive = 0
        sliprates = 0.d0
        do I=1,5
            tmp1 = Msimplex(I+1,1)
            tmp2 = IPOSV(I)
            if (tmp1 > tol) then
                Nactive = Nactive + 1
                activesID(Nactive) = tmp2
                sliprates(tmp2) = tmp1
            end if          
        end do

    end subroutine findslips
    

    subroutine findslips_alamel(Pg1, Pg2, CRSS_g1, CRSS_g2, relax_penalties, D, slipratesG1, slipratesG2,           &
                                 activesIDG1, activesIDG2, aIDg1, aIDg2, relaxation_slips, active_relax)
        
        ! input
        real(kind=8)    ::  CRSS_G1(:), CRSS_G2(:), D(3,3), Pg1(:,:,:), Pg2(:,:,:), relax_penalties(3)
        ! output
        integer         ::  aIDg1, aIDg2, active_relax(3), activesIDG1(size(CRSS_G1)), activesIDG2(size(CRSS_G1))
        real(kind=8)    ::  slipratesG1(size(CRSS_G1)), slipratesG2(size(CRSS_G1)), relaxation_slips(3)
        !
        integer         ::  STAT, IZROV(52), IPOSV(10), CCOL, M, RC13NUM, RC23NUM, I, J, CID, ACTNUMG1, ACTNUMG2
        real(kind=8)    ::  FF, X(52), CRSS_R13, CRSS_R23, RC13, RC23, SIMPMAT(12, 53)
        !
        intent(in)      ::  CRSS_G1, CRSS_G2, D, Pg1, Pg2, relax_penalties
        intent(out)     ::  slipratesG1, slipratesG2,  activesIDG1, activesIDG2, aIDg1, aIDg2,                      &
                            relaxation_slips, active_relax
        
   
        ! initialization
        SIMPMAT  = 0.D0
        ACTNUMG1 = 0
        ACTNUMG2 = 0
        RC13NUM  = 0
        RC23NUM  = 0
        activesIDG1 = 0
        activesIDG2 = 0
        slipratesG1   = 0.D0
        slipratesG2   = 0.D0
        CRSS_R23 = relax_penalties(1)
        CRSS_R13 = relax_penalties(2)
        relaxation_slips = 0.d0
        active_relax = 0
            
        FF = 0.0
        CCOL = 0
        RC13 = 0.d0
        RC23 = 0.d0
        
        !***  PREPARE FOR SIMPLX METHOD INPUT ***
        ! --- first column
        !CONSTANT COLUMM
        SIMPMAT(1,1)= 0.0   
        !Grain-1           
        SIMPMAT(2,1) = D(1,1)
        SIMPMAT(3,1) = D(2,2)
        SIMPMAT(4,1) = D(2,3)
        SIMPMAT(5,1) = D(1,3)
        SIMPMAT(6,1) = D(1,2)
        !Grain-2
        SIMPMAT(7,1) = D(1,1)
        SIMPMAT(8,1) = D(2,2)
        SIMPMAT(9,1) = D(2,3)
        SIMPMAT(10,1)= D(1,3)
        SIMPMAT(11,1)= D(1,2)

        !SECOND TO SIX ROW, dE=sum(dr)	!The first grain, from colum 2-25
        DO I=1,12		
            CCOL=I+1
            !First 12 colums
            SIMPMAT(1,CCOL) = -CRSS_G1(I)
            SIMPMAT(2,CCOL) = -Pg1(1,1,I)
            SIMPMAT(3,CCOL) = -Pg1(2,2,I)
            SIMPMAT(4,CCOL) = -Pg1(2,3,I)
            SIMPMAT(5,CCOL) = -Pg1(1,3,I)
            SIMPMAT(6,CCOL) = -Pg1(1,2,I)
            !From 12 to 24 colums
            CCOL=I+13
            SIMPMAT(1,CCOL) = -CRSS_G1(I+12)
            SIMPMAT(2,CCOL) = Pg1(1,1,I)
            SIMPMAT(3,CCOL) = Pg1(2,2,I)
            SIMPMAT(4,CCOL) = Pg1(2,3,I)
            SIMPMAT(5,CCOL) = Pg1(1,3,I)
            SIMPMAT(6,CCOL) = Pg1(1,2,I)		
        ENDDO
      
        !The first grain, from colum 26-49
        do I=1, 12		           
            CCOL=I+25
            !First 12 colums
            SIMPMAT(1,CCOL) = -CRSS_G2(I)
            SIMPMAT(7,CCOL) = -Pg2(1,1,I)
            SIMPMAT(8,CCOL) = -Pg2(2,2,I)
            SIMPMAT(9,CCOL) = -Pg2(2,3,I)
            SIMPMAT(10,CCOL)= -Pg2(1,3,I)
            SIMPMAT(11,CCOL)= -Pg2(1,2,I)
            !From 12 to 24 colums
            CCOL=I+37
            SIMPMAT(1,CCOL) = -CRSS_G2(I+12)
            SIMPMAT(7,CCOL) = Pg2(1,1,I)
            SIMPMAT(8,CCOL) = Pg2(2,2,I)
            SIMPMAT(9,CCOL) = Pg2(2,3,I)
            SIMPMAT(10,CCOL)= Pg2(1,3,I)
            SIMPMAT(11,CCOL)= Pg2(1,2,I)		
        ENDDO

        !For colum 50-51 * This is very important lines, since quite easy to make error.
        !What I use is, Grain1: D = d + RC ; Grain2: D = d - RC
        !RC = (R+RT)/2
        !R13 = 0 0 1      R23 = 0 0 0
        !      0 0 0            0 0 1
        !      0 0 0            0 0 0
        !The last two colums are for Relaxation,
        SIMPMAT(4,50) = -0.5D0     !RC23 for Grain1
        SIMPMAT(5,51) = -0.5D0     !RC13 for Grain1
        SIMPMAT(4,52) =  0.5D0     !RC23 for Grain1, negative value
        SIMPMAT(5,53) =  0.5D0     !RC13 for Grain1, negative value

        SIMPMAT(9,50) =  0.5D0     !RC23 for Grain2
        SIMPMAT(10,51)=  0.5D0     !RC13 for Grain2
        SIMPMAT(9,52) = -0.5D0     !RC23 for Grain2, negative value
        SIMPMAT(10,53)= -0.5D0	   !RC13 for Grain2, negative value

        SIMPMAT(1,50) = -CRSS_R23    
        SIMPMAT(1,51) = -CRSS_R13    
        SIMPMAT(1,52) = -CRSS_R23    
        SIMPMAT(1,53) = -CRSS_R13    

        !If the first colume are negative, change the sign of this line
        DO I = 2, 11
            IF (SIMPMAT(I,1) < 0.0D0) THEN
	            DO J = 1, 53
	                FF = SIMPMAT(I,J)
	                SIMPMAT(I,J) = -FF
	            ENDDO
            ENDIF
        ENDDO

        M = 10
        CALL simplx(SIMPMAT,M,52,12,53,0,0,10,STAT,IZROV,IPOSV)
        if(STAT/=0) THEN
            WRITE(*,*) 'ERROR: No solution found from simplex!' 
            read(*,*)
            STOP
        end if

        !----------------------------------------------------------------------------
        !Retreive information
        aIDg1 = 0
        aIDg2 = 0
        DO I = 1, M
            CID    = IPOSV(I)
            X(CID) = SIMPMAT(I+1,1)
            
	        IF ((CID < 25) .AND. (CID > 0)) THEN
                !for grain 1
		        IF ((ABS(SIMPMAT(I+1,1)) > tol) .AND. (ABS(SIMPMAT(I+1,1)) < tolinv)) THEN
		            slipratesG1(CID) = SIMPMAT(I+1,1)
                    aIDg1 = aIDg1 + 1
		            activesIDG1(aIDg1) = CID
                ENDIF
	        ELSE IF((CID < 49) .AND. (CID > 24)) THEN
                !for grain 2
		        CID = CID - 24
		        IF ((ABS(SIMPMAT(I+1,1)) > tol) .AND. (ABS(SIMPMAT(I+1,1)) < tolinv)) THEN
		            slipratesG2(CID) = SIMPMAT(I+1,1)
		            aIDg2 = aIDg2 + 1
		            activesIDG2(aIDg2) = CID
		        ENDIF

	       ELSE IF (CID == 49) THEN
	            !for RC23
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC23 = SIMPMAT(I+1, 1)
		            RC23NUM = RC23NUM + 1
		        ENDIF

	       ELSE IF (CID == 50) THEN
	            !for RC13
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC13 = SIMPMAT(I+1, 1)
		            RC13NUM = RC13NUM + 1
		        ENDIF
	       ELSE IF (CID == 51) THEN
	            !for RC23
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC23 = -SIMPMAT(I+1,1)
		            RC23NUM = RC23NUM + 1
		        ENDIF
	       ELSE IF (CID == 52) THEN
	            !for RC13
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
                    RC13 = -SIMPMAT(I+1,1)
                    RC13NUM = RC13NUM + 1
		        ENDIF

	       ELSE
	            WRITE(*,*) 'ERROR: Alamel solution error!' 
				write(*,*) 'Sliprates1 g1'
				write(*,*) slipratesG1
				write(*,*) 'Sliprates1 g2'
				write(*,*) slipratesG2
				write(*,*) 'Dp'
				write(*,*) D
	            read(*,*)
                STOP	   
	       ENDIF
	   	   	  
        ENDDO 

        relaxation_slips = (/ RC13, RC23, 0.d0 /) 
        active_relax = (/ RC13NUM, RC23NUM, 0 /)

		!if ((aIDg1+aIDg2) /= 8) then
		!	write(*,*) 'Number of active slip systems in cluster not equal 8!'
		!	write(*,*) 'Grain 1'
		!	write(*,*) slipratesG1(activesIDG1(1:aIDg1))
		!	write(*,*) 'Grain 2'
		!	write(*,*) slipratesG2(activesIDG2(1:aIDg2))
  !      end if
    
    end subroutine findslips_alamel

          

    subroutine findslips_alamel3(Pg1, Pg2, CRSS_g1, CRSS_g2, relax_penalties, D, slipratesG1, slipratesG2,          &
                                 activesIDG1, activesIDG2, aIDg1, aIDg2, relaxation_slips, active_relax)
        
        ! input
        real(kind=8)    ::  Pg1(:,:,:), Pg2(:,:,:), D(3,3), CRSS_G1(:), CRSS_G2(:), relax_penalties(3)
        ! output
        integer         ::  aIDg1, aIDg2, active_relax(3), activesIDG1(size(CRSS_g1)), activesIDG2(size(CRSS_g1))
        real(kind=8)    ::  slipratesG1(size(CRSS_g1)), slipratesG2(size(CRSS_g1)), relaxation_slips(3)
        !
        integer         ::  STAT, IZROV(54), IPOSV(10), CCOL, M, RC13NUM, RC23NUM, RC12NUM, I, J, CID,              &
                            ACTNUMG1, ACTNUMG2
        real(kind=8)    ::  FF, X(54), CRSS_R13, CRSS_R23, CRSS_R12, RC13, RC23, RC12, SIMPMAT(12, 55)

        intent(in)      ::  CRSS_G1, CRSS_G2, D, Pg1, Pg2, relax_penalties
        intent(out)     ::  slipratesG1, slipratesG2,  activesIDG1, activesIDG2, aIDg1, aIDg2,                      &
                            relaxation_slips, active_relax
        

        ! initialization
        SIMPMAT  = 0.D0
        ACTNUMG1 = 0
        ACTNUMG2 = 0
        RC13NUM  = 0
        RC23NUM  = 0
		RC12NUM  = 0
        slipratesG1   = 0.D0
        slipratesG2   = 0.D0
        activesIDG1   = 0
        activesIDG2   = 0
        CRSS_R23 = relax_penalties(1)
        CRSS_R13 = relax_penalties(2)
		CRSS_R12 = relax_penalties(3)
        relaxation_slips = 0.d0
        active_relax = 0
            
        FF = 0.0
        CCOL = 0
        RC13 = 0.d0
        RC23 = 0.d0
		RC12 = 0.d0
        
        !***  PREPARE FOR SIMPLX METHOD INPUT ***
        ! --- first column
        !CONSTANT COLUMM
        SIMPMAT(1,1)= 0.0   
        !Grain-1           
        SIMPMAT(2,1) = D(1,1)
        SIMPMAT(3,1) = D(2,2)
        SIMPMAT(4,1) = D(2,3)
        SIMPMAT(5,1) = D(1,3)
        SIMPMAT(6,1) = D(1,2)
        !Grain-2
        SIMPMAT(7,1) = D(1,1)
        SIMPMAT(8,1) = D(2,2)
        SIMPMAT(9,1) = D(2,3)
        SIMPMAT(10,1)= D(1,3)
        SIMPMAT(11,1)= D(1,2)

        !SECOND TO SIX ROW, dE=sum(dr)	!The first grain, from colum 2-25
        DO I=1,12		
            CCOL=I+1
            !First 12 colums
            SIMPMAT(1,CCOL) = -CRSS_G1(I)
            SIMPMAT(2,CCOL) = -Pg1(1,1,I)
            SIMPMAT(3,CCOL) = -Pg1(2,2,I)
            SIMPMAT(4,CCOL) = -Pg1(2,3,I)
            SIMPMAT(5,CCOL) = -Pg1(1,3,I)
            SIMPMAT(6,CCOL) = -Pg1(1,2,I)
            !From 12 to 24 colums
            CCOL=I+13
            SIMPMAT(1,CCOL) = -CRSS_G1(I+12)
            SIMPMAT(2,CCOL) = Pg1(1,1,I)
            SIMPMAT(3,CCOL) = Pg1(2,2,I)
            SIMPMAT(4,CCOL) = Pg1(2,3,I)
            SIMPMAT(5,CCOL) = Pg1(1,3,I)
            SIMPMAT(6,CCOL) = Pg1(1,2,I)		
        ENDDO
      
        !The first grain, from colum 26-49
        do I=1, 12		           
            CCOL=I+25
            !First 12 colums
            SIMPMAT(1,CCOL) = -CRSS_G2(I)
            SIMPMAT(7,CCOL) = -Pg2(1,1,I)
            SIMPMAT(8,CCOL) = -Pg2(2,2,I)
            SIMPMAT(9,CCOL) = -Pg2(2,3,I)
            SIMPMAT(10,CCOL)= -Pg2(1,3,I)
            SIMPMAT(11,CCOL)= -Pg2(1,2,I)
            !From 12 to 24 colums
            CCOL=I+37
            SIMPMAT(1,CCOL) = -CRSS_G2(I+12)
            SIMPMAT(7,CCOL) = Pg2(1,1,I)
            SIMPMAT(8,CCOL) = Pg2(2,2,I)
            SIMPMAT(9,CCOL) = Pg2(2,3,I)
            SIMPMAT(10,CCOL)= Pg2(1,3,I)
            SIMPMAT(11,CCOL)= Pg2(1,2,I)		
        ENDDO

        !For colum 50-51 * This is very important lines, since quite easy to make error.
        !What I use is, Grain1: D = d + RC ; Grain2: D = d - RC
        !RC = (R+RT)/2
        !R13 = 0 0 1      R23 = 0 0 0		R12 = 0 1 0
        !      0 0 0            0 0 1			  0 0 0
        !      0 0 0            0 0 0			  0 0 0
        !The last two colums are for Relaxation,
        SIMPMAT(4,50) = -0.5D0     !RC23 for Grain1
        SIMPMAT(5,51) = -0.5D0     !RC13 for Grain1
		SIMPMAT(6,52) = -0.5D0     !RC12 for Grain1
        SIMPMAT(4,53) =  0.5D0     !RC23 for Grain1, negative value
        SIMPMAT(5,54) =  0.5D0     !RC13 for Grain1, negative value
		SIMPMAT(6,55) =  0.5D0     !RC12 for Grain1, negative value

        SIMPMAT(9,50) =  0.5D0     !RC23 for Grain2
        SIMPMAT(10,51)=  0.5D0     !RC13 for Grain2
		SIMPMAT(11,52)= -0.5D0     !RC12 for Grain2
        SIMPMAT(9,53) = -0.5D0     !RC23 for Grain2, negative value
        SIMPMAT(10,54)= -0.5D0	   !RC13 for Grain2, negative value
		SIMPMAT(11,55)=  0.5D0     !RC12 for Grain2, negative value

        SIMPMAT(1,50) = -CRSS_R23    
        SIMPMAT(1,51) = -CRSS_R13    
		SIMPMAT(1,52) = -CRSS_R12  
        SIMPMAT(1,53) = -CRSS_R23    
        SIMPMAT(1,54) = -CRSS_R13
		SIMPMAT(1,55) = -CRSS_R12  		

        !If the first colume are negative, change the sign of this line
        DO I = 2, 11
            IF (SIMPMAT(I,1) < 0.0D0) THEN
	            DO J = 1, 53
	                FF = SIMPMAT(I,J)
	                SIMPMAT(I,J) = -FF
	            ENDDO
            ENDIF
        ENDDO

        M = 10
        CALL simplx(SIMPMAT,M,54,12,55,0,0,10,STAT,IZROV,IPOSV)
        if(STAT/=0) THEN
            WRITE(*,*) 'ERROR: No solution found from simplex!' 
            read(*,*)
            STOP
        end if

        !----------------------------------------------------------------------------
        !Retrive information
        aIDg1 = 0
        aIDg2 = 0
        DO I = 1, M
            CID    = IPOSV(I)
            X(CID) = SIMPMAT(I+1,1)
            
	        IF ((CID < 25) .AND. (CID > 0)) THEN
                !for grain 1
		        IF ((ABS(SIMPMAT(I+1,1)) > tol) .AND. (ABS(SIMPMAT(I+1,1)) < tolinv)) THEN
		            slipratesG1(CID) = SIMPMAT(I+1,1)
                    aIDg1 = aIDg1 + 1
		            activesIDG1(aIDg1) = CID
                ENDIF
	        ELSE IF((CID < 49) .AND. (CID > 24)) THEN
                !for grain 2
		        CID = CID - 24
		        IF ((ABS(SIMPMAT(I+1,1)) > tol) .AND. (ABS(SIMPMAT(I+1,1)) < tolinv)) THEN
		            slipratesG2(CID) = SIMPMAT(I+1,1)
		            aIDg2 = aIDg2 + 1
		            activesIDG2(aIDg2) = CID
		        ENDIF

	       ELSE IF (CID == 49) THEN
	            !for RC23
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC23 = SIMPMAT(I+1, 1)
		            RC23NUM = RC23NUM + 1
		        ENDIF

	       ELSE IF (CID == 50) THEN
	            !for RC13
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC13 = SIMPMAT(I+1, 1)
		            RC13NUM = RC13NUM + 1
		        ENDIF
		   ELSE IF (CID == 51) THEN
	            !for RC12
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC12 = SIMPMAT(I+1,1)
		            RC12NUM = RC12NUM + 1
		        ENDIF
	       ELSE IF (CID == 52) THEN
	            !for RC23
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
		            RC23 = -SIMPMAT(I+1,1)
		            RC23NUM = RC23NUM + 1
		        ENDIF
	       ELSE IF (CID == 53) THEN
	            !for RC13
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
                    RC13 = -SIMPMAT(I+1,1)
                    RC13NUM = RC13NUM + 1
		        ENDIF
		   ELSE IF (CID == 54) THEN
	            !for RC12
		        IF ((SIMPMAT(I+1,1) > tol) .AND. (SIMPMAT(I+1,1) < tolinv)) THEN
                    RC12 = -SIMPMAT(I+1,1)
                    RC12NUM = RC12NUM + 1
		        ENDIF

	       ELSE
                WRITE(*,*) 'ERROR: Alamel3 solution error!' 
				write(*,*) 'Sliprates1 g1'
				write(*,*) slipratesG1
				write(*,*) 'Sliprates1 g2'
				write(*,*) slipratesG2
				write(*,*) 'Dp'
				write(*,*) D
	            read(*,*)
                STOP	   
	       ENDIF
	   	   	  
        ENDDO 

        relaxation_slips = (/ RC13, RC23, RC12 /) 
        active_relax = (/ RC13NUM, RC23NUM, RC12NUM /)
    

!		if ((aIDg1+aIDg2+RC23NUM+RC13NUM+RC12NUM) /= 10) then
!			write(*,1000) aIDg1+aIDg2
!1000        format ('Number of active slip systems in cluster: ', i3)
!            write(*,1001) RC23NUM+RC13NUM+RC12NUM
!1001        format ('Number of active relaxations: ', i3)
!			write(*,*) 'Grain 1'
!			write(*,*) slipratesG1(activesIDG1(1:aIDg1))
!			write(*,*) 'Grain 2'
!			write(*,*) slipratesG2(activesIDG2(1:aIDg2))
!        end if
    
    end subroutine findslips_alamel3
                                 
         
                                 
                                 
    subroutine getRSS(Ns, S, rss)
    
        integer         ::  i, Ns
        real(kind=8)    ::  S(3,3), rss(Ns), x(5), B(Ns, 5)

        intent(in)      ::  S, Ns
        intent(out)     ::  rss
 
        x = m2voigt(S)
        
        B = 0.d0
        do i=1, Ns
            B(i,:) = (/ 2.d0*P(1,1,i)+P(2,2,i), 2.d0*P(2,2,i)+P(1,1,i), 2.d0*P(2,3,i), 2.d0*P(1,3,i), 2.d0*P(1,2,i) /)
        end do
        
        rss = mdot(B, x)
    
    end subroutine getRSS
                                 
                                 
                           
                                 
    subroutine solveambSVD(Pi, D, stress, crss, gm, found_SVD_solution)
    
        integer         ::  i, j, Nact, pactID(24), pactID2(24), pID, statu, INFO, negative_ind
        real(kind=8)    ::  Pi(3,3,24), crss(24), stress(3,3), D(3,3), rss(24), WORK(100), gm(24), &
                            tmpD(3,3), tmpdvec(5), maxdiff, Dvec(5)
        logical         ::  found_SVD_solution 
        
        real(kind=8), allocatable   :: A(:,:), S(:), SM(:,:), U(:,:), VT(:,:), Apinv(:,:), gmPA(:)
        
        intent(in)      ::  Pi, D, stress, crss
        intent(out)     ::  gm, found_SVD_solution
        
        Dvec = m2voigt(D)
        call getRSS(Nslips, stress, rss)
        
        pactID = 0
        j = 0
        do i=1, 24
            if (rss(i)/crss(i) > (1.d0-tol)) then
                 j = j + 1
                pactID(j) = i
            end if
        end do
        Nact = j
        
        
        A = 0.d0
        found_SVD_solution = .false.
        do while (.not.(found_SVD_solution))
            allocate( A(Nact,5), Apinv(Nact,5), S(min(Nact,5)), SM(Nact,5), U(Nact,Nact), VT(5,5), gmPA(Nact), STAT=statu)
            if (statu /= 0) then
                write(*,*) 'Allocation in SVD unsuccesfull.'
                read(*,*)
                stop
            end if
            
            do i=1, Nact
                pID = pactID(i)
                A(i,:) = (/ Pi(1,1,pID), Pi(2,2,pID), Pi(2,3,pID), Pi(1,3,pID), Pi(1,2,pID) /)
            end do
        
            call dgesvd('A','A',Nact,5,A,Nact,S,U,Nact,VT,5,WORK,100,INFO)
            if (INFO /= 0) then
                write(*,*) 'SVD failed to converge.'
                read(*,*)
                stop
            end if
            
            SM = 0.d0
            do i=1, min(5,Nact)
                if (S(i) > 1.d-20) then
                    SM(i,i) = 1.d0/S(i)
                else
                    SM(i,i) = S(i)
                end if
            end do
            
            Apinv = matmul(U, matmul(SM, VT))
            
            gmPA = mdot(Apinv, Dvec)

            gm = 0.d0
            
            do i=1, Nact
                pID = pactID(i)
                gm(pID) = gmPA(i)
            end do
                         
            ! to be further testet whether it is a valid solution
            ! check if gm is a valid solution
            tmpD = 0.d0
            do i=1, Nact
                pID = pactID(i)
                tmpD = tmpD + Pi(:,:,pID)*gm(pID)
            end do
            
            tmpdvec = m2voigt(tmpD)
            maxdiff = maxval(abs(Dvec - tmpdvec))
            if (maxdiff > 1e-3) then
                write(*,*) 'S matrix in SVD singular'
                exit
            end if
                
            ! if Nact falls below 5, the subroutine exit with found_SVD_solution = .false. this is rarely
            if (any(gm < -tol)) then
                negative_ind = minloc(gm, DIM=1)
                ! remove the index of the most negative gm from the set of indeces
                ! of potentially active slip systems
                j = 1
                pactID2 = 0
                do i=1, Nact
                    if (pactID(i) /= negative_ind) then
                        pactID2(j) = pactID(i)
                        j = j + 1
                    end if
                end do
                Nact = Nact - 1
                pactID = pactID2
            else
                found_SVD_solution = .true.
            end if
            
            deallocate(A, Apinv, S, SM, U, VT, gmPA, STAT=statu)
            if (statu /= 0) then
                write(*,*) 'Deallocation in SVD unsuccesfull.'
                read(*,*)
                stop
            end if
            
            if (Nact < 5) exit
            
        end do 
        
        
    end subroutine solveambSVD
                                 
            

    subroutine getstress_taylor(sliprates, activesID, Nact, CRSS, S)
        
        ! input
        integer         ::  activesID(:), Nact
        real(kind=8)    ::  CRSS(:), sliprates(:)
        ! output
        real(kind=8)    ::  S(3,3)
        !
        integer         ::  i, j, ipiv(5), info, FID, Nscor, ind(5), ALIVE, pID
        real(kind=8)    ::  B(5,5), x(5), S_corners(3,3,12), rss(size(CRSS)), tmp1(5), tmp2(size(CRSS))
        logical         ::  flag, already_in
    
        intent(in)      ::  sliprates, activesID, Nact, CRSS
        intent(out)     ::  S

        if (Nact == 5) then
            B = 0.d0
            do i=1, Nact
                pID = activesID(i)
                B(i,:) = (/ 2.d0*P(1,1,pID)+P(2,2,pID), 2.d0*P(2,2,pID)+P(1,1,pID), 2.d0*P(2,3,pID), 2.d0*P(1,3,pID), 2.d0*P(1,2,pID) /)
            end do
            x = (/ (CRSS(activesID(j)), j=1,5) /)
            call dgesv(5, 1, B, 5, ipiv, x, 5, info)      
            S = voigt2m(x)
        elseif (Nact < 5) then
            if (not(all(sslookup(1,:) == (/ 1, 2, 4, 5, 7 /) ))) then
                FID = 1112
                open(FID, FILE='sslookup.dat', STATUS='OLD', IOSTAT=ALIVE)          
                do i=1, 12288
                    read(FID,*) (sslookup(i,j), j=1,5)
                end do
                close(FID)
            end if
        
            S_corners = 0.d0
            Nscor = 0
            do i=1, 12288
                flag = .true.
                ind = sslookup(i,:)
                do j=1, Nact
                    if (.not.(any(activesID(j) == ind))) then
                        flag = .false.
                        exit
                    end if
                end do
                if (flag == .true.) then
                    B = 0.d0
                    do j=1, 5
                        pID = ind(j)
                        B(j,:) = (/ 2.d0*P(1,1,pID)+P(2,2,pID), 2.d0*P(2,2,pID)+P(1,1,pID), 2.d0*P(2,3,pID), 2.d0*P(1,3,pID), 2.d0*P(1,2,pID) /)
                    end do
                    x = (/ (CRSS(ind(j)), j=1,5) /)
                    call dgesv(5, 1, B, 5, ipiv, x, 5, info)
                    S = voigt2m(x)
                    call getRSS(Nslips, S, rss)
                    ! check whether tau and gm have the same signs for all ss
                    ! and check whether the yield criterion is obeyed
                    tmp1 = (/ (rss(ind(j)), j=1,5) /)
                    tmp2 = (/ (abs(rss(j))/CRSS(j), j=1,Nslips) /)
                    if (all(tmp1 >= 0.d0) .and. (all(tmp2 <= (1.d0+tol)))) then
                        already_in = .false.
                        if (Nscor == 0) then
                            S_corners(:,:,1) = S
                            Nscor = 1
                        else
                            do j=1, Nscor
                                if (maxval(abs(S_corners(:,:,j) - S)) < tol) then
                                    already_in = .true.
                                    exit
                                end if
                            end do
                            if (already_in == .false.) then
                                Nscor = Nscor + 1
                                S_corners(:,:,Nscor) = S
                            end if
                        end if
                    end if
                end if
            end do
            S = 0.d0
            do i=1, Nscor
                S = S + S_corners(:,:,i)
            end do
            S = S/Nscor 
        else
            write(*,*) 'Nact for a basis solution larger than 5'
            read(*,*)
            stop
        end if
    
    end subroutine getstress_taylor
    

    subroutine getstress_alamel(Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,               &
								NactG1, NactG2, QbQg1T, QbQg2T, Qb, active_relax, SlocG1, SlocG2)
    
        integer         ::  activesIDG1(:), activesIDG2(:), NactG1, NactG2, i, ipiv(10), info,      &
                            active_relax(3)
        real(kind=8)    ::  Pg1(:,:,:), Pg2(:,:,:), CRSS_g1(:), CRSS_g2(:), B(10,10),               &
                            crss_cluster(10), Binv(10,10), SlocG1(3,3), SlocG2(3,3),                &
                            Sg1(3,3), Sg2(3,3), QbQg1T(3,3), QbQg2T(3,3), Qb(3,3), x(10)

        intent(in)      ::  Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,                   &
                            NactG1, NactG2, QbQg1T, QbQg2T, Qb, active_relax
        intent(out)     ::  SlocG1, SlocG2
    
    
        B = 0.d0
        crss_cluster = 0.d0
        
        ! grain1
        do i = 1, NactG1
            B(i, 1:5) = (/                                                                          &
                         2.d0*Pg1(1,1,activesIDG1(i)) + Pg1(2,2,activesIDG1(i)),                    &
                         2.d0*Pg1(2,2,activesIDG1(i)) + Pg1(1,1,activesIDG1(i)),                    &
                         2.d0*Pg1(2,3,activesIDG1(i)),                                              &
                         2.d0*Pg1(1,3,activesIDG1(i)),                                              &
                         2.d0*Pg1(1,2,activesIDG1(i))                                               &
                       /) 
        enddo
        ! grain2
        
        do i = 1, NactG2
            B(NactG1+i, 6:10) = (/  2.d0*Pg2(1,1,activesIDG2(i)) + Pg2(2,2,activesIDG2(i)),         &
                                    2.d0*Pg2(2,2,activesIDG2(i)) + Pg2(1,1,activesIDG2(i)),         &
                                    2.d0*Pg2(2,3,activesIDG2(i)),                                   &
                                    2.d0*Pg2(1,3,activesIDG2(i)),                                   &
                                    2.d0*Pg2(1,2,activesIDG2(i)) /) 
        end do
        ! relaxations
        if (active_relax(1) > 0) then
            B(NactG1+NactG2+1,:) =                                                                  &
                               (/   symR1(1,1) - symR1(3,3),                                        &
                                    symR1(2,2) - symR1(3,3),                                        & 
                                    2.d0*symR1(2,3),                                                &
                                    2.d0*symR1(1,3),                                                &
                                    2.d0*symR1(1,2),                                                &
                                    -symR1(1,1) + symR1(3,3),                                       &
                                    -symR1(2,2)+symR1(3,3),                                         &
                                    -2.d0*symR1(2,3),                                               &
                                    -2.d0*symR1(1,3),                                               &
                                    -2.d0*symR1(1,2) /)
        end if
        if (active_relax(2) > 0) then                                                               
            B(NactG1+NactG2+active_relax(1)+1,:) =                                                  &
                               (/   symR2(1,1) - symR2(3,3),                                        &
                                    symR2(2,2) - symR2(3,3),                                        &
                                    2.d0*symR2(2,3),                                                &
                                    2.d0*symR2(1,3),                                                &
                                    2.d0*symR2(1,2),                                                &
                                    -symR2(1,1) + symR2(3,3),                                       &
                                    -symR2(2,2)+symR2(3,3),                                         &
                                    -2.d0*symR2(2,3),                                               &
                                    -2.d0*symR2(1,3),                                               &
                                    -2.d0*symR2(1,2) /)     
        end if

        CRSS_cluster(:NactG1)   = CRSS_g1(activesIDG1)
        CRSS_cluster((NactG1+1):(NactG1+NactG2)) = CRSS_g2(activesIDG2)
        
        call dgesv(10, 1, B, 10, ipiv, CRSS_cluster, 10, info)
        x = CRSS_cluster

        ! express stress into the each grain's coord sys and global sys
        Sg1 = reshape( (/ x(1), x(5), x(4), x(5), x(2), x(3), x(4), x(3), -x(1)-x(2) /), (/3,3/) )     
        Sg2 = reshape( (/ x(6), x(10), x(9), x(10), x(7), x(8), x(9), x(8), -x(6)-x(7) /), (/3,3/) )             
        SlocG1  = transform(Sg1, transp3(QbQg1T))
        SlocG2  = transform(Sg2, transp3(QbQg2T)) 
    
        ! this is only a check of stress components relations due to relaxations, can be DELETED
        if ((active_relax(1) > 0) .AND. (abs(Sg1(1,3) - Sg2(1,3)) > tol)) then
            write(*,*) 'Stress problem relax 1!'
        end if
        if ((active_relax(2) > 0) .AND. (abs(Sg1(2,3) - Sg2(2,3)) > tol)) then
            write(*,*) 'Stress problem relax 2!'
        end if

    end subroutine getstress_alamel
    

    subroutine getstress_alamel3(Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,              & 
							     NactG1, NactG2, QbQg1T, QbQg2T, Qb, active_relax,                  &
                                 SlocG1, SlocG2)
    
        integer         ::  activesIDG1(:), activesIDG2(:), NactG1, NactG2, i, ipiv(10), info,      &
                            active_relax(3)
        real(kind=8)    ::  Pg1(:,:,:), Pg2(:,:,:), CRSS_g1(:), CRSS_g2(:), B(10,10),               &
                            crss_cluster(10), Binv(10,10), SlocG1(3,3), SlocG2(3,3),                &
                            Sg1(3,3), Sg2(3,3), QbQg1T(3,3), QbQg2T(3,3), Qb(3,3), x(10)

        intent(in)      ::  Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,                   &
                            NactG1, NactG2, QbQg1T, QbQg2T, Qb, active_relax
        intent(out)     ::  SlocG1, SlocG2
    
    
        B = 0.d0
        crss_cluster = 0.d0
        
        ! grain1
        do i = 1, NactG1
            B(i, 1:5) = (/                                                                          &
                         2.d0*Pg1(1,1,activesIDG1(i)) + Pg1(2,2,activesIDG1(i)),                    &
                         2.d0*Pg1(2,2,activesIDG1(i)) + Pg1(1,1,activesIDG1(i)),                    &
                         2.d0*Pg1(2,3,activesIDG1(i)),                                              &
                         2.d0*Pg1(1,3,activesIDG1(i)),                                              &
                         2.d0*Pg1(1,2,activesIDG1(i))                                               &
                       /) 
        enddo
        ! grain2
        
        do i = 1, NactG2
            B(NactG1+i, 6:10) = (/  2.d0*Pg2(1,1,activesIDG2(i)) + Pg2(2,2,activesIDG2(i)),         &
                                    2.d0*Pg2(2,2,activesIDG2(i)) + Pg2(1,1,activesIDG2(i)),         &
                                    2.d0*Pg2(2,3,activesIDG2(i)),                                   &
                                    2.d0*Pg2(1,3,activesIDG2(i)),                                   &
                                    2.d0*Pg2(1,2,activesIDG2(i)) /) 
        end do
        ! relaxations
        if (active_relax(1) > 0) then
            B(NactG1+NactG2+1,:) = (/   symR1(1,1) - symR1(3,3),                                    &
                                        symR1(2,2) - symR1(3,3),                                    & 
                                        2.d0*symR1(2,3),                                            &
                                        2.d0*symR1(1,3),                                            &
                                        2.d0*symR1(1,2),                                            &
                                        -symR1(1,1) + symR1(3,3),                                   &
                                        -symR1(2,2) + symR1(3,3),                                   &
                                        -2.d0*symR1(2,3),                                           &
                                        -2.d0*symR1(1,3),                                           &
                                        -2.d0*symR1(1,2) /)
        end if
        if (active_relax(2) > 0) then
            B(NactG1+NactG2+active_relax(1)+1,:) =                                                  &
                                   (/   symR2(1,1) - symR2(3,3),                                    &
                                        symR2(2,2) - symR2(3,3),                                    &
                                        2.d0*symR2(2,3),                                            &
                                        2.d0*symR2(1,3),                                            &
                                        2.d0*symR2(1,2),                                            &
                                        -symR2(1,1) + symR2(3,3),                                   &
                                        -symR2(2,2) + symR2(3,3),                                   &
                                        -2.d0*symR2(2,3),                                           &
                                        -2.d0*symR2(1,3),                                           &
                                        -2.d0*symR2(1,2) /)   
        end if
        if (active_relax(3) > 0) then
		    B(NactG1+NactG2+active_relax(1)+active_relax(2)+1,:) =                                  &
                                   (/   symR3(1,1) - symR3(3,3),                                    &
                                        symR3(2,2) - symR3(3,3),                                    &
                                        2.d0*symR3(2,3),                                            &
                                        2.d0*symR3(1,3),                                            &
                                        2.d0*symR3(1,2),                                            &
                                        symR3(1,1) - symR3(3,3),                                    &
                                        symR3(2,2) - symR3(3,3),                                    &
                                        2.d0*symR3(2,3),                                            &
                                        2.d0*symR3(1,3),                                            &
                                        2.d0*symR3(1,2) /) 
        end if

        CRSS_cluster(:NactG1)   = CRSS_g1(activesIDG1)
        CRSS_cluster((NactG1+1):(NactG1+NactG2)) = CRSS_g2(activesIDG2)
        
        call dgesv(10, 1, B, 10, ipiv, CRSS_cluster, 10, info)
        x = CRSS_cluster

        ! express stress into the each grain's coord sys and global sys
        Sg1 = reshape( (/ x(1), x(5), x(4), x(5), x(2), x(3), x(4), x(3), -x(1)-x(2) /), (/3,3/) )     
        Sg2 = reshape( (/ x(6), x(10), x(9), x(10), x(7), x(8), x(9), x(8), -x(6)-x(7) /), (/3,3/) )             
        SlocG1  = transform(Sg1, transp3(QbQg1T))
        SlocG2  = transform(Sg2, transp3(QbQg2T))                       
        
        ! this is only a check of stress components relations due to relaxations, can be DELETED
        if ((active_relax(1) > 0) .AND. (abs(Sg1(1,3) - Sg2(1,3)) > tol)) then
            write(*,*) 'Stress problem relax 1!'
        end if
        if ((active_relax(2) > 0) .AND. (abs(Sg1(2,3) - Sg2(2,3)) > tol)) then
            write(*,*) 'Stress problem relax 2!'
        end if
        if ((active_relax(3) > 0) .AND. (abs(Sg1(1,2) + Sg2(1,2)) > tol)) then
            write(*,*) 'Stress problem relax 3!'
        end if

    end subroutine getstress_alamel3                         

                         
    
                                 
                                 
!   ***************** TALYLOR MODEL MAIN ROUTINES ****************                         
    
    subroutine taylor(ndim, dguess, Dsolved, Ssolved, info, resid, nfev)
    
        integer         ::  j, INC, statu, s(2), tmp(2), ndim, info, nfev
        real(kind=8)    ::  Sdir(3,3), Lp(3,3), dguess(ndim), Dsolved(3,3), Ssolved(3,3), &
							resid(5)

        intent(inout)   ::  ndim, dguess
        intent(out)     ::  Dsolved, Ssolved, info, resid, nfev
		
		
        allocate(Q0g1_list, SOURCE=Qg1_list, STAT=s(1))
        allocate(Rg1(3,3,Ngrains), STAT=s(2))
        
        if (any(s /= 0)) then
            write(*,*) 'Allocation of initial orientation matrices was unsuccesfull.'
            read(*,*)
            stop
        end if
		
        ! initialize rotation matrix to unity
        do j = 1, Ngrains
            Rg1(:,:,j) = eye(3)
        end do
            
        totalslip = 0.d0
        outINC    = 0
        Fglob     = eye(3)
		
		call alloc_outputvars_taylor()
		
        if (sum(dguess) .LT. 1.d-10) then
            ndim = 0
        end if
	
		! LOOP OVER TIME STEPS
		write(*,*) 'Time integration START'
        do INC = 1, Nsteps
            write(*,*) INC
            
            if (any(INC == output_steps)) then
                is_output_step = .TRUE.
                outINC = outINC + 1 
            else
                is_output_step = .FALSE.
            end if

            if (ndim /= 0) then
    		    call solve_mixBC(ndim, dguess, Dsolved, Ssolved, info, resid, nfev)
				write(*,*) resid
                Lp = Dsolved ! TODO OFF-DIAGONAL INPUT
			else
				Dsolved = voigt2m(d_prescribed)
				Lp = Dsolved
            end if
			
            Lp = Lp + W_prescribed
            Fglob = mmult(eye(3) + Lp*dt, Fglob)
            
            ! update the orientation matrices, all hardening variables and print output results
            stress_only = .FALSE.
            call taylor_step(Lp, Ssolved)
        end do
              
		deallocate(Q0g1_list, Rg1, STAT=s(1))
        
        if (s(1) /= 0) then
            write(*,*) 'Deallocation of orientation matrices was unsuccesfull.'
            read(*,*)
            stop
        end if
			  
    end subroutine taylor

    
                                 
                                 
    subroutine taylor_step(Lp_glob, Ssolved)
    
        integer         ::  INC, i, j, g, active_relax(3), Nact, statu, s(4), tmp(2)
        real(kind=8)    ::  Lp_glob(3,3), Lp(3,3), Dp(3,3), stress_loc(3,3), Q(3,3), Wp(3,3), W_slip(3,3), W_lattice(3,3),  &
                            tmp1(3,3), tmp2(3,3), Ssolved(3,3), coef, stress_glob(3,3), slipinc, slipratesSVD(24)
        logical         ::  success
        
        integer, allocatable        :: activesID(:)
        real(kind=8), allocatable   :: hvar(:), CRSS(:), sliprates(:)

        intent(in)      ::  Lp_glob
        intent(out)     ::  Ssolved

        
        ! consistency check of the hardening model 
        tmp = shape(hvarG1_list)
        if (tmp(1) /= Nhv) then
            write(*,*) 'Inconsistency between hardening model and number of hardening parameters.'
			write(*,*) tmp, Nhv
            read(*,*)
            stop
        end if        

        allocate(hvar(Nhv), activesID(Nslips), CRSS(Nslips), sliprates(Nslips), STAT=statu)           
        if (statu /= 0) write(*,*) 'Allocation in taylor_step unsuccesfull.'

        
		crss_mean = 0.d0
		Ssolved   = 0.d0
            
		! LOOP OVER THE GRAINS
        do j = 1, Ngrains
            Q = Qg1_list(:,:,j)
            ! update the critical resolved shear stress
            select case (hardening_law)
            case(1)
                CRSS = hvarG1_list(1,j)
            case(2)
                CRSS = hvarG1_list(74:,j)
            case default
                write(*,*) 'Hardening model not recognized.'
                read(*,*)
                stop
            end select
                        
			
            Lp = transform(Lp_glob, Q)
            Dp = getsym(Lp)

            select case (SCYLon)
            case (0)
                call findslips(CRSS, Dp, sliprates, activesID, Nact)
                call getstress_taylor(sliprates, activesID, Nact, CRSS, stress_loc)
            case (1)
            ! use of SCYL
                call SCYL_Taylor(Dp, CRSS, stress_loc, sliprates, activesID, Nact)
            end select
            
            stress_glob = transform(stress_loc, transp3(Q))
            
            if (.not. stress_only) then                 
                Wp = getskw(Lp)
                    
                ! solve Taylor ambiguity by SVD
                if ((SCYLon == 0) .and. (Nact >= 5)) then
                    call solveambSVD(P, Dp, stress_loc, CRSS, slipratesSVD, success)
                    if (success) then 
						sliprates = slipratesSVD
					else
						write(*,*) 'SVD problem'
					end if
                end if
                
                slipinc = sum(sliprates)*dt
                totalslip = totalslip + slipinc
                hvar = hvarG1_list(:,j)
                            
                if (hardening_law /= 1) then
                ! update of the hardening state variables for both grains
                    call hardening(totalslip, slipinc, activesID, Nact, hparam, hvar)
                    hvarG1_list(:,j) = hvar
                end if
                ! calculate spin tensor from slip activity
                W_slip = 0.d0
                do i=1, Nslips
                    W_slip = W_slip + Omega(:,:,i)*sliprates(i)
                end do

                W_lattice = Wp - W_slip
                tmp1 = mmult(eye(3) + W_lattice*dt/2.d0, Rg1(:,:,j))
                tmp2 = minv3(eye(3) - W_lattice*dt/2.d0)
                Rg1(:,:,j) = mmult(tmp2, tmp1)
                Q = mmult(transp3(Rg1(:,:,j)), Q0g1_list(:,:,j))
                Qg1_list(:,:,j) = Q
            end if
                
			crss_mean = crss_mean + sum(CRSS)/size(CRSS)/Ngrains    
			Ssolved   = Ssolved + stress_glob/Ngrains
        
            if ((is_output_step == .TRUE.) .and. (stress_only == .FALSE.)) &
                call save_results_taylor(outINC, j, stress_loc, stress_glob, sliprates, Q, CRSS, activesID, Nact, hvar)
            
        end do
              
    end subroutine taylor_step
    
    

    
    
     
    subroutine save_results_taylor(outINC, g, stress_loc, stress_glob, sliprates, Q, CRSS, activesID, Nact, hvar)

! dictionary of possible outputs specified by "outputvars"
!    1, 'average_stress'
!    2, 'average_slip'
!    3, 'stress_loc'
!    4, 'stress_glob'
!    5, 'euler_angles'
!    6, 'crss'
!    7, 'sliprates'
!    8, 'activesID'
!    9, 'num_actives'
!    10,'statevar'
!    11,'relaxation_slips'



        integer         ::  g, k, outINC, activesID(:), Nact
        real(kind=8)    ::  sliprates(:), stress_glob(3,3), Q(3,3), CRSS(:), hvar(:), stress_loc(3,3)
        
        intent(in)      ::  outINC, g, stress_glob, sliprates, Q, CRSS, activesID, Nact, hvar, stress_loc

        do k = 1, size(outputvars)
            select case(outputvars(k))
				case(-1)
					exit
                case(1)
                    out_average_stress(:,outINC) = out_average_stress(:,outINC) + m2voigt(stress_glob)/Ngrains
                case(2)
                    out_average_slip(outINC) = out_average_slip(outINC) + sum(sliprates)*dt/Ngrains ! average slip per grain
                case(3)
                    out_stress_locG1(:,g,outINC) = m2voigt(stress_loc)
                case(4)
                    out_stress_globG1(:,g,outINC) = m2voigt(stress_glob)
                case(5)
                    out_euler_anglesG1(:,g,outINC) = rotm2euler(Q)
                case(6)
                    out_crssG1(:,g,outINC) = CRSS
                case(7)
                    out_slipratesG1(:,g,outINC) = sliprates
                case(8)
                    out_activesIDG1(:,g,outINC) = activesID
                case(9)
                    out_num_activesG1(g,outINC) = Nact
                case(10)
                    out_statevarG1(:,g,outINC) = hvar
            end select
        end do
                        
    end subroutine save_results_taylor
                                 
                                 
                                 
                                 
!   ***************** ALAMEL MODEL MAIN ROUTINES ****************                         
    
    subroutine alamel(ndim, dguess, Dsolved, Ssolved, info, resid, nfev)
    
        integer         ::  j, INC, statu, s(4), tmp(2), ndim, info, nfev
        real(kind=8)    ::  Sdir(3,3), Lp(3,3), dguess(ndim), Dsolved(3,3), Ssolved(3,3),     &
							resid(5)
        real(kind=8), allocatable   ::  hvarG1(:), hvarG2(:)

        intent(inout)   ::  ndim, dguess
        intent(out)     ::  Dsolved, Ssolved, info, resid, nfev
		
        allocate(Q0g1_list, SOURCE=Qg1_list, STAT=s(1))
        allocate(Q0g2_list, SOURCE=Qg2_list, STAT=s(2))
        allocate(Q0b_list,  SOURCE=Qb_list,  STAT=s(3))
        allocate(Rg1(3,3,Nclus), Rg2(3,3,Nclus), STAT=s(4))
        
        if (any(s /= 0)) then
            write(*,*) 'Allocation of initial orientation matrices was unsuccesfull.'
            read(*,*)
            stop
        end if
		
        ! initialize rotation matrix to unity
        do j = 1, Nclus
            Rg1(:,:,j) = eye(3)
            Rg2(:,:,j) = eye(3)
        end do
            
        totalslipG1 = 0.d0
        totalslipG2 = 0.d0
        outINC      = 0
        Fglob       = eye(3)
		
        if (sum(dguess) .LT. 1.d-10) then
            ndim = 0
        end if
        
		call alloc_outputvars_alamel()
	
		! LOOP OVER TIME STEPS
		write(*,*) 'Time integration START'
        do INC = 1, Nsteps
            write(*,*) INC
            
            if (any(INC == output_steps)) then
                is_output_step = .TRUE.
                outINC = outINC + 1 
            else
                is_output_step = .FALSE.
            end if
            
            if (ndim /= 0) then
    		    call solve_mixBC(ndim, dguess, Dsolved, Ssolved, info, resid, nfev)
				write(*,*) resid
                Lp = Dsolved ! TODO OFF-DIAGONAL INPUT
			else
				Dsolved = voigt2m(d_prescribed)
				Lp = Dsolved
            end if
			
            Lp = Lp + W_prescribed
            Fglob = mmult(eye(3) + Lp*dt, Fglob)
            
            ! update the orientation matrices, all hardening variables and print output results
            stress_only = .FALSE.
            call alamel_step(Lp, Ssolved)
        end do
              
		deallocate(Q0g1_list, Q0g2_list, Q0b_list, Rg1, Rg2, STAT=s(1))
        
        if (s(1) /= 0) then
            write(*,*) 'Deallocation of orientation matrices was unsuccesfull.'
            read(*,*)
            stop
        end if
			  
    end subroutine alamel

    

                                 
                                 
    subroutine alamel_step(Lp, Ssolved)
    
        integer         ::  INC, i, j, g, active_relax(3), NactG1, NactG2, statu, s(4), tmp(2)
        real(kind=8)    ::  Lp(3,3), Lp_gb(3,3), Dp_gb(3,3), stress_locG1(3,3), stress_locG2(3,3), relaxation_slips(3),         & 
                            Qb(3,3), Qg1(3,3), Qg2(3,3), QbQg1T(3,3), QbQg2T(3,3), CRSS_R23, CRSS_R13, Wp(3,3), W_slip(3,3),    &
                            W_lattice(3,3), tmp1(3,3), tmp2(3,3),  Lp_g_relaxed(3,3), Ssolved(3,3), coef, Lp_gb_relaxed(3,3),   &
                            stress_globG1(3,3), stress_globG2(3,3), slipincG1, slipincG2, Dp_gb_relaxed(3,3), slipratesSVD(24)
        logical         ::  success
        
        integer, allocatable        ::  activesIDG1(:), activesIDG2(:)
        real(kind=8), allocatable   ::  hvarG1(:), hvarG2(:), Pg1(:,:,:), Pg2(:,:,:), CRSS_g1(:), CRSS_g2(:), slipratesG1(:), slipratesG2(:)

        intent(in)      ::  Lp
        intent(out)     ::  Ssolved

        
        ! consistency check of the hardening model 
        tmp = shape(hvarG1_list)
        if (tmp(1) /= Nhv) then
            write(*,*) 'Inconsistency between hardening model and number of hardening parameters.'
			write(*,*) tmp, Nhv
            read(*,*)
            stop
        end if        

        allocate( hvarG1(Nhv), hvarG2(Nhv), activesIDG1(Nslips), activesIDG2(Nslips), Pg1(3,3,Nslips), Pg2(3,3,Nslips),         & 
                    CRSS_g1(Nslips), CRSS_g2(Nslips), slipratesG1(Nslips), slipratesG2(Nslips), STAT=statu)
            
        if (statu /= 0) write(*,*) 'Allocation in alamel_step unsuccesfull.'
        
		crss_mean = 0.d0
		Ssolved   = 0.d0
            
		! LOOP OVER THE GRAINS
        do j = 1, Nclus
            Qg1 = Qg1_list(:,:,j)
            Qg2 = Qg2_list(:,:,j)
            Qb  = Qb_list(:,:,j)
                
            ! update the critical resolved shear stress
            select case (hardening_law)
            case(1)
                CRSS_g1 = hvarG1_list(1,j)
                CRSS_g2 = hvarG2_list(1,j)
            case(2)
                CRSS_g1 = hvarG1_list(74:,j)
                CRSS_g2 = hvarG2_list(74:,j)
            case default
                write(*,*) 'Hardening model not recognized.'
                read(*,*)
                stop
            end select
                        
			
            Lp_gb = transform(Lp, Qb)
            Dp_gb = getsym(Lp_gb)
                
            QbQg1T = mmult(Qb, transp3(Qg1))
            QbQg2T = mmult(Qb, transp3(Qg2))
        
            do i=1, 12
                Pg1(:,:,i) = transform(P(:,:,i), QbQg1T)
                Pg2(:,:,i) = transform(P(:,:,i), QbQg2T)
            end do
            Pg1(:,:,13:) = -Pg1(:,:,:12)
            Pg2(:,:,13:) = -Pg2(:,:,:12)

            
            select case (SCYLon)
            case(0)
                if (Nslips /= 24) write(*,*) 'Nslips must equal 24 if Single Crystal Yield Locus method is OFF'
                
			    select case (grain_interaction)
			    case(2)
				    call findslips_alamel(Pg1, Pg2, CRSS_g1, CRSS_g2, relax_penalties, Dp_gb, slipratesG1, slipratesG2,         &
									      activesIDG1, activesIDG2, NactG1, NactG2, relaxation_slips, active_relax)
				    call getstress_alamel(Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2, NactG1, NactG2,                 &
									      QbQg1T, QbQg2T, Qb, active_relax, stress_locG1, stress_locG2)
                case(3)
				    call findslips_alamel3(Pg1, Pg2, CRSS_g1, CRSS_g2, relax_penalties, Dp_gb, slipratesG1, slipratesG2,        &
									      activesIDG1, activesIDG2, NactG1, NactG2, relaxation_slips, active_relax)
				    call getstress_alamel3(Pg1, Pg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2, NactG1, NactG2,                &
									       QbQg1T, QbQg2T, Qb, active_relax, stress_locG1, stress_locG2)
                end select
                
            case(1)
                call SCYL_Alamel(CRSS_g1, CRSS_g2, relax_penalties, Dp_gb, slipratesG1, slipratesG2, activesIDG1,       &
                                 activesIDG2, NactG1, NactG2, relaxation_slips, QbQg1T, QbQg2T, stress_locG1, stress_locG2)
            end select
            
            stress_globG1 = transform(stress_locG1, transp3(Qg1))
            stress_globG2 = transform(stress_locG2, transp3(Qg2))
            
            ! velocity gradients after relaxation for each grain  (expressed in grain boundary sys)
            if (.not. stress_only) then
                do g=1, 2
                    Lp_gb_relaxed = Lp_gb + R1*relaxation_slips(1)*(-1.d0)**g + R2*relaxation_slips(2)*(-1.d0)**g
                    Dp_gb_relaxed = Dp_gb + symR1*relaxation_slips(1)*(-1.d0)**g + symR2*relaxation_slips(2)*(-1.d0)**g
                    ! in case of ALAMEL3
					if (grain_interaction==3) then
                        Lp_gb_relaxed = Lp_gb_relaxed - R3*relaxation_slips(3)
                        Dp_gb_relaxed = Dp_gb_relaxed - symR3*relaxation_slips(3)
					end if
                        
                    ! do update for both grains in the cluster
                    select case (g)
                        case (1)
                            ! into grain's coord sys
                            Lp_g_relaxed = transform(Lp_gb_relaxed, transp3(QbQg1T))     
                            Wp = getskw(Lp_g_relaxed)
                            
                            ! solve Taylor ambiguity by SVD for each grain separately
                            if ((NactG1 >= 5) .and. (SCYLon == 0)) then
                                call solveambSVD(Pg1, Dp_gb_relaxed, stress_locG1, CRSS_g1, slipratesSVD, success)
                                if (success) slipratesG1 = slipratesSVD
                            end if
                                
                            slipincG1 = sum(slipratesG1)*dt
                            totalslipG1 = totalslipG1 + slipincG1
                            hvarG1 = hvarG1_list(:,j)
                            
                            if (hardening_law /= 1) then
                            ! update of the hardening state variables for both grains
                                call hardening(totalslipG1, slipincG1, activesIDG1, NactG1, hparam, hvarG1)
                                hvarG1_list(:,j) = hvarG1
                            end if
                            ! calculate spin tensor from slip activity
                            W_slip = 0.d0
                            do i=1, Nslips
                                W_slip = W_slip + Omega(:,:,i)*slipratesG1(i)
                            end do

                            W_lattice = Wp - W_slip
                            W_lattice = transform(W_lattice, Rg1(:,:,j))
                            tmp1 = mmult(eye(3) + W_lattice*dt/2.d0, Rg1(:,:,j))
                            tmp2 = minv3(eye(3) - W_lattice*dt/2.d0)
                            Rg1(:,:,j) = mmult(tmp2, tmp1)
                            Qg1 = mmult(transp3(Rg1(:,:,j)), Q0g1_list(:,:,j))
                            Qg1_list(:,:,j) = Qg1
                        case (2)
                            ! into grain's coord sys
                            Lp_g_relaxed = transform(Lp_gb_relaxed, transp3(QbQg2T))     
                            Wp = getskw(Lp_g_relaxed)
                            
                            ! solve Taylor ambiguity by SVD for each grain separately
                            if ((NactG2 >= 5) .and. (SCYLon == 0)) then
                                call solveambSVD(Pg2, Dp_gb_relaxed, stress_locG2, CRSS_g2, slipratesSVD, success)
                                if (success) slipratesG2 = slipratesSVD
                            end if
                                
                            slipincG2 = sum(slipratesG2)*dt
                            totalslipG2 = totalslipG2 + slipincG2
                            hvarG2 = hvarG2_list(:,j)
                            
                            if (hardening_law /= 1) then
                            ! update of the hardening state variables for both grains
                                call hardening(totalslipG2, slipincG2, activesIDG2, NactG2, hparam, hvarG2)
                                hvarG2_list(:,j) = hvarG2
                            end if
                            ! calculate spin tensor from slip activity
                            W_slip = 0.d0
                            do i=1, Nslips
                                W_slip = W_slip + Omega(:,:,i)*slipratesG2(i)
                            end do

                            W_lattice = Wp - W_slip
                            W_lattice = transform(W_lattice, Rg2(:,:,j))
                            tmp1 = mmult(eye(3) + W_lattice*dt/2.d0, Rg2(:,:,j))
                            tmp2 = minv3(eye(3) - W_lattice*dt/2.d0)
                            Rg2(:,:,j) = mmult(tmp2, tmp1)
                            Qg2 = mmult(transp3(Rg2(:,:,j)), Q0g2_list(:,:,j))
                            Qg2_list(:,:,j) = Qg2
                    end select
                end do
                
				if (rotate_boundary == 1) then
                    call updateOrientation(Q0b_list(:,:,j), Qb, coef)
                    Qb_list(:,:,j) = Qb
                end if

            end if
                
			crss_mean = crss_mean + (sum(CRSS_g1)/size(CRSS_g1) + sum(CRSS_g2)/size(CRSS_g2))/(2*Nclus)    
			Ssolved   = Ssolved + (stress_globG1 + stress_globG2)/(2*Nclus)
        
            if ((is_output_step == .TRUE.) .and. (stress_only == .FALSE.)) then
                call save_results(outINC, j, stress_locG1, stress_locG2, stress_globG1, stress_globG2,                              &
                        slipratesG1, slipratesG2, relaxation_slips, Qg1, Qg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,           &
                        NactG1, NactG2, hvarG1, hvarG2)
            end if
            
        end do
              
    end subroutine alamel_step
    
    

    
    
     
    subroutine save_results(outINC, c, stress_locG1, stress_locG2, stress_globG1, stress_globG2,                                    &
                            slipratesG1, slipratesG2, relaxation_slips, Qg1, Qg2, CRSS_g1, CRSS_g2, activesIDG1, activesIDG2,       &
                            NactG1, NactG2, hvarG1, hvarG2)

!    1,'average_stress'
!    2,'average_slip'
!    3,'stress_loc'
!    4,'stress_glob'
!    5,'euler_angles'
!    6,'crss'
!    7,'sliprates'
!    8,'activesID'
!    9,'num_actives'
!   10,'statevar'
!   11,'relaxsliprates'



        integer         ::  c, k, outINC, activesIDG1(:), activesIDG2(:), NactG1, NactG2
        real(kind=8)    ::  slipratesG1(:), slipratesG2(:), stress_globG1(3,3), stress_globG2(3,3), Qg1(3,3), Qg2(3,3),             & 
                            CRSS_g1(:), CRSS_g2(:), hvarG1(:), hvarG2(:), stress_locG1(3,3),                                        &
                            stress_locG2(3,3), relaxation_slips(2)
        
        intent(in)      ::  outINC, c, stress_globG1, stress_globG2, slipratesG1, slipratesG2, Qg1, Qg2,                            &
                            CRSS_g1, CRSS_g2, activesIDG1, activesIDG2, NactG1, NactG2, hvarG1, hvarG2, relaxation_slips,           &
                            stress_locG1, stress_locG2

        do k = 1, size(outputvars)
            select case(outputvars(k))
				case(-1)
					exit
                case(1)
                    out_average_stress(:,outINC) = out_average_stress(:,outINC) + (m2voigt(stress_globG1) + m2voigt(stress_globG2))/(2*Nclus)
                case(2)
                    out_average_slip(outINC) = out_average_slip(outINC) + (sum(slipratesG1) + sum(slipratesG2))*dt/(2*Nclus) ! average slip per grain
                case(3)
                    out_stress_locG1(:,c,outINC) = m2voigt(stress_locG1)
                    out_stress_locG2(:,c,outINC) = m2voigt(stress_locG2)
                case(4)
                    out_stress_globG1(:,c,outINC) = m2voigt(stress_globG1)
                    out_stress_globG2(:,c,outINC) = m2voigt(stress_globG2)
                case(5)
                    out_euler_anglesG1(:,c,outINC) = rotm2euler(Qg1)
                    out_euler_anglesG2(:,c,outINC) = rotm2euler(Qg2)
                case(6)
                    out_crssG1(:,c,outINC) = CRSS_g1
                    out_crssG2(:,c,outINC) = CRSS_g2
                case(7)
                    out_slipratesG1(:,c,outINC) = slipratesG1
                    out_slipratesG2(:,c,outINC) = slipratesG2
                case(8)
                    out_activesIDG1(:,c,outINC) = activesIDG1
                    out_activesIDG2(:,c,outINC) = activesIDG2
                case(9)
                    out_num_activesG1(c,outINC) = NactG1
                    out_num_activesG2(c,outINC) = NactG2
                case(10)
                    out_statevarG1(:,c,outINC) = hvarG1
                    out_statevarG2(:,c,outINC) = hvarG2
                case(11)
                    out_relaxsliprates(:,c,outINC) = relaxation_slips
            end select
        end do
                        
    end subroutine save_results

                            
                            
    ! Single Crystal Yield Locus
     
    subroutine SCYL_Taylor(Dp, CRSS, S, sliprates, activesID, Nact)
        
        ! input
        real(kind=8)    ::  Dp(3,3), CRSS(:)
        ! output
        integer         ::  Nact, activesID(size(CRSS))
        real(kind=8)    ::  sliprates(size(CRSS)), S(3,3)
        !
        integer         ::  i, j, mask(5)
        real(kind=8)    ::  phi, lamdot, xi(size(CRSS)), res(3,3), dphi(3,3), errsig, SCYLexp_orig
        
        intent(in)      ::  Dp, CRSS
        intent(out)     ::  S, sliprates, activesID, Nact
    
        
        S = Dp
        SCYL_Dp = m2voigt(Dp)
        call SCYLphi(S, CRSS, phi)
        S = 1.d0/(phi + 1.d0)*S
        
        errsig = 1.d-9
        if (SCYLexp > 20) then
            SCYLexp_orig = SCYLexp
            SCYLexp = 20.d0
            call findstress(S, Dp, CRSS, errsig)
            SCYLexp = SCYLexp_orig
            errsig = 1.d-9
        end if
        call findstress(S, Dp, CRSS, errsig)

        lamdot = inner(Dp, S)       
        call SCYLgmdot(S, CRSS, lamdot, sliprates, activesID, Nact)
              
    end subroutine SCYL_Taylor
    
    
    subroutine SCYLphi(S, CRSS, phi)
    
        ! input
        real(kind=8)    ::  CRSS(:)
        ! output
        real(kind=8)    ::  phi
        !
        integer         ::  i
        real(kind=8)    ::  S(3,3), Sn(3,3), xi(size(CRSS)), VM, tmp(size(CRSS)), maxtmp
        
        intent(in)      ::  CRSS
        intent(out)     ::  phi
        
        xi = 1.d0
        phi = 0.d0
        ! scaling by the maximum element
        tmp = 0.d0
        do i=1, Nslips
            tmp(i) = max(0.d0, inner(S, P(:,:,i))/CRSS(i))
        end do
        maxtmp = maxval(tmp, dim=1)
        tmp = tmp/maxtmp
        
        do i=1, Nslips
            phi = phi + xi(i)* tmp(i)**SCYLexp
        end do  
        phi = phi**(1.d0/SCYLexp)
        phi = maxtmp*phi - 1.d0
    
    end subroutine SCYLphi
    
    
    subroutine SCYLdphi(S, CRSS, dphi)
    
        ! input
        real(kind=8)    ::  CRSS(:), S(3,3)
        ! output
        real(kind=8)    ::  dphi(3,3)
        !
        integer         ::  i, j
        real(kind=8)    ::  Sn(3,3), xi(size(CRSS)), phi, VM, tmp(size(CRSS)), maxtmp
        
        intent(in)      ::  S, CRSS
        intent(out)     ::  dphi
        
        call SCYLphi(S, CRSS, phi)
        xi = 1.d0
        dphi = 0.d0
        ! scaling by the maximum element
        tmp = 0.d0
        do i=1, Nslips
            tmp(i) = max(0.d0, inner(S, P(:,:,i))/((phi+1.d0)*CRSS(i)))
        end do
        maxtmp = maxval(tmp, dim=1)
        tmp = tmp/maxtmp
        
        do i=1, Nslips
            dphi = dphi + xi(i)* tmp(i)**(SCYLexp-1) *P(:,:,i)/CRSS(i)
        end do 
        
        dphi = dphi* maxtmp**(SCYLexp-1)
    
    end subroutine SCYLdphi
    
    
    subroutine SCYLddphi(S, CRSS, ddphi)
    
        ! input
        real(kind=8)    ::  CRSS(:), S(3,3)
        ! output
        real(kind=8)    ::  ddphi(5,5)
        !
        integer         ::  i, j
        real(kind=8)    ::  Sn(3,3), dphi(3,3), dphivec(5), xi(size(CRSS)), phi, Pvec(5), VM, tmp(size(CRSS)), maxtmp
        
        intent(in)      ::  S, CRSS
        intent(out)     ::  ddphi
        
        call SCYLphi(S, CRSS, phi)
        call SCYLdphi(S, CRSS, dphi)
        dphivec = m2voigt(dphi)
        
        xi = 1.d0
        ddphi = 0.d0
        ! scaling by the maximum element
        tmp = 0.d0
        do i=1, Nslips
            tmp(i) = max(0.d0, inner(S, P(:,:,i))/((phi+1.d0)*CRSS(i)))
        end do
        maxtmp = maxval(tmp, dim=1)
        tmp = tmp/maxtmp
        
        do i=1, Nslips
            Pvec = m2voigt(P(:,:,i))
            ddphi = ddphi + xi(i)*(tmp(i)**(SCYLexp-2)) * (1.d0/CRSS(i))**2 * tensordot(Pvec, Pvec)
        end do         
        ddphi = ddphi* maxtmp**(SCYLexp-2)
        ddphi = (SCYLexp-1)/(phi + 1.d0) * (ddphi - tensordot(dphivec, dphivec))
    
    end subroutine SCYLddphi
    
    
    subroutine SCYLgmdot(S, CRSS, lamdot, gm, activesID, Nact)
    
        ! input
        real(kind=8)    ::  CRSS(:), S(3,3), lamdot
        ! output
        integer         ::  activesID(size(CRSS)), Nact
        real(kind=8)    ::  gm(size(CRSS)) 
        !
        integer         ::  i, j
        real(kind=8)    ::  xi(size(CRSS)), tmp(size(CRSS)), maxtmp
    
        intent(in)      ::  S, CRSS, lamdot
        intent(out)     ::  gm, activesID, Nact
        
        gm = 0.d0
        xi = 1.d0
        tmp = 0.d0
        do i=1, Nslips
            tmp(i) = max(0.d0, inner(S, P(:,:,i))/CRSS(i))
        end do
        maxtmp = maxval(tmp, dim=1)
        tmp = tmp/maxtmp
        
        do i=1, Nslips
            gm(i) = xi(i)*tmp(i)**(SCYLexp-1) /CRSS(i)
        end do
        gm = lamdot*gm*maxtmp**(SCYLexp-1)
        j = 0
        do i = 1, Nslips 
            if (gm(i) > 1.d-6) then
                j = j + 1
                activesID(j) = i
            end if
        end do
        Nact = j
    
    
    end subroutine SCYLgmdot
    
                       
    subroutine SCYL_Alamel(CRSS_g1, CRSS_g2, relax_penalties, Dp_gb, slipratesG1, slipratesG2, activesIDG1, activesIDG2,    &
                           NactG1, NactG2, relaxation_slips, QbQg1T, QbQg2T, SlocG1, SlocG2)
    
        ! input
        real(kind=8)            ::  CRSS_g1(:), CRSS_g2(:), relax_penalties(3), Dp_gb(3,3), QbQg1T(3,3), QbQg2T(3,3)
        ! output
        integer                 ::  activesIDG1(size(CRSS_g1)), activesIDG2(size(CRSS_g1)), NactG1, NactG2
        real(kind=8)            ::  slipratesG1(size(CRSS_g1)), slipratesG2(size(CRSS_g1)), SlocG1(3,3), SlocG2(3,3),       &
                                    relaxation_slips(3)
        !
        real(kind=8)            ::  Dp_gb_relaxed(3,3), lamdot, Dg1(3,3), Dg2(3,3)
        ! decalration of variables for minpack
        integer                 ::  info, nfev
        integer, parameter      ::  n = 2, ml = 2, mu = 2, ldfjac = 2, lr = 3, maxfev = 200, mode = 1, nprint = 0
        double precision        ::  factor = 100., xtol = 1.d-6, epsfcn
        double precision        ::  x(n), fvec(n), diag(n), fjac(n,n), r(lr),                                               &
                                    qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
        
        intent(in)              ::  CRSS_g1, CRSS_g2, relax_penalties, Dp_gb, QbQg1T, QbQg2T
        intent(out)             ::  slipratesG1, slipratesG2, SlocG1, SlocG2, activesIDG1, activesIDG2, NactG1, NactG2,     &
                                    relaxation_slips

        
        ! make those global                  
        SCYL_QbQg1T = QbQg1T
        SCYL_QbQg2T = QbQg2T  
        SCYL_Dp     = m2voigt(Dp_gb)     
        SCYL_CRSSg1(:Nslips) = CRSS_g1
        SCYL_CRSSg2(:Nslips) = CRSS_g2
        
        x = 0.d0
        epsfcn = eps4jacobian
        crss_mean = (sum(CRSS_g1)+sum(CRSS_g2))/(size(CRSS_g1)+size(CRSS_g2))

        call hybrd(fcn_scyl_ala,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,                                                     &
                            mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,                                                  &
                            qtf,wa1,wa2,wa3,wa4)
        
        relaxation_slips = 0.d0
        relaxation_slips(:2) = x*2.d0
        ! grain 1
        SlocG1 = SCYL_S1solved
        Dp_gb_relaxed = Dp_gb - symR1*relaxation_slips(1) - symR2*relaxation_slips(2)
        Dg1  = transform(Dp_gb_relaxed, transp3(QbQg1T))
        lamdot = inner(Dg1, SlocG1)
        call SCYLgmdot(SlocG1, CRSS_g1, lamdot, slipratesG1, activesIDG1, NactG1)
        
        ! grain 2
        SlocG2 = SCYL_S2solved
        Dp_gb_relaxed = Dp_gb + symR1*relaxation_slips(1) + symR2*relaxation_slips(2)
        Dg2  = transform(Dp_gb_relaxed, transp3(QbQg2T))
        lamdot = inner(Dg2, SlocG2)
        call SCYLgmdot(SlocG2, CRSS_g2, lamdot, slipratesG2, activesIDG2, NactG2)
        
    end subroutine SCYL_Alamel
    
    
    
                            
    subroutine solve_mixBC(n, dguess, Dsolved, Ssolved, info, resid, nfev)
    
        integer         ::  i, j, stt, n, mask(5), info, nfev
        real(kind=8)    ::  xtol, d(5), dguess(n), Dsolved(3,3), sdirn(5), vm, rhs5(5), Ssolved(3,3), resid(5),     &
                            x(n), fvec(n)

        intent(in)      ::  n, dguess
        intent(out)     ::  Dsolved, Ssolved, info, resid, nfev


		! orientations and hardening variables must be not updated during iterations on BC
        stress_only = .TRUE.
    
!       assembling the rhs vector for the residuals
        allocate(rhs(n), STAT=stt)          
        if (stt /= 0) write(*,*) 'Allocation of rhs failed.'
            
!       1. first VonMises-normalization of elements defining stress direction sdir -> sdirn
        call vonMises_dvec(sdir_prescribed, 'stress', ind_sdir, vm, sdirn)       
        rhs5 = sdirn
!       2. add the absolute stress terms
        forall (i=1:5, ind_sabs(i) == 1) rhs5(i) = sabs_prescribed(i)
!       3. slice only those n with either ind_sdir or ind_sabs == TRUE     
        j = 0
        rhs = 0.d0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                rhs(j) = rhs5(i)
            end if
        end do
		
		
        if (j /= n) then
            write(*,*) 'Problem dimension mismatch!'
            read(*,*)
            stop
        end if
        
        x   = dguess
        call solver(fcn_nojac, n, x, fvec, nfev, info)
		
		resid = -1000000.d0
		resid(1:n) = fvec
            
        d = 0.d0
        forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
        j = 0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                d(i) = x(j)
            end if
        end do
            
        d = rescaleD(d, ind_d)
        Dsolved = voigt2m(d)
        
        if (grain_interaction == 1) then
            call taylor_step(Dsolved, Ssolved)
        elseif ((grain_interaction == 2) .or. (grain_interaction == 3)) then
            call alamel_step(Dsolved, Ssolved)
        end if
        
        deallocate(rhs, STAT=stt)
        if (stt /= 0) write(*,*) 'Deallocation of rhs failed.'
        
    end subroutine solve_mixBC
                            
    
    
    subroutine solver(fcn, n, x, fvec, nfev, info)
		
        ! decalration of variables for minpack
        integer n, ml, mu, info, nfev, ldfjac, lr
        integer, parameter      ::  maxfev = 200, mode = 1, nprint = 0
        double precision        ::  factor = 100., xtol = 1.d-6, epsfcn
        double precision x(n), fvec(n), diag(n), fjac(n,n),                                                                 &
                         qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
        double precision, allocatable  :: r(:)
						
		integer			::  niter,i,j,maxresloc, mask(n), stt
		real(kind=8)	::	dx, d(5), maxres, sol_list(n,20), maxres_list(20), fvec_list(n,20), x0(n), xnew(n), vm,         &
                            SCYLexp_orig, SCYLexp_safe, xsafe(n), k
        logical         ::  success
		
        external fcn
        interface
            subroutine fcn(n, x, fvec, iflag)
            integer n, iflag
            double precision x(n), fvec(n)
            end subroutine fcn
        end interface
        
        lr     = (n*(n+1))/2
        ldfjac = n
        ml     = n
        mu     = n
        epsfcn = eps4jacobian
!       assembling the rhs vector for the residuals
        allocate(r(lr), STAT=stt)         
        if (stt /= 0) write(*,*) 'Allocation of r failed.'
        
        call hybrd(fcn_nojac,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,                                                        &
                            mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,                                                  &
                            qtf,wa1,wa2,wa3,wa4)
        maxres = maxval(abs(fvec))
        niter = 1
        
        select case (SCYLon)
        case (1)
            if (SCYLexp > 6) then
                SCYLexp_orig = SCYLexp
                SCYLexp = 6              
                call hybrd(fcn_nojac,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,                                                &
                                mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,                                              &
                                qtf,wa1,wa2,wa3,wa4)
                maxres = maxval(abs(fvec))
                if (maxres < 1.d-3) then
                    xsafe = x
                    k = 6.d0
                    SCYLexp_safe = SCYLexp
                    success = .TRUE.
                    do while (k > 0)
                        if (success) then
                            SCYLexp = min(SCYLexp + k, SCYLexp_orig)
                        else
                            if (k >= 2.d0) then
                                k = k - 1
                            else
                                k = k/2.d0
                            end if
                            SCYLexp = min(SCYLexp_safe + k, SCYLexp_orig)
                            x = xsafe
                        end if
                            
                        call hybrd(fcn_nojac,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,                                        &
                                    mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,                                          &
                                    qtf,wa1,wa2,wa3,wa4)
                        maxres = maxval(abs(fvec))
                        if (maxres < 1.d-2) then
                            d = 0.d0
				            forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
				            j = 0
				            do i=1, 5
					            if (ind_d(i) == 0) then
						            j = j + 1
						            d(i) = x(j)
					            end if
				            end do
				            d = rescaleD(d, ind_d)
                            x = d
                            xsafe = x
                            SCYLexp_safe = SCYLexp
                            success = .TRUE.
                        else
                            success = .FALSE.
                        end if
                    end do
                else
                    write(*,*) 'SCYL poor for exp 6.'
                end if
            end if
    	        	
        case (0) 
            ! NOT the case of Single Crystal Yield Locus
		    maxres_list = 0.d0
		    sol_list    = 0.d0
		    maxres_list(niter) = maxres
		    sol_list(:,niter)  = x
            fvec_list(:,niter) = fvec
		    x0 = x
            ! do extra iterations on randomly perturbated initial guess
            do while ((maxres > 1.d-2) .and. (niter < 10))
				! prepare a randomly perturbated guess
				d = 0.d0
				forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
				j = 0
				do i=1, 5
					if (ind_d(i) == 0) then
						j = j + 1
						d(i) = x0(j)
					end if
				end do
				d = rescaleD(d, ind_d)
				call random_number(dx)
				x = d + 0.2*dx - 0.1d0
				niter = niter + 1
				! run iteration with a new guess
                call hybrd(fcn_nojac,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,                                                    &
                                mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,                                                  &
                                qtf,wa1,wa2,wa3,wa4)     
                maxres = maxval(abs(fvec))
				maxres_list(niter) = maxres
				sol_list(:,niter)  = x
                fvec_list(:,niter) = fvec
            end do
			
			if ((niter == 10) .and. (maxres > 1.d-2)) then
				x = sol_list(:,minloc(maxres_list, dim=1))
                fvec = fvec_list(:,minloc(maxres_list, dim=1))
				write(*,1001) minval(maxres_list, dim=1)
	1001		format('Requirement on residual not satified, min max residual is:', F8.5)
			end if
        end select
        
        deallocate(r, STAT=stt)         
        if (stt /= 0) write(*,*) 'Deallocation of r failed.'
                            
    end subroutine solver                     
    
        
    subroutine updateOrientation(QB0, QB, coef)
    
        real(kind=8)    ::  C0(3,3), C(3,3), QB0(3,3), QB(3,3), invF(3,3), vec1(3), vec2(3), vec3(3), &
							E(2,2), eigvec(2,2), eigval(3), vv1(3), vv2(3), Se, S, coef, a1, b1, c1
		real(kind=8), parameter :: p = 1.6075d0, three1(3) = (/1.d0, 1.d0, 1.d0/)
		integer			::  nrot
		
		intent(in)		:: QB0
		intent(out)		:: QB, coef
        
        ! initializing the ellipsoid as a sphere
        C0 = eye(3)
        invF = minv3(Fglob)
        ! get deformed ellipsoid
        C = mmult(transp3(invF), mmult(C0, invF))
        ! get axes x y z of coord system given by QB matrix
        vec1 = QB0(1,:)
        vec2 = QB0(2,:)
        vec3 = QB0(3,:)

        ! apply deformation on a planar interface given by x and y axes, and
        ! obtain normal axis z and after renormalizing and
        ! reorthonormalizing obtain new QB matrix
        vec1 = mdot(Fglob, vec1)
        vec1 = vec1/norm(vec1)
        vec2 = mdot(Fglob, vec2)
        vec3 = cross(vec1, vec2)
        vec3 = vec3/norm(vec3)
        vec2 = cross(vec3,vec1)
        vec2 = vec2/norm(vec2)
        QB(1,:) = vec1
        QB(2,:) = vec2
        QB(3,:) = vec3
		
		! for the case of ALAMEL3
		if (grain_interaction == 3) then 
			! rotate ellipsoid to this coord system
			C = mmult(QB, mmult(C, transp3(QB)))
			
			! get ellipse as a intersection of plane given by x and y from QB
			! and ellipsoid (coord z = 0)
			E = C(1:2,1:2)
			call jacobi(E,2,2,eigval,eigvec,nrot)
			! calc area of this ellipse 
			Se = pi/sqrt(eigval(1)*eigval(2))
			
			! rotate QB so that its x and y axes are aligned with principal
			! direction of ellipse representing interface boundary
			vv1(1:2) = eigvec(1:2,1)
			vv1(3) = 0.d0
			vv1 = mdot(transp3(QB), vv1)
			vv2(1:2) = eigvec(1:2,2)
			vv2(3) = 0.d0
			vv2 = mdot(transp3(QB), vv2)
			
			if ( all(sign(three1, cross(vv1, vv2)) == sign(three1, QB(3,:))) ) then
				QB(1,:) = vv1
				QB(2,:) = vv2
			else
				QB(1,:) = vv2
				QB(2,:) = vv1
			end if
		end if
		
        
		! calc surface of the ellipsoid
        call jacobi(C,3,3,eigval,eigvec,nrot)
        a1 = 1.d0/sqrt(eigval(1))
        b1 = 1.d0/sqrt(eigval(2))
        c1 = 1.d0/sqrt(eigval(3))
        S = 4.d0*pi*(( (a1*b1)**p + (a1*c1)**p + (b1*c1)**p )/3)**(1.d0/p)
        ! this ratio is coef for orientation in output file
        coef = Se/S
		
    end subroutine updateOrientation
                            
                            
    subroutine hardening(totalslip, slipinc, activesID, Nact, hparam, hvar)

        ! input
        integer             ::  activesID(:), Nact
        real(kind=8)        ::  hparam(:), totalslip, slipinc
        ! output
        real(kind=8)        ::  hvar(:)
        !
        integer             ::  reverseID(size(activesID)), Npas, ID, passiveID(size(activesID))
        real(kind=8)        ::  tau0, th2, th3, th4, gm2, gm3, qP, qL, qR, qLR,                                             &
                                gmP, gmL, gmR, tauI, tauL(size(activesID)), tauP(size(activesID)), tauR(size(activesID)),   &
                                tauLsat(size(activesID)), tauPsat(size(activesID)), tauRsat(size(activesID)),               &
                                hL(size(activesID)), hP(size(activesID)), hR(size(activesID))
        character(20)       ::  alg
        
        intent(in)          ::  totalslip, slipinc, activesID, Nact, hparam
        intent(inout)       ::  hvar
    
        select case (hardening_law)
        case(1)
            tau0 = hparam(1)
            th2  = hparam(2)
            th3  = hparam(3) 
            th4  = hparam(4)
            gm2  = hparam(5)
            gm3  = hparam(6)
            qP   = hparam(7)
            qL   = hparam(8)
            qR   = hparam(9)
            qLR  = hparam(10)
            gmP  = hparam(11)
            gmL  = hparam(12)
            gmR  = hparam(13)
            
            tauI = hvar(1)
            tauL = hvar(2:1+Nslips)
            tauP = hvar(2+Nslips:1+2*Nslips)
            tauR = hvar(2+2*Nslips:1+3*Nslips)
            
            tauLsat = 0.d0
            tauPsat = 0.d0
            tauRsat = 0.d0
            
            reverseID(:Nact)   = activesID(:Nact) + (Nslips/2)*(1-2*(activesID(:Nact)/(Nslips/2)))
            
            Npas = 0
            do ID = 1, Nslips
                if (any(ID == activesID(:Nact)) .or. any(ID == reverseID(:Nact))) then
                    cycle
                else
                    Npas = Npas + 1
                    passiveID(Npas) = ID
                end if
            end do

            ! saturation stresses have no history and are calculated based on the current slip activity status
            
            ! for the active slip systems only
            tauPsat(activesID(:Nact)) = qP*tauI
            ! for reverse of active slip systems only
            tauLsat(reverseID(:Nact)) = qLR*tauI
            tauRsat(reverseID(:Nact)) = -qR*tauI
            ! for passive slip systems only
            tauLsat(passiveID(:Npas)) = qL*tauI
            
            alg = 'expl'
            if (alg == 'expl') then
                hL = (tauLsat - tauL)/gmL
                hP = (tauPsat - tauP)/gmP
                hR = (tauRsat - tauR)/gmR
                tauL = tauL + hL*slipinc   ! latent hardening
                tauP = tauP + hP*slipinc   ! polarisation
                tauR = tauR + hR*slipinc   ! reverse hardening - Bauschinger effect
            else if (alg == 'impl') then
                tauL = (1.d0/(gmL+slipinc))*(tauLsat*slipinc + tauL*gmL)
                tauP = (1.d0/(gmP+slipinc))*(tauPsat*slipinc + tauP*gmP)
                tauR = (1.d0/(gmR+slipinc))*(tauRsat*slipinc + tauR*gmR)
            end if

            ! isotropic hardening (extended Voce law)
            tauI = tau0 + th2*gm2*(1.d0-exp(-totalslip/gm2)) + th3*gm3*(1.d0-exp(-totalslip/gm3)) + th4*totalslip

            ! update the hardening variables
            hvar(1)                     = tauI 
            hvar(2:1+Nslips)            = tauL 
            hvar(2+Nslips:1+2*Nslips)   = tauP 
            hvar(2+2*Nslips:1+3*Nslips) = tauR
            hvar(2+3*Nslips:)           = tauI + tauL + tauP + tauR
!        case default
!            write(*,*) 'Hardening model not recognized.'
!            stop
        end select
        
    end subroutine hardening
                            
                            
    subroutine YLfull(N, T_list, input_type, tol, YL, j)
      
        integer             ::  i, j, N, statu, tmp(2), info, nfev, Nout, mask(5)
        real(kind=8)        ::  T11, T22, T33, T12, T13, T23, T(3,3), T_list(5,N), YL(5,N), Ssolved(3,3), Dsolved(3,3), &
                                resid(5), dguess(5), tol, VM
        character(len=6)    ::  input_type
        
        intent(in)          ::  N, T_list, input_type, tol
        intent(out)         ::  YL, j
        
        stress_only = .TRUE.
        
        j = 0
        do i=1, N  
            T11 = 1.d0/6.d0*(sqrt(3.d0)+3.d0)*T_list(1,i) + 1.d0/6.d0*(sqrt(3.d0)-3.d0)*T_list(2,i)
            T22 = 1.d0/6.d0*(sqrt(3.d0)-3.d0)*T_list(1,i) + 1.d0/6.d0*(sqrt(3.d0)+3.d0)*T_list(2,i)
            T33 = -T11-T22
            T12 = sqrt(2.d0)/2.d0*T_list(5,i)
            T13 = sqrt(2.d0)/2.d0*T_list(4,i)
            T23 = sqrt(2.d0)/2.d0*T_list(3,i)
            T   = reshape( (/ T11, T12, T13, T12, T22, T23, T13, T23, T33 /), (/3,3/) )
            
            if (input_type == 'strain') then
                if (grain_interaction == 1) then
                    call taylor_step(T, Ssolved)
                else
                    call alamel_step(T, Ssolved)
                end if
                resid = 1.d-100
            elseif (input_type == 'stress') then
                
                dguess = m2voigt(T)
                mask = 1
                call vonMises_dvec(dguess, 'stress', mask, VM, sdir_prescribed)
                ind_sdir = 1
                ind_sabs = 0
                ind_d    = 0
                sabs_prescribed = 0.d0
                d_prescribed    = 0.d0
                call solve_mixBC(5, dguess, Dsolved, Ssolved, info, resid, nfev)
                write(*,*) i
            end if
            
            if (maxval(resid, dim=1) < tol) then
                j = j + 1
                YL(:,j) = m2voigt(Ssolved)
            end if
			
            if (modulo(i,1000) == 0) then
                write(*,*) i
            end if
        end do
    end subroutine ylfull
    
    
    

    subroutine fcn(n, x, fvec, fjac, ldfjac, iflag)
        integer n, ldfjac, iflag
        double precision x(n), fvec(n), fjac(ldfjac,n)

    !----------
    !if iflag = 1 calculate the functions at x and
    !return this vector in fvec. do not alter fjac.
    !if iflag = 2 calculate the jacobian at x and
    !return this matrix in fjac. do not alter fvec.
    !---------
        integer             ::  i, j, k
        real(kind=8)        ::  Dp(3,3), Ssolved(3,3), xi(n), tmp(5), y0(n), d(5), vm, y(n)

        stress_only = .TRUE.
        
        d = 0.d0
        forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
        j = 0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                d(i) = x(j)
            end if
        end do
        
        d = rescaleD(d, ind_d)
        Dp = voigt2m(d)
        Dp = Dp/VonMises(Dp,'strain')
        d = m2voigt(Dp)
        
        select case(iflag)
        case(1)
            if (grain_interaction == 1) then
                call taylor_step(Dp, Ssolved)
            elseif ((grain_interaction == 2) .or. (grain_interaction == 3)) then
                call alamel_step(Dp, Ssolved)
            end if
            ! renormalize stress elemetns but only those related to stress direction
            call vonMises_dvec(m2voigt(Ssolved), 'stress', ind_sdir, vm, tmp)
            ! slice only those n with either ind_sdir or ind_sabs == TRUE     
            j = 0
            y = 0.d0
            do i=1, 5
                if (ind_d(i) == 0) then
                    j = j + 1
                    y(j) = tmp(i)
                end if
            end do
            fvec = y - rhs
            
            j = 0
            do i=1, 5
                if (ind_d(i) == 0) then
                    j = j + 1
                    if (ind_sabs(i) == 1) then
                        fvec(j) = fvec(j)/crss_mean
                    end if
                end if
            end do 
            
        case(2)
            if (grain_interaction == 1) then
                call taylor_step(Dp, Ssolved)
            elseif ((grain_interaction == 2) .or. (grain_interaction == 3)) then
                call alamel_step(Dp, Ssolved)
            end if
            ! renormalize stress elemetns but only those related to stress direction
            call vonMises_dvec(m2voigt(Ssolved), 'stress', ind_sdir, vm, tmp)
            ! slice only those n with either ind_sdir or ind_sabs == TRUE     
            j = 0
            y0 = 0.d0
            do i=1, 5
                if (ind_d(i) == 0) then
                    j = j + 1
                    y0(j) = tmp(i)
                end if
            end do
            
            j = 0
            x = 0.d0
            do i=1, 5
                if (ind_d(i) == 0) then
                    j = j + 1
                    x(j) = d(i)
                end if
            end do
            
            do i = 1, n
                xi = x
                xi(i) = x(i) + eps4jacobian/2.d0
                d = 0.d0
                forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
                j = 0
                do k=1, 5
                    if (ind_d(k) == 0) then
                        j = j + 1
                        d(k) = xi(j)
                    end if
                end do
                d = rescaleD(d, ind_d)
                Dp = voigt2m(d)
                Dp = Dp/VonMises(Dp,'strain')
                
                if (grain_interaction == 1) then
                    call taylor_step(Dp, Ssolved)
                elseif ((grain_interaction == 2) .or. (grain_interaction == 3)) then
                    call alamel_step(Dp, Ssolved)
                end if
                ! renormalize stress elemetns but only those related to stress direction
                call vonMises_dvec(m2voigt(Ssolved), 'stress', ind_sdir, vm, tmp)
                ! slice only those n with either ind_sdir or ind_sabs == TRUE     
                j = 0
                y = 0.d0
                do k=1, 5
                    if (ind_d(k) == 0) then
                        j = j + 1
                        y(j) = tmp(k)
                    end if
                end do
                fjac(:,i) = (y - y0)/eps4jacobian
            end do
        end select

    end subroutine fcn
    
    
    
    subroutine fcn_nojac(n, x, fvec, iflag)
        integer n, iflag
        double precision x(n), fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
        integer             ::  i, j
        real(kind=8)        ::  Dp(3,3), Ssolved(3,3), xi(n), tmp(5), y(n), d(5), vm

        stress_only = .TRUE.
        
        d = 0.d0
        forall (i=1:5, ind_d(i) == 1)  d(i) = d_prescribed(i)
        j = 0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                d(i) = x(j)
            end if
        end do
        
        d = rescaleD(d, ind_d)
        
        Dp = voigt2m(d)
        Dp = Dp/VonMises(Dp,'strain')
        
        if (grain_interaction == 1) then
            call taylor_step(Dp, Ssolved)
        elseif ((grain_interaction == 2) .or. (grain_interaction == 3)) then
            call alamel_step(Dp, Ssolved)
        end if
        ! renormalize stress elemetns but only those related to stress direction
        call vonMises_dvec(m2voigt(Ssolved), 'stress', ind_sdir, vm, tmp)
        ! slice only those n with either ind_sdir or ind_sabs == TRUE     
        j = 0
        y = 0.d0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                y(j) = tmp(i)
            end if
        end do
        fvec = y - rhs
        
        j = 0
        do i=1, 5
            if (ind_d(i) == 0) then
                j = j + 1
                if (ind_sabs(i) == 1) then
                    fvec(j) = fvec(j)/crss_mean
                end if
            end if
        end do
		
    end subroutine fcn_nojac
    
    
    
    subroutine fcn_scyl_tay(n, x, fvec, fjac, ldfjac, iflag)
        integer n, ldfjac, iflag
        double precision x(n), fvec(n), fjac(ldfjac,n)

    !----------
    !if iflag = 1 calculate the functions at x and
    !return this vector in fvec. do not alter fjac.
    !if iflag = 2 calculate the jacobian at x and
    !return this matrix in fjac. do not alter fvec.
    !---------
    
        integer             ::  tmpSCYLexp
        real(kind=8)        ::  ddphi(5,5), dphi(3,3), dphivec(5), S(3,3), phi, lamdot
        
        S = voigt2m(x(:5))
        lamdot = x(6)
            
        call SCYLphi(S, SCYL_CRSS, phi)
        call SCYLdphi(S, SCYL_CRSS, dphi)
        
        dphivec = m2voigt(dphi)
        
        select case(iflag)
        case(1)
            fvec(:5) = lamdot*dphivec - SCYL_Dp
            fvec(6)  = phi
        case(2)
            call SCYLddphi(S, SCYL_CRSS, ddphi)
            fjac(:5,:5) = ddphi
            fjac(6,:5)  = dphivec
            fjac(:5,6)  = dphivec
            fjac(6,6)   = 0.d0
        end select

    end subroutine fcn_scyl_tay
    
        
    subroutine fcn_scyl_ala(n, x, fvec, iflag)
        integer n, iflag
        double precision x(n), fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
        real(kind=8)            ::  Dp(3,3), Dpvec(5), S1(3,3), S2(3,3), phi,           &
                                    SCYLexp_orig, errsig, DpG1(3,3), DpG2(3,3)
        
        ! grain 1
        Dpvec = SCYL_Dp
        Dpvec(4) = Dpvec(4) - x(1)
        Dpvec(3) = Dpvec(3) - x(2)
        !
        Dp = voigt2m(Dpvec) 
        DpG1  = transform(Dp, transp3(SCYL_QbQg1T))
        
        S1 = DpG1 ! Von Mises initial guess
        call SCYLphi(S1, SCYL_CRSSg1, phi)
        S1 = 1.d0/(phi + 1.d0)*S1
        
        errsig = 1.d-9        
        if (SCYLexp > 20) then
            SCYLexp_orig = SCYLexp
            SCYLexp = 20.d0
            call findstress(S1, DpG1, SCYL_CRSSg1, errsig)
            SCYLexp = SCYLexp_orig
            errsig = 1.d-9
        end if
        call findstress(S1, DpG1, SCYL_CRSSg1, errsig)
        
        ! grain 2
        Dpvec = SCYL_Dp
        Dpvec(4) = Dpvec(4) + x(1)
        Dpvec(3) = Dpvec(3) + x(2)
        !
        Dp = voigt2m(Dpvec)
        DpG2  = transform(Dp, transp3(SCYL_QbQg2T))
        
        S2 = DpG2 ! Von Mises initial guess
        call SCYLphi(S2, SCYL_CRSSg2, phi)
        S2 = 1.d0/(phi + 1.d0)*S2
        
        errsig = 1.d-9        
        if (SCYLexp > 20) then
            SCYLexp_orig = SCYLexp
            SCYLexp = 20.d0
            call findstress(S2, DpG2, SCYL_CRSSg2, errsig)
            SCYLexp = SCYLexp_orig
            errsig = 1.d-9
        end if
        call findstress(S2, DpG2, SCYL_CRSSg2, errsig)
        
        ! stress in grain's CS
        SCYL_S1solved = S1
        SCYL_S2solved = S2
        
        ! stress in grain boundary's CS
        S1  = transform(S1, SCYL_QbQg1T)
        S2  = transform(S2, SCYL_QbQg2T)
        
        fvec = 0.d0
        fvec(1) = (S1(2,3) - S2(2,3))/crss_mean
        fvec(2) = (S1(1,3) - S2(1,3))/crss_mean

        
    end subroutine fcn_scyl_ala 
        
    
    ! subroutine readori_taylor(euangf, ori, GNUM)
    
        ! integer             ::  FID, GNUM, I, J, ALIVE
        ! real(kind=8)        ::  TT, FF, EANGAR(5,Ngrainmax)
        ! real(kind=8), allocatable   ::  ori(:,:,:)
        ! character(20)       ::  STATUS, dummy, euangf
        
        ! intent(in)          ::  euangf
        ! intent(out)         ::  ori
        
            ! !--**	Read euler angles    **-- 
        ! FID=1112
        ! OPEN(FID,FILE=euangf,STATUS='OLD',IOSTAT=ALIVE)

        ! READ(FID,*) dummy
        ! READ(FID,*) dummy
        ! READ(FID,*) GNUM, TT, FF    
        
        ! allocate(ori(3,3,GNUM))
        
        ! DO I=1,GNUM
            ! READ(FID,*) (EANGAR(J,I),J=1,5)
	        ! IF(DABS(TT)>1.0D-12) THEN
	            ! EANGAR(4,I) = TT
	        ! ENDIF
        ! ENDDO
        ! CLOSE(FID)
                   
        ! DO I=1,GNUM
            ! ori(:,:,I) = euler2rotm(EANGAR(1,I), EANGAR(2,I), EANGAR(3,I))
        ! ENDDO

    ! end subroutine readori_taylor
    
    
    
    ! subroutine readori_alamel(euangf, orig1, orig2, origb, GNUM)
    
        ! integer             ::  FID, GNUM, I, J, ALIVE
        ! real(kind=8)        ::  TT, FF, EANGAR(5,Ngrainmax)
        ! real(kind=8), allocatable   ::  orig1(:,:,:), orig2(:,:,:), origb(:,:,:)
        
        ! character(20)       ::  STATUS, dummy, euangf
        
        ! intent(in)          ::  euangf
        ! intent(out)         ::  orig1, orig2, origb
        
            ! !--**	Read euler angles    **-- 
        ! FID=1112
        ! OPEN(FID,FILE=euangf,STATUS='OLD',IOSTAT=ALIVE)

        ! READ(FID,*) dummy
        ! READ(FID,*) dummy
        ! READ(FID,*) GNUM, TT, FF    
        
        ! allocate(orig1(3,3,GNUM/2), orig2(3,3,GNUM/2), origb(3,3,GNUM/2))
        
        ! DO I=1,GNUM
            ! READ(FID,*) (EANGAR(J,I),J=1,5)
	        ! IF(DABS(TT)>1.0D-12) THEN
	            ! EANGAR(4,I) = TT
	        ! ENDIF
        ! ENDDO
        ! CLOSE(FID)
                   
        ! DO I=1,GNUM/2
            ! orig1(:,:,I) = euler2rotm(EANGAR(1,2*I-1), EANGAR(2,2*I-1), EANGAR(3,2*I-1))
            ! orig2(:,:,I) = euler2rotm(EANGAR(1,2*I), EANGAR(2,2*I), EANGAR(3,2*I))
            ! origb(:,:,I)  = eye(3)
        ! ENDDO
        ! GNUM = GNUM/2
    ! end subroutine readori_alamel

    
end module crystal_plasticity