module py2f_yldfun
    implicit none
    
    integer :: NS
    real(8), allocatable :: S_list(:,:,:)
    real(8), allocatable :: Snat_list(:,:)
    
    contains
    
    
    subroutine py2f_stresses(pyS_all)
        integer :: s(2)
        real(8) :: pyS_all(:,:,:)
    
        intent(in) :: pyS_all
    
        s = 0.d0
        if (allocated(S_list)) deallocate(S_list, STAT=s(1))
        allocate(S_list, SOURCE=pyS_all, STAT=s(2))
    
        if (any(s /= 0)) then
            write(*,*) 'Deallocation/allocation not successful.'
            stop
        else
            NS = size(S_list, DIM=3)
        end if
    end subroutine py2f_stresses
    
    subroutine py2f_stresses_nat(pyS_all)
        integer :: s(2)
        real(8) :: pyS_all(:,:)
    
        intent(in) :: pyS_all
    
        s = 0.d0
        if (allocated(Snat_list)) deallocate(Snat_list, STAT=s(1))
        allocate(Snat_list, SOURCE=pyS_all, STAT=s(2))
    
        if (any(s /= 0)) then
            write(*,*) 'Deallocation/allocation not successful.'
            stop
        else
            NS = size(Snat_list, DIM=2)
        end if
    end subroutine py2f_stresses_nat
    
    
    subroutine get_residual(coef, model, a, NT, sigY, R)   
        integer :: i, NT
        real(8) :: svec(6), coef(200), a, sigY, R, f
        character(15) :: model
        
        intent(out) :: R
    
!       c12  = coef(1)
!       c13  = coef(2)
!       c21  = coef(3)
!       c23  = coef(4)
!       c31  = coef(5)
!       c32  = coef(6)
!       c44  = coef(7)
!       c55  = coef(8)
!       c66  = coef(9)

        R = 0.d0
        do i=1, NS
            svec = [S_list(1,1,i), &
                    S_list(2,2,i), &
                    S_list(3,3,i), &
                    S_list(2,3,i), &
                    S_list(1,3,i), &
                    S_list(1,2,i)]
                    
            select case (model)
			case('KB')
               call KB(svec, coef(:6), a, f)
            case('YLD2004_09p')
               call Yld2004_9p(svec, coef(:9), a, f)
            case('YLD2004_18p')
                call Yld2004_18p(svec, coef(:18), a, f)
            case('YLD2004_27p')
                call Yld2004_27p(svec, coef(:27), a, f)
            case('YLD2011_18p')
                call Yld2011_18p(svec, coef(:18), a, f)
            case('YLD2011_27p')
                call Yld2011_27p(svec, coef(:27), a, f)
            case('YLD2013_27p')
                call Yld2013_27p(svec, coef(:27), a, f)
            case('YLD2013_Xp')
                call Yld2013_Xp(svec, coef, a, NT, f)
			case('YLD2004NAT_18p')
			    call Yld2004nat_18p(svec, coef(:18), a, f)
            case('YLD2004NORT_21p')
                call Yld2004nort_21p(svec, coef(:21), a, f)
            end select        
            R = R + (f-sigY)**2
            
        end do 
    
    end subroutine get_residual

    
    subroutine get_residual_extra_weight(coef, model, a, NT, sigY, weight, R)   
        integer :: i, NT
        real(8) :: svec(6), coef(200), a, sigY, R, f, w, weight
        character(11) :: model
        
        intent(out) :: R
    
!       c12  = coef(1)
!       c13  = coef(2)
!       c21  = coef(3)
!       c23  = coef(4)
!       c31  = coef(5)
!       c32  = coef(6)
!       c44  = coef(7)
!       c55  = coef(8)
!       c66  = coef(9)

        R = 0.d0
        do i=1, NS
            svec = [S_list(1,1,i), &
                    S_list(2,2,i), &
                    S_list(3,3,i), &
                    S_list(2,3,i), &
                    S_list(1,3,i), &
                    S_list(1,2,i)]
                    
            select case (model)
            case('YLD2004_09p')
               call Yld2004_9p(svec, coef(:9), a, f)
            case('YLD2004_18p')
                call Yld2004_18p(svec, coef(:18), a, f)
            case('YLD2004_27p')
                call Yld2004_27p(svec, coef(:27), a, f)
            case('YLD2011_18p')
                call Yld2011_18p(svec, coef(:18), a, f)
            case('YLD2011_27p')
                call Yld2011_27p(svec, coef(:27), a, f)
            case('YLD2013_27p')
                call Yld2013_27p(svec, coef(:27), a, f)
            case('YLD2013_Xp')
                call Yld2013_Xp(svec, coef, a, NT, f)
            case('YLD2004NORT_21p')
                call Yld2004nort_21p(svec, coef(:21), a, f)
            end select    
            
            w = 1.d0
            if (i .EQ. 1) w = weight
            R = R + w*(f-sigY)**2
            
        end do 
    
    end subroutine get_residual_extra_weight
    
!    subroutine get_residual_FACET(facets, Nfacets, exps, N, s0s, pref, q, R)   
!        integer :: i, Nfacets(1000), N
!        real(8) :: svec(6), facets(6,1000), exps(1000), s0s(1000), pref(1000), q, R, f
!        
!        intent(out) :: R
!
!        R = 0.d0
!        do i=1, NS
!            svec = [S_list(1,1,i), &
!                    S_list(2,2,i), &
!                    S_list(3,3,i), &
!                    S_list(2,3,i), &
!                    S_list(1,3,i), &
!                    S_list(1,2,i)]
!                    
!            call facet(svec, facets, Nfacets, exps, N, s0s, pref, q, f)   
!            R = R + (f-1.d0)**2     
!        end do 
!    
!    end subroutine get_residual_FACET
!    
!    subroutine get_residual_FACETR(facets0, Nfacets0, exps0, N0, s0s0, &
!                                   pref, q, facets, exps, s0s, N, ortho, R)
!        logical :: ortho
!        integer :: i, N, Nfacets0(2), N0
!        real(8) :: svec(6), facets0(6,6), exps0(2), s0s0(2), pref(1000), &
!                   facets(6,1000), exps(1000), s0s(1000), q, R, f, nn(6)
!        
!        intent(out) :: R
!        
!        R = 0.d0
!        do i=1, NS
!            svec = [S_list(1,1,i), &
!                    S_list(2,2,i), &
!                    S_list(3,3,i), &
!                    S_list(2,3,i), &
!                    S_list(1,3,i), &
!                    S_list(1,2,i)]
!            
!            call facetR(svec, facets0, Nfacets0, exps0, N0, s0s0, pref, q, &
!                            facets, exps, s0s, N, ortho, f)
!            f = f - 1.d0
!            if (f .LE. 0.d0) then
!                R = R + f**2
!            else
!                R = R + 100.d0*f**1.5d0
!!                R = R + f**2
!            end if
!        end do 
!    
!    end subroutine get_residual_FACETR

                                   
    subroutine getJac(coef, Ncoef, fixcoefID, Nfixcoef, model, a, NT, S_list, Nstresses, D)
    
        ! input/output
        integer :: NT, Nstresses, Ncoef, fixcoefID(200), Nfixcoef
        real(8) :: coef(200), fixcoef(200), a, D(Ncoef-Nfixcoef, Nstresses), S_list(3,3,10000)
        character(11) :: model
        
        integer :: i, j, m, n
        real(8) :: f, k, svec(6), dc, c(200)
        logical :: isfixed
        
        intent(out) :: D
        
        dc = 1.d-6
        
        m = 0
        do i=1, Ncoef
            isfixed = .False.
            do n=1, Nfixcoef
                if (fixcoefID(n) .EQ. i) then
                    isfixed = .True.
                    exit
                end if
            end do
            if (isfixed) then
                cycle
            else
                m = m + 1
            end if
            
            c = coef
            c(i) = c(i) + dc
            do j=1, Nstresses
                svec = [S_list(1,1,j), &
                        S_list(2,2,j), &
                        S_list(3,3,j), &
                        S_list(2,3,j), &
                        S_list(1,3,j), &
                        S_list(1,2,j)]
                select case (model)
                case('YLD2004_09p')
                   Ncoef = 9
                   call Yld2004_9p(svec, c(:Ncoef), a, f)
                case('YLD2004_18p')
                    Ncoef = 18
                    call Yld2004_18p(svec, c(:Ncoef), a, f)
                case('YLD2004_27p')
                    Ncoef = 27
                    call Yld2004_27p(svec, c(:Ncoef), a, f)
                case('YLD2011_18p')
                    Ncoef = 18
                    call Yld2011_18p(svec, c(:Ncoef), a, f)
                case('YLD2011_27p')
                    Ncoef = 27
                    call Yld2011_27p(svec, c(:Ncoef), a, f)
                case('YLD2013_27p')
                    Ncoef = 27
                    call Yld2013_27p(svec, c(:Ncoef), a, f)
                case('YLD2013_Xp')
                    Ncoef = 9*NT
                    call Yld2013_Xp(svec, c, a, NT, f)
                case('YLD2004NORT_21p')
                    Ncoef = 21
                    call Yld2004nort_21p(svec, c(:Ncoef), a, f)
                end select
                
                k = 1.d0/f
                svec = k*svec
                
                select case (model)
                case('YLD2004_09p')
                   call Yld2004_9p(svec, coef(:Ncoef), a, f)
                case('YLD2004_18p')
                    call Yld2004_18p(svec, coef(:Ncoef), a, f)
                case('YLD2004_27p')
                    call Yld2004_27p(svec, coef(:Ncoef), a, f)
                case('YLD2011_18p')
                    call Yld2011_18p(svec, coef(:Ncoef), a, f)
                case('YLD2011_27p')
                    call Yld2011_27p(svec, coef(:Ncoef), a, f)
                case('YLD2013_27p')
                    call Yld2013_27p(svec, coef(:Ncoef), a, f)
                case('YLD2013_Xp')
                    call Yld2013_Xp(svec, coef, a, NT, f)
                case('YLD2004NORT_21p')
                    call Yld2004nort_21p(svec, coef(:Ncoef), a, f)
                end select
                
                D(m,j) = (f - 1.d0)/dc
            end do          
        end do
        
    end subroutine getJac
    
    
end module py2f_yldfun