
! ForCLaP, the Fortran Collinear Lagrange Points Calculator.
!
! This program calculates the positions of the collinear Lagrange points in the circular restricted three-body problem.
!
! * Radiation pressure is included.
! * Newton-Raphson convergence components are displayed.
! * Selected grouping of terms by sign and magnitude was done to help avoid round-off errors.
!
! ForCLaP contains Modern Fortran updates to the Fortran 95 program described in:
!
! Stenborg, TN 2008, "Collinear Lagrange Point Solutions in the Circular Restricted Three-Body Problem with Radiation
! Pressure using Fortran", in RW Argyle, PS Bunclark & JR Lewis (eds), Astronomical Data Analysis Software and Systems
! XVII, Astronomical Society of the Pacific Conference Series, vol 394, pp. 734-737.


module module_forclap

    ! Define how double precision (dp) is implemented here.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    ! Enforce explicit variable declaration.
    implicit none

contains
    logical function divergence_detected(fr, fr_deriv, point_name, r)

        ! Test for runaway divergence, by checking for NaN, a value not equal to itself.

        ! Input parameters.
        real(dp), intent(in) :: r, fr, fr_deriv
        character(len = 2), intent(in) :: point_name

        if ((r /= r) .OR. (fr /= fr) .OR. (fr_deriv /= fr_deriv)) then
            write(*, *) "Runaway Newton-Raphson divergence detected for these parameters."
            write(*, *) point_name, " calculations cancelled."
            divergence_detected = .TRUE.
        else
            divergence_detected = .FALSE.
        end if

    end function

    logical function precision_reached(fr_deriv, fr_deriv_old, r, r_old)

        ! Test if r and fr_deriv are both precise to five decimal places.

        ! Input parameters.
        real(dp), intent(in) :: r, r_old, fr_deriv, fr_deriv_old

        precision_reached = (abs(1.0_dp - abs(r / r_old)) <= 0.0001_dp) &
                             .AND. (abs(1.0_dp - abs(fr_deriv / fr_deriv_old)) <= 0.0001_dp)

    end function

    subroutine newton_raphson_step(fr, fr_deriv, r, r_old, i)

        ! Input parameters.
        real(dp), intent(in) :: fr, fr_deriv

        ! Input/Output parameters.
        real(dp), intent(inout) :: r
        integer, intent(inout) :: i

        ! Output parameters.
        real(dp), intent(out) :: r_old

        ! Compute a new r at i + 1.
        r_old = r
        r = r - fr / fr_deriv
        i = i + 1

    end subroutine

end module


program forclap

    ! Include the local module.
    use module_forclap

    ! Enforce explicit variable declaration.
    implicit none

    ! Set the Routh stability value, ~0.03852089650455.
    real(dp), parameter :: routh_stability = (9.0_dp - sqrt(69.0_dp)) / 18.0_dp

    ! Variable declarations.
    real(dp) :: beta, mass_ratio, r, r_old
    real(dp) :: fr, fr_deriv, fr_deriv_old
    real(dp) :: fr_deriv_numerator, fr_deriv_denominator
    integer :: i, ios

    write(*, *) "ForCLaP, the Fortran Collinear Lagrange Points Calculator."

    ! Get a mass ratio value, ensuring it's numeric.
    write(*, *) "Enter desired secondary/primary mass ratio (e.g., 0.000955 for the Jupiter/Sun system):"
    loop_mass_ratio: do
        read(*, *, iostat = ios) mass_ratio
        if (ios /= 0) then
            write(*, *) "Please enter a numeric value:"
        else if (mass_ratio < 0.0_dp) then
            write(*, *) "Please enter a positive numeric value:"
        else if (mass_ratio >= routh_stability) then
            write(*, *) "Please enter a ratio less than the Routh stability value ((9 - sqrt(69)) / 18):"
        else
            exit loop_mass_ratio
        end if
    end do loop_mass_ratio

    ! Get a beta value, ensuring it's numeric.
    write(*, *) "Enter desired ratio of solar radiation pressure to gravitational force, beta (e.g., 0.2):"
    loop_beta: do
        read(*, *, iostat = ios) beta
        if (ios /= 0) then
            write(*, *) "Please enter a numeric value:"
        else if (beta < 0.0_dp) then
            write(*, *) "Please enter a positive numeric value:"
        else
            exit loop_beta
        end if
    end do loop_beta


    write(*, "(/,A)", advance = "no") "Newton-Raphson convergence on r2 for L1, [ mass ratio : beta ] = [ "
    write(*, "(E19.13, A, E19.13, A)") mass_ratio, " : ", beta, " ]"
    write(*, "(A6, 3X, A20, 3X, A20, 3X, A20)") "i", "r2", "f(r2)", "f'(r2)"

    ! Initialise L1 Newton-Raphson parameters.
    i = 0
    r = 0.07_dp
    fr_deriv = 0.0_dp

    loop_newton_L1: do

        ! Compute a new f(r).
        fr = ((r ** 5 + 3.0_dp * r ** 3 - (3.0_dp * r ** 4 + beta * r ** 2)) &
              / (2.0_dp * r ** 4 + r ** 2 + 1.0_dp - (r ** 5 + r ** 3 + 2.0_dp * r))) - mass_ratio

        ! Compute a new f'(r).
        fr_deriv_numerator = 4.0_dp * r ** 7 + (26.0_dp - beta) * r ** 4 + (9.0_dp - 2.0_dp * beta) * r ** 2 &
                             - (r ** 8 + 3.0_dp * beta * r ** 6 + (14.0_dp - 4.0_dp * beta) * r ** 5 &
                                + 24.0_dp * r ** 3 + 2.0_dp * beta * r)
        fr_deriv_denominator = r ** 10 + 6.0_dp * r ** 8 + 9.0_dp * r ** 6 + 9.0_dp * r ** 4 + 6.0_dp * r ** 2 &
                               + 1.0_dp - (4.0_dp * r ** 9 + 6.0_dp * r ** 7 + 12.0_dp * r ** 5 + 6.0_dp * r ** 3 &
                                           + 4.0_dp * r)
        fr_deriv_old = fr_deriv
        fr_deriv = fr_deriv_numerator / fr_deriv_denominator

        write(*, "(I6, 3X, E20.13, 3X, E20.13, 3X, E20.13)") i, r, fr, fr_deriv

        call newton_raphson_step ( fr, fr_deriv, r, r_old, i )

        if (divergence_detected(fr, fr_deriv, "L1", r)) exit loop_newton_L1
        if (precision_reached(fr_deriv, fr_deriv_old, r, r_old)) exit loop_newton_L1

    end do loop_newton_L1


    write(*, "(/,A)", advance = "no") "Newton-Raphson convergence on r2 for L2, [ mass ratio : beta ] = [ "
    write(*, "(E19.13, A, E19.13, A)") mass_ratio, " : ", beta, " ]"
    write(*, "(A6, 3X, A20, 3X, A20, 3X, A20)") "i", "r2", "f(r2)", "f'(r2)"

    ! Initialise L2 Newton-Raphson parameters.
    i = 0
    r = 0.07_dp
    fr_deriv = 0.0_dp

    loop_newton_L2: do

        ! Compute a new f(r).
        fr = ((r ** 5 + 3.0_dp * r ** 4 + 3.0_dp * r ** 3 + beta * r ** 2) &
              / (r ** 2 + 2.0_dp * r + 1.0_dp - (r ** 5 + 2.0_dp * r ** 4 + r ** 3))) - mass_ratio

        ! Compute a new f'(r).
        fr_deriv_numerator = r ** 8 + 4.0_dp * r ** 7 + (6.0_dp + 3.0_dp * beta) * r ** 6 &
                             + (14.0_dp + 4.0_dp * beta) * r ** 5 + (26.0_dp + beta) * r ** 4 &
                             + 24.0_dp * r ** 3 + (9.0_dp + 2.0_dp * beta) * r ** 2 + 2.0_dp * beta * r
        fr_deriv_denominator = r ** 10 + 4.0_dp * r ** 9 + 6.0_dp * r ** 8 + 2.0_dp * r ** 7 + 2.0_dp * r ** 3 &
                               + 6.0_dp * r ** 2 + 4.0_dp * r + 1.0_dp &
                               - (7.0_dp * r ** 6 + 12.0_dp * r ** 5 + 7.0_dp * r ** 4)
        fr_deriv_old = fr_deriv
        fr_deriv = fr_deriv_numerator / fr_deriv_denominator

        write(*, "(I6, 3X, E20.13, 3X, E20.13, 3X, E20.13)") i, r, fr, fr_deriv

        call newton_raphson_step ( fr, fr_deriv, r, r_old, i )

        if (divergence_detected(fr, fr_deriv, "L2", r)) exit loop_newton_L2
        if (precision_reached(fr_deriv, fr_deriv_old, r, r_old)) exit loop_newton_L2

    end do loop_newton_L2


    write(*, "(/,A)", advance = "no") "Newton-Raphson convergence on r2 for L3, [ mass ratio : beta ] = [ "
    write(*, "(E19.13, A, E19.13, A)") mass_ratio, " : ", beta, " ]"
    write(*, "(A6, 3X, A20, 3X, A20, 3X, A20)") "i", "r1", "f(r1)", "f'(r1)"

    ! Initialise L3 Newton-Raphson parameters.
    i = 0
    r = 1.0_dp
    fr_deriv = 0.0_dp

    loop_newton_L3: do

        ! Compute a new f(r).
        fr = (((1.0_dp - beta) * r ** 2 + (1.0_dp - beta) * 2.0_dp * r + 1.0_dp &
              - (r ** 5 + 2.0_dp * r ** 4 + r ** 3 + beta)) &
              / (r ** 5 + 3.0_dp * r ** 4 + 3.0_dp * r ** 3)) - mass_ratio

        ! Compute a new f'(r).
        fr_deriv_numerator = -(r ** 6 + 4.0_dp * r ** 5 + (2.0_dp - beta) * 3.0_dp * r ** 4 &
                               + (1.0_dp + beta) * 14.0_dp * r ** 3 &
                               + (1.0_dp - beta) * (26.0_dp * r ** 2 + 24.0_dp * r + 9.0_dp))
        fr_deriv_denominator = r ** 8 + 6.0_dp * r ** 7 + 15.0_dp * r ** 6 + 18.0_dp * r ** 5 + 9.0_dp * r ** 4
        fr_deriv_old = fr_deriv
        fr_deriv = fr_deriv_numerator / fr_deriv_denominator

        write(*, "(I6, 3X, E20.13, 3X, E20.13, 3X, E20.13)") i, r, fr, fr_deriv

        call newton_raphson_step ( fr, fr_deriv, r, r_old, i )

        if (divergence_detected(fr, fr_deriv, "L3", r)) exit loop_newton_L3
        if (precision_reached(fr_deriv, fr_deriv_old, r, r_old)) exit loop_newton_L3

    end do loop_newton_L3

    write(*, *) ""

end program forclap
