! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

     module run_star_extras

       use star_lib
       use star_def
       use const_def
       use math_lib


       implicit none

       !real(dp) :: t_spindown
       real(dp), save :: ages(5), Teffs(5), Prots(5), total_AMs(5)
       real(dp) :: dJdt_hist
       integer :: target_age_ID

       ! these routines are called by the standard run_star check_model
       contains

         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! this is the place to set any procedure pointers you want to change
            ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


            ! the extras functions in this file will not be called
            ! unless you set their function pointers as done below.
            ! otherwise we use a null_ version which does nothing (except warn).

            s% extras_startup => extras_startup
            s% extras_start_step => extras_start_step
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns

            s% how_many_extra_history_header_items => how_many_extra_history_header_items
            s% data_for_extra_history_header_items => data_for_extra_history_header_items
            s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
            s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

            s% other_torque => torque_magnetic_braking

            target_age_ID = 1

         end subroutine extras_controls


         subroutine extras_startup(id, restart, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ages      = 0._dp
            Teffs     = 0._dp
            Prots     = 0._dp
            total_AMs = 0._dp


         end subroutine extras_startup


         integer function extras_start_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_start_step = 0


         end function extras_start_step


         ! returns either keep_going, retry, or terminate.
         integer function extras_check_model(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            logical :: do_retry

            real(dp) :: target_age(6), precision

            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_check_model = keep_going
            if (.false. .and. s% star_mass_h1 < 0.35d0) then
               ! stop when star hydrogen mass drops to specified level
               extras_check_model = terminate
               write(*, *) 'have reached desired hydrogen mass'
               return
            end if



            do_retry = .false.
            ! if you want to check multiple conditions, it can be useful
            ! to set a different termination code depending on which
            ! condition was triggered.  MESA provides 9 customizeable
            ! termination codes, named t_xtra1 .. t_xtra9.  You can
            ! customize the messages that will be printed upon exit by
            ! setting the corresponding termination_code_str value.
            ! termination_code_str(t_xtra1) = 'my termination condition'

            precision = 0.01
            target_age(1) = 0.3d9
            target_age(2) = 0.5d9
            target_age(3) = 1.3d9
            target_age(4) = 2.0d9
            target_age(5) = 3.5d9
            target_age(6) = 1d99



            call step_at_age(id, target_age(target_age_ID), precision, do_retry, ierr)
            if (do_retry) extras_check_model = retry


            if (s% omega_avg_surf / s% omega_crit_avg_surf > 1d0) then
               extras_check_model = terminate
               write(*, *) 'have reached critical rotation'
               return
            end if
            ! by default, indicate where (in the code) MESA terminated
            if (extras_check_model == terminate) s% termination_code = t_extras_check_model
         end function extras_check_model


         integer function how_many_extra_history_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_columns = 1
         end function how_many_extra_history_columns


         subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(n)
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           names(1) = 'dJdt'
           vals(1) = dJdt_hist
         end subroutine data_for_extra_history_columns


         integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_columns = 0
         end function how_many_extra_profile_columns


         subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k

            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! note: do NOT add the extra names to profile_columns.list
            ! the profile_columns.list is only for the built-in profile column options.
            ! it must not include the new column names you are adding here.

            ! here is an example for adding a profile column
            !if (n /= 1) stop 'data_for_extra_profile_columns'
            !names(1) = 'beta'
            !do k = 1, nz
            !   vals(k,1) = s% Pgas(k)/s% P(k)
            !end do


         end subroutine data_for_extra_profile_columns


         integer function how_many_extra_history_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_header_items = 0
         end function how_many_extra_history_header_items


         subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra history header item
            ! also set how_many_extra_history_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha

         end subroutine data_for_extra_history_header_items


         integer function how_many_extra_profile_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_header_items = 0
         end function how_many_extra_profile_header_items


         subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra profile header item
            ! also set how_many_extra_profile_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha

         end subroutine data_for_extra_profile_header_items


         ! returns either keep_going or terminate.
         ! note: cannot request retry; extras_check_model can do that.
         integer function extras_finish_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_finish_step = keep_going

            ! to save a profile,
               ! s% need_to_save_profiles_now = .true.
            ! to update the star log,
               ! s% need_to_update_history_now = .true.


            ! see extras_check_model for information about custom termination codes
            ! by default, indicate where (in the code) MESA terminated

            !if ((t_spindown / s% dt) .lt. 10) then
            !    s% dt_next = s% dt * 0.5d0
            !    write(*,*) "Warning: Torque too large. Decreasing timestep to ",s% dt_next
            !end if

            if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step



         end function extras_finish_step


         subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k
            ierr = 0
            call star_ptr(id, s, ierr)

            if (ierr /= 0) return

            write(*,*) '***********************************'
            do k =1,5
              write(*,*) 'age [Gyr]', ages(k), 'Teff [K]', Teffs(k), 'P_rot [d]', Prots(k), 'log(J_tot)', total_AMs(k)
            end do
            write(*,*) '***********************************'

         end subroutine extras_after_evolve


         subroutine torque_magnetic_braking(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr

            real(dp) :: jdot_total, wsum, omega_e, omega_crit
            real(dp), allocatable :: w(:)
            real(dp) :: P_phot, tcz_bot_Hp, fk, KM, omega_sat, omega_sun, omega_threshold
            real(dp) :: Ro, Ro_crit, P_phot_sun
            real(dp) :: tcz_bot_Hp_sun, KM_sun
            real(dp) :: mixing_legnth_ocz_bot_Hp, conv_vel_ocz_bot_Hp, r_ocz_bot, m_ocz_bot, r_ocz_bot_Hp
            real(dp) :: subcell_ocz_bot, subcell_ocz_bot_Hp, subcell_ocz_top, subcell_phot
            real(dp) :: sch_k0, sch_k1, sch_unstable(10000), turnover_time_k
            integer :: k, k_ocz_bot, k_ocz_bot_Hp, k_ocz_top, k_phot
            integer :: zonefound, nunstable

            type (star_info), pointer :: s
            include 'formats'
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            if (s% doing_relax) return

            ! ! note that can set extra_omegadot instead of extra_jdot if that is more convenient
            ! ! set one or the other, not both.  set the one you are not using to 0 as in the following line.
            ! s% extra_jdot(:) = 0
            ! s% extra_omegadot(:) = 0
            omega_crit = star_surface_omega_crit(id, ierr) ! this forces a call to set_surf_avg_rotation_info to ensure things are up
                                                           ! to date with the state

            ! Initialize variables
            s% extra_jdot(:) = 0d0
            jdot_total = 0d0
            omega_e = s% omega(1)

            if (omega_e < 0._dp) write(*,*) 'Warning: Negative surface rotation.'

            tcz_bot_Hp_sun = 717053.416 ! seconds, pre-calibrated
            ! KM_sun = 3.9715223932d35 ! pre-calibrated
            KM_sun = 4.9715223932d35 ! pre-calibrated
            P_phot_sun = 100179.864 ! pre-calibrated
            fk = 5.655 ! pre-calibrated
            omega_sat = 3.863d-5 ! rad/s
            omega_sun = 2.863d-6 ! rad/s

            ! Calculate jdot_total using a simple magnetic braking law
            ! e.g. Kawaler (1988) with n=1.5 and a solar-calibrated K.  See also Matt et al. (2015).
            ! I will first try van Saders & Pinsonneault (2013) which is a modified Kawaler law

            ! find location where tau=2/3
            do k = 2, s% nz
              if ((s% tau(k) > 2._dp/3._dp)) then
                k_phot = k
                subcell_phot = ((2._dp/3._dp) - s% tau(k_phot-1)) / (s% tau(k_phot) - s% tau(k_phot-1))
                exit
              end if
            end do

            ! calculate photospheric pressure at tau=2/3
            P_phot = s% Peos(k_phot-1) + (s% Peos(k_phot-1) - s% Peos(k_phot)) * subcell_phot

            ! find Schwarzschild unstable zones
            nunstable = 0
            do k = 1, s% nz - 1
              if ((s% grada_face(k) - s% gradr(k)) < 0) then
                nunstable = nunstable + 1
                sch_unstable(nunstable) = k
              end if
            end do

            ! find the bottom of the outer convection zone
            ! if the star is fully convective, choose the core
            zonefound = 0
            k_ocz_bot = 2
            if (s% mass_conv_core == s% star_mass) then
                k_ocz_bot = s% nz - 1
            ! otherwise, find the last unstable shell followed by a gap of >100 shells
            else
               do k = 1, nunstable
                  if ((sch_unstable(k) > 10) .and. ((sch_unstable(k+1) - sch_unstable(k)) > 100) .and. (zonefound == 0)) then
                     k_ocz_bot = sch_unstable(k)
                     zonefound = 1
                  end if
               end do

               ! if no gap (e.g. no convective core) take the last unstable shell
               if (zonefound == 0) then
                  k_ocz_bot = MAXVAL(sch_unstable)
               end if

               ! correct for subshell
               sch_k0 = s% grada_face(k_ocz_bot) - s% gradr(k_ocz_bot)
               sch_k1 = s% grada_face(k_ocz_bot+1) - s% gradr(k_ocz_bot+1)
               subcell_ocz_bot = - sch_k0 / (sch_k1 - sch_k0)

            end if

            ! find the top of the outer convection zone
            k_ocz_top = 2
            subcell_ocz_top = 0.
            do k = k_ocz_bot, 1, -1
               if ((s% grada_face(k) - s% gradr(k)) > 0) then
                  k_ocz_top = k

                  ! correct for subshell
                  ! sch_k0 = s% grada_face(k_ocz_top-1) - s% gradr(k_ocz_top-1)
                  ! sch_k1 = s% grada_face(k_ocz_top) - s% gradr(k_ocz_top)
                  ! subcell_ocz_top = sch_k1 / (sch_k1 - sch_k0)
                  sch_k0 = s% grada_face(k_ocz_top) - s% gradr(k_ocz_top)
                  sch_k1 = s% grada_face(k_ocz_top+1) - s% gradr(k_ocz_top+1)
                  subcell_ocz_top = - sch_k0 / (sch_k1 - sch_k0)

                  exit
               endif
            enddo

            ! find the location r_bcz_Hp := 1 pressure scale height above the base of convective zone
            tcz_bot_Hp = 0.
            r_ocz_bot = s% r(k_ocz_bot) + (s% r(k_ocz_bot+1) - s% r(k_ocz_bot)) * subcell_ocz_bot
            m_ocz_bot = s% m(k_ocz_bot) + (s% m(k_ocz_bot+1) - s% m(k_ocz_bot)) * subcell_ocz_bot
            r_ocz_bot_Hp = s% r(k_ocz_bot) + (s% r(k_ocz_bot+1) - s% r(k_ocz_bot)) * subcell_ocz_bot + &
                     (s% mlt_mixing_length(k_ocz_bot) + (s% mlt_mixing_length(k_ocz_bot+1) - s% mlt_mixing_length(k_ocz_bot)) * subcell_ocz_bot )
            do k = k_ocz_bot, k_ocz_top, -1
               if (s% r(k) > r_ocz_bot_Hp) then
                     k_ocz_bot_hp = k
                     subcell_ocz_bot_Hp = (r_ocz_bot_Hp - s% rmid(k)) / (s% rmid(k+1) - s% rmid(k) )

                     conv_vel_ocz_bot_Hp = s% conv_vel(k) + (s% conv_vel(k+1) - s% conv_vel(k)) * subcell_ocz_bot_Hp
                     mixing_legnth_ocz_bot_Hp = s% mlt_mixing_length(k) + (s% mlt_mixing_length(k+1) - s% mlt_mixing_length(k)) * subcell_ocz_bot_Hp

                     tcz_bot_Hp = mixing_legnth_ocz_bot_Hp / conv_vel_ocz_bot_Hp
                     exit
               end if
            end do

            KM = KM_sun * &
               pow(s% photosphere_r, 3.1d0)* &
               pow(s% star_mass, -0.22d0) * &
               pow(s% photosphere_L, 0.56d0) * &
               pow(P_phot / P_phot_sun, 0.44d0)

            ! Rotation period and Rossby number
            Ro = (2._dp*pi / omega_e) / tcz_bot_Hp
            Ro_crit = (2._dp*pi / omega_sun) / tcz_bot_Hp_sun * 1._dp

            ! Saturation boundary:
            ! omega_sat <= omega_e * tcz_bot_Hp / tcz_bot_Hp_sun  => saturated branch
            if ((Ro <= Ro_crit) .and. (omega_sat <= omega_e * tcz_bot_Hp / tcz_bot_Hp_sun)) then
               write(*, *) "Saturation"
               jdot_total = - fK * KM * omega_e * pow2(omega_sat / omega_sun)
            end if
            ! Unsaturated branch
            if ((Ro <= Ro_crit) .and. (omega_sat > omega_e * tcz_bot_Hp / tcz_bot_Hp_sun)) then
               write(*, *) "Unsaturated"
               jdot_total = - fK * KM * omega_e * pow2((omega_e * tcz_bot_Hp) / (omega_sun * tcz_bot_Hp_sun))
            end if
            ! No braking above critical Rossby number
            if (Ro > Ro_crit) then
               write(*, *) "No braking"
               jdot_total = 0d0
            end if

            write(*, *) "The total jdot is ", jdot_total, " cgs units. Fractional AM loss w.r.t total AM is", (jdot_total*s% dt)/s% total_angular_momentum

            dJdt_hist = jdot_total
            !--------------------------------------------
            ! Moment-of-inertia-like weights
            ! w(k) ~ dm * r^2
            !--------------------------------------------
            allocate(w(s% nz))
            w(:) = 0d0
            do k = 1, s% nz
               w(k) = s% dm_bar(k) * s% rmid(k) * s% rmid(k)
            end do

            wsum = sum(w(:))
            if (wsum <= 0d0) then
               ierr = 1
               deallocate(w)
               return
            end if

            !--------------------------------------------
            ! Distribute total torque to maintain rigid rotation
            ! extra_jdot(k) = shell torque / shell mass
            !--------------------------------------------
            do k = 1, s% nz
               s% extra_jdot(k) = (jdot_total * w(k) / wsum) / s% dm_bar(k)
            end do

         end subroutine torque_magnetic_braking

         subroutine step_at_age(id, target_age, precision, do_retry, ierr)
             integer, intent(in) :: id
             real(dp), intent(inout) :: target_age
             real(dp), intent(in) :: precision
             logical, intent(out) :: do_retry
             integer, intent(out) :: ierr
             logical :: ok_time_step
             type(star_info), pointer :: s

             call star_ptr(id, s, ierr)
             if (ierr /= 0) then
                 write(*,*) 'Error: grid_io: set_logD_mix: star_ptr failed'
                 return
             endif

             do_retry  = .false.
             ok_time_step = .false.
             if (ABS(s% star_age - target_age)/target_age < precision) then
                 ok_time_step = .true.
             endif

             ! retry if target Xc was missed
             if (s% star_age > target_age .and. .not. ok_time_step) then
                 s%dt = 0.5*s%dt
                 do_retry = .true.
                 write(*,*) 'Reducing time step for specific age.'
                 return
             endif

             ! Adjust Xc to save if model is going to be saved
             if (ok_time_step .and. .not. do_retry) then
                 target_age_ID = target_age_ID + 1
             endif

             if ((target_age == 3.5d9) .and. (ages(5) == 0._dp) .and. ok_time_step) then
               ages(5)      = s% star_age/1d9
               Teffs(5)     = s% Teff
               Prots(5)     = 1._dp/(s% omega_avg_surf*86400/(2._dp*pi))
               total_AMs(5) = safe_log10(s% total_angular_momentum)

             else if ((target_age == 2.0d9) .and. (ages(4) == 0._dp) .and. ok_time_step) then
               ages(4)      = s% star_age/1d9
               Teffs(4)     = s% Teff
               Prots(4)     = 1._dp/(s% omega_avg_surf*86400/(2._dp*pi))
               total_AMs(4) = safe_log10(s% total_angular_momentum)

             else if ((target_age == 1.3d9) .and. (ages(3) == 0._dp) .and. ok_time_step) then
               ages(3)      = s% star_age/1d9
               Teffs(3)     = s% Teff
               Prots(3)     = 1._dp/(s% omega_avg_surf*86400/(2._dp*pi))
               total_AMs(3) = safe_log10(s% total_angular_momentum)

             else if ((target_age == 0.5d9) .and. (ages(2) == 0._dp) .and. ok_time_step) then
               ages(2)      = s% star_age/1d9
               Teffs(2)     = s% Teff
               Prots(2)     = 1._dp/(s% omega_avg_surf*86400/(2._dp*pi))
               total_AMs(2) = safe_log10(s% total_angular_momentum)

             else if ((target_age == 0.3d9) .and. (ages(1) == 0._dp) .and. ok_time_step) then
               ages(1)      = s% star_age/1d9
               Teffs(1)     = s% Teff
               Prots(1)     = 1._dp/(s% omega_avg_surf*86400/(2._dp*pi))
               total_AMs(1) = safe_log10(s% total_angular_momentum)

             end if

         end subroutine step_at_age

end module run_star_extras
