      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      ! self-defined variables for the whole module to use
      ! character(len=6) :: model_number
      ! character(len=80) :: summary_filename
      ! integer :: summary_unit
      ! logical :: flag_gyre, flag_diffusion, flag_mass_loss, zams_model_saved

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
         
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
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

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         ! if ( (s% center_h1 < 0.01d0) ) then
         !    ! stop when center H1 mass drops to specified level
         !    extras_check_model = terminate
         !    write(*, *) 'have reached terminal main sequence.'
         !    return
         ! end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      ! subroutine default_other_torque(id, ierr)
      !    integer, intent(in) :: id
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    integer :: k
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
      !    ! note that can set extra_omegadot instead of extra_jdot if that is more convenient
      !    ! set one or the other, not both.  set the one you are not using to 0 as in the following line.
      !    s% extra_jdot(:) = 0
      !    s% extra_omegadot(:) = 0
      ! end subroutine default_other_torque


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
           if ((s% tau(k) > 6.66666666666666E-1)) then
             k_phot = k
             subcell_phot = ((6.66666666666666E-1) - s% tau(k_phot-1)) / (s% tau(k_phot) - s% tau(k_phot-1))
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

         KM = KM_sun * 1 * &
            (s% photosphere_r)**3.1d0 * &
            (s% star_mass)**(-0.22d0) * &
            (s% photosphere_L)**0.56d0 * &
            (P_phot / P_phot_sun)**0.44d0

         ! Rotation period and Rossby number
         Ro = (2d0*3.14159 / omega_e) / tcz_bot_Hp
         Ro_crit = (2d0*3.14159 / omega_sun) / tcz_bot_Hp_sun * 1.0

         ! Saturation boundary:
         ! omega_sat <= omega_e * tcz_bot_Hp / tcz_bot_Hp_sun  => saturated branch
         if ((Ro <= Ro_crit) .and. (omega_sat <= omega_e * tcz_bot_Hp / tcz_bot_Hp_sun)) then
            write(*, *) "Saturation"
            jdot_total = - fK * KM * omega_e * (omega_sat / omega_sun)**(2.0d0)
         end if
         ! Unsaturated branch
         if ((Ro <= Ro_crit) .and. (omega_sat > omega_e * tcz_bot_Hp / tcz_bot_Hp_sun)) then
            write(*, *) "Unsaturated"
            jdot_total = - fK * KM * omega_e * ((omega_e * tcz_bot_Hp) / (omega_sun * tcz_bot_Hp_sun))**(2.0d0)
         end if
         ! No braking above critical Rossby number
         if (Ro > Ro_crit) then
            write(*, *) "No braking"
            jdot_total = 0d0
         end if

         write(*, *) "The total jdot is ", jdot_total, " cgs units."

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


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 10
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)

         real(dp) :: vals(n), dr, dc
         logical :: entered_g_mode_cavity

         real(dp) :: delta_nu_int, delta_Pg_int, delta_nu02_int, delta_Pg_outward, nu_max
         real(dp) :: I_tot, I_env, P_phot, ocz_turnover_time_g, ocz_turnover_time_bot_Hp
         real(dp) :: mixing_legnth_ocz_bot_Hp, conv_vel_ocz_bot_Hp, r_ocz_bot, m_ocz_bot, r_ocz_bot_Hp 
         real(dp) :: subcell_ocz_bot, subcell_ocz_bot_Hp, subcell_ocz_top, subcell_phot
         real(dp) :: sch_k0, sch_k1, sch_unstable(10000), turnover_time_k
         integer :: k, k_ocz_bot, k_ocz_bot_Hp, k_ocz_top, k_phot
         integer :: zonefound, nunstable

         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.

         !## proper handling of g mode period spacing integration (by WB)
         delta_Pg_int = 0.
         entered_g_mode_cavity = .false.
         do k = s% nz, 2, -1
            dr = s% rmid(k-1) - s% rmid(k)
            if (s% brunt_N2(k) > 0) then
               entered_g_mode_cavity = .true.
               delta_Pg_int = delta_Pg_int + sqrt(s% brunt_N2(k))/s% r(k)*dr
            else
               if (entered_g_mode_cavity) exit
            end if
         end do
         delta_Pg_int = (3.1415926535**2.) * sqrt(2.0) / delta_Pg_int

         !## integration of Dnu expressed in 1/s (by YL)
         Delta_nu_int = 0.
         do k = 2, s% nz, 1
            dr = s% rmid(k-1) - s% rmid(k)
            Delta_nu_int = Delta_nu_int + dr/s% csound(k)
         end do 
         Delta_nu_int = 1./(2.*Delta_nu_int)

         !## integration of dnu02 expressed in 1/s (by YL)
         delta_nu02_int = 0.
         nu_max = s% nu_max * 1d-6 ! 3090. * s% star_mass * pow(s%photosphere_r, -2.) * pow(s% Teff, -0.5)
         do k = 2, s% nz, 1
            ! dr = s% rmid(k-1) - s% rmid(k)
            dc = s% csound(k-1) - s% csound(k)
            delta_nu02_int = delta_nu02_int + dc / s%rmid(k)
         end do
         delta_nu02_int = delta_nu02_int - s% csound(1)/s% r(1)
         delta_nu02_int = -2 * Delta_nu_int * delta_nu02_int / (3.1415926535**2. * nu_max)


         !## rotation modeling (added by NS & JvS, edited by YL)
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         ! find location where tau=2/3
         do k = 2, s% nz
           if ((s% tau(k) > 6.66666666666666E-1)) then
             k_phot = k
             subcell_phot = ((6.66666666666666E-1) - s% tau(k_phot-1)) / (s% tau(k_phot) - s% tau(k_phot-1))
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
         ocz_turnover_time_bot_Hp = 0.
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

                  ocz_turnover_time_bot_Hp = mixing_legnth_ocz_bot_Hp / conv_vel_ocz_bot_Hp
                  exit
            end if
         end do

         ! !compute the "global" turnover time
         ! ocz_turnover_time_g = 0.
         ! do k = k_ocz_top, k_ocz_bot
         !    turnover_time_k = ((s%rmid(k-1)-s%rmid(k))/s%conv_vel(k))
         !    if (k == k_ocz_top) then
         !       ocz_turnover_time_g = ocz_turnover_time_g + turnover_time_k * subcell_ocz_top
         !    else if (k == k_ocz_bot) then
         !       ! (s% rmid(k+1) - s% rmid(k)) * subcell_ocz_bot / ( s%conv_vel(k) +  )
         !       ocz_turnover_time_g = ocz_turnover_time_g + turnover_time_k * subcell_ocz_bot
         !    else
         !       ocz_turnover_time_g = ocz_turnover_time_g + turnover_time_k
         !    end if
         ! end do    

         ! calculate the moment of inertia - direct integration
         I_tot = 0.
         I_env = 0.
         do k = 1, s% nz
            I_tot = I_tot + (2./3.) * s% rmid(k)**2. * s% dq(k)
            if (k < k_ocz_bot) then
               I_env = I_env + (2./3.) * s% rmid(k)**2. * s% dq(k)
            else if (k == k_ocz_bot) then
               I_env = I_env + (2./3.) * s% rmid(k)**2. * subcell_ocz_bot * s% dq(k)
            end if
         end do
         
         ! multiple by mass to get correct cgs units
         I_tot = I_tot * s% mstar
         I_env = I_env * s% mstar


         ! save new columns
         names(1) = 'delta_Pg_int'
         vals(1) = delta_Pg_int 
         names(2) = 'Delta_nu_int'
         vals(2) = Delta_nu_int
         names(3) = 'delta_nu02_int'
         vals(3) = delta_nu02_int
         names(4) = "I_tot"
         vals(4) = I_tot
         names(5) = "I_env"
         vals(5) = I_env
         names(6) = "P_phot"
         vals(6) = P_phot
         names(7) = "ocz_turnover_time_bot_Hp"
         vals(7) = ocz_turnover_time_bot_Hp
         names(8) = "r_ocz_bot"
         vals(8) = r_ocz_bot
         names(9) = "m_ocz_bot"
         vals(9) = m_ocz_bot
         names(10) = "P_surf"
         vals(10) = s% Peos(1)
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

         character(len=50), dimension(2) :: strings_to_replace(2), replace_with_strings(2)
         character(len=50) :: freq_min, freq_max
         real(dp) :: width !, obs_Teff_min, obs_Teff_max, obs_L_min, obs_L_max
         

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going


         ! ! save a profile, near ZAMS
         ! if (s%center_h1 < 0.70) and (.not. zams_model_saved) then
         !    zams_model_saved = .true.
         !    s% need_to_save_profiles_now = .true.
         ! endif

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      function trapz(y, x) result(integral)
         real(dp), dimension(:), intent(in) :: y, x
         real(dp) :: integral
         integer :: n

         ! integral = 0.0
         ! do i = 1, size(x) - 1
         !    integral = integral + 0.5_dp * (x(i+1) - x(i)) * (y(i+1) + y(i))
         ! end do
         n = size(x)
         integral = sum(0.5_dp*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))

      end function trapz

      function apply_mask(input_array, mask) result(output_array)
         real(dp), dimension(:), intent(in) :: input_array
         logical, dimension(:), intent(in) :: mask
         real(dp), dimension(:), allocatable :: output_array

         integer :: i, counter

         ! Check for size mismatch
         if (size(input_array) /= size(mask)) then
            print *, "Error: Size mismatch between input array and mask."
            return
         end if

         counter = count(mask)
         allocate(output_array(counter))

         counter = 0
         do i = 1, size(input_array)
            if (mask(i)) then
                  counter = counter + 1
                  output_array(counter) = input_array(i)
            end if
         end do
      end function apply_mask

      end module run_star_extras