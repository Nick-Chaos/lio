!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Grid_integration()
   use garcha_mod      , only: integrate_density, NCO, NUNP
   implicit none
   double precision :: density
   integer :: electrons
   electrons=2*NCO+NUNP

   integrate_density=1
   call g2g_solve_groups(1,density,0)
   write(*,*) "Total Integrated Density: ", density, "N. electrons: ", electrons, "error: ", density-dble(electrons)

end subroutine Grid_integration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
