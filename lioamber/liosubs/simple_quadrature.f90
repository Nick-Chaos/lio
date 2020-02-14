subroutine simple_quad()
use rhoint
use constants, only : bohr
implicit none
double precision :: xmin, xmax, ymin, ymax, zmin, zmax
double precision :: dx,dy,dz
double precision :: integral
double precision :: x,y,z

xmin=w_rho_xmin/bohr
xmax=w_rho_xmax/bohr
ymin=w_rho_ymin/bohr
ymax=w_rho_ymax/bohr
zmin=w_rho_zmin/bohr
zmax=w_rho_zmax/bohr
dx=w_rho_dx/bohr
dy=w_rho_dy/bohr
dz=w_rho_dz/bohr

if (write_int_rho == 'xy') then
open(unit=978, file='rho_z.dat')

do z=zmin,zmax,dz
  integral=0.d0
  do x=xmin+dx*0.5d0, xmax, dx
    do y=ymin+dy*0.5d0, ymax, dx
      integral=integral+obtainrho(x,y,z)*dx*dy
      write(55978,*) x,y,z,integral
    end do
  end do
  write(978,*) z*bohr, integral
end do

close(978)
end if
contains

double precision function obtainrho(x,y,z)
use garcha_mod, only: r, Pmat_vec
use basis_data, only: M, ncont, nuc, nshell, a, c
implicit none
double precision, intent(in) :: x,y,z
double precision, dimension(3) :: eval_p
double precision :: p_val, p_func, p_dist
double precision, allocatable :: p_array(:)
double precision, parameter   :: expmax = 10.0D0
integer :: ns, np, nd

integer :: ii, jj, ni, jjj, kkk


ns = nshell(0)
np = nshell(1)
nd = nshell(2)

allocate(p_array(ns + 3*np + 6*nd))


eval_p(1)=x
eval_p(2)=y
eval_p(3)=z

p_val = 0.D0

! Calculate function values at this voxel, store in energy
! weighted Rho. s functions
do ii = 1, ns
   p_dist = 0.D0
   do jj = 1, 3
      p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
   enddo

   p_func = 0.D0
   do ni = 1, ncont(ii)
      if ((a(ii,ni)*p_dist) < expmax) &
         p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
   enddo
   p_array(ii) = p_func
enddo

! p functions
do ii = ns+1, ns+np, 3
   p_dist = 0.D0
   do jj = 1, 3
      p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
   enddo

   p_func = 0.D0
   do ni = 1, ncont(ii)
      if ((a(ii,ni)*p_dist) < expmax) &
         p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
   enddo

   do jj = 1,3
      p_array(ii+jj-1) = p_func * (eval_p(jj) - r(Nuc(ii),jj))
   enddo
enddo

! d functions
do ii = ns+np+1, M, 6
   p_dist = 0.D0
   do jj = 1,3
      p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
   enddo

   p_func = 0.D0
   do ni = 1, ncont(ii)
      if ((a(ii,ni) * p_dist) < expmax) &
         p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
   enddo

   kkk = 0
   do jj = 1, 3
      do jjj = 1, jj
         kkk = kkk + 1
         p_array(ii+kkk-1) = p_func * (eval_p(jj)  - r(Nuc(ii),jj)) *&
                                (eval_p(jjj) - r(Nuc(ii),jjj))
      enddo
   enddo
enddo

! Calculate density
kkk = 0
do ii = 1 , M
   do jj = ii, M
      kkk   = kkk + 1
      p_val = p_val + Pmat_vec(kkk) * p_array(ii) * p_array(jj)
   enddo
enddo

obtainrho = p_val
return
end function obtainrho


end subroutine

