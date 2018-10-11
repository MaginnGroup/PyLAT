subroutine calcdistances(nummol,comx, comy, comz, Lx, Ly, Lz,r)

integer, intent(in)	:: nummol
real, dimension(:), intent(in) :: comx
real, dimension(:), intent(in) :: comy
real, dimension(:), intent(in) :: comz
real, intent(in)	:: Lx
real, intent(in)	:: Ly
real, intent(in)	:: Lz
real	:: dx,dy,dz
real, dimension(0:nummol-1,0:nummol-1),intent(out)	:: r
	do i=1,nummol-1
		do j=i+1, nummol
			dx=comx(i)-comx(j)
			dy=comy(i)-comy(j)
			dz=comz(i)-comz(j)
			dx=dx-Lx*float(nint(dx/Lx))
			dy=dy-Ly*float(nint(dy/Ly))
			dz=dz-Lz*float(nint(dz/Lz))
			r(i-1,j-1)=sqrt(dx**2+dy**2+dz**2)
			r(j-1,i-1)=sqrt(dx**2+dy**2+dz**2)
		end do
	end do
return 
end subroutine calcdistances
