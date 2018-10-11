subroutine calccom(n, nummol, x, y, z,  mol, amass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, comxt, comyt, comzt)

integer, intent(in)	:: n, nummol
real, dimension(:), intent(in)	:: x,y,z
integer, dimension(:), intent(in)	:: mol
real, dimension(:), intent(in)	:: amass
real, dimension(:), intent(in)	:: molmass
real, intent(in)	:: Lx, Ly, Lz, Lx2, Ly2, Lz2
real, dimension(nummol), intent(out)	:: comxt, comyt, comzt
real, dimension(nummol)	:: xt, yt, zt
real, dimension(n)	:: ux, uy, uz
xt(:)=0.
comxt(:)=0.
comyt(:)=0.
comzt(:)=0.
do i=1,n
	if (xt(mol(i)) == 0.) then
		xt(mol(i)) = x(i)
		yt(mol(i)) = y(i)
		zt(mol(i)) = z(i)
		ux(i) = x(i)
		uy(i) = y(i)
		uz(i) = z(i)
	else
		if (x(i)-xt(mol(i)) > Lx2) then
			ux(i) = x(i)-Lx
		elseif (xt(mol(i))- x(i)  > Lx2) then
			ux(i) = x(i)+lx
		else
			ux(i) = x(i)
		endif

		if (y(i)-yt(mol(i)) > Ly2) then
			uy(i) = y(i)-Ly
		elseif (yt(mol(i))- y(i)  > Ly2) then
			uy(i) = y(i)+ly
		else
			uy(i) = y(i)
		endif

		if (z(i)-zt(mol(i)) > Lz2) then
			uz(i) = z(i)-Lz
		elseif (zt(mol(i))- z(i)  > Lz2) then
			uz(i) = z(i)+lz
		else
			uz(i) = z(i)
		endif
	endif

	comxt(mol(i)) = comxt(mol(i)) + ux(i)*amass(i)
	comyt(mol(i)) = comyt(mol(i)) + uy(i)*amass(i)
	comzt(mol(i)) = comzt(mol(i)) + uz(i)*amass(i)
enddo
do i=1,nummol
	comxt(i) = comxt(i)/molmass(i)
	comyt(i) = comyt(i)/molmass(i)
	comzt(i) = comzt(i)/molmass(i)
enddo
return
end subroutine calccom
