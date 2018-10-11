subroutine ipcorr(closest, skip, L1, L2, L3, CL,moltype,correlation)

integer, intent(in)	:: skip, L1, L2, L3, CL
integer, dimension(:,:,:), intent(in)	:: closest
integer, dimension(:), intent(in)	:: moltype
real, intent(out)	:: correlation(1:CL+1,1:L3,1:L3)
integer	:: norm
integer :: C1, C2

correlation = 0.d0

do i=1,L2
	do j=1,L3
		do k=skip+1,skip+CL
			do l=k,k+CL
				if (closest(k,i,j) == closest(l,i,j)) then
					correlation(l-k+1,moltype(i)+1,j) = correlation(l-k+1,moltype(i)+1,j)+1
				endif
			enddo
		enddo
	enddo
enddo

do i=1,L3
	do j=1,L3
		norm = correlation(1,i,j)
		do k=1,CL+1
			correlation(k,i,j) = correlation(k,i,j)/norm
		enddo
	enddo
enddo
end subroutine ipcorr
