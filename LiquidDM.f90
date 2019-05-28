program LiquidDM
implicit none

                        !binding energy arrays for neutron and proton arrays

real:: prevSepEn, sen, prevSepEz, sez                 !separation energy values for the neutron drip and proton drip binding energies
integer:: z, n, a
integer:: nfinn, zfinn, nfinz, zfinz

real,    allocatable :: ben(:), bez(:)
integer, allocatable:: ndripzdep(:), ndripndep(:)         !array of z values for neutron drip the first neutron dependent the second proton dependent
integer, allocatable:: zdripndep(:), zdripzdep(:)
integer, allocatable::NZline(:)

open(unit=32, file="ProtonDripline.dat", status="unknown")
open(unit=33, file="NeutronDripline.dat", status="unknown")
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
write(*,*) "Outer loop of the proton drip line must be smaller than the inner loop"
write(*, "(A)", advance="NO") "desired n range for protron drip: nfinz="      !outer loop limit value should be smaller than the inner loop limit value 
read *, nfinz                        

write(*,*) "inner loop of the proton drip line must be larger than the outer loop"
write(*, "(A)", advance="NO") "desired z range for protron drip: "          
read *, zfinz 

write(*,*) "Outer loop of the Neutron drip line must be smaller than the inner loop"
write(*, "(A)", advance="NO") "desired z range for neutron drip: "      !outer loop limit value should be smaller than the inner loop limit value      
read *, zfinn 

write(*,*) "inner loop of the Neutron drip line must be larger than the outer loop"
write(*, "(A)", advance="NO") "desired n range for neutron drip: "       
read *, nfinn
                         
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////

allocate(bez(zfinz))
allocate(zdripndep(nfinz))
!allocate(zdripzdep(nfinz))!IRELIVANT

allocate(ben(nfinn))
allocate(ndripzdep(zfinn))
!allocate(ndripndep(nfinn))!IRELIVANT
!allocate(NZline(nfinn))

!Z=1
!do n=1, 120
!ndripndep(n)=0
!NZLINE(N)=z
!Z=Z+1
!end do 
!do z=1, 120
!zdripzdep(z)=0
!end do


!=============================================
!===========SEC 1 Proton drip line============
!=============================================
	nlpA: do n=1, nfinz
		
		zlpB: do z=1, zfinz
		call BindingE(z,n,bez(z))
	    	end do zlpB
	
	   !finding the separation energy: the difference between binding energies incremented z 
	   sepEzB: do z=1, zfinz
	   sez=bez(z+1)-bez(z)

	        if (sez<0.) then 
	        !}}}}}}}}}}}}}}}
	        zdripndep(n)=z!}  !this array is the proton drip line z values dependent on n
	        !}}}}}}}}}}}}}}}
	        goto 111       
	        end if
	        
	        prevSepEz=sez   
	   end do sepEzB
	   
	111 continue 	
	end do nlpA
!=============================================
!=========SEC 3 Neutron drip line=============
!=============================================    
	zloopA: do z=1, zfinn

		nloopB: do n=1, nfinn
		call BindingE(z, n, ben(n))
		end do nloopB
        
        !finding the separation energy: the difference between binding energies incremented n
    	sepEnB: do n=1, nfinn
    	sen=ben(n+1)-ben(n)

    	    if (sen<0.) then 
    	       !}}}}}}}}}}}}}}}}
	        ndripzdep(z)=n!}  this array is the NEUTRON DRIP LINE N VALUES DEPENDENT ON Z 
	       !}}}}}}}}}}}}}}}} 
	        goto 333
	        end if
	        
	    prevSepEn=sen
	    end do sepEnB
	    
	333 continue
	end do zloopA

!=============================================
!========SEC 5 printing/writing===============
!=============================================
    prn1A: do n=1, nfinz !SEC 1 
            !||||||||||||||||||||||||||||||||
            print *, "n=", n, "z=", zdripndep(n), "proton drip line value"         !printing the proton drip line z values dependent on n
            write(32,*) n, zdripndep(n)
	    !||||||||||||||||||||||||||||||||
    end do prn1A

    prn4A: do z=1, zfinn  !SEC 3
            !||||||||||||||||||||||||||||||||
            print *, "z=", z, "n=", ndripzdep(z), "Neutron drip line value"        !printing the Neutron drip line n values dependent on z
            write(33,*) ndripzdep(z), z
	    !||||||||||||||||||||||||||||||||
    end do prn4A
  
end program LiquidDM
!-------------------------------------
!-------------------------------------
!-------------------------------------
!-------------------------------------
!-------------------------------------
!subroutine to find the binding energy
subroutine BindingE(z,n,be)
implicit none
real:: be, p1, p2, p3, p4, p5
integer:: z, n, a
a=z+n
	!-----------------------------
	!----------LQDM---------------
	!-----------------------------
	p1 =15.5*(A)
	p2=  -16.8*(A**(2./3.))
	p3=-(0.72*(Z*(Z-1.))*(A**(-1./3.)))
	p4= -(23.*((Z-N)**2.))/A
    !||||||||||||||||||||||||||||||||||
	if(    (mod(z,2)==0)  ) then
	if(    (mod(n,2)==0)  ) then
	p5=34.*(a**(-3./4.))
	!print *, "p5 z,n= even"
	end if
	if(    (mod(n,2)==1)  ) then
	p5=0
	!a is odd
	end if
	end if !main if
	if(    (mod(z,2)==1)  ) then
	if(    (mod(n,2)==1)  ) then
	p5=-34.*(a**(-3./4.))
	!print *, "p5 z,n= odd"
	end if
	if(    (mod(n,2)==0)  ) then
	p5=0
	!a is odd
	end if
	end if !main if
	!p5=0 !ignoring P5
	!-----------------------------
	BE= p1+p2+p3+p4+p5 !creating array of binding energy dependent on N
	!-----------------------------
	!-----------------------------
	!-----------------------------
end subroutine bindingE
