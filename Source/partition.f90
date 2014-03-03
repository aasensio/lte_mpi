module atomic_partition_mod
implicit none
contains
!-----------------------------------------------------------------
! This subroutine calculates the atomic partition function for a given element
!	nel: atomic number
!	t: array with the temperature
!	u1: partition function for the neutral element (also an array of the same size as t)
!	u2: partition function for the 1st ionization (also an array of the same size as t)
!	u3: partition function for the 2nd ionization (also an array of the same size as t)	
!     A. D. Wittmann, Goettingen (1975).(1974, solar phys. 35, 11)
!-----------------------------------------------------------------
	subroutine partition_atomic(t, nel, u1, u2, u3)
	real(kind=8) :: t(:), u1(size(t)), u2(size(t)), u3(size(t))
	real(kind=8), allocatable :: x(:), y(:), coeff(:), zz(:)
	integer :: nel, n, i

		n = size(t)
		allocate(x(n))
		allocate(y(n))
		
		x = alog(5040.e0/real(t))
		y = 1.d-3*t
		select case(nel)
! H
			case(1)
				allocate(coeff(6))
				allocate(zz(n))
				zz = log10(t)
				coeff(1) = -2.61655891d2
				coeff(2) = 1.63428326d2
				coeff(3) = -4.06133526d1
				coeff(4) = 5.03282928d0
				coeff(5) = -3.10998364d-1
				coeff(6) = 7.66654594d-3
				zz=dlog(t)
				u1=(((((coeff(6)*zz+coeff(5))*zz)+coeff(4))*zz+coeff(3))*zz+coeff(2))*zz+coeff(1)
				u1=exp(u1)
				where(t > 16000.d0)
					u1 = dexp(7.19420668d-1)
				endwhere
				u2 = 1.d0
				u3 = 0.d0
				deallocate(coeff)
				deallocate(zz)
! He
			case(2)
				u1 = 1.d0 + 0.d0*t
				u2 = 2.d0 + 0.d0*t
				u3 = 1.d0 + 0.d0*t

! Li
			case(3)
				u1 = 2.081d0 - y*(6.8926d-2-y*1.4081d-2)
				where(t > 6.d3)
					u1 = 3.4864 + t*(-7.3292d-4+t*8.5586d-8)
				endwhere
				u2 = 1.d0 + 0.d0*t
      			u3 = 2.d0 + 0.d0*t

! Be
			case(4)
				do i = 1, n
					u1(i) = max(1.d0+0.d0*t(i),0.631d0+7.032d-5*t(i))
				enddo
				u2 = 2.d0 + 0.d0*t
				u3 = 1.d0 + 0.d0*t

! B
			case(5)
				u1 = 5.9351d0 + 1.0438d-2*y
				u2 = 1.d0 + 0.d0*t
				u3 = 2.d0 + 0.d0*t
			
! C
			case(6)
				u1 = 8.6985d0 + y * (2.0485d-2 + y * (1.7629d-2 - 3.9091d-4 * y))
				where(t > 1.2d4)
					u1 = 13.97d0 + t*(-1.3907d-3 + 9.0844d-8 * t)
				endwhere
				u2 = 5.838d0 + 1.6833d-5 * t
				where(t > 2.4d4)
					u2 = 10.989d0 + t * (-6.9347d-4 + t * 2.0861d-8)
				endwhere
				u3 = 1.d0 + 0.d0*t
				where(t > 1.95d4)
					u2 = -0.555d0 + 8.d-5 * t
				endwhere

! N
			case(7)
				u1 = 3.9914d0 + y * (1.7491d-2 - y * (1.0148d-2 - y * 1.7138d-3))
				where(t > 8800.d0 .and. t <= 1.8d4)
					u1 = 2.171d0 + 2.54d-4*t
				endwhere
				where(t > 1.8d4)
					u1 = 11.396d0 + t * (-1.7139d-3 + t * 8.633d-8)
				endwhere
				u2 = 8.060d0 + 1.420d-4*t
				where(t > 3.3d4)
					u2 = 26.793d0 + t * (-1.8931d-3 + t * 4.4612d-8)
				endwhere
				u3 = 5.9835d0 + t * (-2.6651d-5 + t * 1.8228d-9)
				where(t > 7310.5d0)
					u3 = 5.89d0
				endwhere
! O
			case(8)
				u1 = 8.29d0 + 1.10d-4*t
				where(t > 1.9d4)
					u1 = 66.81 + t * (-6.019d-3 + t * 1.657d-7)
				endwhere
				do i = 1, n
					u2(i) = max(4.d0+0.d0*t(i),3.51d0+8.d-5*t(i))
				enddo
				where(t > 2.d4)
					u2 = 68.7+t*(-4.216d-3+t*6.885d-8)
				endwhere
				u3 = 7.865d0 + 1.1348d-4*t
				
! F
			case(9)
				u1 = 4.5832d0+y*(0.77683d0+y*(-.20884d0+y*(2.6771d-2-1.3035d-3*y)))
				where(t > 8750.d0 .and. t <= 2.d4)
					u1 = 5.9d0
				endwhere
				where(t > 2.d4)
					u1 = 15.16+t*(-9.229d-4+t*2.312d-8)
				endwhere
				u2 = 8.15d0 + 8.9d-5*t
      			u3 = 2.315d0 + 1.38d-4*t
! Ne
			case(10)
				u1 = 1.d0
				where(t > 2.69d4)
					u1 = 26.3d0+t*(-2.113d-3+t*4.359d-8)
				endwhere
				u2 = 5.4d0 + 4.d-5 * t
				u3 = 7.973d0 + 7.956d-5 * t
! Na
			case(11)
				do i = 1, n
					u1(i) = max(2.+0.*t(i),1.72+9.3d-5*t(i))
				enddo
				where(t > 5400.d0 .and. t <= 8.5d3)
					u1 = -0.83+5.66d-4*t
				endwhere
				where(t > 8.5d3)
					u1 = 4.5568d0+t*(-1.2415d-3+t*1.3861d-7)
				endwhere
				u2 = 1.d0
				u3 = 5.69d0 + 5.69d-6 * t
! Mg
			case(12)
				u1 = 1.d0+exp(-4.027262-x*(6.173172+x*(2.889176+x*(2.393895+.784131*x))))
				where(t > 8.d3)
					u1 = 2.757+t*(-7.8909e-4+t*7.4531e-8)
				endwhere
				u2 = 2.d0+exp(-7.721172-x*(7.600678+x*(1.966097+.212417*x)))
				where(t > 2.d4)
					u2 = 7.1041+t*(-1.0817e-3+t*4.7841d-8)
				endwhere
				u3 = 1.d0

! Al
			case(13)
				u1 = 5.2955+y*(.27833-y*(4.7529d-2-y*3.0199d-3))
				do i = 1, n
					u2(i) = max(1.d0, 0.725+3.245d-5*t(i))
				enddo
				where(t > 2.24d4)
					u2 = 61.06+t*(-5.987e-3+t*1.485e-7)
				endwhere
				do i = 1, n
					u3(i) = max(2.d0,1.976d0+3.43d-6*t(i))
				enddo
				where(t > 1.814d4)
					u3 = 3.522+t*(-1.59d-4+t*4.382d-9)
				endwhere
! Si
			case(14)
				u1 = 6.7868+y*(.86319+y*(-.11622+y*(.013109-6.2013d-4*y)))
				where(t > 1.04d4)
					u1 = 86.01+t*(-1.465d-2+t*7.282d-7)
				endwhere
				u2 = 5.470+4.e-5*t
				where(t > 1.8d4)
					u2 = 26.44+t*(-2.22d-3+t*6.188d-8)
				endwhere
				do i = 1, n
					u3(i) = max(1.d0,0.911+1.1d-5*t(i))
				enddo
				where(t > 3.33d4)
					u3 = 19.14+t*(-1.408d-3+t*2.617d-8)
				endwhere
! P
			case(15)
				u1 = 4.2251+y*(-.22476+y*(.057306-y*1.0381e-3))
				where (t > 6.d3)
					u1 = 1.56+5.2e-4*t
				endwhere
				u2 = 4.4151+y*(2.2494+y*(-.55371+y*(.071913-y*3.5156e-3)))
				where(t > 7250.d0)
					u2 = 4.62+5.38d-4*t
				endwhere
				u3 = 5.595+3.4e-5*t
! S
			case(16)
				u1 = 7.5+2.15e-4*t
				where(t > 1.16d4)
					u1 = 38.76+t*(-4.906e-3+t*2.125e-7)
				endwhere
				u2 = 2.845+2.43e-4*t
				where(t > 1.05d4)
					u2 = 6.406+t*(-1.68e-4+t*1.323e-8)	
				endwhere
				u3 = 7.38d0 + 1.88d-4*t
! Cl
			case(17)
				u1 = 5.2+6.e-5*t
				where(t > 1.84d4)
					u1 = -81.6+4.8e-3*t
				endwhere
				u2 = 7.0+2.43e-4*t
      			u3 = 2.2+2.62e-4*t
! Ar
			case(18)
				u1 = 1.d0
				u2 = 5.20+3.8e-5*t
				u3 = 7.474+1.554e-4*t
! K
			case(19)
				u1 = 1.9909+y*(.023169-y*(.017432-y*4.0938e-3))
				where(t > 5800.d0)
					u1 = -9.93+2.124e-3*t
				endwhere
				u2 = 1.d0
				u3 = 5.304d0 + 1.93d-5*t
! Ca
			case(20)
				u1 = 1.+exp(-1.731273-x*(5.004556+x*(1.645456+x*(1.326861+.508553*x))))
				u2 = 2.+exp(-1.582112-x*(3.996089+x*(1.890737+.539672*x)))
      			u3 = 1.d0+0.*t
! Sc
			case(21)
				u1 = 4.+exp(2.071563+x*(-1.2392+x*(1.173504+.517796*x))) 
      			u2 = 3.+exp(2.988362+x*(-.596238+.054658*x))
      			u3 = 10.+0.*t
! Ti
			case(22)
				u1 = 5.+exp(3.200453+x*(-1.227798+x*(.799613+.278963*x)))
				where(t > 5.5d3)
					u1 = 16.37+t*(-2.838e-4+t*5.819e-7)
				endwhere
				u2 = 4.+exp(3.94529+x*(-.551431+.115693*x))
      			u3 = 16.4+8.5e-4*t			
! V
			case(23)
				u1 = 4.+exp(3.769611+x*(-.906352+x*(.724694+.1622*x)))
				u2 = 1.+exp(3.755917+x*(-.757371+.21043*x))
      			u3 = -18.+1.03e-2*t
				where(t < 2.25d3)
					u3 = 2.4d-3 * t
				endwhere
! Cr
			case(24)
				u1 = 7.+exp(1.225042+x*(-2.923459+x*(.154709+.09527*x))) 
      			u2 = 6.+exp(.128752-x*(4.143973+x*(1.096548+.230073*x)))
      			u3 = 10.4+2.1e-3*t
! Mn
			case(25)
				u1 = 6.+exp(-.86963-x*(5.531252+x*(2.13632+x*(1.061055+.265557*x)))) 
      			u2 = 7.+exp(-.282961-x*(3.77279+x*(.814675+.159822*x)))
      			u3 = 10.+t*0.
! Fe
			case(26)
				u1 = 9.+exp(2.930047+x*(-.979745+x*(.76027+.118218*x)))
				where(t < 4.d3)
					u1 = 15.85+t*(1.306e-3+t*2.04e-7)
				endwhere
				where(t > 9.d3)
					u1 = 39.149+t*(-9.5922e-3+t*1.2477e-6)
				endwhere
				u2 = 10.+exp(3.501597+x*(-.612094+.280982*x))
				where(t > 1.8d4)
					u2 = 68.356+t*(-6.1104e-3+t*5.1567e-7)
				endwhere
				u3 = 17.336+t*(5.5048e-4+t*5.7514e-8)
! Co
			case(27)
				u1 = 8.65+4.9e-3*t 
      			u2 = 11.2+3.58e-3*t
      			u3 = 15.0+1.42e-3*t
! Ni
			case(28)
				u1 = 9.+exp(3.084552+x*(-.401323+x*(.077498-.278468*x))) 
      			u2 = 6.+exp(1.593047-x*(1.528966+.115654*x))
      			u3 = 13.3+6.9e-4*t
! Cu
			case(29)
				do i = 1, n
					u1(i) = max(2.d0,1.50+1.51d-4*t(i))
				enddo
				where(t > 6250.d0)
					u1 = -.3+4.58e-4*t
				endwhere
				do i = 1, n
					u2(i) = max(1.d0,0.22+1.49d-4*t(i))
				enddo
				u3 = 8.025+9.4e-5*t
! Zn
			case(30)
				do i = 1, n
					u1(i) = max(1.d0,0.632d0+5.11d-5*t(i))
				enddo
				u2 = 2.d0
				u3 = 1.d0
! Ga
			case(31)
				u1 = 1.7931+y*(1.9338+y*(-.4643+y*(.054876-y*2.5054e-3)))
				where(t > 6.d3)
					u1 = 4.18+2.03e-4*t
				endwhere
				u2 = 1.d0
				u3 = 2.d0
! Ge
			case(32)
				u1 = 6.12+4.08e-4*t 
      			u2 = 3.445+1.78e-4*t
      			u3 = 1.1+0.*t
! As
			case(33)
				u1 = 2.65+3.65e-4*t
				u2 = -.25384+y*(2.284+y*(-.33383+y*(.030408-y*1.1609e-3)))
				where(t > 1.2d4)
					u2 = 8.d0
				endwhere
				u3 = 8.d0
! Se
			case(34)
				u1 = 6.34+1.71e-4*t 
      			u2 = 4.1786+y*(-.15392+3.2053e-2*y)
      			u3 = 8.+0.*t
! Br
			case(35)
				u1 = 4.12+1.12e-4*t 
      			u2 = 5.22+3.08e-4*t
      			u3 = 2.3+2.86e-4*t
! Kr
			case(36)
				u1 = 1.d0
				u2 = 4.11+7.4e-5*t
      			u3 = 5.35+2.23e-4*t
! Rb
			case(37)
				do i = 1, n
					u1(i) = max(2.d0,1.38+1.94d-4*t(i))
				enddo
				where (t> 6250.d0)
					u1 = -14.9+2.79e-3*t
				endwhere
				u2 = 1.d0
				u3 = 4.207+4.85e-5*t
! Sr
			case(38)
				u1 = .87127+y*(.20148+y*(-.10746+y*(.021424-y*1.0231e-3)))
				where(t > 6500.d0)
					u1 = -6.12+1.224e-3*t
				endwhere
				do i = 1, n
					u2(i) = max(2.d0,0.84+2.6d-4*t(i))
				enddo
				u3 = 1.d0
! Y
			case(39)
				u1 = .2+2.58e-3*t 
      			u2 = 7.15+1.855e-3*t
      			u3 = 9.71+9.9e-5*t
! Zr
			case(40)
				u1 = 76.31+t*(-1.866e-2+t*2.199e-6)
				where(t < 6236d0)
					u1 = 6.8+t*(2.806e-3+t*5.386e-7)
				endwhere
				u2 = 4.+exp(3.721329-.906502*x)
      			u3 = 12.3+1.385e-3*t
! Nb
			case(41)
				do i = 1, n
					u1(i) = max(1.d0,-0.19+1.43d-2*t(i))
				enddo
				u2 = -4.+1.015e-2*t
      			u3 = 25.+0.*t
! Mo
			case(42)
				do i = 1, n
					u1(i) = max(7.d0,2.1+1.5d-3*t(i))
				enddo
				where(t < 7.d3)
					u1 = -38.1+7.28e-3*t
				endwhere
				u2 = 1.25+1.17e-3*t
				where(t > 6900.d0)
					u2 = -28.5+5.48e-3*t
				endwhere
				u3 = 24.04+1.464e-4*t
! Tc
			case(43)
				u1 = 4.439+y*(.30648+y*(1.6525+y*(-.4078+y*(.048401-y*2.1538e-3))))
				where(t > 6.d3)
					u1 = 24.d0
				endwhere
				u2 = 8.1096+y*(-2.963+y*(2.369+y*(-.502+y*(.049656-y*1.9087e-3))))
				where(t > 6.d3)
					u2 = 17.d0
				endwhere
				u3 = 220.d0
! Ru
			case(44)
				u1 = -3.+7.17e-3*t 
      			u2 = 3.+4.26e-3*t
      			u3 = 22.+0.*t
! Rh
			case(45)
				u1 = 6.9164+y*(3.8468+y*(.043125-y*(8.7907e-3-y*5.9589e-4))) 
      			u2 = 7.2902+y*(1.7476+y*(-.038257+y*(2.014e-3+y*2.1218e-4)))
      			u3 = 30.+t*0.
! Pd
			case(46)
				do i = 1, n
					u1(i) = max(1.d0,-1.75+9.86d-4*t(i))
				enddo
				u2 = 5.60+3.62e-4*t
      			u3 = 20.+t*0.
! Ag
			case(47)
				do i = 1, n
					u1(i) = max(2.d0,1.537+7.88d-5*t(i))
					u2(i) = max(1.d0,0.73d0+3.4d-5*t(i))
				enddo
				u3 = 6.773+1.248e-4*t
! Cd
			case(48)
				do i = 1, n
					u1(i) = max(1.d0,0.43+7.6d-5*t(i))
				enddo
				u2 = 2.d0
				u3 = 1.d0
! In
			case(49)
				u1 = 2.16+3.92e-4*t 
      			u2 = 1.0+t*0.
      			u3 = 2.00+t*0.
! Sn
			case(50)
				u1 = 2.14+6.16e-4*t 
      			u2 = 2.06+2.27e-4*t
      			u3 = 1.05+0.*t
! Sb
			case(51)
				u1 = 2.34+4.86e-4*t 
      			u2 = .69+5.36e-4*t
      			u3 = 3.5+0.*t
! Te
			case(52)
				u1 = 3.948+4.56e-4*t 
      			u2 = 4.2555+y*(-.25894+y*(.06939-y*2.4271e-3))
				where(t > 1.2d4)
					u2 = 7.d0
				endwhere
				u3 = 5.d0
! I
			case(53)
				do i = 1, n
					u1(i) = max(4.d0,3.8d0+9.5d-5*t(i))
				enddo
				u2 = 4.12+3.e-4*t
      			u3 = 7.+0.*t
! Xe
			case(54)
				u1 = 1.d0
				u2 = 3.75+6.876e-5*t
      			u3 = 4.121+2.323e-4*t
! Cs
			case(55)
				do i = 1, n
					u1(i) = max(2.d0,1.56d0+1.67d-4*t(i))
				enddo
				where(t > 4750.d0)
					u1 = -2.680+1.04e-3*t
				endwhere
				u2 = 1.d0
				u3 = 3.769+4.971e-5*t
! Ba
			case(56)
				do i = 1, n
					u1(i) = max(1.d0,-1.8+9.85e-4*t(i))
				enddo
				where(t > 6850.d0)
					u1 = -16.2+3.08e-3*t
				endwhere
				u2 = 1.11+5.94e-4*t
      			u3 = 1.00+0.*t
! La
			case(57)
				u1 = 15.42+9.5e-4*t
				where(t > 5060.d0)
					u1 = 1.+3.8e-3*t
				endwhere
				u2 = 13.2+3.56e-3*t
      			u3 = 12.+0.*t
! Ce
			case(58)
				u1 = 9.+exp(5.202903+x*(-1.98399+x*(.119673+.179675*x))) 
      			u2 = 8.+exp(5.634882-x*(1.459196+x*(.310515+.052221*x)))
      			u3 = 9.+exp(3.629123-x*(1.340945+x*(.372409+x*(.03186-.014676*x))))
! Pr
			case(59)
				u2 = 9.+exp(4.32396-x*(1.191467+x*(.149498+.028999*x))) 
      			u1 = u2 
      			u3 = 10.+exp(3.206855+x*(-1.614554+x*(.489574+.277916*x)))
! Nb
			case(60)
				u1 = 9.+exp(4.456882+x*(-2.779176+x*(.082258+x*(.50666+.127326*x)))) 
      			u2 = 8.+exp(4.689643+x*(-2.039946+x*(.17193+x*(.26392+.038225*x))))
      			u3 = u2
! Pm
			case(61)
				u1 = 20.d0
				u2 = 25.d0
				u3 = 100.d0
! Sm
			case(62)
				u1 = 1.+exp(3.549595+x*(-1.851549+x*(.9964+.566263*x))) 
      			u2 = 2.+exp(4.052404+x*(-1.418222+x*(.358695+.161944*x)))
      			u3 = 1.+exp(3.222807-x*(.699473+x*(-.056205+x*(.533833+.251011*x))))
! Eu
			case(63)
				u1 = 8.+exp(1.024374-x*(4.533653+x*(1.540805+x*(.827789+.286737*x))))
      			u2 = 9.+exp(1.92776+x*(-1.50646+x*(.379584+.05684*x)))
      			u3 = 8.+0.*t
! Gd
			case(64)
				u1=5.+exp(4.009587+x*(-1.583513+x*(.800411+.388845*x))) 
      			u2=6.+exp(4.362107-x*(1.208124+x*(-.074813+x*(.076453+.055475*x))))
      			u3=5.+exp(3.412951-x*(.50271+x*(.042489-4.017e-3*x)))
! Tb
			case(65)
				u1 = 16.+exp(4.791661+x*(-1.249355+x*(.570094+.240203*x))) 
      			u2 = 15.+exp(4.472549-x*(.295965+x*(5.88e-3+.131631*x)))
      			u3 = u2
! Dy
			case(66)
				u1 = 17.+exp(3.029646-x*(3.121036+x*(.086671-.216214*x))) 
      			u2 = 18.+exp(3.465323-x*(1.27062+x*(-.382265+x*(.431447+.303575*x))))
      			u3 = u2
! Ho
			case(67)
				u3 = 16.+exp(1.610084-x*(2.373926+x*(.133139-.071196*x))) 
      			u1 = u3
      			u2 = u3
! Er
			case(68)
				u1 = 13.+exp(2.895648-x*(2.968603+x*(.561515+x*(.215267+.095813*x)))) 
      			u2 = 14.+exp(3.202542-x*(.852209+x*(-.226622+x*(.343738+.186042*x))))
      			u3 = u2
! Tm
			case(69)
				u1 = 8.+exp(1.021172-x*(4.94757+x*(1.081603+.034811*x))) 
      			u2 = 9.+exp(2.173152+x*(-1.295327+x*(1.940395+.813303*x)))
      			u3 = 8.+exp(-.567398+x*(-3.383369+x*(.799911+.554397*x)))
! Yb
			case(70)
				u1 = 1.+exp(-2.350549-x*(6.688837+x*(1.93869+.269237*x))) 
      			u2 = 2.+exp(-3.047465-x*(7.390444+x*(2.355267+.44757*x)))
      			u3 = 1.+exp(-6.192056-x*(10.560552+x*(4.579385+.940171*x)))
! Lu
			case(71)
				u1 = 4.+exp(1.537094+x*(-1.140264+x*(.608536+.193362*x)))
				do i = 1, n
					u2(i) = max(1.d0,0.66+1.52d-4*t(i))
				enddo
				where(t > 5250.d0)
					u2 = -1.09+4.86e-4*t
				endwhere
				u3 = 5.d0
! Hf
			case(72)
				u1 = 4.1758+y*(.407+y*(.57862-y*(.072887-y*3.6848e-3))) 
      			u2 = -2.979+3.095e-3*t
      			u3 = 30. +0.*t
! Ta
			case(73)
				u1 = 3.0679+y*(.81776+y*(.34936+y*(7.4861e-3+y*3.0739e-4))) 
      			u2 = 1.6834+y*(2.0103+y*(.56443-y*(.031036-y*8.9565e-4)))
      			u3 = 15.+0.*t
! W
			case(74)
				u1 = .3951+y*(-.25057+y*(1.4433+y*(-.34373+y*(.041924-y*1.84e-3))))
				where(t > 1.2d4)
					u1 = 23.d0
				endwhere
				u2 = 1.055+y*(1.0396+y*(.3303-y*(8.4971e-3-y*5.5794e-4)))
      			u3 = 20.+0.*t
! Re
			case(75)
				u1 = 5.5671+y*(.72721+y*(-.42096+y*(.09075-y*3.9331e-3)))
				where(t > 1.2d4)
					u1 = 29.d0
				endwhere
				u2 = 6.5699+y*(.59999+y*(-.28532+y*(.050724-y*1.8544e-3)))
				where(t > 1.2d4)
					u2 = 22.d0
				endwhere
				u3 = 20.d0
! Os
			case(76)
				u1 = 8.6643+y*(-.32516+y*(.68181-y*(.044252-y*1.9975e-3))) 
      			u2 = 9.7086+y*(-.3814+y*(.65292-y*(.064984-y*2.8792e-3)))
      			u3 = 10.+0.*t
! Ir
			case(77)
				u1 = 11.07+y*(-2.412+y*(1.9388+y*(-.34389+y*(.033511-1.3376e-3*y))))
				where(t > 1.2d4)
					u1 = 30.d0
				endwhere
				u2 = 15.d0
				u3 = 20.d0
! Pt
			case(78)
				u1 = 16.4+1.27e-3*t 
      			u2 = 6.5712+y*(-1.0363+y*(.57234-y*(.061219-2.6878e-3*y)))
      			u3 = 15.d0
! Au
			case(79)
				u1 = 1.24+2.79e-4*t 
      			u2 = 1.0546+y*(-.040809+y*(2.8439e-3+y*1.6586e-3))
      			u3 = 7.+0.*t
! Hg
			case(80)
				u1 = 1.0+0.*t 
      			u2 = 2.+0.*t
				do i = 1, n
					u3(i) = max(1.d0,.669+3.976e-5*t(i))
				enddo
! Tl
			case(81)
				do i = 1, n
					u1(i) = max(2.d0,0.63+3.35e-4*t(i))
				enddo
				u2 = 1.d0
				u3 = 2.d0
! Pb
			case(82)
				do i = 1, n
					u1(i) = max(1.d0,0.42+2.35e-4*t(i))					
				enddo
				where(t > 6125.d0)
					u1 = -1.2+5.e-4*t
				endwhere
				do i = 1, n
					u2(i) = max(2.d0,1.72+7.9e-5*t(i))
				enddo
				u3 = 1.d0
! Bi
			case(83)
				u1 = 2.78+2.87e-4*t
				do i = 1, n
					u2(i) = max(1.d0,.37+1.41e-4*t(i))
				enddo
				u3 = 2.5d0
! Po
			case(84)
				u1 = 5.+0.*t 
      			u2 = 5.+0.*t
      			u3 = 4.+0.*t
! At
			case(85)
				u1 = 4.+0.*t  
      			u2 = 6.+0.*t
      			u3 = 6.+0.*t
! Rn	
			case(86)
				u1 = 1. +0.*t
      			u2 = 4.+0.*t
      			u3 = 6.+0.*t
! Fr
			case(87)
				u1 = 2. +0.*t
      			u2 = 1.+0.*t
      			u3 = 4.5+0.*t
! Ra
			case(88)
				u1 = 1. +0.*t
      			u2 = 2.+0.*t
      			u3 = 1.+0.*t
! Ac
			case(89)
				u1 = 6. +0.*t
      			u2 = 3.+0.*t
      			u3 = 7.+0.*t
! Th
			case(90)
				u1 = 8. +0.*t
      			u2 = 8.+0.*t
      			u3 = 8.+0.*t		
! Pa
			case(91)
				u1 = 50. +0.*t
      			u2 = 50.+0.*t
      			u3 = 50.+0.*t
! U
			case(92)
				u1 = 25. +0.*t
      			u2 = 25.+0.*t
      			u3 = 25.+0.*t
			case DEFAULT
				print *, 'Partition function not available...'
		end select
		deallocate(x)
		deallocate(y)
	end subroutine partition_atomic
end module atomic_partition_mod
