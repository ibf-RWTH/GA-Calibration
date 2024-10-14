      subroutine caleulang(Mx,vv,ising)     
          implicit none
          integer irdummy,ising,ising_1
          real(8) det
          real(8) vv(3)
          real(8) Mx(3,3),rm(3,3),um(3,3)
          real(8) pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180/pi
          g2r=pi/180
          ising_1=0
!         call PDECOMPOSITION(Mx,Um,Rm,ising_1)
          call polar_decomp(Mx,Um,Rm,ising_1)
          if(ising_1/=0)then
              write(6,*) 'pdecompsition of Mx is failure'
!             call pm(Mx,3,3) !Commented out to make the subroutine work with windows compiler
              ising=1
              return
          endif
          call icams_Q2Eang(Rm,vv(1),vv(2),vv(3))
          return
      end subroutine caleulang

      subroutine icams_misori(v1,v2,ang)
          implicit none
          real(8) v1(3),v2(3),ang,x1,x2
          real(8) QM1(3,3),QM2(3,3),dQM(3,3)
          real(8) pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180/pi
          g2r=pi/180
          call icams_Eang2Q(v1(1),v1(2),v1(3),QM1)
          call icams_Eang2Q(v2(1),v2(2),v2(3),QM2)
          dQM=matmul(QM2,transpose(QM1))
          x1=dQM(1,1)+dQM(2,2)+dQM(3,3)
          x2=(x1-1.d0)/2
          if(dabs(x2)>1.d0)x2=1.d0*sign(1.d0,x2)
          ang=dabs(pi/2-dasin(x2)) *r2g
          return
      end subroutine icams_misori

!****************************************************************
      subroutine pdecomposition(Mx,UMx,RMx,ising)
          implicit none
          integer ising
          real(8) Mx(3,3),ce(3,3),RMx(3,3),UMx(3,3),IUMx(3,3)
          real(8) eb1(3,3),eb2(3,3),eb3(3,3)
          real(8) ev1,ev2,ev3,det
          ising=0
          ce=matmul(transpose(Mx),Mx)
!         call spectral(ce,ev1,ev2,ev3,eb1,eb2,eb3,ising) !Commented out to make the subroutine work with windows compiler
          if(ev1<=0 .or. ev2<=0 .or. ev3<=0 .or. ising/=0)then
              write(6,*) 'eigen value of ce <0'
              print*, ev1,ev2,ev3
              ising=1
              return
          endif
          UMx=dsqrt(ev1)*eb1+dsqrt(ev2)*eb2+dsqrt(ev3)*eb3
          IUMx=1/dsqrt(ev1)*eb1+1/dsqrt(ev2)*eb2+1/dsqrt(ev3)*eb3
          RMx=matmul(Mx,IUMx)
          return 
      end subroutine pdecomposition

!cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Eang2Q_OLD(p1,p,p2,QM)
!cccccccccccccccccccccccccccccccccccccccccccccccc
!         rotate from X[100],Y[010],Z[001] to v1,v2,v3
!cccccccccccccccccccccccccccccccccccccccccccccccc
          implicit none
          real(8) QM(3,3)
          real(8) p1,p,p2,xp1,xp,xp2
          real(8) c1,c,c2,s1,s,s2
          real(8) pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180/pi
          g2r=pi/180
          xp1=p1*g2r
          xp =p *g2r
          xp2=p2*g2r
          c1=dcos(xp1)
          s1=dsin(xp1)
          c =dcos(xp )
          s =dsin(xp )
          s2=dsin(xp2)
          c2=dcos(xp2)
          QM(1,1)=+c1*c2-s1*s2*c
          QM(1,2)=+s1*c2+c1*s2*c
          QM(1,3)=+s2*s
          QM(2,1)=-c1*s2-s1*c2*c
          QM(2,2)=-s1*s2+c1*c2*c
          QM(2,3)=+c2*s
          QM(3,1)=+s1*s
          QM(3,2)=-c1*s
          QM(3,3)=+c
          return
      end subroutine icams_Eang2Q_OLD

!cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine icams_Eang2QAnxin(p1,p,p2,QM)
!cccccccccccccccccccccccccccccccccccccccccccccccc
!     rotate from X[100],Y[010],Z[001] to v1,v2,v3
!cccccccccccccccccccccccccccccccccccccccccccccccc
          implicit none
          real(8) QM(3,3)
          real(8) p1,p,p2,xp1,xp,xp2
          real(8) c1,c,c2,s1,s,s2
          real(8) pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180/pi
          g2r=pi/180
          xp1=p1*g2r
          xp =p *g2r
          xp2=p2*g2r
          c1=dcos(xp1)
          s1=dsin(xp1)
          c =dcos(xp )
          s =dsin(xp )
          s2=dsin(xp2)
          c2=dcos(xp2)
          QM(1,1)=+c1*c2-s1*s2*c
          QM(1,2)=+s1*c2+c1*s2*c
          QM(1,3)=+s2*s
          QM(2,1)=-c1*s2-s1*c2*c
          QM(2,2)=-s1*s2+c1*c2*c
          QM(2,3)=+c2*s
          QM(3,1)=+s1*s
          QM(3,2)=-c1*s
          QM(3,3)=+c
          return
      end subroutine icams_Eang2QAnxin

!********************************************************************** 
      subroutine HI(M,HI1M,HI2M,HI3M)
!---------HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M
          IMPLICIT NONE
          real(8) M(3,3),HI1M,HI2M,HI3M 
            HI1M = M(1,1)+M(2,2)+M(3,3)
            HI2M =(M(1,1)+M(2,2)+M(3,3))**2/2.d0-M(1,1)**2/2.d0-M(2,2)**2/2.d0-
     &            M(3,3)**2/2.d0-M(1,2)*M(2,1)-M(1,3)*M(3,1)-M(2,3)*M(3,2)
            HI3M =+M(1,1)*M(2,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*
     &            M(1,2)*M(2,3)-M(1,1)*M(2,3)*M(3,2)-
     &            M(2,1)*M(1,2)*M(3,3)-M(3,1)*M(1,3)*M(2,2)
          RETURN  
      end subroutine HI


!**********************************************************************
      subroutine pv(M,ni)
          implicit none
          integer ni,i
          real(8) M(ni)
          write(*,100) (M(i),i=1,ni)          
!         write(6,*) 
!100      format(48e20.8)       
100       format(48e12.4)       
          return
      end subroutine pv


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> Calculates inverse of Matrix A(n,n)
!> @param[in]   A                Matrix A
!> @param[in]   n                dimension of matrix
!> @param[out]  B                Inverse of matrix A
!> @param[out]  ising            Error flag
!ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussj(A,n,B,ising)
!ccccccccccccccccccccccccccccccccccccccccccccccc
          implicit none
          integer   n,ising
          integer   i,icol,irow,j,k,l,ll,i1,i2
          integer   indxc(n),indxr(n),ipiv(n)
          real(8)   A(n,n),B(n,n),c(n,n),vx(n)
          real(8)   big,dum,pivinv,x1,x2,x3
!-------------------------------------------
!         write(6,*) 'coming to jd'
!         call pm(A,3,3)
!         call flush(6)
          ising=0
          C=A
          ipiv=0
!---------loop for the pivot procedures, from 1 to n
          do i=1,n
              big=0.d0
              do j=1,n
                  do k=1,n
                      if(ipiv(j).ne.1 .and. ipiv(k)==0)then
                          if(dabs(a(j,k)).ge.big)then
                              big=dabs(a(j,k))
                              irow=j
                              icol=k
                          endif
                      endif
                      if(ipiv(j).ne.1 .and. ipiv(k).gt.1)then
!                         print*,'sigular matrix in gauss_jordan'
!                         write(6,*)'sigular matrix in gauss_jordan'
!                         call flush(6)
                          ising=1
                          return
                      endif
                  enddo
              enddo
!------------check whether the pivot element is zero or not
              if(a(irow,icol)==0.d0)then
!                 print*,'sigular matrix in gauss_jordan'
!                 print*,'indices is:', irow,icol
!                 write(6,*)'sigular matrix in gauss_jordan'
!                 write(6,*)'indices is:', irow,icol
!                 call flush(6)
                  ising=1
                  return
              endif
!-----------------------------------------------------------
!             if one component is selected as pivot element
!             the second indice is important, so it is marked
!             from 0 to 1 in ipiv(:) array
!             after convert, it is the row number
!----------------------------------------------------------
              ipiv(icol)=ipiv(icol)+1
!-------------record the row and collum number for ith pivot element
              indxr(i)=irow
              indxc(i)=icol
!-------------change pivot element to diagonal position
              if(irow.ne.icol)then
                  vx=a(irow,:)
                  a(irow,:)=a(icol,:)
                  a(icol,:)=vx
              endif
!-------------eliminate the elements besides a(icol,icol)
              pivinv=1.d0/a(icol,icol)
              a(icol,icol)=1.d0
              a(icol,:)=a(icol,:)*pivinv 
              do i2=1,n
                  if(i2.ne.icol)then
                      dum=a(i2,icol)
                      a(i2,icol)=0.d0
                      a(i2,:)=a(i2,:)-a(icol,:)*dum
                  endif
              enddo
          enddo
!---------after maximum pivot strategy elimination
!---------rearrage the left matrix
          do l=n,1,-1
             if(indxr(l).ne.indxc(l))then
                vx=a(:,indxr(l))
                a(:,indxr(l))=a(:,indxc(l))
                a(:,indxc(l))=vx
             endif
          enddo
          B=A
          A=C
          return
      end subroutine gaussj

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @author 
   !> Anxin Ma
   ! DESCRIPTION: 
   !> c Performs the polar decomposition of the	
   !>tensor F=RU using the Cayley-Hamilton theorem.
   !> @param[in]   F                Deformation gradient
   !> @param[out]  cross product    returns vector
!--------------------------------------------------------------------------
      subroutine polar_decomp(F,U,R,ising)
!--------------------------------------------------------------------------
          implicit none
          integer i,j,k,iflag1,iflag2,ising
          real(8) F(3,3),U(3,3),R(3,3),RT(3,3),C(3,3),CS(3,3),UINV(3,3)
          real(8) x1,x2,x3,f1,f2,c1,c2,c3,u1,u2,u3,b,b1,b2,b3,b4
          C=matmul(transpose(F),F)
          CS=matmul(C,C)
!---------1st, 2st, 3rd invariant for tensor C
          C1=C(1,1)+C(2,2)+C(3,3)
          C2=(C1**2.d0-(CS(1,1)+CS(2,2)+CS(3,3)))/2
          C3=+C(1,1)*(C(2,2)*C(3,3)-C(2,3)*C(3,2))-C(1,2)*(C(2,1)*C(3,3)-C(2,3)*C(3,1))+C(1,3)*(C(2,1)*C(3,2)-C(2,2)*C(3,1))
!---------3rd invariant for tensor U
          U3=dsqrt(C3)
          X1= 2.0**5.0 /27.0 * (2.0*C1**3.0-9.0*C1*C2+27.0*C3)
          X2= 2.**10.0 /27.0 * (4.0*C2**3.0 - C1**2.0*C2**2.0 + 4.0*C1**3.0*C3 - 18.0 * C1*C2*C3 + 27.0 * C3**2.0)
          
          IF(X2<0)X2=0
              F1=X1+dsqrt(X2)
              IFLAG1=0
          IF(F1<0)IFLAG1=1
              F2=X1-dsqrt(X2)
              IFLAG2=0
          IF(F2<0) IFLAG2=1
          IF(IFLAG1==1) F1=-F1
          IF(IFLAG2==1) F2=-F2
          X3= -2.0/3.0*C1 + F1**(1.0/3.0) + F2**(1.0/3.0)
          IF(IFLAG1==1) X3= -2.0/3.0*C1 - F1**(1.0/3.0) + F2**(1.0/3.0)
          IF(IFLAG2==1) X3= -2.0/3.0*C1 + F1**(1.0/3.0) - F2**(1.0/3.0)
!---------1st, 2nd invariant for tensor U
          B=-2.0*C1
          if(X3==B)then
              U1=dsqrt(C1+2.0*dsqrt(C2))
          else
              x1=dsqrt(2.0*C1+X3)
              if(x1==0)then
                  ising=1
                  return
              endif      
              U1= 0.5 * ( x1 + dsqrt(2.0*C1 - X3 + 16.0*dsqrt(C3)/x1) )
          endif
          U2=dsqrt(C2+2.0*U3*U1)
          B1= U3**2.0 * (U3+U1*C1) + U1**2.0 * (U1*C3+U3*C2)
          if(B1==0)then
              ising=1
              return
          endif
          B2= U1*(U1*U2-U3) / B1
          B3=-(U1*U2-U3) * (U3+U1*C1) / B1
          B4= (U2*U3*(U3+U1*C1) + U1**2.0 * (U2*C2+C3))/B1
          UINV=B2*CS + B3*C
          Uinv(1,1)=Uinv(1,1)+B4
          Uinv(2,2)=Uinv(2,2)+B4
          Uinv(3,3)=Uinv(3,3)+B4
          R=matmul(F,UINV)
          U=matmul(transpose(R),F)
          return
      end subroutine polar_decomp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> c     Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
!> c     by n, stored in a physical np by np array. On output, elements of a above the diagonal are
!> c     destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
!> c     logical and physical dimensions as a, whose columns contain, on output, the normalized
!> c     eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
!> @param[in]   a,b              3x1 vectors
!> @param[out]  cross product    returns vector
!--------------------------------------------------------------------------
      subroutine jacobi(a,n,np,d,v,nrot,ising)
          implicit none
!------------------------------------------------------------------------------		      
          integer :: n,np,nrot,nmax
          real(8) :: a(np,np),d(np),v(np,np)
          parameter (nmax=500)   !???
          integer i,ip,iq,j,ising
          real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
          ising=0
          do ip=1,n  !initialize to the identity matrix.
              do iq=1,n
                  v(ip,iq)=0.
              enddo
              v(ip,ip)=1.
          enddo 
          do ip=1,n
              b(ip)=a(ip,ip) !initialize b and d to the diagonal of a.
              d(ip)=b(ip)
              z(ip)=0.  !this vector will accumulate terms of the form tapq
          enddo        !as in equation (11.1.14).
          nrot=0
          do i=1,50
              sm=0.
              do ip=1,n-1  !sum off-diagonal elements.
                  do iq=ip+1,n
                      sm=sm+dabs(a(ip,iq))
                  enddo 
              enddo 
              if(sm==0.)return   
              !===================!
              ! sucessful return  !
              !===================!
              if(i.lt.4)then       
                   tresh=0.2*sm/n**2  !...on the first three sweeps.
              else
                   tresh=0.           !...thereafter.
              endif
              do ip=1,n-1
                  do iq=ip+1,n
                      g=100.*dabs(a(ip,iq))
                      !after four sweeps, skip the rotation if the off-diagonal element is small.
                      if((i.gt.4).and.(dabs(d(ip))+g==dabs(d(ip))).and.(dabs(d(iq))+g==dabs(d(iq))))then
                          a(ip,iq)=0.
                      else if(dabs(a(ip,iq)).gt.tresh)then
                          h=d(iq)-d(ip)
                          if(dabs(h)+g==dabs(h))then
                              t=a(ip,iq)/h          !t = 1/(2*theta)
                          else
                              theta=0.5*h/a(ip,iq)  !equation (11.1.10).
                              t=1./(dabs(theta)+dsqrt(1.+theta**2))
                              if(theta.lt.0.)t=-t
                          endif
                          c=1./dsqrt(1+t**2)
                          s=t*c
                          tau=s/(1.+c)
                          h=t*a(ip,iq)
                          z(ip)=z(ip)-h
                          z(iq)=z(iq)+h
                          d(ip)=d(ip)-h
                          d(iq)=d(iq)+h
                          a(ip,iq)=0.
                          do j=1,ip-1  !case of rotations 1 = j < p.
                              g=a(j,ip)
                              h=a(j,iq)
                              a(j,ip)=g-s*(h+g*tau)
                              a(j,iq)=h+s*(g-h*tau)
                          enddo 
                          do j=ip+1,iq-1 !case of rotations p < j < q.
                              g=a(ip,j)
                              h=a(j,iq)
                              a(ip,j)=g-s*(h+g*tau)
                              a(j,iq)=h+s*(g-h*tau)
                          enddo 
                          do j=iq+1,n !case of rotations q < j = n.
                              g=a(ip,j)
                              h=a(iq,j)
                              a(ip,j)=g-s*(h+g*tau)
                              a(iq,j)=h+s*(g-h*tau)
                          enddo
                          do j=1,n
                              g=v(j,ip)
                              h=v(j,iq)
                              v(j,ip)=g-s*(h+g*tau)
                              v(j,iq)=h+s*(g-h*tau)
                          enddo 
                          nrot=nrot+1
                      endif
                  enddo 
              enddo 
              do ip=1,n
                  b(ip)=b(ip)+z(ip)
                  d(ip)=b(ip) !update d with the sum of tapq ,
                  z(ip)=0.    !and reinitialize z.
              enddo
          enddo 
          write(*,*) 'too many iterations in jacobi, give up'
          ising=1
          return
      end subroutine jacobi

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff
! DESCRIPTION: 
!> Creates Identity Matrix of dimension n
!> @param[in]   n                Dimension of Matrix
!> @param[out]  MI               Identity matrix
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine IdentityMatrix(n,MI)
!--------------------------------------------------------------------------
          implicit none
          integer,intent(in) :: n
          real(8),intent(out) :: MI(n,n)
          !Local variables
          integer :: i            !Loop integer
          !Start assembly
          MI=0.0
          do i=1,n
              MI(i,i)=1.0                        
          enddo
          return
      end subroutine IdentityMatrix

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Returns cross product
!> @param[in]   a,b              3x1 vectors
!> @param[out]  crossproduct    returns vector
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine crossproduct(abcross, a, b)
          implicit none
!--------------------------------------------------------------------------
          real(8), dimension(3), intent(out) :: abcross
          real(8), dimension(3), intent(in) :: a, b
          abcross(1) = a(2) * b(3) - a(3) * b(2)
          abcross(2) = a(3) * b(1) - a(1) * b(3)
          abcross(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine crossproduct

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Returns voigt index
!> @param[in]   o,p             Input indices
!> @param[in]   dof             Degree of freedom
!> @param[out]  voigt index
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function voigt(o,p,dof)
!--------------------------------------------------------------------------		      
      implicit none
          integer :: o, p, dof
          
          IF (o==1 .AND. p==1) THEN 
              voigt=1
          END IF
          IF (o==2 .AND. p==2) THEN 
              voigt=2
          END IF
          IF (o==3 .AND. p==3) THEN 
              voigt=3
          END IF
!         3D 
          IF (dof==3) THEN
              IF ((o==1 .AND. p==2) .OR. (o==2 .AND. p==1)) THEN 
                  voigt=6
              END IF
              IF ((o==2 .AND. p==3) .OR. (o==3 .AND. p==2)) THEN 
                  voigt=4
              END IF
              IF ((o==1 .AND. p==3) .OR. (o==3 .AND. p==1)) THEN 
                  voigt=5
              END IF
          END IF
!         2D 
          IF (dof==2) THEN
              IF ((o==1 .AND. p==2) .OR. (o==2 .AND. p==1)) THEN 
                  voigt=4
              END IF
              IF ((o==2 .AND. p==3) .OR. (o==3 .AND. p==2)) THEN 
                  voigt=666
              END IF
              IF ((o==1 .AND. p==3) .OR. (o==3 .AND. p==1)) THEN 
                  voigt=666
              END IF
          END IF
      end function voigt

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Returns of entry (o,p) delta matrix
!> @param[in]   a,b              3x1 vectors
!> @param[out]  delta            one/zero
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      real function delta(o,p)
!--------------------------------------------------------------------------		      
          implicit none
          integer, intent(in) :: o,p
          if (o==p) then 
              delta=1.0
          else
              delta=0.0
          end if
      end function delta

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Returns outer product
!> @param[in]   a,b              
!> @param[out]  cross product    returns vector
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!      real function outerprod(a,b)
!          Funktion äußeres Produkt
!          IMPLICIT none
!          REAL(8), DIMENSION(:), INTENT(IN) ::a,b
!          REAL(8), DIMENSION(SIZE(A),SIZE(B)) :: temp

!          temp(:,:) = SPREAD(a,dim=2,ncopies=SIZE(B))*SPREAD(b,dim=1,ncopies=SIZE(A))
!          outerprod = temp
!          return	  
!      end function outerprod
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> Returns determinante of 3x3 Matrix
!> @param[in]  matrix33	 3x3 Matrix              
!> @param[out]  det		 determinante
!--------------------------------------------------------------------------
      subroutine get_determinant(matrix33,det)
!--------------------------------------------------------------------------
          implicit none
          real(8), dimension(3,3), intent(in) :: matrix33
          real(8) :: v1, v2, v3
          real(8),intent(out) :: det
          v1 = matrix33(1,1)*(matrix33(2,2)*matrix33(3,3)-matrix33(2,3)*matrix33(3,2))
          v2 = matrix33(1,2)*(matrix33(2,1)*matrix33(3,3)-matrix33(2,3)*matrix33(3,1))
          v3 = matrix33(1,3)*(matrix33(2,1)*matrix33(3,2)-matrix33(2,2)*matrix33(3,1))
          det = v1-v2+v3
          return
      end subroutine get_determinant
   
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> Converts 3x3 matrix (symetric) to vector 
!> @param[in]   Am	 3x3 matrix              
!> @param[out]  Av	 Vector
!--------------------------------------------------------------------------
      subroutine icams_conv33to6(Am,Av)
!!--------------------------------------------------------------------------
          implicit none
          integer :: i
          integer, dimension(9) :: ib1
          integer, dimension(9) :: ib2
          real(8), dimension(3,3), intent(in) :: Am
          real(8), dimension(6), intent(out) :: Av

!         Abaqus stress vector sequence
          ib1=[1,2,3,1,1,2,2,3,3]
          ib2=[1,2,3,2,3,3,1,1,2]
          do i=1,6
              Av(i)=Am(ib1(i),ib2(i))
          enddo
          return
      end subroutine icams_conv33to6

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> Converts 3x3 matrix to vector (asymetric)
!> @param[in]  Am	 3x3 matrix              
!> @param[out]  Av	 Vector
!--------------------------------------------------------------------------
      subroutine icams_conv33to9(Am,Av)
!--------------------------------------------------------------------------
          implicit none
          integer :: i
          integer, dimension(9) :: ib1
          integer, dimension(9) :: ib2
          real(8), dimension(3,3), intent(in) :: Am
          real(8), dimension(9), intent(out) :: Av
          ib1=[1,2,3,1,1,2,2,3,3]
          ib2=[1,2,3,2,3,3,1,1,2]
          do i=1,9
              Av(i) = Am(ib1(i),ib2(i))
          enddo
          return
      end subroutine icams_conv33to9

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
! DESCRIPTION: 
!> Converts vector to 3x3 matrix (symmetric) 
!> @param[out]  Am	 3x3 matrix              
!> @param[in]  Av	 Vector
!--------------------------------------------------------------------------
      subroutine icams_conv6to33(Av,Am)
!--------------------------------------------------------------------------
          implicit none
          integer :: i
          integer, dimension(9) :: ib1
          integer, dimension(9) :: ib2
          real(8), dimension(3,3), intent(out) :: Am
          real(8), dimension(6), intent(in) :: Av
          ib1=[1,2,3,1,1,2,2,3,3] 
          ib2=[1,2,3,2,3,3,1,1,2]
          do i=1,6
              Am(ib1(i),ib2(i))=Av(i)
          enddo
          do i=7,9
              Am(ib1(i),ib2(i))=Av(i-3)
          enddo
          return
      end subroutine icams_conv6to33

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @author 
   !> Anxin Ma
   !> \f$ Z_1 X_2 Z_3 \f$ convention from "http://en.wikipedia.org/wiki/Euler_angles" \n
   !> \f$ c_i=cos(...) \f$ and \f$ s_i=sin(...) \f$
   ! DESCRIPTION: 
   !> Calculates based on the Euler angles the Rotation Matrix
   !> @param[in]  phi1               first euler angle \f$ \Psi \f$
   !> @param[in]  PHI                euler angle \f$ \Theta \f$
   !> @param[in]  phi2               second euler angle \f$ \Phi \f$
   !> @param[out] QM                 Rotation matrix
   !> \image html EulerAngles.png
   !> \image html EulerAnglesImage.png
!--------------------------------------------------------------------------
      subroutine icams_Eang2Q(phi1,PHI,phi2,QM)
!--------------------------------------------------------------------------
!         rotate from X[100],Y[010],Z[001] to v1,v2,v3
          
          implicit none
          real(8), intent(in):: phi1
          real(8), intent(in):: PHI
          real(8), intent(in):: phi2
          real(8), intent(out):: QM(3,3)
          !Local parameters
          real(8) xp1,xp,xp2
          real(8) c1,c2,c3,s1,s2,s3
          real(8) pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180.0/pi
          g2r=pi/180.0
          xp1=phi1*g2r
          xp =PHI *g2r
          xp2=phi2*g2r
          c1=dcos(xp1)
          s1=dsin(xp1)
          c2=dcos(xp )
          s2=dsin(xp )
          s3=dsin(xp2)
          c3=dcos(xp2)
          QM(1,1)=+c1*c3-s1*s3*c2
          QM(1,2)=+s1*c3+c1*s3*c2
          QM(1,3)=+s3*s2
          QM(2,1)=-c1*s3-s1*c3*c2
          QM(2,2)=-s1*s3+c1*c3*c2
          QM(2,3)=+c3*s2
          QM(3,1)=+s1*s2
          QM(3,2)=-c1*s2
          QM(3,3)=+c2
          return
      end subroutine icams_Eang2Q

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
!> \f$ Z_1 X_2 Z_3 \f$ convention from "http://en.wikipedia.org/wiki/Euler_angles" \n
!> \f$ c_i=cos(...) \f$ and \f$ s_i=sin(...) \f$
! DESCRIPTION: 
!> ???
!> @param[in]  phi1               first euler angle \f$ \Psi \f$
!> @param[in]  PHI                euler angle \f$ \Theta \f$
!> @param[in]  phi2               second euler angle \f$ \Phi \f$
!> @param[out] QM                 Rotation matrix
!> \image html EulerAngles.png
!> \image html EulerAnglesImage.png
!--------------------------------------------------------------------------
      subroutine icams_angax2QM(ang,u,v,w,QM)
!--------------------------------------------------------------------------
          implicit none
          real(8) :: QM(3,3)
          real(8) :: s,c,u2,v2,w2,ang,u,v,w,x1
          real(8) :: pi,r2g,g2r
          pi=dacos(-1.d0)
          r2g=180.0/pi
          g2r=pi/180.0
          x1=dsqrt(u**2+v**2+w**2)
          u2=u/x1
          v2=v/x1
          w2=w/x1
          s=dsin(ang)
          c=dcos(ang)
          QM(1,1)=(1.0-u2**2)*c+u2**2
          QM(2,2)=(1.0-v2**2)*c+v2**2
          QM(3,3)=(1.0-w2**2)*c+w2**2
          QM(1,2)=u2*v2*(1.0-c)+w2*s
          QM(2,1)=u2*v2*(1.0-c)-w2*s
          QM(1,3)=u2*w2*(1.0-c)-v2*s
          QM(3,1)=u2*w2*(1.0-c)+v2*s
          QM(2,3)=v2*w2*(1.0-c)+u2*s
          QM(3,2)=v2*w2*(1.0-c)-u2*s
          return
      end subroutine icams_angax2QM

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Anxin Ma
!> convention from "http://en.wikipedia.org/wiki/Euler_angles" \n
! DESCRIPTION: 
!> Calculates based on the Rotation Matrix the Euler angles
!> @param[in]   QM               Rotation matrix
!> @param[out]  p1               first euler angle \f$ \Psi \f$
!> @param[out]  p                euler angle \f$ \Theta \f$
!> @param[out]  p2               second euler angle \f$ \Phi \f$
!--------------------------------------------------------------------------
      subroutine icams_Q2Eang(QM,phi1,PHI,phi2)
!--------------------------------------------------------------------------
          implicit none
          real(8), intent(in) :: QM(3,3)
          real(8), intent(out):: phi1,PHI,phi2
!---------Local parameters
          real(8) sqhkl,squvw,sqhk,val
          real(8) pi,r2g,g2r,Tol
          Tol=1.d-15
          pi=dacos(-1.d0)
          r2g=180.0/pi
          g2r=pi/180.0
!---------------------------------------------
!             v1   v2    v3
!
!           | u   v2_1    h  |
!     QM =  | v   v2_2    k  |  with v2 = v3 x v1   
!           | w   v2_3    l  | 
!
!             100   010   001
!          X|  u   v2_1    h  |
!     QM = Y|  v   v2_2    k  |    
!          Z|  w   v2_3    l  | 
!---------------------------------------------
!    if the roation tensor is defined as following
!    QM_ij:= (dx/dX)_ij = dx_i/dX_j
!
!               X     Y     X
!          v1|  u   v2_1    h  |
!     QM'= v2|  v   v2_2    k  |    
!          v3|  w   v2_3    l  | 
!    the ratation must be take 
!
!    QM=transpose(QM')
!---------------------------------------------
          squvw=dsqrt(QM(1,1)**2+QM(2,1)**2+QM(3,1)**2)
          sqhkl=dsqrt(QM(1,3)**2+QM(2,3)**2+QM(3,3)**2)
          sqhk =dsqrt(QM(1,3)**2+QM(2,3)**2           )
          
          val=QM(3,3)/sqhkl
          if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
          PHI=dacos(val)
          
          if(PHI < TOL) then
              phi2=0.0
              val=QM(1,1)/squvw
              if(QM(2,1) <= 0.d0) then
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
                  phi1=dacos(val)
              else
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
!                 phi1=2*pi-dacos(val)
                  phi1=-dacos(val)
              endif
          else
!             val=QM(2,3)/sqhk
              val=QM(2,3)/dsin(PHI)
              if(QM(1,3) >= 0.d0) then
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
                  phi2=dacos(val)
              else
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
!                 phi2=2*pi-dacos(val)
                  phi2=-dacos(val)
              endif
              val=-QM(3,2)/dsin(PHI)
              if(QM(3,1) >= 0.d0) then
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
                  phi1=dacos(val)
              else
                  if(dabs(val)>1.d0)val=1.d0*sign(1.d0,val)
!                 phi1=2*pi-dacos(val)
                  phi1=-dacos(val)
              endif
          endif
          phi1=phi1*r2g
          PHI=PHI*r2g
          phi2=phi2*r2g
          return
      end subroutine icams_Q2Eang

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> returns entry of permutation matrix (a.k.a. levi-civita operator)
!> @param[in]   i,j,k            Indices
!> @param[out]  levicivita       
!--------------------------------------------------------------------------
      real function levicivita(i,j,k)
!--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: i,j,k
          if (i == 1 .and. j == 2 .and. k == 3) then
              levicivita = 1.0
          elseif (i == 3 .and. j == 1 .and.k == 2) then
              levicivita = 1.0
          elseif (i == 2 .and. j == 3 .and. k == 1) then
              levicivita = 1.0
          elseif (i == 1 .and. j == 3 .and. k == 2) then
              levicivita = -1.0
          elseif (i == 3 .and. j == 2 .and. k == 1) then
              levicivita = -1.0
          elseif (i == 2 .and. j == 1 .and. k == 3) then
              levicivita = -1.0
          else 
              levicivita = 0.0
          end if
      end function levicivita

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels, Anxin Ma
! DESCRIPTION: 
!> returns Frobenius norm of matrix 
!> @param[in]   matrix    Matrix
!> @param[in]   rows
!> @param[in]   columns
!> @param[out]  norm      Frobenius norm    
!--------------------------------------------------------------------------
      subroutine icams_matrixnorm(matrix,rows,columns,norm)
!--------------------------------------------------------------------------
          implicit none
          INTEGER :: I,J,rows, columns
          real(8), intent(in) :: matrix(rows,columns)
          real(8), intent(out) :: norm  
          norm = 0.d0
          do i=1,rows
              do j=1,columns 
                  norm=norm+matrix(i,j)**2.d0 
              end do 
          end dO
          norm=dsqrt(norm)
          return 
      end subroutine icams_matrixnorm
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
