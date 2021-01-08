!     
!     3lpt initial condition with transverse mode
!     by Takayuki Tatekawa
!

      Implicit None
      Integer Nx, Nparticles
      parameter(Nx=256)
      parameter(Nparticles=Nx**3)

      REAL*8 dx,epsilon,astart,omegam,omegav,h0,dt,etat
      REAL*8 ekin,egrav,egint
      REAL*8 a1
      REAL xp(Nparticles),yp(Nparticles),zp(Nparticles)
      REAL xx(3,Nparticles),vv(3,Nparticles)
      REAL vxp(Nparticles),vyp(Nparticles),vzp(Nparticles)
      REAL*8 S1(3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S1d(3,3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S2(3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S2d(3,3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S3A(3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S3B(3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 Psi2(0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 Psi3A(0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 Psi3B(0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 S3T1,S3T2
      REAL*8 S3T(3,0:Nx+1,0:Nx+1,0:Nx+1)
      REAL*8 dat(0:2*Nparticles-1)
      REAL*8 dr(0:Nx-1,0:Nx-1,0:Nx-1), di(0:Nx-1, 0:Nx-1, 0:Nx-1)
      REAL*8 CC,DD,EE,EEt,FFa,FFb,FFat,FFbt
      REAL*8 FFtb, FFtbt
      REAL*8 x1,x2,x3,divS1,divS2
      REAL*8 Pi2,kkf
      character*100 filein,fileout
      INTEGER np1,np2,np3,nstep,ii,j,j1,j2,j3,k,aa,FL,nn(3),mode1
      integer Npart(6),FlagSfr,FlagFeedback,FlagCooling
      integer FlagStellarAge,FlagMetals,NallHW(6),flagentrics,Nall(6)
      integer NpF,NumFiles
      real unused(15),hfactor,velfactor
      real*8 Massarr(6),afactor,redshift
      real*8 Omega0,OmegaL0,hlittle
      real*8 Om, OL, D0
      real*8 RBox

      filein='ZA-init.dat'

      fileout='3lpt-init.dat'

      write (*,*) "Please choose mode. 1: 2LPT,
     &  2: 3LPT (longitudinal mode only), 3: 3LPT"
      read (*,*) mode1

      open(unit=4,file=filein,status='unknown',form='unformatted')
      read(4)Npart,Massarr,afactor,redshift,FlagSfr,FlagFeedback,
     &        Nall,FlagCooling,NumFiles,RBox,Omega0,OmegaL0,
     &        hlittle,FlagStellarAge,FlagMetals,NallHW,flagentrics,
     &        unused

      read(4)(xp(k),yp(k),zp(k),k=1,Npart(2))
      read(4)(vxp(k),vyp(k),vzp(k),k=1,Npart(2))
      read(4)(k,k=1,Npart(2))
      close(4)

      np1=INT((Npart(2)+1)**(1.0d0/3.0d0))
      np2=np1
      np3=np1
      nn(1)=Nx
      nn(2)=Nx
      nn(3)=Nx
      dx=dble(Rbox)/dble(np1)

      Pi2=ATAN(1.0d0)*8.0d0

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  j=j+1
                  S1(1,j1,j2,j3)=(xp(j)/dx-DBLE(j1-1))
                  S1(2,j1,j2,j3)=(yp(j)/dx-DBLE(j2-1))
                  S1(3,j1,j2,j3)=(zp(j)/dx-DBLE(j3-1))

                  IF (S1(1,j1,j2,j3).gt.1.0d1) THEN
                     S1(1,j1,j2,j3)=S1(1,j1,j2,j3)-DBLE(np1)
                  ELSE IF (S1(1,j1,j2,j3).lt.(-1.0d1)) THEN
                     S1(1,j1,j2,j3)=S1(1,j1,j2,j3)+DBLE(np1)
                  ENDIF

                  IF (S1(2,j1,j2,j3).gt.1.0d1) THEN
                     S1(2,j1,j2,j3)=S1(2,j1,j2,j3)-DBLE(np2)
                  ELSE IF (S1(2,j1,j2,j3).lt.(-1.0d1)) THEN
                     S1(2,j1,j2,j3)=S1(2,j1,j2,j3)+DBLE(np2)
                  ENDIF

                  IF (S1(3,j1,j2,j3).gt.1.0d1) THEN
                     S1(3,j1,j2,j3)=S1(3,j1,j2,j3)-DBLE(np3)
                  ELSE IF (S1(3,j1,j2,j3).lt.(-1.0d1)) THEN
                     S1(3,j1,j2,j3)=S1(3,j1,j2,j3)+DBLE(np3)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         CC=vxp(1)/S1(1,1,1,1)

!
!     2nd order setup
!

         DO j2=1,np2
            DO j3=1,np3
               S1(1,0,j2,j3)=S1(1,np1,j2,j3)
               S1(2,0,j2,j3)=S1(2,np1,j2,j3)
               S1(3,0,j2,j3)=S1(3,np1,j2,j3)
               S1(1,np1+1,j2,j3)=S1(1,1,j2,j3)
               S1(2,np1+1,j2,j3)=S1(2,1,j2,j3)
               S1(3,np1+1,j2,j3)=S1(3,1,j2,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
               S1(1,j1,0,j3)=S1(1,j1,np2,j3)
               S1(2,j1,0,j3)=S1(2,j1,np2,j3)
               S1(3,j1,0,j3)=S1(3,j1,np2,j3)
               S1(1,j1,np2+1,j3)=S1(1,j1,1,j3)
               S1(2,j1,np2+1,j3)=S1(2,j1,1,j3)
               S1(3,j1,np2+1,j3)=S1(3,j1,1,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
               S1(1,j1,j2,0)=S1(1,j1,j2,np3)
               S1(2,j1,j2,0)=S1(2,j1,j2,np3)
               S1(3,j1,j2,0)=S1(3,j1,j2,np3)
               S1(1,j1,j2,np3+1)=S1(1,j1,j2,1)
               S1(2,j1,j2,np3+1)=S1(2,j1,j2,1)
               S1(3,j1,j2,np3+1)=S1(3,j1,j2,1)
            ENDDO
         ENDDO

!
!     2nd order
!

         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  divS1=5.0d-1*
     &                 (S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3)
     &                 +S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3)
     &                 +S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))

                  Psi2(j1,j2,j3)=5.0d-1*(divS1*divS1
     &                 -( (S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3))
     &                   * (S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3))
     &                   + (S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3))
     &                   * (S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3))
     &                   + (S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))
     &                   * (S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))
     &                   -2.0d0*( (S1(1,j1,j2+1,j3)-S1(1,j1,j2-1,j3))
     &                   * (S1(1,j1,j2+1,j3)-S1(1,j1,j2-1,j3))
     &                   + (S1(2,j1,j2,j3+1)-S1(2,j1,j2,j3-1))
     &                   * (S1(2,j1,j2,j3+1)-S1(2,j1,j2,j3-1))
     &                   + (S1(3,j1+1,j2,j3)-S1(3,j1-1,j2,j3))
     &                   * (S1(3,j1+1,j2,j3)-S1(3,j1-1,j2,j3)))) )
               ENDDO
            ENDDO
         ENDDO

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  dat(2*j)=Psi2(j1,j2,j3)
                  dat(2*j+1)=0.0d0
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, 1)

         j=0
         DO j3=0,np3-1
            DO j2=0,np2-1
               DO j1=0,np1-1

                  kkf=-2.0d0*(cos(Pi2*DBLE(j1)/DBLE(np1))
     &                 +cos(Pi2*DBLE(j2)/DBLE(np2))
     &                 +cos(Pi2*DBLE(j3)/DBLE(np3))-3.0d0)

                  IF (kkf.ne.0.0d0) THEN
                     dat(2*j)=dat(2*j)/kkf
                     dat(2*j+1)=dat(2*j+1)/kkf
                  ENDIF         

                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, -1)

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  Psi2(j1,j2,j3)=dat(j*2) *(1.0d0/DBLE(Nx))**3
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         DO j2=1,np2
            DO j3=1,np3
               Psi2(0,j2,j3)=Psi2(np1,j2,j3)
               Psi2(np1+1,j2,j3)=Psi2(1,j2,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
               Psi2(j1,0,j3)=Psi2(j1,np2,j3)
               Psi2(j1,np2+1,j3)=Psi2(j1,1,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
               Psi2(j1,j2,0)=Psi2(j1,j2,np3)
               Psi2(j1,j2,np3+1)=Psi2(j1,j2,1)
            ENDDO
         ENDDO

         DO j1=1,np1
            DO j2=1,np1
               DO j3=1,np1
                  S2(1,j1,j2,j3)=
     &                 (Psi2(j1+1,j2,j3)-Psi2(j1-1,j2,j3))*5.0d-1
                  S2(2,j1,j2,j3)=
     &                 (Psi2(j1,j2+1,j3)-Psi2(j1,j2-1,j3))*5.0d-1
                  S2(3,j1,j2,j3)=
     &                 (Psi2(j1,j2,j3+1)-Psi2(j1,j2,j3-1))*5.0d-1
               ENDDO
            ENDDO
         ENDDO

!
!     3rd order setup
!

         DO j2=1,np2
            DO j3=1,np3
               S2(1,0,j2,j3)=S2(1,np1,j2,j3)
               S2(2,0,j2,j3)=S2(2,np1,j2,j3)
               S2(3,0,j2,j3)=S2(3,np1,j2,j3)
               S2(1,np1+1,j2,j3)=S2(1,1,j2,j3)
               S2(2,np1+1,j2,j3)=S2(2,1,j2,j3)
               S2(3,np1+1,j2,j3)=S2(3,1,j2,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
               S2(1,j1,0,j3)=S2(1,j1,np2,j3)
               S2(2,j1,0,j3)=S2(2,j1,np2,j3)
               S2(3,j1,0,j3)=S2(3,j1,np2,j3)
               S2(1,j1,np2+1,j3)=S2(1,j1,1,j3)
               S2(2,j1,np2+1,j3)=S2(2,j1,1,j3)
               S2(3,j1,np2+1,j3)=S2(3,j1,1,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
               S2(1,j1,j2,0)=S2(1,j1,j2,np3)
               S2(2,j1,j2,0)=S2(2,j1,j2,np3)
               S2(3,j1,j2,0)=S2(3,j1,j2,np3)
               S2(1,j1,j2,np3+1)=S2(1,j1,j2,1)
               S2(2,j1,j2,np3+1)=S2(2,j1,j2,1)
               S2(3,j1,j2,np3+1)=S2(3,j1,j2,1)
            ENDDO
         ENDDO
         
!
!     3rd order A
!

         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1

                  Psi3A(j1,j2,j3)=
     &                 ((S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3))*
     &                 (S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))
     &                 -(S1(2,j1,j2,j3+1)-S1(2,j1,j2,j3-1))*
     &                 (S1(3,j1,j2+1,j3)-S1(3,j1,j2-1,j3)))
     &                 *(S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3)) *1.25d-1
     &                 +((S1(2,j1,j2,j3+1)-S1(2,j1,j2,j3-1))*
     &                 (S1(3,j1+1,j2,j3)-S1(3,j1-1,j2,j3))
     &                 -(S1(2,j1+1,j2,j3)-S1(2,j1-1,j2,j3))*
     &                 (S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1)))
     &                 *(S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3)) *1.25d-1
     &                 +((S1(2,j1+1,j2,j3)-S1(2,j1-1,j2,j3))*
     &                 (S1(3,j1,j2+1,j3)-S1(3,j1,j2-1,j3))
     &                 -(S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3))*
     &                 (S1(3,j1+1,j2,j3)-S1(3,j1-1,j2,j3)))
     &                 *(S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1)) *1.25d-1

               ENDDO
            ENDDO
         ENDDO

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  dat(2*j)=Psi3A(j1,j2,j3)
                  dat(2*j+1)=0.0d0
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, 1)

         j=0
         DO j3=0,np3-1
            DO j2=0,np2-1
               DO j1=0,np1-1

                  kkf=-2.0d0*(cos(Pi2*DBLE(j1)/DBLE(np1))
     &                 +cos(Pi2*DBLE(j2)/DBLE(np2))
     &                 +cos(Pi2*DBLE(j3)/DBLE(np3))-3.0d0)

                  IF (kkf.ne.0.0d0) THEN
                     dat(2*j)=dat(2*j)/kkf
                     dat(2*j+1)=dat(2*j+1)/kkf
                  ENDIF         

                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, -1)

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  Psi3A(j1,j2,j3)=dat(j*2)*(1.0d0/DBLE(Nx))**3
                 j=j+1
               ENDDO
            ENDDO
         ENDDO

         DO j2=1,np2
            DO j3=1,np3
               Psi3A(0,j2,j3)=Psi3A(np1,j2,j3)
               Psi3A(np1+1,j2,j3)=Psi3A(1,j2,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
               Psi3A(j1,0,j3)=Psi3A(j1,np2,j3)
               Psi3A(j1,np2+1,j3)=Psi3A(j1,1,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
               Psi3A(j1,j2,0)=Psi3A(j1,j2,np3)
               Psi3A(j1,j2,np3+1)=Psi3A(j1,j2,1)
            ENDDO
         ENDDO

         DO j1=1,np1
            DO j2=1,np1
               DO j3=1,np1
                  S3A(1,j1,j2,j3)=
     &                 (Psi3A(j1+1,j2,j3)-Psi3A(j1-1,j2,j3))*5.0d-1
                  S3A(2,j1,j2,j3)=
     &                 (Psi3A(j1,j2+1,j3)-Psi3A(j1,j2-1,j3))*5.0d-1
                  S3A(3,j1,j2,j3)=
     &                 (Psi3A(j1,j2,j3+1)-Psi3A(j1,j2,j3-1))*5.0d-1
               ENDDO
            ENDDO
         ENDDO

!
!     3rd order B
!

         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  divS1=5.0d-1*
     &                 (S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3)
     &                 +S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3)
     &                 +S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))
                  divS2=5.0d-1*
     &                 (S2(1,j1+1,j2,j3)-S2(1,j1-1,j2,j3)
     &                 +S2(2,j1,j2+1,j3)-S2(2,j1,j2-1,j3)
     &                 +S2(3,j1,j2,j3+1)-S2(3,j1,j2,j3-1))

                  Psi3B(j1,j2,j3)=5.0d-1*(divS1*divS2
     &                 -( (S1(1,j1+1,j2,j3)-S1(1,j1-1,j2,j3))
     &                   * (S2(1,j1+1,j2,j3)-S2(1,j1-1,j2,j3))
     &                   + (S1(2,j1,j2+1,j3)-S1(2,j1,j2-1,j3))
     &                   * (S2(2,j1,j2+1,j3)-S2(2,j1,j2-1,j3))
     &                   + (S1(3,j1,j2,j3+1)-S1(3,j1,j2,j3-1))
     &                   * (S2(3,j1,j2,j3+1)-S2(3,j1,j2,j3-1))
     &                   -2.0d0*((S1(1,j1,j2+1,j3)-S1(1,j1,j2-1,j3))
     &                   * (S2(1,j1,j2+1,j3)-S2(1,j1,j2-1,j3))
     &                   + (S1(2,j1,j2,j3+1)-S1(2,j1,j2,j3-1))
     &                   * (S2(2,j1,j2,j3+1)-S2(2,j1,j2,j3-1))
     &                   + (S1(3,j1+1,j2,j3)-S1(3,j1-1,j2,j3))
     &                   * (S2(3,j1+1,j2,j3)-S2(3,j1-1,j2,j3)))) )

               ENDDO
            ENDDO
         ENDDO

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  dat(2*j)=Psi3B(j1,j2,j3)
                  dat(2*j+1)=0.0d0
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, 1)

         j=0
         DO j3=0,np3-1
            DO j2=0,np2-1
               DO j1=0,np1-1

                  kkf=-2.0d0*(cos(Pi2*DBLE(j1)/DBLE(np1))
     &                 +cos(Pi2*DBLE(j2)/DBLE(np2))
     &                 +cos(Pi2*DBLE(j3)/DBLE(np3))-3.0d0)

                  IF (kkf.ne.0.0d0) THEN
                     dat(2*j)=dat(2*j)/kkf
                     dat(2*j+1)=dat(2*j+1)/kkf
                  ENDIF         

                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, -1)

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  Psi3B(j1,j2,j3)=dat(j*2)*(1.0d0/DBLE(Nx))**3
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         DO j2=1,np2
            DO j3=1,np3
               Psi3B(0,j2,j3)=Psi3B(np1,j2,j3)
               Psi3B(np1+1,j2,j3)=Psi3B(1,j2,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
               Psi3B(j1,0,j3)=Psi3B(j1,np2,j3)
               Psi3B(j1,np2+1,j3)=Psi3B(j1,1,j3)
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
               Psi3B(j1,j2,0)=Psi3B(j1,j2,np3)
               Psi3B(j1,j2,np3+1)=Psi3B(j1,j2,1)
            ENDDO
         ENDDO

         DO j1=1,np1
            DO j2=1,np1
               DO j3=1,np1
!                  S3B(1,j1,j2,j3)=
!     &                 (Psi3B(j1+1,j2,j3)-Psi3B(j1-1,j2,j3))/dx*5.0d-1
!                  S3B(2,j1,j2,j3)=
!     &                 (Psi3B(j1,j2+1,j3)-Psi3B(j1,j2-1,j3))/dx*5.0d-1
!                  S3B(3,j1,j2,j3)=
!     &                 (Psi3B(j1,j2,j3+1)-Psi3B(j1,j2,j3-1))/dx*5.0d-1
                  S3B(1,j1,j2,j3)=
     &                 (Psi3B(j1+1,j2,j3)-Psi3B(j1-1,j2,j3))/dx*5.0d-1
                  S3B(2,j1,j2,j3)=
     &                 (Psi3B(j1,j2+1,j3)-Psi3B(j1,j2-1,j3))/dx*5.0d-1
                  S3B(3,j1,j2,j3)=
     &                 (Psi3B(j1,j2,j3+1)-Psi3B(j1,j2,j3-1))/dx*5.0d-1
               ENDDO
            ENDDO
         ENDDO

!
!       3rd order transverse mode setup
!

         WRITE (*,*) '3rd T'

         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  DO ii=1,3

                  S1d(ii,1,j1,j2,j3)=
     &        (S1(ii,j1+1,j2,j3)-S1(ii,j1-1,j2,j3))*5.0d-1
                  S1d(ii,2,j1,j2,j3)=
     &        (S1(ii,j1,j2+1,j3)-S1(ii,j1,j2-1,j3))*5.0d-1
                  S1d(ii,3,j1,j2,j3)=
     &        (S1(ii,j1,j2,j3+1)-S1(ii,j1,j2,j3-1))*5.0d-1

                  S2d(ii,1,j1,j2,j3)=
     &        (S2(ii,j1+1,j2,j3)-S2(ii,j1-1,j2,j3))*5.0d-1
                  S2d(ii,2,j1,j2,j3)=
     &        (S2(ii,j1,j2+1,j3)-S2(ii,j1,j2-1,j3))*5.0d-1
                  S2d(ii,3,j1,j2,j3)=
     &        (S2(ii,j1,j2,j3+1)-S2(ii,j1,j2,j3-1))*5.0d-1

               ENDDO

               ENDDO
            ENDDO
         ENDDO

         DO j2=1,np2
            DO j3=1,np3
                  DO ii=1,3
               S1d(ii,1,0,j2,j3)=S1d(ii,1,np1,j2,j3)
               S1d(ii,2,0,j2,j3)=S1d(ii,2,np1,j2,j3)
               S1d(ii,3,0,j2,j3)=S1d(ii,3,np1,j2,j3)
               S1d(ii,1,np1+1,j2,j3)=S1d(ii,1,1,j2,j3)
               S1d(ii,2,np1+1,j2,j3)=S1d(ii,2,1,j2,j3)
               S1d(ii,3,np1+1,j2,j3)=S1d(ii,3,1,j2,j3)
               S2d(ii,1,0,j2,j3)=S2d(ii,1,np1,j2,j3)
               S2d(ii,2,0,j2,j3)=S2d(ii,2,np1,j2,j3)
               S2d(ii,3,0,j2,j3)=S2d(ii,3,np1,j2,j3)
               S2d(ii,1,np1+1,j2,j3)=S2d(ii,1,1,j2,j3)
               S2d(ii,2,np1+1,j2,j3)=S2d(ii,2,1,j2,j3)
               S2d(ii,3,np1+1,j2,j3)=S2d(ii,3,1,j2,j3)
               ENDDO
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j3=0,np3
                  DO ii=1,3
               S1d(ii,1,j1,0,j3)=S1d(ii,1,j1,np2,j3)
               S1d(ii,2,j1,0,j3)=S1d(ii,2,j1,np2,j3)
               S1d(ii,3,j1,0,j3)=S1d(ii,3,j1,np2,j3)
               S1d(ii,1,j1,np2+1,j3)=S1d(ii,1,j1,1,j3)
               S1d(ii,2,j1,np2+1,j3)=S1d(ii,2,j1,1,j3)
               S1d(ii,3,j1,np2+1,j3)=S1d(ii,3,j1,1,j3)
               S2d(ii,1,j1,0,j3)=S2d(ii,1,j1,np2,j3)
               S2d(ii,2,j1,0,j3)=S2d(ii,2,j1,np2,j3)
               S2d(ii,3,j1,0,j3)=S2d(ii,3,j1,np2,j3)
               S2d(ii,1,j1,np2+1,j3)=S2d(ii,1,j1,1,j3)
               S2d(ii,2,j1,np2+1,j3)=S2d(ii,2,j1,1,j3)
               S2d(ii,3,j1,np2+1,j3)=S2d(ii,3,j1,1,j3)
               ENDDO
            ENDDO
         ENDDO

         DO j1=0,np1+1
            DO j2=0,np2+1
                  DO ii=1,3
                     S1d(ii,1,j1,j2,0)=S1d(ii,1,j1,j2,np3)
                     S1d(ii,2,j1,j2,0)=S1d(ii,2,j1,j2,np3)
                     S1d(ii,3,j1,j2,0)=S1d(ii,3,j1,j2,np3)
                     S1d(ii,1,j1,j2,np3+1)=S1d(ii,1,j1,j2,1)
                     S1d(ii,2,j1,j2,np3+1)=S1d(ii,2,j1,j2,1)
                     S1d(ii,3,j1,j2,np3+1)=S1d(ii,3,j1,j2,1)
                     S2d(ii,1,j1,j2,0)=S2d(ii,1,j1,j2,np3)
                     S2d(ii,2,j1,j2,0)=S2d(ii,2,j1,j2,np3)
                     S2d(ii,3,j1,j2,0)=S2d(ii,3,j1,j2,np3)
                     S2d(ii,1,j1,j2,np3+1)=S2d(ii,1,j1,j2,1)
                     S2d(ii,2,j1,j2,np3+1)=S2d(ii,2,j1,j2,1)
                     S2d(ii,3,j1,j2,np3+1)=S2d(ii,3,j1,j2,1)
                  ENDDO
            ENDDO
         ENDDO

!
!       3rd order transverse mode
!

         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  DO ii=1,3

                     S3T1
     &   = ((S1d(ii,1,j1+1,j2,j3)*S2d(1,1,j1+1,j2,j3))
     &     -(S1d(ii,1,j1-1,j2,j3)*S2d(1,1,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,2,j1+1,j2,j3)*S2d(1,2,j1+1,j2,j3))
     &     -(S1d(ii,2,j1-1,j2,j3)*S2d(1,2,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,3,j1+1,j2,j3)*S2d(1,3,j1+1,j2,j3))
     &     -(S1d(ii,3,j1-1,j2,j3)*S2d(1,3,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,1,j1,j2+1,j3)*S2d(2,1,j1,j2+1,j3))
     &     -(S1d(ii,1,j1,j2-1,j3)*S2d(2,1,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,2,j1,j2+1,j3)*S2d(2,2,j1,j2+1,j3))
     &     -(S1d(ii,2,j1,j2-1,j3)*S2d(2,2,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,3,j1,j2+1,j3)*S2d(2,3,j1,j2+1,j3))
     &     -(S1d(ii,3,j1,j2-1,j3)*S2d(2,3,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(ii,1,j1,j2,j3+1)*S2d(3,1,j1,j2,j3+1))
     &     -(S1d(ii,1,j1,j2,j3-1)*S2d(3,1,j1,j2,j3-1)))
     &      *5.0d-1
     &   + ((S1d(ii,2,j1,j2,j3+1)*S2d(3,2,j1,j2,j3+1))
     &     -(S1d(ii,2,j1,j2,j3-1)*S2d(3,2,j1,j2,j3-1)))
     &      *5.0d-1
     &   + ((S1d(ii,3,j1,j2,j3+1)*S2d(3,3,j1,j2,j3+1))
     &     -(S1d(ii,3,j1,j2,j3-1)*S2d(3,3,j1,j2,j3-1)))
     &      *5.0d-1

                     S3T2
     &   = ((S1d(1,1,j1+1,j2,j3)*S2d(ii,1,j1+1,j2,j3))
     &     -(S1d(1,1,j1-1,j2,j3)*S2d(ii,1,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(1,2,j1+1,j2,j3)*S2d(ii,2,j1+1,j2,j3))
     &     -(S1d(1,2,j1-1,j2,j3)*S2d(ii,2,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(1,3,j1+1,j2,j3)*S2d(ii,3,j1+1,j2,j3))
     &     -(S1d(1,3,j1-1,j2,j3)*S2d(ii,3,j1-1,j2,j3)))
     &      *5.0d-1
     &   + ((S1d(2,1,j1,j2+1,j3)*S2d(ii,1,j1,j2+1,j3))
     &     -(S1d(2,1,j1,j2-1,j3)*S2d(ii,1,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(2,2,j1,j2+1,j3)*S2d(ii,2,j1,j2+1,j3))
     &     -(S1d(2,2,j1,j2-1,j3)*S2d(ii,2,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(2,3,j1,j2+1,j3)*S2d(ii,3,j1,j2+1,j3))
     &     -(S1d(2,3,j1,j2-1,j3)*S2d(ii,3,j1,j2-1,j3)))
     &      *5.0d-1
     &   + ((S1d(3,1,j1,j2+1,j3)*S2d(ii,1,j1,j2+1,j3))
     &     -(S1d(3,1,j1,j2,j3-1)*S2d(ii,1,j1,j2,j3-1)))
     &      *5.0d-1
     &   + ((S1d(3,2,j1,j2+1,j3)*S2d(ii,2,j1,j2+1,j3))
     &     -(S1d(3,2,j1,j2,j3-1)*S2d(ii,2,j1,j2,j3-1)))
     &      *5.0d-1
     &   + ((S1d(3,3,j1,j2+1,j3)*S2d(ii,3,j1,j2+1,j3))
     &     -(S1d(3,3,j1,j2,j3-1)*S2d(ii,3,j1,j2,j3-1)))
     &      *5.0d-1

                  S3T(ii,j1,j2,j3)=S3T1-S3T2

               ENDDO
               ENDDO
            ENDDO
         ENDDO

         DO ii=1,3

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  dat(2*j)=S3T(ii,j1,j2,j3)
                  dat(2*j+1)=0.0d0
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, 1)

         j=0
         DO j3=0,np3-1
            DO j2=0,np2-1
               DO j1=0,np1-1

                  kkf=-2.0d0*(cos(Pi2*DBLE(j1)/DBLE(np1))
     &                 +cos(Pi2*DBLE(j2)/DBLE(np2))
     &                 +cos(Pi2*DBLE(j3)/DBLE(np3))-3.0d0)
!     &                 *((DBLE(np1)/Pi2)**2)

!                  kkf=-((DBLE(np1)/Pi2)**2)
!     &                 *(DBLE(j1*j1+j2*j2+j3*j3))

                  IF (kkf.ne.0.0d0) THEN
                     dat(2*j)=dat(2*j)/kkf
                     dat(2*j+1)=dat(2*j+1)/kkf
                  ENDIF         

                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         CALL fourn (dat, nn, 3, -1)

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  S3T(ii,j1,j2,j3)=dat(j*2)*(1.0d0/DBLE(Nx))**3
!                  S3T(ii,j1,j2,j3)=dat(j*2)
                  j=j+1
               ENDDO
            ENDDO
         ENDDO

         ENDDO

!
!       3rd order output
!     
         


         EE=-3.0d0/7.0d0
         EEt=-6.0d0/7.0d0

         IF (mode1.ge.2) THEN
            FFa=-1.0d0/3.0d0
            FFat=-1.0d0
            FFb=1.0d1/2.1d1
            FFbt=1.0d1/7.0d0
         ELSE
            FFa=0.0d0
            FFat=0.0d0
            FFb=0.0d0
            FFbt=0.0d0
         ENDIF

         IF (mode1.ge.3) THEN
            FFtb=-3.0d0/9.8d1
            FFtbt=-9.0d0/9.8d1
         ELSE
            FFtb=0.0d0
            FFtbt=0.0d0
         ENDIF

         j=0
         DO j3=1,np3
            DO j2=1,np2
               DO j1=1,np1
                  j=j+1

                  xx(1,j)=REAL(xp(j)
     &                 +(EE*S2(1,j1,j2,j3)
     &                 +FFa*S3a(1,j1,j2,j3)
     &                 +FFb*S3b(1,j1,j2,j3)
     &                 +FFtb*S3T(1,j1,j2,j3))*dx)

                  IF (xx(1,j).lt.0.0) THEN
                     xx(1,j)=xx(1,j)+RBox
                  ELSE IF (xx(1,j).ge.RBox) THEN
                     xx(1,j)=xx(1,j)-RBox
                  ENDIF

                  xx(2,j)=REAL(yp(j)
     &                 +(EE*S2(2,j1,j2,j3)
     &                 +FFa*S3a(2,j1,j2,j3)
     &                 +FFb*S3b(2,j1,j2,j3)
     &                 +FFtb*S3T(2,j1,j2,j3))*dx)

                  IF (xx(2,j).lt.0.0) THEN
                     xx(2,j)=xx(2,j)+RBox
                  ELSE IF (xx(2,j).ge.RBox) THEN
                     xx(2,j)=xx(2,j)-RBox
                  ENDIF

                  xx(3,j)=REAL(zp(j)
     &                 +(EE*S2(3,j1,j2,j3)
     &                 +FFa*S3a(3,j1,j2,j3)
     &                 +FFb*S3b(3,j1,j2,j3)
     &                 +FFtb*S3T(3,j1,j2,j3))*dx)

                  IF (xx(3,j).lt.0.0) THEN
                     xx(3,j)=xx(3,j)+RBox
                  ELSE IF (xx(3,j).ge.RBox) THEN
                     xx(3,j)=xx(3,j)-RBox
                  ENDIF

                  vv(1,j)=REAL(vxp(j)
     &                 +CC*(EEt*S2(1,j1,j2,j3)
     &                 +FFat*S3a(1,j1,j2,j3)
     &                 +FFbt*S3b(1,j1,j2,j3)
     &                 +FFtbt*S3T(1,j1,j2,j3)))

                  vv(2,j)=REAL(vyp(j)
     &                 +CC*(EEt*S2(2,j1,j2,j3)
     &                 +FFat*S3a(2,j1,j2,j3)
     &                 +FFbt*S3b(2,j1,j2,j3)
     &                 +FFtbt*S3T(2,j1,j2,j3)))

                  vv(3,j)=REAL(vzp(j)
     &                 +CC*(EEt*S2(3,j1,j2,j3)
     &                 +FFat*S3a(3,j1,j2,j3)
     &                 +FFbt*S3b(3,j1,j2,j3)
     &                 +FFtbt*S3T(3,j1,j2,j3)))

               ENDDO
            ENDDO
         ENDDO

         DO k=1,Npart(2)
            xp(k)=xx(1,k)
            yp(k)=xx(2,k)
            zp(k)=xx(3,k)
            vxp(k)=vv(1,k)
            vyp(k)=vv(2,k)
            vzp(k)=vv(3,k)
         ENDDO

      open(unit=31,file=fileout(FL),status='unknown',
     &        form='unformatted')

      write(31)Npart,Massarr,afactor,redshift,FlagSfr,FlagFeedback,
     &        Nall,FlagCooling,NumFiles,dble(RBox),Omega0,OmegaL0,
     &        hlittle,FlagStellarAge,FlagMetals,NallHW,flagentrics,
     &        unused

      write(31)(xp(k),yp(k),zp(k),k=1,Npart(2))
      write(31)(vxp(k),vyp(k),vzp(k),k=1,Npart(2))
      write(31)(k,k=1, Npart(2))

      close(31)

      END

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL*8 data(*)

      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,
     &     ip3,k1,k2,n,nprev,nrem,ntot
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 
      ntot=1
      do idim=1,ndim
         ntot=ntot*nn(idim)
      enddo
      nprev=1
      do idim=1,ndim
         n=nn(idim)
         nrem=ntot/(n*nprev)
         ip1=2*nprev
         ip2=ip1*n
         ip3=ip2*nrem
         i2rev=1
         do i2=1,ip2,ip1
            if(i2.lt.i2rev)then
               do i1=i2,i2+ip1-2,2
                  do i3=i1,ip3,ip2
                     i3rev=i2rev+i3-i2
                     tempr=data(i3)
                     tempi=data(i3+1)
                     data(i3)=data(i3rev)
                     data(i3+1)=data(i3rev+1)
                     data(i3rev)=tempr
                     data(i3rev+1)=tempi
                  enddo
               enddo
            endif
            ibit=ip2/2
 1          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
               i2rev=i2rev-ibit
               ibit=ibit/2
               goto 1
            endif
            i2rev=i2rev+ibit
         enddo
         ifp1=ip1
 2       if(ifp1.lt.ip2)then
            ifp2=2*ifp1
            theta=isign*6.28318530717959d0/(ifp2/ip1)
            wpr=-2.d0*sin(0.5d0*theta)**2
            wpi=sin(theta)
            wr=1.d0
            wi=0.d0
            do i3=1,ifp1,ip1
               do i1=i3,i3+ip1-2,2
                  do i2=i1,ip3,ifp2
                     k1=i2 
                     k2=k1+ifp1
                     tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                     tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                     data(k2)=data(k1)-tempr
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr
                     data(k1+1)=data(k1+1)+tempi
                  enddo
               enddo
               wtemp=wr 
               wr=wr*wpr-wi*wpi+wr
               wi=wi*wpr+wtemp*wpi+wi
            enddo
            ifp1=ifp2
            goto 2
         endif
         nprev=n*nprev
      enddo
      return
      END
