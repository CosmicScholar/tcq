! Last update: Feb 27, 2003, by Cheinway Hwang, Cheng-Gi Wang, Yu-Shen Hsiao
!
! NAME
! tcq - FORTRAN90 program to compute terrain correction for gravity reduction using Gaussian quadrature
! 
! SYNOPSIS
! tcq -Cloc.xyh -IDTM1.grd3 -ODTM2.grd3 -Dinner_radius -Router_radius -Gout.txt [-Ttc_type]
!
! DESCRIPTION
! -C: file of longitude, latitude and height where terrain effect is computed
! -I: file of DTM for computing the inner zone effect 
! -O: file of DTM for computing the outer zone effect
! -D: radius (in km) of inner zone
! -R: radius (in km) of outer zone
! -G: output file of longitude, latitude and correction
!
! OPTIONS
! -T: 1 for terrain correction at terrain surface, 0 for terrain correction at sea level
!      [default: 1]
!
! EXAMPLES
!
! d:\researching\tc\tcq -Cinput_1019p.txt -Idtm3s.grd3 -Odtm30s.grd3 -D20 -R200 -Gtcq_2km_200km.txt -T1
! See the batch job exp.bat(total five cases) and test it at 
! http://space.cv.nctu.edu.tw/terrain.htm
! 
!
! NOTES
! 1. The input file format is longitude,latitude(in degree),orthometric height(in meter)
!          
! 2. The output file format is longitude,latitude(in degree),orthometric height(in meter) ,
!    innermost effect,inner zone terrain correction,outer zone terrain correction,
!    total terrain correction(in mgal)    
!
! 3. The DTM format is in .grd3 (binary)
!    record 1: nx,ny,west,east,south,north,dx,dy
!    record 2: ((grid(i,j),j=1,ny),i=1,nx)
!
!    where nx,ny = number of grids in longitude and latitude,
!          dx,dy = grid intervals (in degree)  in longitude and latitude
!    grid(i,j) = grid value at: long = (i-1)*dx, lat = (j-1)*dy
!
! REFERENCES
! 
! Authors: 
! Cheinway Hwang, Cheng-Gi Wang, Yu-Shen Hsiao
! Dept of Civil Engineering
! National Chiao Tung University
! 1001 Ta Hsueh Road, Hsinchu
! Taiwan
!
! All rights reserved
!
!Revised by zhikui guo for marine gravity survey, 2016/11/03
!================================进度条===============================================================
MODULE CMD_Progress
Implicit None
private
Logical , parameter , public :: CMD_PROGRESS_ABSOLUTE = .true.
Type , public :: CLS_CMD_Progress
Integer , private :: N , lens , i
Character :: M = "*" , O = "."
Character(len=64) :: Prefix
Contains
Procedure :: Set
Procedure :: Put
End Type CLS_CMD_Progress

contains

Subroutine Set( this , N , L )
Class( CLS_CMD_Progress ) :: this
Integer , Intent( IN ) :: N , L
this % N = N
this % lens = L
this % i = 0
this % Prefix = " Progress: " !//
End Subroutine Set

Subroutine Put( this , K , bAbsol )
Class( CLS_CMD_Progress ) :: this
Integer , Intent( IN ) :: K
Logical , optional :: bAbsol
Character(len=1) :: br
integer :: jm
this % i = this % i + K
if ( present( bAbsol ) ) then
if ( bAbsol ) this % i = K
end if
if ( this % i > this % n ) this % i = this % n
jm = Nint( real( this%i * this%lens ) / real( this%N ) )
if ( this%i < this%n ) then
br = char(13)
else
br = char(10)
end if
!write( * , '(5a,f6.2,2a)',advance="no") trim(this%Prefix) , '[' , &
write( * , '(5a,f6.2,2a\)') trim(this%Prefix) , '[' , & !// 如您的编译器不支持，请用上方语句代替
repeat(this%M , jm ) , repeat( this%O , this%lens-jm ) , '] ' , this%i*100.0/this%N , "%" , br
End Subroutine Put

END MODULE CMD_Progress
	
	include "link_fnl_shared.h"
	  USE NUMERICAL_LIBRARIES !invoke IMAL library
	  !=====使用进度条模块=====
		USE CMD_Progress
		
	  implicit none
	  integer*4 i,j,nx,ny,npt,npt_all,left,right,up,down,rx,ry,&
			   ii,jj,m,n,ldata,nxi,nyi,nxo,nyo,tc_type,&
				ileft,iright,jup,jdown,kx,ky,k	

 
	  real*8 lon,lat,ht,d2r,west,east,south,north,tc1,tc2,tc,fact,&
			 xi1,xi2,yi1,yi2,dxi,dyi,xo1,xo2,yo1,yo2,dxo,dyo,pi,time,&
			 inner_radius,outer_radius,radius1,radius2,lat_min,lon_min
	  real*8 rk,kc,hx,hy,aa,s0,ka,G,density,tc0,dtor,factor,tmp,cel
	  real*8 dx,dy,xmin,xmax,ymin,ymax,partA,partB,fun_cel
	   
	  real*8,allocatable,dimension(:,:)::data 
	  real*8,allocatable,dimension(:)::xa,ya
      !real*8,allocatable,dimension(:)::xa_debug,ya_debug !using to debug output
 
	  real*4,allocatable,dimension(:,:)::ain,ao,a
	  logical check     
      integer,external :: GetLineNum
! Variables for input arguments. These variables should not be used in
! the main program. 
	  integer  par,nargs,iargc 
	  character*80 ifile1,ifile2,ifile3,ofile1,cha,tbuf
	  logical io(10),ltc
	  data par/6/
	type( CLS_CMD_Progress ) ::Progress !DYI progressbar
	  
! get command-line arguments 
	  nargs=iargc()
	  if(nargs.eq.0) then
	  call write_error
	  stop
	  end if
	  tc_type=1
	  ltc=.true.
	  do i=1,par
	  io(i)=.true.
	  end do

	   i=0
	  do ii=1,nargs
	  call getarg(ii,cha)
	  if(cha(1:1).eq.'-') then
		if(cha(2:2).eq.'C' .or. cha(2:2).eq.'c') then
		ifile1=cha(3:)
				i=i+1
		io(i)=.false.

		elseif(cha(2:2).eq.'I' .or. cha(2:2).eq.'i') then
		ifile2=cha(3:)
				i=i+1
		io(i)=.false.

		elseif(cha(2:2).eq.'O' .or. cha(2:2).eq.'o') then
		ifile3=cha(3:)
				i=i+1
		io(i)=.false.

		elseif(cha(2:2).eq.'D' .or. cha(2:2).eq.'d') then
		tbuf=cha(3:)
		read(tbuf,*)inner_radius
				i=i+1
		io(i)=.false.

		elseif(cha(2:2).eq.'R' .or. cha(2:2).eq.'r') then
		tbuf=cha(3:)
		read(tbuf,*)outer_radius
				i=i+1
		io(i)=.false.

		elseif(cha(2:2).eq.'G' .or. cha(2:2).eq.'g') then
		ofile1=cha(3:)
				i=i+1
		io(i)=.false.
		
		!elseif(cha(2:2).eq.'N' .or. cha(2:2).eq.'n') then
		!tbuf=cha(3:)
		!read(tbuf,*)npt_all
		!		i=i+1
		!io(i)=.false.
! Optional parameters
		elseif(cha(2:2).eq.'T' .or. cha(2:2).eq.'t') then
		tbuf=cha(3:)
		read(tbuf,*)tc_type
		if(tc_type.eq.0) ltc=.false.


		else
			call write_error
		end if
	  else
	  call write_error
	  end if
	  end do

	  do i=1,par
	  if(io(i)) then
	  call write_error
	  end if
	end do
	  
! open computed point file
	open(10,file=ifile1)
    !get numbers of data
    !write(*,*)'总数据个数:', GetLineNum(ifile1)
    npt_all=GetLineNum(ifile1)
! open inner zone DTM
	open(11,file=ifile2,form='unformatted')
! open outer zoner DTM
	open(12,file=ifile3,form='unformatted')
! open output file
	open(60,file=ofile1)
! read DTM for inner zone       
	read(11)nxi,nyi,xi1,xi2,yi1,yi2,dxi,dyi
	allocate(ain(nxi,nyi))
	read(11)((ain(i,j),j=1,nyi),i=1,nxi)
	

! read DTM for outer zone
	read(12)nxo,nyo,xo1,xo2,yo1,yo2,dxo,dyo
	allocate(ao(nxo,nyo))
	read(12)((ao(i,j),j=1,nyo),i=1,nxo)

	pi=4.d0*datan(1.d0)
	dtor=datan(1.d0)/45.d0 ! radian per degree
	density=2.67d3 ! Crust density: kg/m**3
	G=6.67d-11 !Gravitational constant: m**3/(kg*sec**2)
	

	radius1=inner_radius/(2*pi*6371/360.d0) ! radius1 in degree
	radius2=outer_radius/(2*pi*6371/360.d0) ! radius2 in degree
 ! read computed points file    
!====================调用进度条========================================
call Progress % Set( npt_all ,25 )!// 1700次，显示长度25
Progress % Prefix = "Calculating: " !// 前方提示文字，不是必须
Progress % M = "#" !// 已完成部分的字符，不是必须
Progress % O = "." !// 未完成部分的字符，不是必须

	write(*,*)'Number of calculate points: ', npt_all
	npt=0
!/////////////////////////////////////////////////////////////////////////////////////
1   read(10,*,end=2)lon,lat,ht
	  if(npt.ge.npt_all)goto 2
	  !=========下面这条放在循环里------
	call Progress % Put(npt , CMD_PROGRESS_ABSOLUTE ) !// 调用进度条绝对方式
 ! Calculate the innermost zone effect
	   
	i=(lon-xi1)/dxi+1.01
	j=(lat-yi1)/dyi+1.01
	ileft=i-10
	  if(ileft.lt.1)ileft=1
	  iright=i+10
	  if(iright.gt.nxi)iright=nxi
	  jdown=j-10
	  if(jdown.lt.1)jdown=1
	  jup=j+10
	  if(jup.gt.nyi)jup=nyi
	  kx= iright-ileft+1
	  ky=jup-jdown+1
	  allocate(xa(kx),ya(ky),data(kx,ky))
      !allocate(xa_debug(kx),ya_debug(ky))
	  k=0
	  factor=dcos(lat*dtor)*6371000.d0
	  do ii=ileft, iright
	  k=k+1
	  tmp=xi1+(ii-1)*dxi
	  xa(k)=(tmp-lon)*dtor*factor
      !xa_debug(k)=(tmp)
	  end do
	  k=0
	  do jj=jdown, jup
	  k=k+1
	  tmp=yi1+(jj-1)*dyi
	  ya(k)=(tmp-lat)*dtor*6371000.d0
      !ya_debug(k)=(tmp)
	  end do
!open(110,file='debuginfo_innermost.txt')
	  do ii=ileft, iright
	  do jj=jdown, jup
	  data(ii-ileft+1, jj-jdown+1)=ain(ii,jj)
	  !write(110,'(2f12.6,5f10.3)')xa_debug(ii-ileft+1),ya_debug(jj-jdown+1),ain(ii,jj)
      !write(*,*) ii,jj
	  end do
    end do
    
	  check=.true.
	  ldata=max(kx,ky) 
!C Calculate gradient along y and x using quadratic polynomials.
!c See IMSL manual, p.505
	  hx= dqd2dr(1, 0, 0.d0, 0.d0, kx, xa, ky,ya, data, ldata, check)
	  hy= dqd2dr(0, 1, 0.d0, 0.d0, kx, xa, ky,ya, data, ldata, check)
	  s0=dsqrt(dxi*dtor*factor*dyi*dtor*6371000.d0/pi)
	  partA=2*pi*density*G*s0*1.d5
	  aa= dsqrt(hx**2+hy**2)
	  rk=aa/dsqrt(1.d0+aa*aa)
	  kc=dsqrt(1.d0-rk*rk)
	
	  fun_cel=cel(kc,1.d0,1.d0,1.d0)
	  partB=fun_cel*4.d0*density*G*s0*1.d5/dsqrt(1+aa*aa)
	
	  tc0=partA-partB   !innermost effect
	  deallocate(xa,ya,data)
! Calculate the contribution from the inner zone       
	  rx=radius1/(dxi*dcos(lat*d2r))
	  ry=radius1/dyi 
	  left=max0(1,i-rx)
	right=min0(i+rx,nxi)
	down=max0(1,j-ry)
	up=min0(j+ry,nyi)
	nx=right-left+1
	ny=up-down+1

	allocate(a(nx,ny))
	lon_min=xi1+dfloat(left-1)*dxi
	lat_min=yi1+dfloat(down-1)*dyi
!open(111,file='debuginfo_inner.txt')
	do m=down,up
	do n=left,right
	a(n-left+1,m-down+1)=ain(n,m)
    !write(111,'(2f12.6,5f10.3)')lon_min+(n-left)*dxi,lat_min+(m-down)*dyi,ain(n,m)
	  end do
	end do
	  west=lon_min
	east=xi1+dfloat(right-1)*dxi
	south=lat_min
	north=yi1+dfloat(up-1)*dyi
	  
	call tcin(a,nx,ny,lon_min,lat_min,dxi,dyi,ltc,lon,lat,ht,tc1)

! Calculate the contribution from the outer zone  
	i=(lon-xo1)/dxo+1.01
	j=(lat-yo1)/dyo+1.01

	rx=radius2/(dxo*dcos(lat*d2r))
	  ry=radius2/dyo 
	  left=max0(1,i-rx)
	right=min0(i+rx,nxo)
	down=max0(1,j-ry)
	up=min0(j+ry,nyo)
	deallocate(a)
	nx=right-left+1
	ny=up-down+1
	allocate(a(nx,ny))
 
	lon_min=xo1+(left-1.d0)*dxo
	lat_min=yo1+(down-1.d0)*dyo
!open(112,file='debuginfo_outer.txt')
	do m=down,up
	  do n=left,right
		a(n-left+1,m-down+1)=ao(n,m)
        !write(112,'(2f12.6,5f10.3)')lon_min+(n-left)*dxo,lat_min+(m-down)*dyo,ao(n,m)
	  end do
	end do
	call tcout(a,nx,ny,lon_min,lat_min,dxo,dyo,ltc,&
			   west,east,south,north,lon,lat,ht,tc2)
	deallocate(a)
	tc=tc1+tc2+tc0     !Total terrain correction
	if(ht/=9999.0) then
!  output format
	write(60,'(2f12.6,5f10.3)')lon,lat,ht,tc0,tc1,tc2,tc
	end if
	npt=npt+1
	go to 1
2     write(0,*)'Number of points:', npt
	deallocate(ain,ao)
!c     output the time of the program  executeing    
	call second(time)
	write(*,'(a6,f10.3,a4)')'Time:',time,'sec'
	!write(60,'(a6,f10.3,a4)')'Time:',time,'sec'
	stop
    end

    !get linenumbers of file
    FUNCTION GetLineNum(filename)
    implicit none
    integer GetLineNum
    character*80 filename
    integer::val=0
    real line
    GetLineNum=0
    OPEN(UNIT=77,FILE=filename,STATUS='OLD')
    do while(val==0)
      read(77,*,iostat=val) line
      GetLineNum=GetLineNum+1
    end do
    end
    
	!calculate K(), which is the complete elliptic integral of the first kind defined in Eq.(5)
	FUNCTION CEL(QQC,PP,AA,BB)
	implicit real*8(a-h,o-z)
	  PARAMETER (CA=.0003d0, PIO2=1.5707963268d0)
	  IF(QQC.EQ.0.)PAUSE 'failure in CEL'
	  QC=DABS(QQC)
	  A=AA
	  B=BB
	  P=PP
	  E=QC
	  EM=1.
	  IF(P.GT.0.)THEN
		P=DSQRT(P)
		B=B/P
	  ELSE
		F=QC*QC
		Q=1.-F
		G=1.-P
		F=F-P
		Q=Q*(B-A*P)
		P=DSQRT(F/G)
		A=(A-B)/G
		B=-Q/(G*G*P)+A*P
	  ENDIF
1     F=A
	  A=A+B/P
	  G=E/P
	  B=B+F*G
	  B=B+B
	  P=G+P
	  G=EM
	  EM=QC+EM
	  IF(DABS(G-QC).GT.G*CA)THEN
		QC=DSQRT(E)
		QC=QC+QC
		E=QC*EM
		GO TO 1
	  ENDIF
	  CEL=PIO2*(B+A*EM)/(EM*(EM+P))
	  RETURN
	  END
 


 
	  subroutine tcin(arr,nx,ny,xmin,ymin,dx,dy,ltc,lon,lat,hp,gra)
!c Last update: March 17, 2001
!c 
!c Compute terrain correction due to inner zone
!c 
!c Input: 
!c arr: elevation, in meter. SINT
!c nx,ny: grid sizes of arr
!c xmin,ymin: west, south of borders
!c dx,dy: grid intervals, in degree
! tc: true for terrain correction at surface, false for terrain correction at sea level
!c 
! lon,lat, hp: longitude, latitude, height of computed point, in 
!             degree,degree,meter
!c
!c Output:
!c gra: terrain correction, in mgal
!c
!c
!c Author: Cheinway Hwang, NCTU, Hsinchu, Taiwan
 
	  implicit none
	  integer i,j,nx,ny,degree,ix,iy
 
	  parameter(degree=6)
!c arr contains height
	  real*4 arr(nx,ny)
	  real*8,allocatable,dimension(:)::wxp,wyp,xvec,yvec, &
		workx,worky,xnode,ynode
 
	  real*8 sum,pi,d2r,r,x0,y0,x1,x2,y1,y2,gra,val,kern,umin,&
		  xmin,xmax,ymin,ymax,dx,dy,gconst,density,e2,psi,fact,&
		   latp,lonp,hp,&
		  u,tmp1,tmp2,d2,d,sp,cp,dd,lat,lon,dint,phip,lamp,gx,gy
	  logical ltc
	  data gconst/6.6732d-11/ ! m**3/kg/s**2
!c Averaged density of rock
	  data density/2.67d3/ ! kg/m**3 
	  allocate(wxp(nx),wyp(ny),xvec(nx),yvec(ny),workx(nx),&
			   worky(ny),xnode(nx),ynode(ny)) 
	  pi=4.d0*datan(1.d0)
	  d2r=pi/180.d0
	  e2=0.00669438002290d0 !e2=2*f-f**2
	  fact=dcos(lat*d2r)          
!c Use T/P ellipsoid to compute a mean local radius of the earth
	  r=6378136.3d0*dsqrt(1.d0-e2)/(1.d0-e2*(dsin(lat*d2r))**2)

	  x1=(xmin-lon)*d2r*fact*r
	  x2=(xmin+dfloat(nx-1)*dx-lon)*d2r*fact*r
	  y1=(ymin-lat)*d2r*r
	  y2=(ymin+dfloat(ny-1)*dy-lat)*d2r*r
  
	  gx=r*dx*d2r*fact !meter
	  gy=r*dy*d2r !meter
 
   
!c Compute abscissas and weights in longitude and latitude

	  call GAULEG(x1,x2,xnode,wxp,nx)
 
	  call GAULEG(y1,y2,ynode,wyp,ny)
	 
 ! Compute planar coordinates. the origin is at lon, lat
	  do i=1,nx
	  xvec(i)=(xmin+dfloat(i-1)*dx-lon)*d2r*r*fact
	  end do
	  do j=1,ny
	  yvec(j)=(ymin+dfloat(j-1)*dy-lat)*d2r*r
	  end do
!c Sum over x      
	  do j=1,ny
	  do i=1,nx
	d=xvec(i)**2+yvec(j)**2
	if(d.le.gy) then
	workx(i)=0.d0
	else
! terrain correction at surface
	  if(ltc) then
	  workx(i)=1.d0/dsqrt(d)-1.d0/dsqrt(d+(arr(i,j)-hp)**2)
! terrain correction at sea level
	  else
	  workx(i)=1.d0/dsqrt(d+hp**2)-1.d0/dsqrt(d+arr(i,j)**2)
	end if
	  end if
	  end do

	  sum=0.d0

	  do i=1,nx
	  call interp1d(xvec,workx,gx,nx,degree,xnode(i),val)
	  sum=sum+wxp(i)*val
	  end do

	  worky(j)=sum
	  end do

!c Sum over y
	  sum=0
	  do j=1,ny
	  call interp1d(yvec,worky,gy,ny,degree,ynode(j),val)
	  sum=sum+wyp(j)*val  
	  end do
!c Multiply the sum by G*rho*1.d5 gra is in mgal
	  gra=gconst*density*sum*1.d5
	  return
	  deallocate(wxp,wyp,xvec,yvec,workx,worky,xnode,ynode)
	  end subroutine 


	   subroutine tcout(arr,nx,ny,xmin,ymin,dx,dy,ltc,&
						west,east,south,north,lon,lat,hp,gra)
!c Last update: March 17, 2001
!c 
!c Compute terrain correction due to outer zone
!c 
!c Input: 
!c arr: elevation, in meter. SINT
!c nx,ny: grid sizes of arr
!c xmin,ymin: west, south of borders
!c dx,dy: grid intervals, in degree
! tc: true for terrain correction at surface, false for terrain correction at sea level
! w,e,s,n: borders of the inner zone
!c 
!c lon,lat, hp: longitude, latitude, height of computed point, in 
!             degree,degree,meter
!c
!c Output:
!c gra: terrain correction, in mgal
!c
!c
!c Author: Cheinway Hwang, NCTU, Hsinchu, Taiwan
 
	  implicit none
	  integer i,j,nx,ny,degree,ix,iy
 
	  parameter(degree=6)
!c arr contains height
	  real*4 arr(nx,ny)
	  real*8,allocatable,dimension(:)::wxp,wyp,xvec,yvec, &
			 workx,worky,xnode,ynode
 
	  real*8 sum,pi,d2r,r,x0,y0,x1,x2,y1,y2,gra,val,kern,umin,&
		xmin,xmax,ymin,ymax,dx,dy,gconst,density,e2,psi,fact,&
		u,tmp1,tmp2,d2,d,sp,cp,dd,lat,lon,dint,phip,lamp,gx,gy,&
		west,east,south,north,hp
	  logical ltc
	  data gconst/6.6732d-11/ ! m**3/kg/s**2
!c Averaged density of rock
	  data density/2.67d3/ ! kg/m**3 
	  allocate(wxp(nx),wyp(ny),xvec(nx),yvec(ny),workx(nx),&
			   worky(ny),xnode(nx),ynode(ny))
 
	  pi=4.d0*datan(1.d0)
	  d2r=pi/180.d0
	  e2=0.00669438002290d0 !e2=2*f-f**2
	  fact=dcos(lat*d2r)          
!c Use T/P ellipsoid to compute a mean local radius of the earth
	  r=6378136.3d0*dsqrt(1.d0-e2)/(1.d0-e2*(dsin(lat*d2r))**2)

	  x1=(xmin-lon)*d2r*fact*r
	  x2=(xmin+dfloat(nx-1)*dx-lon)*d2r*fact*r
	  y1=(ymin-lat)*d2r*r
	  y2=(ymin+dfloat(ny-1)*dy-lat)*d2r*r

	west=(west-lon)*d2r*fact*r
	east=(east-lon)*d2r*fact*r
	south=(south-lat)*d2r*r
	north=(north-lat)*d2r*r
  
	  gx=r*dx*d2r*fact !meter
	  gy=r*dy*d2r !meter
 
   
!c Compute abscissas and weights in longitude and latitude

	  call GAULEG(x1,x2,xnode,wxp,nx)
 
	  call GAULEG(y1,y2,ynode,wyp,ny)
	 
 ! Compute planar coordinates. the origin is at lon, lat
	  do i=1,nx
	  xvec(i)=(xmin+dfloat(i-1)*dx-lon)*d2r*r*fact
	  end do
	  do j=1,ny
	  yvec(j)=(ymin+dfloat(j-1)*dy-lat)*d2r*r
	  end do
!c Sum over x      
	  do j=1,ny
	  do   i=1,nx
  
	  if(xvec(i).ge.west .and. xvec(i).le.east .and.& 
		yvec(j).ge.south .and. yvec(j).le.north) then
		workx(i)=0.d0
	  else

		d=xvec(i)**2+yvec(j)**2
		if(d.le.gy) then
		workx(i)=0.d0
		else
! terrain correction at surface
			if(ltc) then
			workx(i)=1.d0/dsqrt(d)-1.d0/dsqrt(d+(arr(i,j)-hp)**2)
! terrain correction at sea level
			 else
			 workx(i)=1.d0/dsqrt(d+hp**2)-1.d0/dsqrt(d+arr(i,j)**2)
			 end if
		end if

	  end if
	   end do

	  sum=0.d0

	  do i=1,nx
	  call interp1d(xvec,workx,gx,nx,degree,xnode(i),val)
	  sum=sum+wxp(i)*val
	  end do

	  worky(j)=sum
	  end do

!c Sum over y
	  sum=0
	  do j=1,ny
	  call interp1d(yvec,worky,gy,ny,degree,ynode(j),val)
	  sum=sum+wyp(j)*val  
	  end do
!c Multiply the sum by G*rho*1.d5 gra is in mgal
	  gra=gconst*density*sum*1.d5
	  return
	  deallocate(wxp,wyp,xvec,yvec,workx,worky,xnode,ynode)
	  end subroutine

	  SUBROUTINE GAULEG(X1,X2,X,W,N)
!C Program from Numerical Recipe
	  IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*8 X1,X2,X(N),W(N)
	  PARAMETER (EPS=3.D-14)
	  M=(N+1)/2
	  XM=0.5D0*(X2+X1)
	  XL=0.5D0*(X2-X1)
	  DO 12 I=1,M
		Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
		  P1=1.D0
		  P2=0.D0
		  DO 11 J=1,N
			P3=P2
			P2=P1
			P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
		  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
		  Z1=Z
		  Z=Z1-P1/PP
		IF(ABS(Z-Z1).GT.EPS)GO TO 1
		X(I)=XM-XL*Z
		X(N+1-I)=XM+XL*Z
		W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
		W(N+1-I)=W(I)
12    CONTINUE
	  RETURN
	  END
			
	  subroutine interp1d(time,data,h,n,deg,x,out)
!c Program to do 1d polynomial interpolation
!c deg: number of data points used in interpolaiton  
	  implicit none
	  integer n,degree,deg,k,k1,k2,np,i
	  real*8 data(*),out,time(*),ya(50),x,x0,h,tmp
	  if(degree.gt.50)stop'increase dim of ya'
!c Use even number of points
	  degree=(deg/2)*2
	  k=(x-time(1))/h+1.01
	  k1=k-degree/2+1
	  if(k1.lt.1)k1=1
	  k2=k+degree/2
	  if(k2.gt.n)k2=n
	  np=k2-k1+1
	  x0=time(k1)
	  k1=k1-1
	  do i=1,np
	  ya(i)=data(k1+i)
	  end do
	  call divint(x0,ya,np,h,x,tmp)
	  out=tmp
	  return
	  end
	  SUBROUTINE DIVINT(X0,YA,N,H,X,Y)
!C==============================================================
!C   Polynomial interpolation using divided Difference
!C==============================================================
!C   Veriables :
!C       X0    ==> coordinate of the first point ( Input )
!c                 ie, its value is ya(1)
!C       YA(N) ==> Data             ( Input )
!C       N     ==> Number of points ( Input )
!C       H     ==> Stepsize         ( Input )
!C       X     ==> coordinate where int is needed ( Input )
!C       Y     ==> Result           ( Output )
!C       DF    ==> Work array
	  IMPLICIT REAL*8(A-H,O-Z)
	  DIMENSION YA(N),DF(50,50)
	  if(n.gt.50) stop'increase dim of DF'
	  S=(X-X0)/H
!c If x happens to be on the knot, then s is an integer
!c
	  if( abs(int(s)-s).lt.1.0d-7) then
	  y=ya(int(s)+1)
	  return
	  end if

	  S1=1.D0
	  S2=1.D0
	  DO J=1,N
		 DF(1,J)=YA(J)
	  ENDDO
	  Y=YA(1)
	  DO  I=2,N
	  DO  J=1,N-I+1
		 DF(I,J)=DF(I-1,J+1)-DF(I-1,J)
	  ENDDO
	  di=dfloat(i)
		 S1=S1*(S-di+2.d0)
		 S2=S2*(di-1.d0)
		 DY=S1/S2*DF(I,1)
		 Y=Y+DY
	  ENDDO
	  RETURN
	  END


	  SUBROUTINE second(TIME)
!c Program to sum the system and user times 
	  REAL*8 TIME
	  REAL*4 T(2)
	  DATA TOT/0.D0/
	  TIME=DTIME(T)+TOT
	  TOT=TIME
	  RETURN
	  END

	subroutine write_error
	write(0,*)'Usage: tcq -Cloc.xyh -Idem1.grd3 -Odem2.grd3& 
	  -Dinner_radius -Router_radius -Gout.txt -Ttc_type'
	write(0,*)'-C: file of longitude, latitude and height where terrain effect is computed'
	write(0,*)'-I: file of DTM for computing the inner zone effect'
	write(0,*)'-O: file of DTM for computing the outer zone effect'
	write(0,*)'-D: radius (in km) of inner zone'
	write(0,*)'-R: radius (in km) of outer zone'
	write(0,*)'-G: output file of longitude, latitude and correction'
	write(0,*)'-T: 1 for terrain correction at terrain surface, 0 for terrain correction at sea level'
	stop
	return
	end
 
