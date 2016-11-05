c usage:
c       grd2xyz f.grd -Z|z2grd3 -Rw/e/s/n -Idx/dy -Gf.grd3
c
c In some cases you have produced a .grd file using gmt commands
c and you wish to read it in a fortran or a C++ program


c Assume f.grd is created by GMT in your machine and f.grd3 is what you
c need.

c compile z2grd3.f in your machine:

c f77 z2grd3.f -o z2grd3
c Convert using both z2grd3 and grd2xyz. west/east/south/north and dx/dy
c must be exactly the same as those in f.grd when using "z2grd3" 

c grd2xyz f.grd -Z|z2grd3 -Rwest/east/south/north -Idx/dy -Gf.grd3

c Now you can use f.grd3 for further applications such as input
c to a fortran program. 
C
C To open an .grd3 file in fortran, please see the
c following code.
      implicit none
      integer i,j,dim,nx,ny,lent
      parameter(dim=5000)
      real*4 arr(dim,dim)
      real*8 xmin,xmax,ymin,ymax,dx,dy,rms
CCCCCCCCCCCCCCCCCCCCCCCCCCC
c variables for getting arguments
      integer*4 ii,nargs,iargc,i1,i2,i3,i4
      logical lg1, lg2, lg3
      character*80 cha, tbuf,ofile1
ccccccccccccccccccccccccccc

C get arguments *************************************************
      nargs=iargc()
      lg1=.true.
      lg2=.true.
      lg3=.true.
      do ii=1,nargs
      call getarg(ii,cha)
      if(cha(1:1).eq.'-') then

		if(cha(2:2).eq.'G' .or. cha(2:2).eq.'g') then
		ofile1=cha(3:)
                lg1=.false.

		elseif(cha(2:2).eq.'I' .or. cha(2:2).eq.'i') then
                i1=index(cha,'/')
                i2=index(cha,' ')
                tbuf=cha(3:i1-1)
 		read(tbuf,*)dx
                tbuf=cha(i1+1:i2-1)
		read(tbuf,*)dy
                lg2=.false.

 		elseif(cha(2:2).eq.'R' .or. cha(2:2).eq.'r') then
                i1=index(cha,'/')
		tbuf=cha(3:i1-1)
		read(tbuf,*)xmin
		i2=index(cha(i1+1:),'/')
                tbuf=cha(i1+1:i1+i2-1)
		read(tbuf,*)xmax
 		i3=index(cha(i1+i2+1:),'/')
		tbuf=cha(i1+i2+1:i1+i2+i3-1)
 		read(tbuf,*)ymin
		i4=index(cha(i1+i2+i3+1:),' ')
                tbuf=cha(i1+i2+i3+1:i1+i2+i3+i4-1)
 		read(tbuf,*)ymax
                lg3=.false.
		else
	        call write_error
                stop
		end if
      else
      call write_error
      stop
      end if
      end do

      if(lg1 .or. lg2 .or. lg3) then
      call write_error
      stop
      end if

      nx=(xmax-xmin)/dx+1.001
      ny=(ymax-ymin)/dy+1.001
      if(nx.gt.dim.or.ny.gt.dim)stop'increase dim'
c Read z table
      read(5,*)((arr(i,j),i=1,nx),j=ny,1,-1)
 
      open(60,file=ofile1,form='unformatted')
! grd1 format
!     open(60,file=ofile1,access='direct',form='binary',recl=64+4*nx*ny)
!         write(60,rec=1)nx,ny,xmin,xmax,ymin,ymax,dx,dy,
!    &((arr(i,j),j=1,ny),i=1,nx)
c write .grd3
      write(60)nx,ny,xmin,xmax,ymin,ymax,dx,dy
      write(60)((arr(i,j),j=1,ny),i=1,nx)
      end


      subroutine write_error
      write(0,*)' Usage:'
      write(0,*)''
      write(0,*)'z2grd3 -Rwest/east/south/north -Idx/dy<z_table(ascii)]'
      write(0,*)'-R: for west, east, south, north'
      write(0,*)'-I: grid sizes along x(longitude) and y(latitude)'
      return
      end
