module integrals
!----------
!   author:
!       Sequencer
!   description:
!       计算物理插值模块(HUST PHY 2013 第四次作业)
!       MIT协议
!   contains:
!       simpsonIntegral(xmin,xmax,N)
!       trapezoidIntegral(xmin,xmax,N)
!   exception:
!       None
!----------
contains
function func(x)
    implicit none
    real*8::x,func
    func = sin(x)
end function func
    
function simpsonIntegral(xmin,xmax,N)
    implicit none
    real*8::xmin,xmax,dx,simpsonIntegral
    integer::N,i
    real*8::X(N+1),Y(N+1)
    dx = (xmax-xmin)/N
    do i=1,N+1
        X(i) = xmin+dx*(i-1)
        Y(i) = func(X(i))
    enddo
    simpsonIntegral = Y(1)+Y(N+1)
    do i=2,N-1
        if (mod(i,2) .eq. 0) then
            simpsonIntegral = simpsonIntegral+4 *func(X(i))
        else
            simpsonIntegral = simpsonIntegral+2*func(X(i))
        endif
    enddo
    simpsonIntegral = simpsonIntegral*dx/3
end function simpsonIntegral

function trapezoidIntegral(xmin,xmax,N)
    implicit none
    real*8::xmin,xmax,dx,trapezoidIntegral
    integer::N,i
    real*8::X(N+1),Y(N+1)
    dx = (xmax-xmin)/N
    do i=1,N+1
        X(i) = xmin+dx*(i-1)
        Y(i) = func(X(i))
    enddo
    trapezoidIntegral = Y(1)+Y(N+1)
    do i=2,N-1
        trapezoidIntegral = trapezoidIntegral + 2*func(X(i))
    enddo
    trapezoidIntegral = trapezoidIntegral*dx/2

end function trapezoidIntegral

end module integrals
program main
    use integrals
    implicit none
    real*8::xmin,xmax,s
    integer::N
    xmin = 1
    xmax = 5
    N = 40
    s = simpsonIntegral(xmin,xmax,N)
    ! s = trapezoidIntegral(xmin,xmax,N)
    write(*,*) s
end program main