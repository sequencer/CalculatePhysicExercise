module linalg
!----------
!   author:
!       Sequencer
!   description:
!       计算物理线性代数辅助模块
!       MIT协议
!   contains:
!       norm(vector)
!       zeros(n,m)
!       identity(n)
!       diag(matrix)
!       column_stack(A,B)
!       row_stack(A,B)
!       switchmatrix(matrix,m,n)
!   exception:
!       wrongShape
!----------
contains

function norm(vector)
    !返回一个向量的长度
    implicit none
    real*8::vector(:,:)
    real*8::norm
    integer::i
    do while (i .le. size(vector))
        norm = norm + vector(i,1)**2
        i = i+1
    enddo
end function norm

function zeros(n,m)
    !返回一个n,m的零矩阵
    implicit none 
    integer n,m
    real*8 zeros(n,m)
    zeros = 0d0
end function zeros

function identity(n)
    !返回一个n阶单位矩阵
    implicit none 
    integer n,m,i
    real*8 identity(n,n)
    identity = 0d0
    forall(i = 1:n) identity(i,i) = 1d0
end function identity

function diag(matrix)
    implicit none
    real*8 :: matrix(:,:)
    real*8 :: diag(size(matrix,dim = 1),size(matrix,dim = 1))
    integer :: shape,i
    shape = size(matrix,dim = 1)
    !形状判断
    if (size(matrix,dim = 1)-size(matrix,dim = 2) .ne. 0) then
        write(*,*) 'wrongShape'
    endif
    diag = matrix
    forall(i = 1:shape) diag(i,i) = 0d0
    diag = matrix - diag
end function diag

function column_stack(A,B)
    implicit none
    real*8,intent(in):: A(:,:),B(:,:)
    real*8 :: column_stack(size(A,dim =1),size(A,dim=2)+size(B,dim =2))
    column_stack = 0
    column_stack(:,1:size(A,dim =2)) = A
    column_stack(:,size(A,dim =2)+1:) = B
end function column_stack

function row_stack(A,B)
    implicit none
    real*8,intent(in):: A(:,:),B(:,:)
    real*8 :: row_stack(size(A,dim=1)+size(B,dim =1),size(A,dim=2))
    row_stack = 0
    row_stack(1:size(A,dim = 1),:) = A
    row_stack(size(A,dim = 1)+1:,:) = B
end function row_stack


subroutine switchmatrix(matrix,m,n)
    !交换矩阵中的某两行
    real*8,intent(inout):: matrix(:,:)
    real*8::temp_mat(1,size(matrix,dim = 2))
    integer shape
    shape = size(matrix,dim = 1)
    temp_mat(1,:) = matrix(m,:)
    matrix(m,:) = matrix(n,:)
    matrix(n,:) = temp_mat(1,:)
end subroutine switchmatrix
end module linalg

module linear_soving
!----------
!   author:
!       Sequencer
!   description:
!       计算物理非线型方程寻根模块(HUST PHY 2013 第二次作业)
!       MIT协议
!   contains:
!       downtri(A,B,X)
!       uptri(A,B,X)
!       gausselimination(A,B,X)
!       lu(matrix,L,U)
!       doolittle(A,B,X)
!       jacobi(A,B,X,f_error,k_max,x_error)
!       successive_relaxation(A,B,X,omega,f_error,k_max,x_error)
!   exception:
!       kOverflow
!       fDivergence
!       wrongShape
!----------
use linalg

contains

subroutine downtri(A,B,X)
    implicit none
    real*8,intent(in) :: A(:,:),B(:,:)
    real*8,intent(out) :: X(size(A,dim = 1),1)
    integer :: shape
    integer :: i
    shape = size(A,dim = 1)
    X = zeros(shape,1)
    i = 1
    do while (i.le.shape)
        X(i,1) = (B(i,1) - dot_product(A(i,:),X(:,1)))/A(i,i)
        i = i+1
    enddo
end subroutine downtri

subroutine uptri(A,B,X)
    implicit none
    real*8,intent(in) :: A(:,:),B(:,:)
    real*8,intent(out) :: X(size(A,dim = 1),1)
    integer :: shape
    integer :: i
    shape = size(A,dim = 1)
    X = zeros(shape,1)
    i = shape
    do while (i.ge.1)
        X(i,1) = (B(i,1) - dot_product(A(i,:),X(:,1)))/A(i,i)
        i = i-1
    enddo
end subroutine uptri

subroutine gausselimination(A,B,X)
    implicit none
    real*8 :: A(:,:),B(:,:)
    real*8,intent(out) :: X(size(A,dim = 1),1)
    real*8 :: matrix(size(A,dim = 1),size(A,dim = 1)+1)
    integer :: i,j,shape,jmaxindex
    if (size(A,dim =1).ne.size(A,dim = 2)) then
        write(*,*) "wrongShape"
    endif
    matrix = column_stack(A,B)
    shape = size(A,dim = 1)
    X = 0
    i = 1
    jmaxindex = 0
    do while (i.le.shape)
        jmaxindex = i-1+ sum(maxloc(matrix(i:shape,i))) ! 用sum提取array中的唯一一个元素
        call switchmatrix(matrix,i,jmaxindex)
        matrix(i,:) = matrix(i,:)/matrix(i,i)
        j = i+1
        do while (j.le.shape)
            matrix(j,:) = matrix(j,:)-matrix(j,i)*matrix(i,:)
            j = j+1
        enddo
        i = i + 1
    enddo
    A(1:shape,1:shape)= matrix(:,1:shape)
    B(1:shape,1)= matrix(:, shape+1)
    call uptri(A,B,X)
end subroutine gausselimination

subroutine lu(matrix,L,U)
    implicit none
    real*8,intent(in) :: matrix(:,:)
    real*8,intent(inout) :: L(size(matrix,dim = 1),size(matrix,dim = 2)),U(size(matrix,dim = 1),size(matrix,dim = 2))
    integer :: shape,i,j
    shape = size(matrix,dim = 1)
    if (size(matrix,dim =1).ne.size(matrix,dim = 2)) then
        write(*,*) "wrongShape"
    endif
    L = identity(shape)
    U = zeros(shape,shape)
    U(1,:) = matrix(1,:)
    L(:,1) =matrix(:,1) / U(1,1)
    i = 2
    do while (i .le. shape)
        j = i
        do while (j.le.shape)
            U(i,j) = matrix(i,i) - dot_product(L(i,:),U(:,j))
            j = j+1
        enddo
        j = i+1
        do while (j.le.shape)
            L(j,i) = (matrix(j,i)-dot_product(L(j,:),U(:,i)))/U(i,i)
            j =j+1
        enddo
        i = i + 1
    enddo
end subroutine lu

subroutine doolittle(A,B,X)
    implicit none
    real*8, intent(in) :: A(:,:),B(:,:)
    real*8, intent(out) :: X(size(B,dim = 1),1)
    real*8 ::L(size(A,dim = 1),size(A,dim = 2)),U(size(A,dim = 1),size(A,dim = 2)),Y(size(B,dim = 1),1)
    call lu(A,L,U)
    call downtri(L,B,Y)
    call uptri(U,Y,X)
end subroutine doolittle

subroutine jacobi(A,B,X,f_error,k_max,x_error)
    implicit none
    real*8, intent(in) :: A(:,:),B(:,:)
    real*8, intent(inout) :: X(:,:),x_error
    real*8 :: Di(size(A,dim = 1),size(A,dim = 2)),G(size(A,dim = 1),size(A,dim = 2)),g0(size(A,dim = 1),1)
    real*8 :: X0(size(X,dim =1),size(X,dim =2))
    real*8 :: f_error
    integer :: k_max,k,shape
    shape = size(A,dim = 1)
    if (size(A,dim =1).ne.size(A,dim = 2)) then
        write(*,*) "wrongShape"
    endif
    G = 1
    Di = diag(G/A)
    G = identity(shape) - matmul(Di,A)
    g0 = matmul(Di,B)
    X = matmul(G,X) + g0
    k = 1
    do while (norm(matmul(A,X)-B)>f_error)
        X0 = X
        X = matmul(G,X)+g0
        k = k+1
        if (k>k_max) then
            write(*,*) "kOverflow"
            call abort()
        endif
    enddo
end subroutine jacobi

subroutine successive_relaxation(A,B,X,omega,f_error,k_max,x_error)
    implicit none
    real*8, intent(in) :: A(:,:),B(:,:),omega,f_error,x_error
    real*8, intent(inout) :: X(:,:)
    real*8 :: X0(size(X,dim =1),size(X,dim =2)),X_copy(size(X,dim =1),size(X,dim =2)),X_i
    integer k_max,shape,k,i

    shape = size(A,dim = 1)
    if (size(A,dim =1).ne.size(A,dim = 2)) then
        write(*,*) "wrongShape"
    endif
    k = 0
    do while (norm(matmul(A,X)-B).gt.f_error)
        i = 1
        do while (i.le.shape)
            X0 = X
            X_copy = X
            x_i = (1-omega)*X_copy(i,1)
            X_copy(i,1) = 0
            X(i,1) = x_i - omega*(sum(matmul(A(i,:),X_copy)-B(i,1))) / A(i,i)
            i = i+1
        enddo
        k = k+1
        if (k>k_max) then
            write(*,*) "kOverflow"
            call abort()
        endif
    enddo
end subroutine

end module linear_soving

module interpolation
!----------
!   author:
!       Sequencer
!   description:
!       计算物理插值模块(HUST PHY 2013 第三次作业)
!       MIT协议
!   contains:
!       larangeInterpolation(X,Y,X0)
!       newtonInterpolation(X,Y,X0)
!   exception:
!       wrongShape
!----------
use linear_soving
implicit none

contains


function sampleFunction(x)
implicit none
real*8 :: x,sampleFunction
sampleFunction = (1 + x ** 2) ** (-1)
end function sampleFunction

function linspace(xmin,xmax,N)
    implicit none
    real*8 :: xmin,xmax,linspace(N)
    integer :: N,i
    linspace = (/((xmin+(xmax-xmin)/(N-1)*i),i = 0,N )/)
end function linspace


subroutine samplingFromFunction(X,Y)
    implicit none
    real*8::Y(:),X(:)
    integer :: i,N
    N = size(X)
    do i=1,N
        Y(i) = sampleFunction(X(i))
    enddo
end subroutine samplingFromFunction

subroutine larangeInterpolation(X,Y,X0,Y0)
    implicit none
    ! X,Y 为采样点的座标集合, X, Y 为拟合后的座标集合
    real*8 :: X(:),Y(:),X0(:)
    real*8 :: Y0(size(X0))
    real*8 :: li
    ! shape为采样点座标的点数, shape0 为拟合后的座标的点数
    integer :: shape,shape0,k,i,j
    ! 判断X,Y的点数是否一致
    if (size(X) .ne. size(Y)) then
        write(*,*) "wrongShape"
        call abort()
    endif
    Y0 = 0
    shape0 = size(X0)
    shape = size(X)
    ! 遍历X0
    do k=1,shape0
        ! 遍历X
        do i=1,shape
            li = 1
            ! 连乘求得li
            do j=1,shape
                if (j == i) then
                    cycle
                endif
                li = li *(X0(k)-X(j))/(X(i)-X(j))
            enddo
            Y0(k) = Y0(k)+Y(i)*li
        enddo
    enddo
end subroutine larangeInterpolation

subroutine newtonInterpolation(X,Y,X0,Y0)
    implicit none
    real*8 :: X(:),Y(:),X0(:)
    real*8 :: Y0(size(X0))
    real*8 :: A(size(X)-1,size(X)-1),B(size(X)-1,1),M(size(X)-1),S(size(X)-1,1)
    real*8 :: xx,yy,li
    integer :: shape,shape0,k,i,j,p
    ! 判断X,Y的点数是否一致
    if (size(X) .ne. size(Y)) then
        write(*,*) "wrongShape"
        call abort()
    endif
    shape = size(X)
    shape0 = size(X0)

    ! 遍历 X0
    do k=1, shape0
        ! 遍历 X
        do i=1,shape-1
            do j=1,i
                ! 一定记得写注释！
                ! A(1,1) = (X2-X1)
                ! A(2,1) = (X3-X1)    A(2,1) = (X3-X1)(X3-X2)
                ! A(3,1) = (X4-X1)    A(3,1) = (X4-X1)(X4-X2)(X4-X3)
                ! ...
                A(i,j) = X(i+1)-X(1)
                    if (j.gt.1) then
                        do p=2,j
                            A(i,j) = A(i,j)*(X(i+1)-X(p))
                        enddo
                    endif
            enddo
        enddo
        B(:,1) = Y(2:)-Y(1)
        S = 0
        call downtri(A,B,S)
        M = 0

        do i=1,shape-1
                M(i) = 1
                do j=1,i
                    M(i) =M(i)*(X0(k)-X(j))
                enddo
        enddo
        Y0(k) = dot_product(M(:),S(:,1)) + Y(1)
    enddo
end subroutine newtonInterpolation

subroutine splineCurve(X,Y,X0,Y0)
    implicit none
    real*8 :: X(:),Y(:),X0(:)
    real*8 :: Y0(size(X0))
    real*8 :: A(size(X)-2,size(X)-2),h(size(X)-1,1),C(size(X)-2,1),M(size(X),1),mu(size(X)-2),la(size(X)-2)
    integer :: i,k,shape,shape0
    ! write(*,"(1f9.5)") Y
    ! write(*,*) Y
    M = 0
    shape = size(X)
    shape0 = size(X0)
    do i=1,shape-1
        h(i,1) = X(i+1)-X(i)
    enddo
    do i=1,shape-2
        A(i,i) = 2
        ! mu(i)
        mu(i) = h(i+1,1)/(h(i,1)+h(i+1,1))
        if (i .ne. shape-1) then
            A(i,i+1) = mu(i)
        endif
        ! la(i)
        la(i) = 1 - mu(i)
        if (i .ne. 1) then
            A(i,i-1) = la(i)
        endif
        ! C(i)
        C(i,1) = 3*((la(i)*(Y(i+1)-Y(i)))/h(i,1) + mu(i)*(Y(i+2)-Y(i+1))/h(i+1,1))
    enddo
    ! Solve M array
    ! M(0) and M(shape) = 0 
    call gausselimination(A,C,M(2:shape-1,1))

    i = 1
    do k=1,shape0
        do while (X(i+1) .lt. X0(k))
            i = i+1
            cycle
        enddo
        if (X(i+1) .lt. X0(k)) then
            i = i+1
        endif
        Y0(k)=1/h(i,1)**3*(h(i,1)+2*(X0(k)-X(i)))*(X0(k)-X(i+1))**2*Y(i)
        Y0(k)=Y0(k)+1/h(i,1)**3*(h(i,1)-2*(X0(k)-X(i+1)))*(X0(k)-X(i))**2*Y(i+1)
        Y0(k)=Y0(k)+1/h(i,1)**2*(X0(k)-X(i))*(X0(k)-X(i+1))**2*M(i,1)
        Y0(k)=Y0(k)+1/h(i,1)**2*(X0(k)-X(i+1))*(X0(k)-X(i))**2*M(i+1,1)
    enddo
    write(*,"(1f9.5)") Y0

end subroutine splineCurve


end module interpolation





program main
    use interpolation
    real*8 :: X(15),Y(15),X0(150),Y0(150)
    X = linspace(-5d0,5d0,15)
    X0 = linspace(-5d0,5d0,150)
    call samplingFromFunction(X,Y)
    call splineCurve(X,Y,X0,Y0)
    write(*,"(1f9.5)") Y0
end program main