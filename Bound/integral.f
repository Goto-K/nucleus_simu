      program integral_comp
      implicit none
c local
      integer   m,n               ! 刻み数n=2**m
      real      h,integral        ! 刻み幅h,数値積分の中間変数
      real      simpson           ! Simpson法の解
      real      x                 ! 積分変数
      real      exact             ! 解析解
      integer   i                 ! ループ用変数
c function:
      real      y                 ! 被積分関数
c begin:
      exact = exp(1.0)-1
      m=6
      n=2**m                      ! 刻み数
      h=1.0/n                     ! 刻み幅
c Simpson法
         integral=y(0.0)+y(1.0)   ! 積分端での値
         do i=1,n/2
            x = h*(2*i-1)
            integral=integral+4.0*y(x)
         end do
         do i=1,n/2-1
            x = h*(2*i)
            integral=integral+2.0*y(x)
         end do
         simpson=integral*h/3.0   ! Simpson法の解
      stop
      end
c
c 被積分関数の定義
c
      real function y(x)
      implicit none
c input:
      real  x
c begin:
      y = exp(x)
      return
      end
