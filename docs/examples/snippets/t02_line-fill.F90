h = 1._R_P/n
do i=1, n                                     ! cell i is [(i-1) h, i h]: its average of sin(2 pi x)
   u(i) = (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h)
enddo
u(1-NG:0)   = u(n+1-NG:n)                     ! periodic ghost cells
u(n+1:n+NG) = u(1:NG)
