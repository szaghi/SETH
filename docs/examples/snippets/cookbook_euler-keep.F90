do i=0, n                                    ! d = 1: the face's normal is x
   call seth_keep_pe6(1, GAMMA, wg(:,i-2:i+3), f(:,i))
enddo
