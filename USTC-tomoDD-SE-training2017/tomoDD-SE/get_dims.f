      subroutine get_dims(x,y,xd,yd) 
      real x(4), y(4), xd, yd 
      integer i 
      real max_x, max_y 
 
      max_x=0.0 
      max_y=0.0 
 
      do i=1,4 
         if(abs(x(i)).gt.max_x) then 
            max_x=abs(x(i)) 
         endif 
         if(abs(y(i)).gt.max_y) then 
            max_y=abs(y(i)) 
         endif 
      enddo 
      xd=max_x 
      yd=max_y 

      end 
