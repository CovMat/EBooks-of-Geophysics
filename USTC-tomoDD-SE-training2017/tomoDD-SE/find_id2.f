	subroutine find_id(eve_sta,staID,evID,k)

	use tomoFDD
	integer eve_sta(MAXEVE,MAXOBS+1)
	integer staID,evID,k,i

	do i=1,eve_sta(evID,1)
	   if(staID.eq.eve_sta(evID,i+1)) then
	      k=i
              goto 200 
           endif
        enddo
        k=0
200     continue
	return
        end
