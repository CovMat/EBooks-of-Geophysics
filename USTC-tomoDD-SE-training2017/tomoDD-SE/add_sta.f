        subroutine add_sta(eve_sta,evID,staID)

 	use tomoFDD
        integer eve_sta(MAXEVE, MAXOBS+1)
        integer evID, staID
        integer i, j

        do i=1, eve_sta(evID,1)
           if(eve_sta(evID,i+1).eq.staID) goto 119
        enddo
        eve_sta(evID,1)=eve_sta(evID,1)+1
        j=eve_sta(evID,1)
        eve_sta(evID,j+1)=staID
119     continue
        return
        end

