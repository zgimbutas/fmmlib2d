cc Copyright (C) 2010: Leslie Greengard, Zydrunas Gimbutas, Vladimir Rokhlin
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the FMM tree plotting routines in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine plot_points2d(iw,z,n)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
c
        do 1200 i=1,n
        write(iw,1000) z(1,i),z(2,i)
 1200   continue
c
 1000   format(6(1x,e11.5))
c
        return
        end
c
c
c
c
c
        subroutine plot_box2d(iw,center,size)
        implicit real *8 (a-h,o-z)
        real *8 center(2)
c
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2
        write(iw,1200)
c
 1000   format(6(1x,e11.5))
 1200   format(80a1)
        return
        end 
c
c
c
c
c
        subroutine plot_label2d(iw,center,size,itag,label)
        implicit real *8 (a-h,o-z)
        real *8 center(2)
c        
        if( label .ge. 1 .and. label .lt. 10 )
     $     write(iw,1010) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 10 .and. label .lt. 100 )
     $     write(iw,1020) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 100 .and. label .lt. 1000 )
     $     write(iw,1030) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 1000 .and. label .lt. 10000 )
     $     write(iw,1040) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 10000 .and. label .lt. 100000 )
     $     write(iw,1050) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 100000 )
     $     write(iw,*) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
 1010   format(a11,i7,a3,i1,a5,e13.7,a3,e13.7,a7)
 1020   format(a11,i7,a3,i2,a5,e13.7,a3,e13.7,a7)
 1030   format(a11,i7,a3,i3,a5,e13.7,a3,e13.7,a7)
 1040   format(a11,i7,a3,i4,a5,e13.7,a3,e13.7,a7)
 1050   format(a11,i7,a3,i5,a5,e13.7,a3,e13.7,a7)
 1060   format(a11,i7,a3,i6,a5,e13.7,a3,e13.7,a7)

        return
        end
