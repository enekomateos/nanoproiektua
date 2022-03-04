module mcf_tipos
! Luzapenean oinarritutako mota estandarrak definitzen ditu (zenbaki osoetarako)
! eta zehaztasunean (zifra esanguratsuen kopurua) errealetarako.
!
! Kalkulu zientifikoa egiteko programa batean, interesgarria da
! behar den zehaztasun adierazi ahal izatea. "selected_int_kind" 
! eta "selected_real_kind" funtzioak erabiltzeak mota egokia 
! aukeratzeko bermea ematen du konputagailuaren arkitektura 
! ezagutu gabe.

!
! Erabilera
!---------------------------------
! program adibidea
! use tipos
!
! real(kind=doble)    :: x
! integer(kind=byte)  :: n
!   ...
! end program adibidea
!----------------------------------
!
!=======================================================================
! Osoak (10en berredurarekin adierazten dira)
!
integer, parameter, public :: int2   = selected_int_kind(2)    
integer, parameter, public :: int4   = selected_int_kind(4)  
integer, parameter, public :: int8   = selected_int_kind(8)  
integer, parameter, public :: int10  = selected_int_kind(10)   

integer, parameter, public :: byte  = int2
integer, parameter, public :: short = int4
integer, parameter, public :: int   = int8
integer, parameter, public :: long  = int10

!
! Tipo errealak
!
integer, parameter, public :: single = selected_real_kind(6) 
integer, parameter, public :: double = selected_real_kind(14)
!
! Beste izen batzuk tipo berdinentzat
!
integer, parameter, public :: sencillo = single
integer, parameter, public :: doble = double
!
integer, parameter, public :: sp = single
integer, parameter, public :: dp = double
integer, parameter, public :: qp=selected_real_kind(33)

end module mcf_tipos
