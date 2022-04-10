program shrodinger
  use mcf_tipos
  use mcf_diagonalizacion

  !- Puntu konpurua:
  integer,parameter  :: n=300
  !- Kontsideratuko dugun tartea [-L,L] da, unitate naturaletan. 
  real(kind=sp), parameter    :: L=10.0_sp
  ! -Hamiltondarrar matrizearen diagonala : D
  ! -Hamiltondarrar goi eta behe diagonalak : E
  ! -puntu diskretuak : xi
  real(kind=sp), allocatable, dimension(:)   :: D, E, xi
  ! -Autobektorean Z(:,:) matrizean jasoko ditugu (zutabeka).
  real, allocatable, dimension(:,:) :: Z
  ! - Aldagai lagungarriak
  integer :: i,j, m, ind, esk_kop
  real(kind=sp) ::  vi, h

  character(len=15) :: fitx_izen1, fitx_izen2, fitx_izen3, n_str

  allocate(D(n), E(n), xi(n))
  allocate(Z(n,n))
  ! Pausua definitu 
  h=L/n

  !Eskailera kopurua eskatu:
  print*, "Sartu eskailera kopurua"
  read*, esk_kop
  
  !Exekuzioari dagozkion fitxategiak sortu:
  write(unit=n_str, fmt="(i10)") esk_kop
  fitx_izen1="pot_"//trim(adjustl(n_str))//".dat"
  fitx_izen2= "E_"//trim(adjustl(n_str))//".dat"
  fitx_izen3= "uf_"//trim(adjustl(n_str))//".dat"

  ! Matrizearen diagonal eta goi eta behe diagonalak ezarri.
  open(unit=110,file=fitx_izen1, status="replace",action="write")
  do i=1, n
     xi(i) = -L + 2.0_sp*L*real(i-1)/real(n-1)
     vi =  eskailera_pot(xi(i),esk_kop)
     D(i) = (vi + 2.0_sp/h**2) 
     E(i) = -1.0_sp/h**2
     write(unit=110,fmt="(100f18.6)")xi(i), vi , d(i), e(i)
  end do
  close(unit=110)

  ! Hasiera batean Z=1 egin behar dugu matrizea 
  ! zuzenean diagonala bada.
  z=0.0_sp
  do i=1,n
     z(i,i)=1.0_sp
  end do
  !- Matrize hiru-diagonala diagonalizatzeko.
  call  TQLI (D,E,Z)
  !- Energiak ordenatzeko txikienetik handienera
  call  EIGSRT(D,Z)
  !- Lortutako energiak idatzi
  open(unit=109,file=fitx_izen2,status="replace",action="write")
  do i=1, n
     write(unit=109,fmt="(100f12.6)")D(i)
  end do
  close(unit=109)

  !ZEINUAK KONPONTZEKO-------------------------------------------------------------------------------
  !Autofuntzio guztiek orientazio berdina izan dezaten "eskuz" egin behar izan dugu hau. Izan ere,
  !hau egin ezean, gora/behera orientatzen ziren guk kontrolatu ahal izan gabe.
  !Hasieran negatiboa bada m=1 jarriko dugu, bestela, m=0. Horrela, (-1)^m faktoreak zeinuarena
  !konponduko digu. "Hasierakoa" zehazteko, finkatutako n-ren %6a hartuko dugu.
  ind=floor(n*0.6)
  if ( Z(ind,1) >= 0.0) then
          m=1
  else if (Z(ind,1) < 0.0) then
          m=-1
  end if
  !---------------------------------------------------------------------------
  open(unit=110,file=fitx_izen3,status="replace",action="write")
  do i=1, n
     write(unit=110,fmt="(100f12.6)") xi(i), ( Z(i,j)*m, j=1, 20)  
  end do
  deallocate(D, E, xi)
  close(unit=110)

contains
  ! potentzial harmonikoa
  function v_pot(x)
    real(kind=sp),intent(in) ::  x
    real(kind=sp)            ::  v_pot
    ! v =1/2 * w^2 x**2 
    ! w = sqrt(2) 
    v_pot = x**2
  end function v_pot

  ! eskailera potentziala
  function eskailera_pot(x,n)
     real(kind=sp), intent(in) :: x
     integer, intent(in)       :: n
     real(kind=sp)             :: eskailera_pot
     integer                   :: i
     real(kind=sp)             :: E, a, k, S, muga

     E=100.0
     a=E/n
     k=2.0

     S=0.0
     do i=1,n
        muga=0.5*k*x**2/a+0.5
       ! print*, muga
        if (muga>=real(i)) then
           S=S+a
        else if (muga<real(i)) then
           S=S+0.0
        end if
     end do

     eskailera_pot=S

  end function eskailera_pot
end program shrodinger

