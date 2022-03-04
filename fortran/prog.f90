program shrodinger
  use mcf_tipos
  use mcf_diagonalizacion
  !- Puntu konpura diferentzia finituak
  integer,parameter  :: n=300
  !- Kontsideratuko dugun tartea [-L,L] da, unitate atomikoak erabiliz a0. 
  real(kind=sp), parameter    :: L=10.0_sp
  ! -Hamiltondarrar matrizearen diagonala : D
  ! -Hamiltondarrar goi eta behe diagonalak : D, E
  ! -puntu diskretuak : xi
  real(kind=sp), allocatable, dimension(:)   :: D, E, xi
  ! -Autobektorean Z(:,:) matrizean jasoko ditugu (zutabeka).
  real, allocatable, dimension(:,:) :: Z
  ! - Aldagai lagungarriak
  integer :: i,j, m, ind !m eta ind hau Beñat Berasategik sartu du
  real(kind=sp) ::  vi, h
  allocate(D(n), E(n), xi(n))
  allocate(Z(n,n))
  ! Pausua definitu 
  h=L/n
  ! Matrizearen diagonal eta goi eta behe diagonalak ezarri.
  open(unit=110,file="pot2.dat",status="replace",action="write")
  do i=1, n
     xi(i) = -L + 2.0_sp*L*real(i-1)/real(n-1)
     vi =  eskailera_pot(xi(i))
     D(i) = (vi + 2.0_sp/h**2) 
     E(i) = -1.0_sp/h**2
     write(unit=110,fmt="(100f18.6)")xi(i), vi , d(i), e(i)
  end do
  close(unit=110)
  ! Hasiera batean Z=1 egin behar dugu matrizea 
  ! zuzenean diagonala bada. Begiratu inplementazio apunteak!
  z=0.0_sp
  do i=1,n
     z(i,i)=1.0_sp
  end do
  !- Matrize hiru-diagonala diagonalizatzeko.
  call  TQLI (D,E,Z)
  !- Energiak ordenatzeko txikienetik handienera
  call  EIGSRT(D,Z)
  !- Lortutako energiak idatzi
  open(unit=109,file="energiak2.dat",status="replace",action="write")
  do i=1, n
     write(unit=109,fmt="(100f12.6)")D(i)
  end do
  close(unit=109)

  !ZEINUAREN KONTUA KONPONDU NAHIAN NIK SARTUTAKOA----------------------------
  !Hasieran negatiboa bada m=1 jarriko dugu, bestela, m=0. Horrela, (-1)^m faktoreak zeinuarena
  !konponduko digu. "Hasierakoa" zehazteko, finkatutako n-ren %10a hartuko dugu, ea funtzionatzen
  !duen. Bestela, aldaketak egin beharko ditugu.
  ind=floor(0.1*n)
  if ( Z(1,ind) > 0.0) then
          m=0
  else if (Z(1,ind) <= 0.0) then
          m=1
  end if
  !---------------------------------------------------------------------------
  open(unit=110,file="uhin_funtzioak2.dat",status="replace",action="write")
  do i=1, n
     write(unit=110,fmt="(100f12.6)") xi(i), ( (-1.0)**m*Z(i,j), j=1, 20)  !(-1)^m guk sartu dugu.
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
  !floor erabilita
  function f_pot(x)
    real(kind=sp),intent(in) :: x
    real(kind=sp)            :: f_pot
    f_pot=floor(x**2)
  end function f_pot

  function eskailera_pot(x)
     real(kind=sp), intent(in) :: x
     real(kind=sp)             :: eskailera_pot
     integer                   :: n, i
     real(kind=sp)             :: E, a, k, S, muga

     n=10
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

    ! print*, S
     eskailera_pot=S

  end function eskailera_pot
end program shrodinger

