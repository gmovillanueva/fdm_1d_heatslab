! ------------------------------------------------------
! Simulator implementation of Non-Steady Stat
! ------------------------------------------------------
program FuelCladdingSimulator
  implicit none
  
  ! Fuel Props
  real :: kFuel, pFuel, cFuel, qFuel, lenFuel 

  ! Cladding Props
  real :: kCladding, pCladding, cCladding, qCladding, lenCladding 

  ! Convection Coefficient
  real :: hCoeff

  ! Reading inputs from users
  read(*,*) kFuel, pFuel, cFuel, qFuel, lenFuel, kCladding, pCladding, cCladding, qCladding, lenCladding, hCoeff

  write(*,*) "Fuel Thermal Conductivity (W/mK) = ", kFuel
  write(*,*) "Fuel Density (kg/m^3) = ", pFuel
  write(*,*) "Fuel Heat Capacity (J/kgK) = ", cFuel
  write(*,*) "Fuel Heat Density (kW/m^3) = ", kFuel
  write(*,*) "Fuel Length = (cm) ", kFuel
  write(*,*) "Cladding Thermal Conductivity (W/mK) = ", kCladding
  write(*,*) "Cladding Density (kg/m^3) = ", pCladding
  write(*,*) "Cladding Heat Capacity (J/kgK) = ", cCladding
  write(*,*) "Cladding Heat Density (kW/m^3) = ", kCladding
  write(*,*) "Cladding Length = (cm) ", kCladding

  ! Int

  ! Output to file
  !open(1, file = '')

end program FuelCladdingSimulator

