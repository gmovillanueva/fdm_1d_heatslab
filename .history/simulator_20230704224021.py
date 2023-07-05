# =============================================================================
# Standard Python Imports
# =============================================================================
import numpy
import math

# =============================================================================
# Analytical Function Solver:
# =============================================================================
def analytical():
    # Would contain the expressions to solve

    # Stepper for solving the values
    # Example of statement below
    # while index <= ints:
        
    return


# =============================================================================
# Numerical Function Solver
# =============================================================================
def numerical():
    global matrix_a
    global vector_b

    # Init. multiple matrices
    matrix_a = numpy.zeros([ints + 1, ints + 1], dtype='float')
    vector_b = numpy.zeros(ints + 1, dtype='float')

    # Loop that builds both matrices
    for i in range(0, ints + 1):
        # Calculating all the i values
        cr = (i + 1)
        cc = i
        cl = (i - 1)

        # Calculating length
        rLen = cr*dx
        cLen = cc*dx
        lLen = cl*dx

        # Calculating Matrix Values

        # Building Matrix A
        if (i > 0) and (i < ints):
            matrix_a[i, i - 1] = -aw
            matrix_a[i, i] = ac
            matrix_a[i, i + 1] = -ae
            vector_b[i] = 0
        elif i == 0:
            ac = ae + st * vp
            matrix_a[i, i] = ac
            matrix_a[i, i + 1] = -ae
            vector_b[i] = q_value
        elif i == ints:
            ac = aw + st * vn + arc * 0.5 / d_valueR
            matrix_a[i, i - 1] = -aw
            matrix_a[i, i] = ac
            vector_b[i] = 0

    # Solving Ax = b
    num_flux = numpy.linalg.solve(matrix_a, vector_b)

    return


# =============================================================================
# Parameters
#
# kFuel         := Fuel Thermal Conductivity (W/mK) = 
# pFuel         := Fuel Density (kg/m^3)
# cFuel         := Fuel Heat Capacity (J/kgK)
# qFuel         := Fuel Heat Density (kW/m^3)
# lenFuel       := Fuel Length cm)
# kCladding     := Cladding Thermal Conductivity (W/mK)
# pCladding     := Cladding Density (kg/m^3)
# cCladding     := Cladding Heat Capacity (J/kgK)
# qCladding     := Cladding Heat Density (kW/m^3)
# lenCladding   := Cladding Length (cm)
# hCoeff        := Convection Coefficient (kW/m^2K)
# timeSteps     := Time Steps
# fileSave      := Save to file
# =============================================================================

kFuel = input("Fuel Thermal Conductivity (W/mK) = ")      
pFuel = input("Fuel Density (kg/m^3) = ")        
cFuel = input("Fuel Heat Capacity (J/kgK) = ")        
qFuel = input("Fuel Heat Density (kW/m^3) = ")        
lenFuel = input("Fuel Length = ")      
kCladding = input("Cladding Thermal Conductivity (W/mK) = ")    
pCladding = input("Cladding Density (kg/m^3) = ")    
cCladding = input("Cladding Heat Capacity (J/kgK) = ")    
qCladding = input("Cladding Heat Density (kW/m^3) = ")    
lenCladding = input("Fuel Thermal Conductivity (W/mK) = ")  
hCoeff = input("Convection Coefficient (kW/m^2K) = ")       
timeSteps = input("Time Steps = ")   
fileSave = input("Export to file?") 

dx = (lenFuel + lenCladding) / timeSteps
