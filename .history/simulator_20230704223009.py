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
        rLen = cr*dr
        cLen = cc*dr
        lLen = cl*dr

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

    return num_flux


# =============================================================================
# Calculate Limited J Value:
#
#
# =============================================================================
def calc_jvalue_limited(sigma, h, phi_old):
    st = sigma
    dr = h

    n = numpy.size(phi_old) - 1
    x_values = numpy.zeros(n)
    j_values = numpy.zeros(n)

    # Calculating Right Flux Limiter
    for i in range(0, n):
        d_valueR = 0
        fl_num = math.pow(phi_old[i + 1] - phi_old[i], 2) * 4
        fl_den = math.pow(phi_old[i + 1] + phi_old[i], 2) * math.pow(dr, 2)
        flimit = fl_num / fl_den
        d_valueR = 1 / math.sqrt(9 * math.pow(st, 2) + flimit)
        #d_valueR = 0

        j_values[i] = -d_valueR * (phi_old[i + 1] - phi_old[i]) / dr
        x_values[i] = 0.5 * h + i * h

    return x_values, j_values


# =============================================================================
# Calculate Non-Limited J Value::
#
# Calculates D Values
#
# =============================================================================
def calc_jvalue_nonlimited(sigma, h, phi_old):
    st = sigma
    dr = h

    n = numpy.size(phi_old) - 1
    x_values = numpy.zeros(n)
    j_values = numpy.zeros(n)

    # Calculating Right Flux Limiter
    for i in range(0, n):
        d_valueR = 1 / (3*st)
        #d_valueR = 0

        j_values[i] = -d_valueR * (phi_old[i + 1] - phi_old[i]) / dr
        x_values[i] = 0.5 * h + i * h

    return x_values, j_values


# =============================================================================
# Calculate Residual:
#
# Balance should be roughly E-14
#
# =============================================================================
def calc_residual(vector_phi):
    total_rez = 0
    for i in range(0, ints + 1):
        rez = vector_b[i] - numpy.dot(matrix_a[i, :], vector_phi[:])
        total_rez += rez
    return total_rez


# =============================================================================
# Calculate Boundaries
# =============================================================================
def calc_balance(phi_old, rez, flag):
    global st
    global q_value
    global dr

    h = dr

    diff_A = 0

    n = numpy.size(phi_old) - 1
    ajL_values = numpy.zeros(n)
    ajR_values = numpy.zeros(n)
    aL_values = numpy.zeros(n)
    aC_values = numpy.zeros(n)
    aR_values = numpy.zeros(n)
    vL_values = numpy.zeros(n)
    vR_values = numpy.zeros(n)

    for i in range(0, n):
        xR_value = i * h + 0.5 * h
        xC_value = i * h
        xL_value = i * h - 0.5 * h

        aR_values[i] = 4 * pi * math.pow(xR_value, 2)
        aC_values[i] = 4 * pi * math.pow(xC_value, 2)
        aL_values[i] = 4 * pi * math.pow(xL_value, 2)

        vR_values[i] = (4/3)*pi*(math.pow(xR_value, 3) - math.pow(xC_value, 3))
        vL_values[i] = (4/3)*pi*(math.pow(xC_value, 3) - math.pow(xL_value, 3))

    if flag == 1:
        for i in range(0, n):
            d_value = 1 / (3 * st)

            if (i > 0) and (i < n):
                jR_value = -d_value * (phi_old[i + 1] - phi_old[i]) / dr
                jL_value = -d_value * (phi_old[i] - phi_old[i - 1]) / dr
                ajL_values[i] = aL_values[i] * jL_value
                ajR_values[i] = aR_values[i] * jR_value
            elif i == 0:
                ajL_values[i] = q_value
            elif i == n:
                ajR_values[i] = aC_values[i] * 0.5 * phi_old[i]

    elif flag == 2:
        for i in range(0, n):
            d_value = 0
            fl_num = math.pow(phi_old[i + 1] - phi_old[i], 2) * 4
            fl_den = math.pow(phi_old[i + 1] + phi_old[i], 2) * math.pow(dr, 2)
            flimit = fl_num / fl_den
            d_value = 1 / math.sqrt(9 * math.pow(st, 2) + flimit)

            if (i > 0) and (i < n):
                jR_value = -d_value * (phi_old[i + 1] - phi_old[i]) / dr
                jL_value = -d_value * (phi_old[i] - phi_old[i - 1]) / dr
                ajL_values[i] = aL_values[i] * jL_value
                ajR_values[i] = aR_values[i] * jR_value
            elif i == 0:
                ajL_values[i] = q_value
            elif i == n:
                ajR_values[i] = aC_values[i] * 0.5 * phi_old[i]

    cell_bal = 0
    bal_den = 0
    for i in range(0, n):
        bal_den = 1 / (ajR_values[i] - ajL_values[i])

    bal = abs(rez / bal_den)
    return bal


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

radius = 1
# radius = 10
ints = 1000
q_value = 1
src_value = 0
sigma_a = 1
sigma_s = 0
#sigma_a = 0
#sigma_s = 10
sigma_t = sigma_a + sigma_s
st = sigma_a + sigma_s
dr = radius/ints
xlist = numpy.arange(0.0, radius, dr)
rlist = numpy.linspace(0.0, radius, ints + 1)
pi = math.pi

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
