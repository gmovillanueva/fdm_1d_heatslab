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

    # Stepper for solving the 
    while index <= ints:
        iradius = rlist[index]
        ipos = math.exp(alpha * iradius)
        ineg = math.exp(alpha * iradius * -1)
        tpos = ipos / iradius * c2
        tneg = ineg / iradius * c1
        analytical_sln[index, 0] = tneg + tpos
        index = index + 1

    return


# =============================================================================
# Numerical Function Solver
#
# Handles the solving part of the numerical solution
#
# Parameters
#
# arp = Refers to Area_(i+1/2)
# arn = Refers to Area_(i-1/2)
# arc = Refers to Centered Area
# vp = Refers to Volume_(i+1/2)
# vn = Refers to Volume_(i-1/2)
# cp = Refers to J_(i+1/2) constants
# cn = Refers to J_(i-1/2) constants
# cl  = Refers to Cell Left
# cc  = Refers to Cell Center
# cr  = Refers to Cell Right
# =============================================================================
def numerical(in_d_values):
    global matrix_a
    global vector_b

    # Defining multiple matrices
    matrix_a = numpy.zeros([ints + 1, ints + 1], dtype='float')
    vector_b = numpy.zeros(ints + 1, dtype='float')

    # Loop that builds both matrices
    for i in range(0, ints + 1):
        # Calculating all the i values
        cr = (i + 0.5)
        cc = i
        cl = (i - 0.5)

        # Calculating Radius
        rr = cr*dr
        rc = cc*dr
        rl = cl*dr

        # Calculating Areas
        arp = 4*pi*math.pow(rr, 2)
        arn = 4*pi*math.pow(rl, 2)
        arc = 4*pi*math.pow(rc, 2)

        # Calculating Volumes
        vp = (4/3)*pi*(math.pow(rr, 3) - math.pow(rc, 3))
        vn = (4/3)*pi*(math.pow(rc, 3) - math.pow(rl, 3))

        # Calculating Constants
        d_valueR = 1
        if i < ints:
            d_valueR = in_d_values[i]
        else:
            d_valueR = in_d_values[i - 1]

        d_valueL = 1
        if i > 0:
            d_valueL = in_d_values[i - 1]
        else:
            d_valueL = in_d_values[i]

        cp = d_valueR/dr
        cn = d_valueL/dr

        # Calculating Matrix Values
        aw = arn * cn
        ae = arp * cp
        ac = ae + aw + st * vp + st * vn

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
# Calculate D_value:
#
# Calculates D Values
#
# =============================================================================
def calc_dvalue(sigma, h, phi_old):
    st = sigma
    dr = h

    n = numpy.size(phi_old) - 1
    d_values = numpy.zeros(n)

    # Calculating Right Flux Limiter
    for i in range(0, n):
        d_valueR = 0
        fl_num = math.pow(phi_old[i + 1] - phi_old[i], 2) * 4
        fl_den = math.pow(phi_old[i + 1] + phi_old[i], 2) * math.pow(dr, 2)
        flimit = fl_num / fl_den
        d_valueR = 1 / math.sqrt(9 * math.pow(st, 2) + flimit)

        d_values[i] = d_valueR

    return d_values


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
# Calculate Balance:
#
# Balance should be roughly E-14
#
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
# Calculate FLD Solution:
#
#
# =============================================================================
def fld_solution():
    flux_guess = 1
    tol = 1E-8
    tol_check = 2
    iter_count = 1
    iter_limit = 100
    phi_new = numpy.zeros(ints + 1, dtype='float')
    phi_old = numpy.zeros([ints + 1, 1], dtype='float')
    phi_old.fill(flux_guess)
    while tol_check > tol:
        d_values = calc_dvalue(st, dr, phi_old)
        phi_new = numerical(d_values)
        tol_check = numpy.linalg.norm(phi_new - phi_old) / numpy.linalg.norm(phi_new)
        phi_old = phi_new.copy()
        if iter_count == iter_limit:
            print("Exited after: {} loops. Last calculated tolerance check: {}.".format(iter_count, tol_check))
            break
        if tol_check <= tol:
            print("Exited after: {} loops. Last calculated tolerance check: {}.".format(iter_count, tol_check))
        iter_count = iter_count + 1
    return phi_new


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
