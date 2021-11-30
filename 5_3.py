
#loads and variables declaration:




#classes:
#material
#propellant
class Material:
    def __init__(self,E,sigmaYield,density,sigmaYieldShear):
        self.density = density
        self.E = E
        self.sigmaYield = sigmaYield
        self.sigmaYieldShear = sigmaYieldShear
        self.v = 1/3

class Propellant:
    def __init__(self, mass,molarMass):
        self.mass = mass
        self.molarMass = molarMass

class Spacecraft:
    def __init__(self,diameter,radius,height):
        self.diameter = diameter
        self.radius = 0.5 * diameter
        self.height = height

class LaunchAccel:
    def __init__(self, long_static, long_dynamic, lat_static, lat_dynamic, SafetyFactor):
        self.long_static = long_static
        self.long_dynamic = long_dynamic
        self.lat_static = lat_static
        self.lat_dynamic = lat_dynamic
        self.lat_max = lat_static + lat_dynamic
        self.long_max = long_static + long_dynamic
        self.SafetyFactor = SafetyFactor





#functions:
#calcEulerBucklingCriticalStress
def calcEulerBucklingCriticalStress(material,R,t,L):
    E = material.E
    I = calcMomentOfInertiaCylinder(R,t)
    A = calcAreaCylinder(R,t)
    sig_crit_BCS= ((math.pi**2) * E * I)/(A*(L**2))
    return sig_crit_BCS

#calcShellBucklingStress (v is poisson ratio)
def calcShellBucklingStress(material, R, t_1, L, p):
    E = material.E
    v = material.v
    Q = calcQ(p, R, E, t_1)
    k = calcK(L, v, R, t_1)


    sig_crit_SBS = (1.983-0.983*math.exp(-23.14*Q))*k*((math.pi**2 * E * t_1**2)/(12*(1-v**2)*L**2))
    return sig_crit_SBS

#calcQ
def calcQ(p, R, E, t_1):
    Q = (p*R**2)/(E*t_1**2)
    return Q

#calcK (v is poisson ratio)
def calcK(L, v, R, t_1):
    k = 0
    for lamb in range(1,10):
        if(k==0):
            k = lamb + (12*(L**4)*(1-(v**2)))/((math.pi**4)*(R**2)*(t_1**2)*lamb)
        tempK = lamb + (12*(L**4)*(1-(v**2)))/((math.pi**4)*(R**2)*(t_1**2)*lamb)
        if(tempK<k):
            k = tempK
    return k


#calcTankMass
def calcTankMass(R,t,l,material):
    #cylinder part:
    volumeCylinder = calcAreaCylinder(R,t) * l

    #hemispherical parts:
    volumeHemispheres = ((4/3 * math.pi * ((R+t)**3))-(4/3 * math.pi * (R**3)))

    volumeTotal = volumeHemispheres + volumeCylinder
    return volumeTotal * material.density

#calcMomentOfInertia
def calcMomentOfInertiaCylinder(R,t):
    return math.pi*(R*3)*t

#calcArea
def calcAreaCylinder(R,t):
    return ((R+t)**2-(R**2))*math.pi

#PLAN
#iterate through R and t
#Calculate the mass of these combinations
# With this mass, work out the forces on the tank
#With the forces calculate the experienced stresses
#Check the critical stresses (sigma>sigma_critical = failure)
# Setup fail pass condition
# Search with binary search for optimized combination


def iterate(material,propellant):
    list_R = []
    list_t = []
    for R in np.arange(0.001, spacecraft.radius, 0.001):
        L = 3
        t = 0
        dt = 0.001
    #check condition to see if the requirements are met
        check = False
        while not check:
            t = t + dt
            mass = calcTankMass(R, t, L, testAluminium) + propellant.mass
            Area = calcAreaCylinder(R, t)
            pressure = calcPressure(R,L,propellant,temp)
            Forces = calcLaunchLoads(mass, launchAccel)
            sigma_Axial = calcsigma_A(Forces[0], Area)


            sigma_cr_Euler = calcEulerBucklingCriticalStress(material,R,t,L)
            sigma_cr_shell = calcShellBucklingStress(material,R,t,L,pressure)



            if (sigma_Axial < sigma_cr_Euler) and (sigma_Axial < sigma_cr_shell) and (sigma_Axial < material.sigmaYield):
                check = True
        print(str(sigma_Axial) + " " + str(sigma_cr_Euler) + " " + str(sigma_cr_shell) + " " + str(mass) + " " + str(R) + " " + str(t) + " "+ str(pressure/100000))
        list_R.append(R)
        list_t.append(t)


    return [list_R,list_t]

#calcLaunchLoads
def calcHoopStress(R,t,int_pressure):

    return(int_pressure*R/t)

def calcAxialStress(R,t,int_pressure):

    return (int_pressure * R / (2*t))


#Class Values
launchAccel = LaunchAccel(4.55, 1, 0.25, 0.8, 2.5)
spacecraft = Spacecraft(3.6, 1.8, 5)
testAluminium = Material(87000000000, 300000000, 2700, 200000000)
propellant = Propellant(817.75,131.293)

#assumptions
temp = 295 #K

print(iterate(testAluminium,propellant)[1])

print(testAluminium.density)