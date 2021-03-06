import math

#Electrostatic force parameters
Vi = 5.0#Work function of wall (Volt)
Vj = 4.9 #Work function of particle (Volt)
eps0 = 8.885e-12 #permitivity of air (C^2/N.m^2)
pois = 0.3
ymod = 1.0e7
dens = 890.0
parDia = 140.0e-6
ha = 6.5e-20

charge = -6.0e-9
ks = 1.0e-4
vel = 2.0
vDash = 0.0
imageConst = 2.0
Zs =300.0e-9 
noOfImpacts = 80
alpha = 0.9
emod = alpha*(1.0-pois*pois)/ymod
k0 = imageConst*Zs/eps0/(4.0*math.pi*pow(parDia*0.5,2))
S = 1.36*pow(emod, 2.0/5.0)*pow(dens,2.0/5.0)*pow(parDia,2)*pow(vel,4.0/5.0)
vDash = k0*charge 
ks = 1.8e-5/(Vi-Vj-vDash)
k1 = ks*Zs/eps0
Qmax = (Vi-Vj)/k0 # (Matsusaka et al, 2000)
minGap = 1.0e-6
print("S ",S)
print("k0 ks k1 ",k0, ks,k1)
print("Qmax ",Qmax)


surf_tens = 0.013
liq_vol = 9.59E-9
cont_ang =0.0

ks = 1.0e-4

f = open("charge-impact.dat","w")
capf = 0.0
gap = 0.0
fv = 0.0
for i in range(noOfImpacts):
    #Electrostatic force
    vDash = k0*charge  
    deltaV = Vi - Vj - vDash
    deltaQ = ks*S*deltaV
    charge += deltaQ

    gap = max(gap, 1e-6)
    esF = charge*charge/(4.0*math.pi*eps0*gap*gap)

    f.write(str(i+1)+" "+str(round(charge*1e9,3))+"\n")
    print(round(charge*1e9,3))


    # Vanderwaal force
    gap = 0
    vGapMn = 1.0e-9
    gap = max(gap, vGapMn)
    fv = -ha*parDia*0.5/(6.*gap*gap)


    # Capillary force
    s_rupture = (1.0+0.5*cont_ang)*pow(liq_vol,1.0/3.0)
    if(gap < s_rupture):
        sepMin = 5.0e-6; #(Hornbaker, Albert et al. 1997, Nase, Vargas et al. 2001)
        separation = max(sepMin, gap)
        rStar = parDia*0.5
        capf = -2.*math.pi*rStar*surf_tens*math.cos(cont_ang)/(1.0+1.0/math.sqrt(1.0+2.*liq_vol/(math.pi*rStar*pow(separation,2))-separation))
minGap = 10.0e-6
print("Maximum Capillary force ", capf)        
print("Maximum Vanderwaal force ", fv)
print("Max Electrostatic Force ", Qmax*Qmax/(4.0*math.pi*eps0*minGap*minGap))
f.close()