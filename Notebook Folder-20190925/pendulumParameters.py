#File pendulumParameters.py

stringDiameter=0.031*0.0254 #meters
stringLength=0.627 #meters, from pivot to top of mass
stringDensity=0.000273/0.738 #kg/m
mMass=0.024245-1.260*stringDensity #kg just the mass (mass+string=24.245 g)
dMass=0.02539 #m diameter of mass, measured with caliper
stringDensity=0.00133/1.31 #kg/m, mass of entire length of string-1.33 g = mass of 1.31 m of string
stringMass=stringDensity*stringLength

#calculate length of pendulum
l=stringLength+dMass/2

#moments of Inertia
IMass=mMass*l*l  #moment of inertia of mass about pivot
IString=1.0/3*stringMass*stringLength**2 #1/3 mL**2 is Moment of intertia of rod about one end
I=IMass+IString

#center of mass
m=mMass+stringMass
lcm=(mMass*l+stringMass*stringLength/2)/m

print('Total Mass ',m)
print('lcm ',lcm,' l ',l)
print('I, IMass, IString:',I, IMass, IString)
