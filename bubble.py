# Axis y-rows, x-columns - origin at left-top corner


# Bubble
# Vitor Vilela
# 14 October 2016

import numpy as np
import math
import matplotlib.pyplot as plt
from skimage import io
import random


# Turn on 
onDrag = 1.
onSaffman = 1.
onVirtualMass = 0.
onWeightBuoyancy = 0.

# Simulation
timestep = 1.e-3
simulation_time = 30.
numberOfparticles = 1


# Pipe radius
R = 1.
# Pipe section length
L = R
# Reference velocity
U = 1.
# Gravity
gravity = -10.

# Bubble modelled as particle
particle_diameter = 1.e-1
particle_density = 1.2
water_viscosity = 2.5e-4
water_density = 998.

# Particle data type
# Time 2 > Time 1
# X and Y correspond to R and Z directions, respectively
class Particle:
  def __init__(self, x1, y1, x2, y2, vx1, vy1, vx2, vy2):
    self.x1 = x1
    self.y1 = y1
    self.x2 = x2
    self.y2 = y2
    self.vx1 = vx1
    self.vy1 = vy1
    self.vx2 = vx2
    self.vy2 = vy2   

# Half-pipe mesh
centered_grid = 200
delta = R/centered_grid

# Centered velocity
vz = np.zeros((centered_grid,centered_grid))
# Centered vorticity
vort = np.zeros((centered_grid,centered_grid))



# Vz analytical function
def Vz(r):
  return U*(1-pow((r/R),2.))

# Vr analytical function
def Vr():
  return 0.

# Vorticity analytical function
def Vorticity(r):
  return 2*U*r/(R*R)



# Discretized - Eulerian velocity field
# Axis y-rows, x-columns - origin at left-top corner

# Flow Centered Velocity
for i in range(centered_grid):
  r = delta*(i+0.5)
  for j in range(centered_grid):  
    vz[j,i] = Vz(r)
 
# Analytical Vorticity
for i in range(centered_grid):
  r = delta*(i+0.5)
  for j in range(centered_grid):  
    vort[j,i] = Vorticity(r)    

# Show Eulerian fields: vz and vort   
#io.imshow(vort)
#io.show()



# Bubble's volume
def Vol(radius):
  return (4/3)*math.pi*pow(radius,3.)

# Bubble's mass
def Mp(radius, density):
  return density*Vol(radius)

# Component of the relative velocity vector VR (i = z or r)
def VRi(v, vp):
  return v-vp

# Module of the relative velocity vector
def VR(vr1, vr2):
  return math.sqrt(vr1*vr1+vr2*vr2)

# Bubbles' Reynolds number
def Re(medium_density, medium_viscosity, VR, Dp):
  #print 'Re: ', medium_density*VR*Dp/medium_viscosity
  return medium_density*VR*Dp/medium_viscosity

# Bubble's diameter
def Dp():
  return particle_diameter

# Drag Coefficient
def Cd(Re):
  if Re > 1.e-3:
    Cd = 24./Re
  else:
    Cd = 0.
  return Cd

# Drag force (z or r direction depending on VRi)
def drag(medium_density, particle_density, Dp, Cd, VR, VRi):
  return 3*medium_density*Cd*VR*VRi/(4*particle_density*Dp)

# Weight-Buoyancy force
def weightBuoyancy(radius, particle_density, medium_density, gravity):
  weight = Mp(radius, particle_density)*gravity
  buoyancy = Vol(radius)*medium_density*gravity
  return buoyancy-weight
  
# Saffman force (z or r direction depending on VRi)
def saffman(Dp, Mp, medium_density, medium_viscosity, WZ, VRi):
  return 1.61*Dp*Dp*math.sqrt(medium_density*medium_viscosity*abs(WZ))*VRi*WZ

# Temporal derivative of the flow velocity (z or r direction depending on v1 and v2)
def DVpDt(v1, v2, timestep):
  return (v2-v1)/timestep

# Virtual mass force (z or r direction depending on DVpDt)
def virtualMass(medium_density, particle_velocity, dVpdt):
  return -0.5*medium_density*particle_velocity*dVpdt

# Euler Temporal Integration - Velocity
def Euler(vp0, timestep, rhs):
  return vp0 + timestep*rhs

# Update particle position
def Move(p0, timestep, v2):
  return p0 + timestep*v2



         
    
# Print particle's position
fig, ax = plt.subplots()
ax.imshow(vz, interpolation='nearest', cmap=plt.cm.gray)  
        
# Particle initialization     
particles = []
random.seed(numberOfparticles)
for p in range(numberOfparticles): 
  x = random.random()
  y = random.random()
  vx = Vr()
  vy = Vz(x)
  particles.append(Particle(x, y, x, y, vx, vy, vx, vy))
  ax.plot(particles[p].x1*centered_grid, particles[p].y1*centered_grid, 'g.')



# Forces acting on particles an time integration

t = 0.
ct = 0

while t < simulation_time:

  t = t + timestep
  print 'time: ', t
  ct = ct + 1
  
  for p in range(numberOfparticles):

    # Move particle's position (x and y)
    particles[p].x2 = Move(particles[p].x1, timestep, particles[p].vx2)
    particles[p].y2 = Move(particles[p].y1, timestep, particles[p].vy2)  

    # Find new particle's velocity integrating ODE (Newton's law) in time
    particle_mass = Mp(0.5*particle_diameter, particle_density)
    VRx = VRi( Vr(), particles[p].vx1 )
    VRy = VRi( Vz(particles[p].x2), particles[p].vy1 )
    relative_velocity = VR(VRx, VRy)
    vorticity = Vorticity(particles[p].x1)  
    Reynolds = Re(water_density, water_viscosity, relative_velocity, particle_diameter)
    drag_coefficient = Cd(Reynolds)    
        
    # DVpDt are based on present and past times  
    RHSx = onDrag*drag(water_density, particle_density, particle_diameter, drag_coefficient, relative_velocity, VRx) + onSaffman*saffman(particle_diameter, particle_mass, water_density, water_viscosity, vorticity, VRy) + onVirtualMass*virtualMass(water_density, particles[p].vx1, DVpDt(particles[p].vx1, particles[p].vx2, timestep))
    RHSy = onWeightBuoyancy*weightBuoyancy(0.5*particle_diameter, particle_density, water_density, gravity) + onDrag*drag(water_density, particle_density, particle_diameter, drag_coefficient, relative_velocity, VRy) + onSaffman*saffman(particle_diameter, particle_mass, water_density, water_viscosity, vorticity, VRx) + onVirtualMass*virtualMass(water_density, particles[p].vy1, DVpDt(particles[p].vy1, particles[p].vy2, timestep))
    
    particles[p].vx2 = Euler(particles[p].vx1, timestep, RHSx)
    particles[p].vy2 = Euler(particles[p].vy1, timestep, RHSy)  
    
    # Particle BC at Pipe Centre (Symmetric)
    if particles[p].x2 < 0.:
      particles[p].x2 = -particles[p].x2
      
    # Particle BC at Pipe Wall (Reflect)
    if particles[p].x2 > R:
      particles[p].x2 = 2*R - particles[p].x2
    
    # Periodicity 
    if particles[p].y2 < 0.:
      particles[p].y2 = L + particles[p].y2 
    if particles[p].y2 > L:
      particles[p].y2 = particles[p].y2 - L        
       
    # Update
    particles[p].x1 = particles[p].x2
    particles[p].y1 = particles[p].y2
    particles[p].vx1 = particles[p].vx2
    particles[p].vy1 = particles[p].vy2    
  
  
  print 'ct: ', ct
  if ct%100 == 0:    
    print 'printing'    
    # Print particles' final position
    for p in range(numberOfparticles): 
      ax.plot(particles[p].x2*centered_grid, particles[p].y2*centered_grid, 'y.')     
    
ax.axis('image')
ax.set_xticks([])
ax.set_yticks([])
plt.show()  















