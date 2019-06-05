# Gas-solid
# Vitor Vilela
# 29 September 2016

import numpy as np
import math
import matplotlib.pyplot as plt
from skimage import io
import random

# Simulate for Stokes number = 0.01, 0.1 and 1.
Stokes = 1.
simulation_time = 1.

L = 1.0
U = 1.0

timestep = 0.001
numberOfparticles = 100  

particle_density = 800.0
air_viscosity = 2.e-5
air_density = 1.

centered_grid = 30
delta = L/centered_grid

staggered_grid = centered_grid+2
node_grid = centered_grid+1

vx = np.zeros((staggered_grid,staggered_grid))
vy = np.zeros((staggered_grid,staggered_grid))
velocity = np.zeros((centered_grid,centered_grid))
vort = np.zeros((node_grid,node_grid))
wz = np.zeros((node_grid,node_grid))



# Particle data type
# Time 2 > Time 1
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


# Discretized - Eulerian velocity and vorticity fields
# Axis y-rows, x-columns - origin at left-top corner

# Vx analytical function
def VX(x, y):
  return U*math.sin(math.pi*x/L)*math.cos(math.pi*y/L)

# Vy analytical function
def VY(x, y):
  return -U*math.cos(math.pi*x/L)*math.sin(math.pi*y/L)

# Wz analytical function
def WZ(x, y):
  return (2*U*math.pi/L)*math.sin(math.pi*x/L)*math.sin(math.pi*y/L)

# Particle mass
def Mp(radius, density):
  return density*(4/3)*math.pi*pow(radius,3.)

# Particle diameter 
def Dp(air_viscosity, U, L, particle_density, Stokes):
  return math.sqrt(18*air_viscosity*L*Stokes/(particle_density*U))

# Component of Relative Velocity vector (Vrx or Vry)
def VRi(v, vp):
  return v-vp

# Module of Relative Velocity vector
def VR(vr1, vr2):
  return math.sqrt(vr1*vr1+vr2*vr2)

# Particle's Reynolds number
def Re(air_density, air_viscosity, VR, Dp):
  return air_density*VR*Dp/air_viscosity

# Drag Coefficient
def Cd(Re):
  return 64./Re

# Drag force (Dx or Dy)
def drag(air_density, particle_density, Dp, Cd, VR, VRi):
  return 3*air_density*Cd*VR*VRi/(4*particle_density*Dp)
  
# Saffman force (Sx or Sy)
def saffman(Dp, Mp, air_density, air_viscosity, WZ, VRi):
  return 1.61*Dp*Dp*math.sqrt(air_density*air_viscosity*abs(WZ))*VRi*WZ

# Euler Temporal Integration - Velocity
def Euler(vp0, timestep, rhs):
  return vp0 + timestep*rhs

# Update particle position
def Move(p0, timestep, v2):
  return p0 + timestep*v2




# Staggered Vx - @face 
for j in range(staggered_grid):
  y = delta*(j-0.5)
  for i in range(staggered_grid):
    x = delta*(i-1)
    vx[j,i] = VX(x, y)
    
# Staggered Vy - @face
for j in range(staggered_grid):
  y = delta*(j-1)
  for i in range(staggered_grid):
    x = delta*(i-0.5)
    vy[j,i] = VY(x, y)   
    
# Centered Velocity Magnitude - @center
for j in range(centered_grid):
  y = delta*(j+0.5)
  for i in range(centered_grid):
    x = delta*(i+0.5)
    velocity[j,i] = math.sqrt(pow(0.5*(vx[j,i+1]+vx[j,i+2]),2.)+pow(0.5*(vy[j+1,i]+vy[j+2,i]),2.))
    
# Numerical and Analytical Vorticity - @node 
for j in range(node_grid):
  y = delta*j
  for i in range(node_grid):
    x = delta*i
    vort[j,i] = (vy[j,i+1]-vy[j,i])/delta - (vx[j+1,i]-vx[j,i])/delta  
    wz[j,i] = WZ(x, y)    
   
# Show Eulerian fields: vx, vy, velocity, vort, wz        
#io.imshow(velocity)
#io.show()
         
    
# Print particle's position
fig, ax = plt.subplots()
ax.imshow(wz, interpolation='nearest', cmap=plt.cm.gray)  
        
# Particle initialization     
particles = []
random.seed(numberOfparticles)
for p in range(numberOfparticles): 
  x = random.random()
  y = random.random()
  vx = VX(x, y)
  vy = VY(x, y)
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

    # Update particle's position (x and y)
    particles[p].x2 = Move(particles[p].x1, timestep, particles[p].vx2)
    particles[p].y2 = Move(particles[p].y1, timestep, particles[p].vy2)  

    # Find new particle's velocity integrating ODE (Newton's law) in time
    particle_diameter = Dp(air_viscosity, U, L, particle_density, Stokes)
    particle_mass = Mp(0.5*particle_diameter, particle_density)
    VRx = VRi( VX(particles[p].x2, particles[p].y2), particles[p].vx1 )
    VRy = VRi( VY(particles[p].x2, particles[p].y2), particles[p].vy1 )
    relative_velocity = VR(VRx, VRy)
    vorticity = WZ(particles[p].x1, particles[p].y1)  
    Reynolds = Re(air_density, air_viscosity, relative_velocity, particle_diameter)
    drag_coefficient = Cd(Reynolds)
        
    RHSx = drag(air_density, particle_density, particle_diameter, drag_coefficient, relative_velocity, VRx) + saffman(particle_diameter, particle_mass, air_density, air_viscosity, vorticity, VRx) 
    RHSy = drag(air_density, particle_density, particle_diameter, drag_coefficient, relative_velocity, VRy) + saffman(particle_diameter, particle_mass, air_density, air_viscosity, vorticity, VRy) 
  
    velx = particles[p].vx2
    vely = particles[p].vy2
    particles[p].vx2 = Euler(particles[p].vx1, timestep, RHSx)
    particles[p].vy2 = Euler(particles[p].vy1, timestep, RHSy)  
    particles[p].x2 =  Move(particles[p].x1, timestep, particles[p].vx2)
    particles[p].y2 =  Move(particles[p].y1, timestep, particles[p].vy2)
    
    # Particle BC - no Reflection
    if particles[p].x2 < 0.:
      particles[p].x2 = 0.
    if particles[p].y2 < 0.:
      particles[p].y2 = 0.
    if particles[p].x2 > L:
      particles[p].x2 = L
    if particles[p].y2 > L:
      particles[p].y2 = L        
       
    # Update
    particles[p].x1 = particles[p].x2
    particles[p].y1 = particles[p].y2
    particles[p].vx1 = particles[p].vx2
    particles[p].vy1 = particles[p].vy2    
  
# Print particles' final position
for p in range(numberOfparticles): 
  ax.plot(particles[p].x2*centered_grid, particles[p].y2*centered_grid, 'r.')     
ax.axis('image')
ax.set_xticks([])
ax.set_yticks([])
plt.show()  















