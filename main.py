from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sph.particle_system import ParticleSystem

from sph.particles import Particle

particles = []

for x in range(5):
    for y in range(5):
        for z in range(5):
            p = Particle(
                5**x * 3**y * 2**z,
                np.array([0.1 + 0.2*x, 0.1 + 0.2*y, 0.1 + 0.2*z]),
                np.array([0.0, 0.0, 0.0]),
                np.array([0.0, 0.0, 0.0]),
                1,
                1,
                np.array([0.0, 0.0, 0.0]),
            )
            particles.append(p)

fig = plt.figure()

ax = fig.add_subplot(projection='3d')

ax.axes.set_xlim3d(left=-0.1, right=1.1)
ax.axes.set_ylim3d(bottom=-0.1, top=1.1)
ax.axes.set_zlim3d(bottom=-0.1, top=1.1)

particle_system = ParticleSystem(particles)

for particle in particle_system.particles:
    ax.scatter(particle.position[0],
               particle.position[1], particle.position[2], color='b')


def animate(i):
    particle_system.tick(0.016)
    ax.clear()
    for particle in particle_system.particles:
        ax.scatter(particle.position[0],
                   particle.position[1], particle.position[2], color='b')
    ax.axes.set_xlim3d(left=-0.1, right=1.1)
    ax.axes.set_ylim3d(bottom=-0.1, top=1.1)
    ax.axes.set_zlim3d(bottom=-0.1, top=1.1)


ani = FuncAnimation(fig, animate, frames=200, interval=500, repeat=False)

plt.show()
