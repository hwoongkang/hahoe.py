import numpy as np
from sph.const import *

from sph.kernels import grad_spiky, laplacian_viscosity, poly_6


class Particle:
    def __init__(self, id, position, velocity, acceleration, mass, radius, color):
        self.id = id
        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.mass = mass
        self.radius = radius
        self.color = color

    def distance_squared(self, other):
        return np.linalg.norm(self.position - other.position) ** 2

    def distance(self, other):
        return np.linalg.norm(self.position - other.position)

    def density(self, others):
        density = 0
        return 1000
        for other in others:
            distance_squared = self.distance_squared(other)
            density += other.mass * poly_6(distance_squared)
        return density * PARTICLE_MASS

    def pressure(self, others):
        density = self.density(others)
        return PRESSURE_KAPPA * density

    def force_gravity(self):
        return np.array([0.0, 0.0, -GRAVITY])

    def force_pressure(self, others):
        force = np.array([0.0, 0.0, 0.0])

        p_i = self.pressure(others)

        for other in others:
            distance = self.distance(other)
            if distance >= DISTANCE_LIMIT:
                continue
            gradient = grad_spiky(self.position, other.position)
            rho_j = other.density(others)
            p_j = other.pressure(others)
            coeff = (p_i + p_j) / rho_j / 2.0
            force += coeff * gradient

        return force * PARTICLE_MASS

    def force_viscosity(self, others):
        force = np.array([0.0, 0.0, 0.0])
        for other in others:
            distance = self.distance(other)
            if distance >= DISTANCE_LIMIT:
                continue
            viscosity_term = laplacian_viscosity(distance)
            rho_j = other.density(others)
            dv = other.velocity - self.velocity
            force += viscosity_term * dv / rho_j
        return force * VISCOSITY_MU * PARTICLE_MASS

    def force_surface_tension(self, others):
        return np.array([0.0, 0.0, 0.0])

    def force(self, others):
        return (
            self.force_gravity()
            + self.force_pressure(others)
            + self.force_viscosity(others)
            + self.force_surface_tension(others)
        )

    def reflect(self):
        if self.position[2] < 0:
            self.position[2] = 0
            self.velocity[2] *= -1

        if self.position[2] > 1:
            self.position[2] = 1
            self.velocity[2] *= -1

        if self.position[0] < 0:
            self.position[0] = 0
            self.velocity[0] *= -1

        if self.position[0] > 1:
            self.position[0] = 1
            self.velocity[0] *= -1

        if self.position[1] < 0:
            self.position[1] = 0
            self.velocity[1] *= -1

        if self.position[1] > 1:
            self.position[1] = 1
            self.velocity[1] *= -1

    def tick(self, force, density, dt):
        acc = force / density
        vel = self.velocity + acc * dt / 2
        pos = self.position + vel * dt
        self.velocity = vel
        self.position = pos
        self.reflect()
