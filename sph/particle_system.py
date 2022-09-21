class ParticleSystem:
    def __init__(self, particles):
        self.particles = particles

    def tick(self, dt):
        forces = {}
        ind = 0
        for particle in self.particles:
            ind += 1
            force = particle.force(self.particles)
            density = particle.density(self.particles)
            forces[particle.id] = (density, force)

        for particle in self.particles:
            density, force = forces[particle.id]

            particle.tick(force, density, dt)
            print(particle.velocity, density)
