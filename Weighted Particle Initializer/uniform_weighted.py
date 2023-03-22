# DEPRECATED: Use faster c++ implementation instead.

import numpy as np
import struct

splight = 2.99792E+10  # speed of light (cm/s)
echarge = 4.80320E-10  # elementary charge (esu)
p_mass  = 1.67262E-24  # mass of a proton (g)
kboltz  = 1.38065E-16  # Boltzmann constant (erg/K)
one_ev  = 1.60218E-12   # 1 eV (erg)

class Particles: 
    def __init__(self, particle_count : int):
        self.nparts = particle_count
        self.x = np.zeros(particle_count)
        self.y = np.zeros(particle_count)
        self.z = np.zeros(particle_count)
        self.vx = np.zeros(particle_count)
        self.vy = np.zeros(particle_count)
        self.vz = np.zeros(particle_count)
        self.w = np.zeros(particle_count)
    
    def write(self, fname : str):
        with open(fname, 'wb') as output:
            output.write(struct.pack('i', self.nparts))
            for i in range(self.nparts):
                output.write(struct.pack('d', self.x[i]))
                output.write(struct.pack('d', self.y[i]))
                output.write(struct.pack('d', self.z[i]))
                output.write(struct.pack('d', self.vx[i]))
                output.write(struct.pack('d', self.vy[i]))
                output.write(struct.pack('d', self.vz[i]))
                output.write(struct.pack('d', self.w[i]))

class Fields:
    def __init__(self, cell_count : int):
        self.ncells = cell_count
        self.bx = np.zeros(cell_count)
        self.by = np.zeros(cell_count)
        self.bz = np.zeros(cell_count)
        self.ex = np.zeros(cell_count)
        self.ey = np.zeros(cell_count)
        self.ez = np.zeros(cell_count)
        self.simtime = 0

    def write(self, fname : str):
        with open(fname, 'wb') as output:
            output.write(struct.pack('i', self.ncells))
            output.write(struct.pack('d', self.simtime))
            output.write(struct.pack('d', 0.0))
            for i in range(self.ncells):
                output.write(struct.pack('d', self.bx[i]))
                output.write(struct.pack('d', self.by[i]))
                output.write(struct.pack('d', self.bz[i]))
            output.write(struct.pack('d', 0.0))
            for i in range(self.ncells):
                output.write(struct.pack('d', self.ex[i]))
                output.write(struct.pack('d', self.ey[i]))
                output.write(struct.pack('d', self.ez[i]))
           
class Parameters:
    def __init__(self):
        with open("params.dat", 'r') as paramsfile: 
            lines = paramsfile.readlines()
        
        self.ncells = int(lines[2].split()[1])
        self.box_length = float(lines[3].split()[1])
        self.B0 = float(lines[6].split()[1])
        self.cos_theta = float(lines[7].split()[1])
        self.e_temp = float(lines[8].split()[1])
        self.ppc = int(lines[14].split()[1])
        self.density = float(lines[15].split()[1])
        self.v_t = np.sqrt(2 * kboltz * self.e_temp / p_mass)
    
    def __repr__(self):
        return f"Using {self.ncells} simulation cells with {self.ppc} particles per cell. \nBox length is {self.box_length} inertial lengths. \nMagnetic field is {self.B0} G. \nMagnetic field angle cosine is {self.cos_theta}. \nElectron temperature is {self.e_temp} K. \nIon number density is {self.density} cm^-3. \nThermal velocity is {self.v_t} cm/s.\n"

# Evaluates the maxwellian distribution given the magnitude of the velocity vector
def evaluate_maxwellian(params : Parameters, abs_v : float) -> float:
    return (params.density / (np.pi**1.5 * params.v_t**3)) * np.exp(-(abs_v * abs_v)/(params.v_t * params.v_t))

def assign_weights(params : Parameters, mode : int, particles : Particles) -> None:
    # Mode 1: Maxwellian
    if (mode == 1): 
        for i in range(particles.nparts):
            particles.w[i] = evaluate_maxwellian(params, np.sqrt((particles.vx[i]*particles.vx[i]) + (particles.vy[i]*particles.vy[i]) + (particles.vz[i]*particles.vz[i])))
    else: 
        print("Improper mode passed to assign_weights; exiting.")
        exit(1)

# Reverses the bits of index in the given radix
# Used to set up a quiet start for x, vx, vy, vz.
def reverse(index : int, radix : int) -> float:
    rev = iquot = irem = 0
    inum = index
    power = 1

    while True: 
        iquot = np.floor(inum/radix)
        irem = inum - radix*iquot
        power /= radix
        rev += irem*power
        inum = iquot
        if not (inum > 0):
            break

    return rev

if __name__ == "__main__":
    params = Parameters()
    print(params)
    
    wpi = 2.0 * np.sqrt(np.pi) * echarge * np.sqrt(params.density / p_mass)
    wci = echarge * params.B0 / (p_mass * splight)
    wpiwci = wpi/wci 
    dx = params.box_length / params.ncells
    lambda_i = splight / wpi

    fields = Fields(params.ncells + 2)
    for i in range(fields.ncells):
        fields.bx[i] = params.cos_theta / wpiwci
        fields.by[i] = np.sqrt(1.0 - (params.cos_theta * params.cos_theta)) / wpiwci
        # All other values are zero; there's nothing to be done since we've initialized them as zeros.
    
    particles = Particles(params.ppc * params.ncells)
    
    for particle_index in range(params.ncells * params.ppc):
        # Quiet start for x within the range of an individual cell
        x_edge = np.floor(particle_index/params.ppc) * dx * lambda_i
        x_variation = reverse(particle_index, 2)
        particles.x[particle_index] = x_edge + x_variation*dx*lambda_i

        # Quiet start for velocities in a cube of -3v_t to +3v_t on each side
        vx_float = reverse(particle_index, 3)
        vx_val = (6*vx_float - 3)*params.v_t
        particles.vx[particle_index] = vx_val
        vy_float = reverse(particle_index, 5)
        vy_val = (6*vy_float - 3)*params.v_t
        particles.vy[particle_index] = vy_val
        vz_float = reverse(particle_index, 7)
        vz_val = (6*vz_float - 3)*params.v_t
        particles.vz[particle_index] = vz_val

        # y and z are intentionally left as zero; w is determined later.    

    assign_weights(params, 1, particles)

    # Convert values back into normalized code units
    # x needs to be in units of lambda_i
    particles.x /= lambda_i
    # y and z are zero anyway
    # velocities should be normalized to c
    particles.vx /= splight
    particles.vy /= splight
    particles.vz /= splight
    # Normalize so that total weight is ncells
    particles.w *= (params.ncells/sum(particles.w))

    particles.write("prt00.out")
    fields.write("fields.out")
