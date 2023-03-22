//==============================================================================
// Start of the file "hybrid_sim.cc"
//
// Version 3b.
//
// Written by Vladimir Florinski
//==============================================================================


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>
#include "geo_memory.hh"
#include "geo_coord.hh"
#include "hybrid_distr.hh"
#include "hybrid_fields.hh"
#include "hybrid_moments.hh"
#include "hybrid_particles.hh"

#define is_boss !(rank % cores_per_node)

// enable generation of a BOV file for visualization of PDF data
// written by Annaleah Ernst, summer 2015
#define PDF_VIS

using namespace std;


// Universal physical constants
const double splight = 2.99792E+10;  // speed of light (cm/s)
const double echarge = 4.80320E-10;  // elementary charge (esu)
const double p_mass  = 1.67262E-24;  // mass of a proton (g)
const double kboltz  = 1.38065E-16;  // Boltzmann constant (erg/K)
const double one_ev  = 1.60218E-12;  // 1 eV (erg)

// size of the species parameter arrays
const int max_species = 10;

// frequency of run time info dumps: use 40 for Fourier analysis
const int tick = 1;

// frequency of data dumps (2 per 1000 orbit run with dt=0.02)
const int save = 110000;

// frequency of intermediate field and moment printing
const int ptick = 314; // 314 for once per gyroperiod

// multipliers for the PDF cube
const int ringResMultiplier = 35;
const int coreResMultiplier = 4;
const int pdfCubeCells = 500;

// File names
const string params_dat = "params.dat";
const string fields_fname = "fields.out";
const string fields_pname = "fields.dat";
const string moments_pname = "moments.dat";
const string spacetime_pname = "spacetime.dat";

// Variable base filename for intermediate field & moment printing
string moments_iname, fields_iname;


//------------------------------------------------------------------------------
// Save particle and field data
//------------------------------------------------------------------------------
// Input  | nsp      | number of species
// Input  | fields   | EM field
// Input  | particles| particle species array
// Input  | simtime  | simulation time
//------------------------------------------------------------------------------

void SaveAll(int nsp, fields_t &fields, particles_t *particles, double simtime)
{
   int rank, cores, spc, npart_saved, batch, skip;
   string prt_fname, prt_pname;
   ofstream prt_file;

   MPI_Comm_size(MPI_COMM_WORLD, &cores);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// This structure will be used for I/O (sequentially). We only activate it on
// the master on parallel runs.
   particles_t particles_disk[nsp];
   if(cores > 1 && is_master) {
      for(spc = 0; spc < nsp; spc++) {
         particles_disk[spc].Activate(particles[spc].GetNpart(), fields.GetLength());
      };
   };

// Save the fields data
   if(is_master) {
      fields.Save(fields_fname, simtime);
      fields.Print(1.0, fields_pname);
   };

// Loop over species - store particles data
   for(spc = 0; spc < nsp; spc++) {

// Open the file on the master
      if(is_master) {
         npart_saved = particles[spc].GetNpart() * cores;
         prt_fname = (spc < 10 ? "prt0" : "prt") + to_string(spc) + ".out";
         prt_file.open(prt_fname, ofstream::out | ofstream::binary);
         prt_file.write((char *)&npart_saved, sizeof(int));

// Write the master's batch and point map
         skip = (particles[spc].GetNpart() / cores) / 100000;
         if(skip < 1) skip = 1;
         particles[spc].Save(prt_file);
         prt_pname = (spc < 10 ? "prt0" : "prt") + to_string(spc) + ".dat";
         particles[spc].Print(3, skip, 1.0, ez, prt_pname);
      };

// Write the particles, in batches (one per core)
      for(batch = 1; batch < cores; batch++) {
         if(is_master) {
            particles_disk[spc].Receive(batch);
            particles_disk[spc].Save(prt_file);
         }
         else {
            if(batch == rank) particles[spc].Send(MASTER);
         };
      };
      if(is_master) prt_file.close();
   };

   if(is_master) cerr << "# Particle data wrote successfully\n";
};


//------------------------------------------------------------------------------
// Load particle and field data from files
//------------------------------------------------------------------------------
// Input  | nsp      | number of species
// InOut  | fields   | EM field
// InOut  | particles| particle species array
// Output | simtime  | simulation time
// Input  | boss_comm| boss communicator
// Input  | node_comm| node communicator
// Return |          | status of read operation
//------------------------------------------------------------------------------

int RestoreAll(int nsp, fields_t &fields, particles_t *particles, double &simtime,
   int boss_comm, int node_comm)
{
   int rank, node_rank, cores, spc, npart_saved, batch, data_bad;
   string prt_fname;
   ifstream prt_file;

   MPI_Comm_size(MPI_COMM_WORLD, &cores);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_rank(node_comm, &node_rank);

// This structure will be used for I/O (sequentially). We only activate it on
// the master on parallel runs.
   particles_t particles_disk[nsp];
   if(cores > 1 && is_master) {
      for(spc = 0; spc < nsp; spc++) {
         particles_disk[spc].Activate(particles[spc].GetNpart(), fields.GetLength());
      };
   };

// Open the fields file and read into the structure on the master
   if(is_master) data_bad = fields.Restore(fields_fname, simtime);

// make sure the number of grid points is the same
   MPI_Bcast(&data_bad, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
   if(data_bad) return 1;

// Broadcast the fields to the bosses and then workers
   if(cores > 1) {
      if(!node_rank) fields.Broadcast(boss_comm, MASTER);
      fields.Broadcast(node_comm, 0);
   };

// Loop over species - read particles data
   for(spc = 0; spc < nsp; spc++) {

// Open the file on the master
      if(is_master) {
         prt_fname = (spc < 10 ? "prt0" : "prt") + to_string(spc) + ".out";
         prt_file.open(prt_fname, ifstream::in | ifstream::binary);
         prt_file.read((char *)&npart_saved, sizeof(int));
         data_bad = (npart_saved != particles[spc].GetNpart() * cores);
      };

// make sure the number of particles is compatible
      MPI_Bcast(&data_bad, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
      if(data_bad) {
         if(is_master) prt_file.close();
         return 1;
      };

// Read in the master's data
      if(is_master) particles[spc].Restore(prt_file);

// Read the particles, in batches (one per core)
      for(batch = 1; batch < cores; batch++) {
         if(is_master) {
            particles_disk[spc].Restore(prt_file);
            particles_disk[spc].Send(batch);
         }
         else {
            if(batch == rank) particles[spc].Receive(MASTER);
         };
      };
      if(is_master) prt_file.close();
   };

   if(is_master) cerr << "# Particle data read successfully\n";
   return 0;
};


//------------------------------------------------------------------------------
// A 1D hybrid plasma simulation - driver program
//------------------------------------------------------------------------------
// Input  | argv[1]  | "-new" (start a new simulation) or "-old" (continue)
// Input  | argv[2]  | simulation duration, in gyro-periods
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

// MPI related variables
   int Ncores, rank;
   int cores_per_node, leftover_cores, cores_used, node;
   int Nbosses, *bosses, *workers, boss, worker;
   MPI_Comm boss_comm, node_comm;
   
// real simulation time
   time_t time_start, time_end;

// code simulation time
   int t = 0, nst = 0;
   double simtime, simtimestart, simtimeend, simtimemax, startft, stopft;

// Size of the simulation
   int Nspecies, Ncells, spc;
   int particles_per_cell[max_species];
   int particles_per_core[max_species];

// Plasma parameters
   int type[max_species];                // distribution type
   double d0[max_species];               // number density (= weight)
   double n0[max_species][4];            // symmetry axis
   double distro_params[max_species][4]; // distribution parameters

// Key parameters
   double xmax, dt;
   double B0[4], b0[4], mu0;
   double de0, pe0, Te0, wpi = 0.0, wci, wpiwci, beta0, betae;

// Computable values related to the magnetosonic wave and driving current
   double v_a, c_s, v_f, l_physical, t_freq;
// Defined parameters related to the magnetosonic wave and driving current
   double j_amplitude, wave_timelength;

// Diagnostic
   double en_ion[max_species], pa_ion[max_species], pa2_ion[max_species];
   double en_field, en_elec,  en_tot, dB2B2, grid_mass, grid_momt[4], grid_enrg;
   ofstream prt_file, st_file;
   ifstream parmfile;

   MPI_Group global_group, boss_group, node_group;

// Get the run time from the command line (specified in gyro-periods)
   if(argc < 3) exit(EXIT_FAILURE);
   simtimemax = atof(argv[2]);
   if(argc == 5) {
      startft = atof(argv[3]);
      stopft = atof(argv[4]);
   }
   else {
      startft = simtimemax + 1.0;
      stopft = -1.0;
   };
   if(!isnormal(simtimemax) || simtimemax < small) simtimemax = 0.0;

// Initialize the MPI subsystem
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &Ncores);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if(is_master) cerr << "#-------------------------------------------------------------------------------\n";

// Find the number of cores per node. To distinguish logical cores we simply
// divide by 2 if "HYPERTHREADS" is set.
   cores_per_node = sysconf(_SC_NPROCESSORS_ONLN);
#ifdef HYPERTHREADS
   cores_per_node /= 2;
#endif

// We don't require the number of core to be divisible by the number of nodes
   if(Ncores == 1) {
      Nbosses = cores_per_node = 1;
      if(is_master) cerr << "# A serial run was requested\n";
   }
   else {
      Nbosses = (Ncores > 1 ? (Ncores - 1) / cores_per_node + 1 : 1);
      if(is_master) {
         cerr << "# A parallel run was requested\n";
         cerr << "# The code will use " << Nbosses << (Nbosses == 1 ? " node\n" : " nodes\n");
         cerr << "# Each node contains " << cores_per_node << " CPU cores\n";
	 cerr << "# Total cores available: " << Nbosses * cores_per_node << "\n";
	 cerr << "# Total cores used: " << Ncores << endl;
      };
   };
    
// the cores that don't divide in nicely if we're running with extra cpus
   leftover_cores = Ncores % cores_per_node;

// initialize the RNG
   srand48(time(NULL) + rank);
//   srand48(rank);

// Build a list of bosses
   bosses = new int[Nbosses];
   for(boss = 0; boss < Nbosses; boss++) bosses[boss] = boss * cores_per_node;

// Set up the boss communicator
   MPI_Comm_group(MPI_COMM_WORLD, &global_group);
   MPI_Group_incl(global_group, Nbosses, bosses, &boss_group);
   MPI_Comm_create(MPI_COMM_WORLD, boss_group, &boss_comm);
   MPI_Group_free(&boss_group);
   delete[] bosses;

// Build a list of workers in our node. This list is different in each node!
   node = rank / cores_per_node;
    
// If we're on the last node and have leftover cores, use fewer cores
   cores_used = (node == Nbosses - 1 && leftover_cores) ? leftover_cores : cores_per_node;
   
   workers = new int[cores_used];
   for(worker = 0; worker < cores_used; worker++) {
      workers[worker] = node * cores_per_node + worker;
   };

// Set up the intra-node communicator
   MPI_Group_incl(global_group, cores_used, workers, &node_group);
   MPI_Comm_create(MPI_COMM_WORLD, node_group, &node_comm);
   MPI_Group_free(&node_group);
   MPI_Group_free(&global_group);
   delete[] workers;

//------------------------------------------------------------------------------
// Read simulation parameters from a parameter file. Each line of the parameter
// file consist of a string - value pair.
   if(is_master) {

      string temps1;
      parmfile.open(params_dat.c_str(), ifstream::in);

// Grid properties
      parmfile >> temps1 >> Nspecies;

      parmfile >> temps1 >> Ncells;
      parmfile >> temps1 >> xmax;
      parmfile >> temps1 >> dt;

      parmfile >> temps1 >> B0[0];
      parmfile >> temps1 >> mu0;
      parmfile >> temps1 >> Te0;

// Driving current properties
      parmfile >> temps1 >> j_amplitude;
      parmfile >> temps1 >> wave_timelength;
      
// Individual distribution properties
      de0 = 0.0;
      for(spc = 0; spc < Nspecies; spc++) {
         parmfile >> temps1 >> type[spc];
         parmfile >> temps1 >> particles_per_cell[spc];
         parmfile >> temps1 >> d0[spc];

// Add up the total number of ions to obtain the electron density
         de0 += d0[spc];

         parmfile >> temps1 >> distro_params[spc][0]; // v_perp,0 / v_0
         parmfile >> temps1 >> distro_params[spc][1]; // v_perp,t / v_t
         parmfile >> temps1 >> distro_params[spc][2]; // v_para,0 / mu_0
         parmfile >> temps1 >> distro_params[spc][3]; // v_para,t / mu_t
      };
      parmfile.close();

// Compute plasma frequencies and beta that will serve as the normalization
// constants. From this point on we switch to the code units:
// Velocity:        c=1
// Number density:  n_e=1
// Frequency:       wpi=1
// Magnetic field:  B=wci/wpi
      wpi = 2.0 * sqrtpi * echarge * sqrt(de0 / p_mass);
      wci = echarge * B0[0] / (p_mass * splight);
      beta0 = fourpi * d0[0] * p_mass * Sqr(distro_params[0][1]) / Sqr(B0[0]);
      betae = 8.0 * M_PI * de0 * kboltz * Te0 / Sqr(B0[0]);
      wpiwci = wpi / wci;
      dt *= wpiwci;
      pe0 = betae / (2.0 * wpiwci * wpiwci);

      cerr << "# " << Nspecies << " ion species\n";
      cerr << "# Plasma to cyclotron frequency ratio is " << wpiwci << endl;
      cerr << "# Core ion beta is " << beta0 << endl;
      cerr << "# Inertial length is " << splight / wpi << " cm\n";

// Magnetic field is in the XY plane
      B0[1] = mu0 / wpiwci;
      B0[2] = sqrt(1.0 - Sqr(mu0)) / wpiwci;
      B0[3] = 0.0;

// Computation required to add magnetosonic wave with driving current
      v_a = B0[0] / sqrt(fourpi * de0 * p_mass); // Alfven speed (cm/s)
      c_s = sqrt((gammaa * kboltz * Te0) / p_mass); // Ion sound speed (cm/s)
      v_f = sqrt((Sqr(v_a) + Sqr(c_s) + sqrt(Sqr(Sqr(v_a) + Sqr(c_s)) - 4.0*Sqr(v_a)*Sqr(c_s)*mu0)) / 2.0); // Fast magnetosonic wave speed (cm/s)
      l_physical = (splight / wpi) * xmax; // Full simulation length (cm)
      t_freq = (twopi / (l_physical / v_f)) / wpi; // Dimensionless time-frequency of the driving current (to multiply by simtime)

// Unit vector along B0
      Copy(B0, b0);
      Normalize(b0);

// Ion distribution properties
      for(spc = 0; spc < Nspecies; spc++) {
         d0[spc] /= de0 * particles_per_cell[spc];
         Copy(b0, n0[spc]);
         cerr << "# Species " << spc << " has " << particles_per_cell[spc] * Ncells << " particles\n";

// Normalize velocity parameters. For spherical distribution the second
// coordinate is mu, which requires no normalization.
         distro_params[spc][0] /= splight;
         distro_params[spc][1] /= splight;
         if(is_cylindrical(type[spc])) {
            distro_params[spc][2] /= splight;
            distro_params[spc][3] /= splight;
         };
      };
   };

// Distribute parameters among all CPUs
   MPI_Bcast(&Nspecies, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(&Ncells, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(&xmax, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(&dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(B0, 4, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(&wpiwci, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(&pe0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   
   MPI_Bcast(type, Nspecies, MPI_INT, MASTER, MPI_COMM_WORLD);
   MPI_Bcast(d0, Nspecies, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   for(spc = 0; spc < Nspecies; spc++) {
      MPI_Bcast(n0[spc], 4, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
      MPI_Bcast(distro_params[spc], 4, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
   };
   MPI_Bcast(particles_per_cell, Nspecies, MPI_INT, MASTER, MPI_COMM_WORLD);

// Make sure the number of particles is divisible by the number of cores
   for(spc = 0; spc < Nspecies; spc++) {
      particles_per_core[spc] = Ncells * particles_per_cell[spc] / Ncores;
      if(Ncells * particles_per_cell[spc] % Ncores) {
         if(is_master) cerr << "# Number of particles is not divisible by cores!\n";
         MPI_Finalize();
         exit(EXIT_FAILURE);
      };
   };

// Allocate particle and field variables
   particles_t particles[Nspecies];
   distribution_t distros[Nspecies];
   for(spc = 0; spc < Nspecies; spc++) {
      particles[spc].Activate(particles_per_core[spc], xmax);
      distros[spc].Activate(type[spc], distro_params[spc]);
   };
   fields_t fields(Ncells, xmax, wpiwci);
   moments_t moments(Ncells, xmax);
   moments_t moments_fs(Ncells, xmax);
   fields.Init(1, B0, pe0);

//------------------------------------------------------------------------------
// Continue from a previous run
   if(!strcmp(argv[1], "-old")) {
      if(RestoreAll(Nspecies, fields, particles, simtimestart, boss_comm, node_comm)) {
         MPI_Finalize();
         exit(EXIT_FAILURE);
      };
   }

// New run - apply initial conditions
   else {
      for(spc = 0; spc < Nspecies; spc++) {
         particles[spc].Load(distros[spc], d0[spc], n0[spc], Ncells);
      };
      simtimestart = 0.0;

// push velocity back for a new run only
      for(spc = 0; spc < Nspecies; spc++) particles[spc].PushV(fields, -dt / 2.0);
      if(is_master) cerr << "# Particles and fields initialized\n";
   };

// Set up the time counters
   simtimemax *= twopi * wpiwci;
   startft *= twopi * wpiwci;
   stopft *= twopi * wpiwci;
   simtimeend = simtimestart + simtimemax - small;
   simtime = simtimestart;
   MPI_Barrier(MPI_COMM_WORLD);
   time_start = time(NULL);

//******************************************************************************

// Simulation loop
   if(is_master) cerr << "#-------------------------------------------------------------------------------\n";
   while(simtime < simtimeend - small) {

// Reset the moments
      moments.Reset();
      moments_fs.Reset();

// Push the particles using fields at t
      for(spc = 0; spc < Nspecies; spc++) {

// Push the particle velocity from n-1/2 to n+1/2
         particles[spc].PushMatrix(fields, dt);
//         particles[spc].PushBoris(fields, dt);

// Push the position to n+1/2. This yields time-centered moments to compute B.
         particles[spc].PushX(dt / 2.0);
         particles[spc].ComputeMoments(moments, true);

// Push the position to n. This yields the free-streaming moments (velocity is
// unchanged) to compute E.
         particles[spc].PushX(dt / 2.0);
         particles[spc].ComputeMoments(moments_fs, true);
      };         

// Collect the moments on the master
      if(Ncores != 1) {
         moments.Collect(node_comm, 0);
         if(is_boss) moments.Collect(boss_comm, MASTER);
         moments_fs.Collect(node_comm, 0);
         if(is_boss) moments_fs.Collect(boss_comm, MASTER);
      };

// Advance the fields using filtered moments
      if(is_master) {
         moments.Filter(0.5);
         fields.AdvanceImplicit(moments, dt);
         if((simtime / wpiwci) / twopi < wave_timelength) { // Apply the driving current if the simulation has not yet reached the switch-off time.
            fields.ApplyDrivingCurrent(simtime, t_freq, j_amplitude);
         };
         moments_fs.Filter(0.5);
         fields.AdvanceElectric(moments_fs, dt / 2.0);
      };

// Broadcast the fields to the bosses and then workers
      if(Ncores != 1) {
         if(is_boss) fields.Broadcast(boss_comm, MASTER);
         fields.Broadcast(node_comm, 0);
      };

//------------------------------------------------------------------------------

// Open the file for FT analysis at the requested time
      if(is_master && !st_file.is_open() && simtime > simtimestart + startft
         && simtime < simtimestart + stopft) {
         st_file.open(spacetime_pname.c_str(), ofstream::out | ofstream::binary);
         cerr << "# Started recording spacetime data\n";
      };

// Print run time info
      if(!(t % tick)) {
// Get one particle from the sample
//         particles[0].GetOne(320, tpart_r, tpart_v, tpart_w);

// Compute and print field and grid energies
         en_field = fields.Energy();
         en_elec = fields.ElectronEnergy();
         en_tot = en_field + en_elec;
         dB2B2 = fields.MagneticVariance();
         if(is_master) {
            moments.GetTotal(grid_mass, grid_momt, grid_enrg);
            cout << setw(14) << setprecision(6) << simtime / wpiwci;
            cout << setw(16) << setprecision(9) << dB2B2;
//                 << setw(12) << setprecision(5) << grid_mass
//                 << setw(12) << setprecision(5) << grid_momt[1]
//                 << setw(12) << setprecision(5) << grid_momt[2]
//                 << setw(12) << setprecision(5) << grid_momt[3]
//                 << setw(12) << setprecision(5) << grid_enrg

//                 << setw(12) << setprecision(5) << tpart_r[2]
//                 << setw(12) << setprecision(5) << tpart_r[3]
//                 << setw(12) << setprecision(5) << tpart_v[1]
//                 << setw(12) << setprecision(5) << tpart_v[2]
//                 << setw(12) << setprecision(5) << tpart_v[3]
//                 << setw(12) << setprecision(5) << tpart_r[1]
//                 << setw(16) << setprecision(9) << en_field
//                 << setw(16) << setprecision(9) << en_elec;
         };

// Compute, collect, and print particle energies and pitch angle spread
         for(spc = 0; spc < Nspecies; spc++) {
            en_ion[spc] = particles[spc].Energy();
            particles[spc].MuAverage(n0[spc], pa_ion[spc], pa2_ion[spc]);
         };
         if(Ncores != 1) {
            MPI_Reduce((is_master ? MPI_IN_PLACE : en_ion), en_ion, Nspecies,
                       MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce((is_master ? MPI_IN_PLACE : pa_ion), pa_ion, Nspecies,
                       MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce((is_master ? MPI_IN_PLACE : pa2_ion), pa2_ion, Nspecies,
                       MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
         };
         for(spc = 0; spc < Nspecies && is_master; spc++) {
            cout << setw(16) << setprecision(9) << en_ion[spc];
            en_tot += en_ion[spc];
            cout << setw(16) << setprecision(9) << pa_ion[spc]
                    / (particles_per_core[spc] * Ncores * d0[spc]);
            cout << setw(16) << setprecision(9)
                 << pa2_ion[spc] / (particles_per_core[spc] * Ncores * d0[spc])
               - Sqr(pa_ion[spc] / (particles_per_core[spc] * Ncores * d0[spc]));
         };
         if(is_master) cout << setw(16) << setprecision(9) << en_tot << endl;

// Dump magnetic field data for later Fourier processing
         if(is_master && st_file.is_open()) {
            fields.Dump(st_file, 2);
            nst++;
         };
      };
      
// Save the data throughout the run
      if(!(t % save) && t) SaveAll(Nspecies, fields, particles, simtime);

// Print intermediate field & moment data
      if(!(t % ptick) && t && is_master) {
         moments_iname = "moments_intermediate_" + to_string(t) + ".dat";
         fields_iname = "fields_intermediate_" + to_string(t) +  ".dat"; 
         moments.Print(1.0, moments_iname);
         fields.Print(1.0, fields_iname);
      };

      t++;
      simtime += dt;

// At the appropriate time close the FT file
      if(is_master && st_file.is_open() && simtime > simtimestart + stopft) {
         st_file.close();
         cerr << "# Stopped recording spacetime data\n";
      };
   };
   if(is_master) cerr << "#-------------------------------------------------------------------------------\n";

//******************************************************************************

// Print the summary info
   MPI_Barrier(MPI_COMM_WORLD);
   time_end = time(NULL);
   if(is_master) {
      cerr << "# Code ran for " << setw(5) << t << " time steps\n";
      cerr << "# Wall execution time was " << setw(5) << (int)(time_end - time_start)
           <<  " s\n";
      cerr << "# Physical simulation time was " << setw(10) << setprecision(3)
           << (simtime - simtimestart) / wpi << " s\n";
      cerr << "# Total simulation time was " << setw(10) << setprecision(3)
           << simtime / wpi << " s\n";
   };

   if(!(!strcmp(argv[1], "-old") && !t)) SaveAll(Nspecies, fields, particles, simtime);

   if(is_master) {
      moments.Print(1.0, moments_pname);
//      fields.Print(wpiwci * magnetic_field, fields_pname);
      fields.Print(1.0, fields_pname);
      
      if(st_file.is_open()) st_file.close();
      cerr << "# Dumped " << nst << " time records to " << spacetime_pname << endl;

#ifdef PDF_VIS // written by Annaleah Ernst, summer 2015

// generate PDF field doc
      for(spc = 0; spc < Nspecies; spc++){

// special handling for size when spc is greater than 0 (ie, not the core)
         particles[spc].GeneratePDF((spc ? ringResMultiplier : coreResMultiplier)
            * beta0 / wpiwci, pdfCubeCells, spc);
      };

// generate python script for visualization of PDF data in this directory
      ofstream py_file;
      py_file.open("RunPDFVis.py", ofstream::out);
      if(py_file) {
         py_file << "# run with python. Assumes install of PyPDFVis.py" << endl;
         py_file << "from os import system, getcwd" << endl;
         py_file << "system('PyPDFVis.py ' + getcwd())" << endl;
         py_file.close();
      };
#endif

   };

// Free up the extra communicators
   if(is_boss) MPI_Comm_free(&boss_comm);
   MPI_Comm_free(&node_comm);

   MPI_Finalize();
   if(is_master) {
      cerr << "# Run completed\n";
      cerr << "#-------------------------------------------------------------------------------\n";
   };
   exit(EXIT_SUCCESS);
};
