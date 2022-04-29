
using CBinding
using MPI
using T8code


MPI.Init()

comm = MPI.COMM_WORLD.val

c"p4est_init"(C_NULL,c"SC_LP_VERBOSE");

conn = c"p4est_connectivity_new_periodic"()

initial_refinement_level = 2

p4est = c"p4est_new_ext"(
            comm,
            conn,
            0, # No minimum initial qudrants per processor
            initial_refinement_level,
            true, # Refine uniformly
            2 * sizeof(Int), # Use Int-Vector of size 2 as quadrant user data
            C_NULL, # No init function
            C_NULL) # No user pointer

c"p4est_destroy"(p4est)
c"p4est_connectivity_destroy"(conn)
