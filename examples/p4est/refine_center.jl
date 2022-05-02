 
# /** \file refine_center.jl
#  *
#  * This 2D example program refines a domain according to a radius.
#  */
# 
# /* p4est has two separate interfaces for 2D and 3D, p4est*.h and p8est*.h.
#  * Most API functions are available for both dimensions.  The header file
#  * p4est_to_p8est.h #define's the 2D names to the 3D names such that most code
#  * only needs to be written once.  In this example, we rely on this. */

using MPI
using T8code
using CBinding

const min_level = 1
const max_level = 7

const center = [0.52,0.6]
const radius = 0.123

const box = [ 0 0 ; 1 0 ; 0 1 ; 1 1 ]

# /** Callback function to decide on refinement.
#  *
#  * Refinement and coarsening is controlled by callback functions.
#  * This function is called for every processor-local quadrant in order; its
#  * return value is understood as a boolean refinement flag.
#  * In this example we use the image file hw32.h to determine the refinement.
#  */
function refine_fn(p4est, which_tree, quadrant)
  # /* The connectivity chosen in main () only consists of one tree. */
  @P4EST_ASSERT(which_tree == 0);

  if quadrant.level[] < min_level
    return Int32(1);
  end

  # /* We do not want to refine deeper than a given maximum level. */
  # if quadrant.level[] > plv
  if quadrant.level[] >= max_level
    return Int32(0);
  end

  coords3d = fill(-42.0,3)
  c"p4est_qcoord_to_vertex"(conn,0,quadrant.x[],quadrant.y[],coords3d)

  w = 1.0 / 2^quadrant.level[]
  coords = coords3d[1:2]

  verts = repeat((coords-center)',4) + w.*box
  radii = sqrt.(sum(verts.^2,dims=2))

  if any(radii .<= radius)
    return Int32(1)
  end

  return Int32(0);
end

# /* Initialize MPI; see sc_mpi.h.
# * If configure --enable-mpi is given these are true MPI calls.
# * Else these are dummy functions that simulate a single-processor run. */
MPI.Init()
mpicomm = MPI.COMM_WORLD.val

NULL = C_NULL

# /* These functions are optional.  If called they store the MPI rank as a
# * static variable so subsequent global p4est log messages are only issued
# * from processor zero.  Here we turn off most of the logging; see sc.h. */
c"sc_init"(mpicomm, 1, 1, NULL, c"SC_LP_ESSENTIAL");
c"p4est_init"(NULL, c"SC_LP_PRODUCTION");
c"P4EST_GLOBAL_PRODUCTIONF"("This is the p4est %dD demo example/steps/%s_step1\n",c"P4EST_DIM", c"P4EST_STRING");

# /* Create a forest that consists of just one quadtree/octree.
# * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
# * checked to execute dimension-dependent code. */
# #ifndef P4_TO_P8
conn = c"p4est_connectivity_new_unitsquare"();
# #else
# conn = c"p8est_connectivity_new_unitcube ();
# #endif

# /* Create a forest that is not refined; it consists of the root octant. */
p4est = c"p4est_new"(mpicomm, conn, 0, NULL, NULL);

# /* Refine the forest recursively in parallel.
# * Since refinement does not change the partition boundary, this call
# * must not create an overly large number of quadrants.  A numerical
# * application would call p4est_refine non-recursively in a loop,
# * repartitioning in each iteration.
# * The P4EST_ASSERT macro only activates with --enable-debug.
# * We check against the data dimensions in example/steps/hw32.h. */
recursive = 1;
c"p4est_refine"(p4est, recursive, refine_fn, NULL);

# /* Partition: The quadrants are redistributed for equal element count.  The
# * partition can optionally be modified such that a family of octants, which
# * are possibly ready for coarsening, are never split between processors. */
partforcoarsen = 0;
c"p4est_partition"(p4est, partforcoarsen, NULL);

# /* If we call the 2:1 balance we ensure that neighbors do not differ in size
# * by more than a factor of 2.  This can optionally include diagonal
# * neighbors across edges or corners as well; see p4est.h. */
balance = true;
if balance
  c"p4est_balance"(p4est, c"P4EST_CONNECT_FACE", NULL);
  c"p4est_partition"(p4est, partforcoarsen, NULL);
end

# /* Write the forest to disk for visualization, one file per processor. */
c"p4est_vtk_write_file"(p4est, NULL, c"P4EST_STRING"*"_step1");

# /* Destroy the p4est and the connectivity structure. */
c"p4est_destroy"(p4est);
c"p4est_connectivity_destroy"(conn);

# /* Verify that allocations internal to p4est and sc do not leak memory.
# * This should be called if sc_init () has been called earlier. */
c"sc_finalize"();

# /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
c"sc_MPI_Finalize"();
