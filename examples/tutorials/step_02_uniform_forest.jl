# vim: tabstop=2:softtabstop=2

using MPI
using T8code
using CBinding


# /* Builds cmesh of 2 prisms that build up a unit cube. 
#  * See step1 for a detailed description.
#  * \param [in] comm   MPI Communicator to use.
#  * \return            The coarse mesh.
#  */
function t8_step2_build_prismcube_coarse_mesh(comm)
    # /* Build a coarse mesh of 2 prism trees that form a cube. */
    cmesh = c"t8_cmesh_new_hypercube"(c"T8_ECLASS_PRISM", comm, 0, 0, 0);
    c"t8_global_productionf"("[step 2] Constructed coarse mesh with 2 prism trees.\n");

    return cmesh;
end

# /* Build a uniform forest on a cmesh 
#  * using the default refinement scheme.
#  * \param [in] comm   MPI Communicator to use.
#  * \param [in] cmesh  The coarse mesh to use.
#  * \param [in] level  The initial uniform refinement level.
#  * \return            A uniform forest with the given refinement level that is
#  *                    partitioned across the processes in \a comm.
#  */
function t8_step2_build_uniform_forest(comm, cmesh, level)
  # /* Create the refinement scheme. */
  scheme = c"t8_scheme_new_default_cxx"();

  # /* Creat the uniform forest. */
  forest = c"t8_forest_new_uniform"(cmesh, scheme, level, 0, comm);

  return forest;
end

# /* Write vtk (or more accurately vtu) files of the forest.
#  * \param [in] forest   A forest.
#  * \param [in] prefix   A string that is used as a prefix of the output files.
#  * 
#  * This will create the file prefix.pvtu
#  * and additionally one file prefix_MPIRANK.vtu per MPI rank.
#  */
function t8_step2_write_forest_vtk(forest, prefix)
    c"t8_forest_write_vtk"(forest, prefix);
end

# /* Destroy a forest. This will free all allocated memory.
#  * \param [in] forest    A forest.
#  * NOTE: This will also free the memory of the scheme and the cmesh, since
#  *       the forest took ownership of them.
#  *       If we do not want this behaviour, but want to reuse for example the cmesh,
#  *       we need to call t8_cmesh_ref (cmesh) before passing it to t8_forest_new_uniform.
#  */
function t8_step2_destroy_forest(forest)
    c"t8_forest_unref"(Ref(forest));
end

# /* The prefix for our output files. */
const prefix = "t8_step2_uniform_forest";

# /* The uniform refinement level of the forest. */
const level = 3;

# /* Initialize MPI. This has to happen before we initialize sc or t8code. */
MPI.Init()
comm = MPI.COMM_WORLD.val

# /* Initialize the sc library, has to happen before we initialize t8code. */
c"sc_init"(comm, 1, 1, C_NULL, c"SC_LP_ESSENTIAL");
# /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
c"t8_init"(c"SC_LP_PRODUCTION");

# /* Print a message on the root process. */
c"t8_global_productionf"(" [step 2] \n");
c"t8_global_productionf"(" [step 2] Hello, this is the step2 example of t8code.\n");
c"t8_global_productionf"(" [step 2] In this example we build our first uniform forest and output it to vtu files.\n");
c"t8_global_productionf"(" [step2] \n");

# /* Create the cmesh from step1 */
cmesh = t8_step2_build_prismcube_coarse_mesh(comm);

# /* Build the uniform forest, it is automatically partitioned among the processes. */
forest = t8_step2_build_uniform_forest(comm, cmesh, level);
# /* Get the local number of elements. */
local_num_elements = c"t8_forest_get_local_num_elements"(forest);
# /* Get the global number of elements. */
global_num_elements = c"t8_forest_get_global_num_elements"(forest);

# /* Print information on the forest. */
c"t8_global_productionf"(" [step 2] Created uniform forest.\n");
c"t8_global_productionf"(" [step 2] Refinement level:\t\t\t%i\n",level);
c"t8_global_productionf"(" [step 2] Local number of elements:\t\t%i\n",local_num_elements);
c"t8_global_productionf"(" [step 2] Global number of elements:\t%li\n",global_num_elements);

# /* Write forest to vtu files. */
t8_step2_write_forest_vtk(forest, prefix);
c"t8_global_productionf"(" [step 2] Wrote forest to vtu files:\t%s*\n",prefix);

# /* Destroy the forest. */
t8_step2_destroy_forest(forest);
c"t8_global_productionf"(" [step 2] Destroyed forest.\n");

c"sc_finalize"();
