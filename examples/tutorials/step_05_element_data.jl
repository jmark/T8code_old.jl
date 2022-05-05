# vim: tabstop=2:softtabstop=2

using MPI
using T8code
using CBinding

NULL = C_NULL

# /* This is step5 of the t8code tutorials.
#  * 
#  * TODO: This file still needs to be documented.
#  * How you can experiment here:
#  *   -
#  *  */

# /* This is our own defined data that we will pass on to the
#  * adaptation callback. */
mutable struct t8_step3_adapt_data_t
  midpoint                    :: Vector{Float64}
  refine_if_inside_radius     :: Float64
  coarsen_if_outside_radius   :: Float64
end

# /* The data that we want to store for each element.
#  * In this example we want to store the element's level and volume. */
mutable struct t8_step5_data_per_element_t
  level   :: Cint
  volume  :: Float64
end

# mutable struct t8_step5_vtk_datafield_t
#   type :: Int32
#   description = 
# end

function t8_step3_adapt_callback(forest,
                         forest_from,
                         which_tree,
                         lelement_id,
                         ts,
                         is_family, 
                         num_elements,
                         elements)

  # /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
  #  * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */
  centroid = Array{Float64}(undef,3);      # /* Will hold the element midpoint. */
  
  # /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
  #  * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
  #  * and we can now access it with t8_forest_get_user_data (forest). */
  adapt_data_ptr = Ptr{t8_step3_adapt_data_t}(c"t8_forest_get_user_data"(forest))
  
  # /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
  #  * If the condition is not true, then the code will abort.
  #  * In this case, we want to make sure that we actually did set a user pointer to forest and thus
  #  * did not get the NULL pointer from t8_forest_get_user_data.
  #  */
  @T8_ASSERT adapt_data_ptr != C_NULL
  
  adapt_data = unsafe_load(adapt_data_ptr)
  
  # /* In order to compute the coordinates of the element, we first need to get the coordinates of
  #  * the tree. t8code does not store the element's coordinates, since they can be interpolated from the
  #  * tree's coordinates. */
  
  tree_vertices = c"t8_forest_get_tree_vertices"(forest_from, which_tree);
  
  # /* Compute the element's centroid coordinates. */
  c"t8_forest_element_centroid"(forest_from, which_tree, elements[], Ref(centroid));
  
  # /* Compute the distance to our sphere midpoint. */
  dist = c"t8_vec_dist"(Ref(centroid), Ref(adapt_data.midpoint));
  if dist < adapt_data.refine_if_inside_radius
    # /* Refine this element. */
    return Int32(1);
  elseif num_elements > 1 && dist > adapt_data.coarsen_if_outside_radius
    # /* Coarsen this family. Note that we check for num_elements > 1 before, since returning < 0
    # * if we do not have a family as input is illegal. */
    return Int32(-1);
  end
  
  # /* Do not change this element. */
  return Int32(0);
end

function t8_step5_build_forest(comm, level)
  cmesh = c"t8_cmesh_new_hypercube_hybrid"(comm, 0, 0);
  scheme = c"t8_scheme_new_default_cxx"();

  # /* Start with a uniform forest. */
  forest = c"t8_forest_new_uniform"(cmesh, scheme, level, 0, comm);

  # /* Adapt, partition, balance and create ghost elements all in the same step. */
  forest_apbg = c"t8_forest_t"();
  forest_apbg_ref = Ref(forest_apbg);
  c"t8_forest_init"(forest_apbg_ref);
  forest_apbg = forest_apbg_ref[];

  adapt_data = t8_step3_adapt_data_t(
    [0.5, 0.5, 1],              # /* Midpoints of the sphere. */
    0.2,                        # /* Refine if inside this radius. */
    0.4                         # /* Coarsen if outside this radius. */
  )

  c"t8_forest_set_user_data"(forest_apbg, Ref(adapt_data));
  c"t8_forest_set_adapt"(forest_apbg, forest, t8_step3_adapt_callback, 0);
  c"t8_forest_set_partition"(forest_apbg, NULL, 0);
  c"t8_forest_set_balance"(forest_apbg, NULL, 0);
  c"t8_forest_set_ghost"(forest_apbg, 1, c"T8_GHOST_FACES");
  c"t8_forest_commit"(forest_apbg);

  return forest_apbg;
end

function t8_step5_create_element_data(forest)
  # t8_locidx_t         num_local_elements;
  # t8_locidx_t         num_ghost_elements;
  # struct t8_step5_data_per_element *element_data;

  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (c"t8_forest_is_committed"(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = c"t8_forest_get_local_num_elements"(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = c"t8_forest_get_num_ghosts"(forest);

  # /* Now we need to build an array of our data that is as long as the number
  #  * of elements plus the number of ghosts. You can use any allocator such as
  #  * new, malloc or the t8code provide allocation macro T8_ALLOC. 
  #  * Note that in the latter case you need
  #  * to use T8_FREE in order to free the memory.
  #  */
  element_data = Vector{t8_step5_data_per_element_t}(undef, num_local_elements + num_ghost_elements)

  # /* Note: We will later need to associate this data with an sc_array in order to exchange the values for
  #  *       the ghost elements, which we can do with sc_array_new_data (see t8_step5_exchange_ghost_data).
  #  *       We could also have directly allocated the data here in an sc_array with
  #  *       sc_array_new_count (sizeof (struct data_per_element), num_local_elements + num_ghost_elements);
  #  */

  # /* Let us now fill the data with something.
  #  * For this, we iterate through all trees and for each tree through all its elements, calling
  #  * t8_forest_get_element_in_tree to get a pointer to the current element.
  #  * This is the recommended and most performant way.
  #  * An alternative is to iterate over the number of local elements and use
  #  * t8_forest_get_element. However, this function needs to perform a binary search
  #  * for the element and the tree it is in, while t8_forest_get_element_in_tree has a
  #  * constant look up time. You should only use t8_forest_get_element if you do not know
  #  * in which tree an element is.
  #  */

  # /* Get the number of trees that have elements of this process. */
  num_local_trees = c"t8_forest_get_num_local_trees"(forest);

  current_index = 0

  for itree = 0:num_local_trees-1
    # /* This loop iterates through all local trees in the forest. */
    # /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
    #  * also a different way to interpret its elements. In order to be able to handle elements
    #  * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */

    tree_class = c"t8_forest_get_tree_class"(forest, itree);
    eclass_scheme = c"t8_forest_get_eclass_scheme"(forest, tree_class);

    tree_vertices = c"t8_forest_get_tree_vertices"(forest, itree);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = c"t8_forest_get_tree_num_elements"(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      current_index += 1
      # /* This loop iterates through all the local elements of the forest in the current tree. */
      # /* We can now write to the position current_index into our array in order to store
      #  * data for this element. */
      # /* Since in this example we want to compute the data based on the element in question,
      #  * we need to get a pointer to this element. */

      element = c"t8_forest_get_element_in_tree"(forest, itree, ielement);

      # /* We want to store the elements level and its volume as data. We compute these
      #  * via the eclass_scheme and the forest_element interface. */
      level = c"t8_element_level"(eclass_scheme, element)
      volume = c"t8_forest_element_volume"(forest, itree, element);

      element_data[current_index] = t8_step5_data_per_element_t(level,volume)
    end # for
  end # for

  return element_data;
end

# /* Each process has computed the data entries for its local elements.
#  * In order to get the values for the ghost elements, we use t8_forest_ghost_exchange_data.
#  * Calling this function will fill all the ghost entries of our element data array with the
#  * value on the process that owns the corresponding element. */
function t8_step5_exchange_ghost_data(forest, data)
  num_elements = c"t8_forest_get_local_num_elements"(forest);
  num_ghosts = c"t8_forest_get_num_ghosts"(forest);

  # /* t8_forest_ghost_exchange_data expects an sc_array (of length num_local_elements + num_ghosts).
  #  * We wrap our data array to an sc_array. */
  sc_array_wrapper = c"sc_array_new_data"(Ref(data), sizeof(t8_step5_data_per_element_t), num_elements + num_ghosts);

  # /* Carry out the data exchange. The entries with indices > num_local_elements will get overwritten. */
  c"t8_forest_ghost_exchange_data"(forest, sc_array_wrapper);

  # /* Destroy the wrapper array. This will not free the data memory since we used sc_array_new_data. */
  c"sc_array_destroy"(sc_array_wrapper);
end

# /* Write the forest as vtu and also write the element's volumes in the file.
#  * 
#  * t8code supports writing element based data to vtu as long as its stored
#  * as doubles. Each of the data fields to write has to be provided in its own
#  * array of length num_local_elements.
#  * We support two types: T8_VTK_SCALAR - One double per element
#  *                  and  T8_VTK_VECTOR - 3 doubles per element
#  */
function t8_step5_output_data_to_vtu(forest, data, prefix)
  num_elements = c"t8_forest_get_local_num_elements"(forest);

  # println("num_elements: ", num_elements)
  # println("length data:  ", length(data))
  # println("volume of 1:  ", data[end].volume)

  # /* Copy the elment's volumes from our data array to the output array. */
  element_volumes = [data[i].volume for i in 1:num_elements]

  # element_volumes = Vector{Float64}(undef,num_elements)
  # for ielem = 1:num_elements
  #   # element_volumes[ielem] = data[ielem].volume;
  #   element_volumes[ielem] = 42.0
  # end

  # # /* The number of user defined data fields to write. */
  num_data = 1;

  # println(c"BUFSIZ")

  # FIXME: c"t8_vtk_data_field_t" chokes on large 'description' field (BUFSIZ > 8,000 !!)

  # # /* For each user defined data field we need one t8_vtk_data_field_t variable */
  vtk_data = c"t8_vtk_data_field_t"(
      # /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
      type = c"T8_VTK_SCALAR",
      # /* The name of the field as should be written to the file. */
      description = "Element Volume",
      data = Ref(element_volumes),
  )

  # /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
  #  * from t8_forest_vtk.h. Despite writing user data, it also offers more control over which 
  #  * properties of the forest to write. */
  write_treeid      = 1;
  write_mpirank     = 1;
  write_level       = 1;
  write_element_id  = 1;
  write_ghosts      = 0;

  # FIXME: 'element_volumes' does not get properly written to file ??

  c"t8_forest_vtk_write_file"(forest, prefix, write_treeid, write_mpirank,
                              write_level, write_element_id, write_ghosts,
                              num_data, Ref(vtk_data));
end

# /* The prefix for our output files. */
const prefix_forest           = "t8_step5_forest";
const prefix_forest_with_data = "t8_step5_forest_with_volume_data";

# /* The uniform refinement level of the forest. */
const level = 3;

# /* Initialize MPI. This has to happen before we initialize sc or t8code. */
MPI.Init()
comm = MPI.COMM_WORLD.val

# /* Initialize the sc library, has to happen before we initialize t8code. */
c"sc_init"(comm, 1, 1, NULL, c"SC_LP_ESSENTIAL");
# /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
# c"t8_init"(c"SC_LP_PRODUCTION");
c"t8_init"(c"SC_LP_VERBOSE");

# /* Print a message on the root process. */
c"t8_global_productionf"(" [step5] \n");
c"t8_global_productionf"(" [step5] Hello, this is the step5 example of t8code.\n");
c"t8_global_productionf"(" [step5] In this example we will store data on our elements and exchange the data of ghost elements.\n");
c"t8_global_productionf"(" [step5] \n");

# /*
# * Setup.
# * Build cmesh and uniform forest.
# */
c"t8_global_productionf"(" [step5] \n");
c"t8_global_productionf"(" [step5] Creating an adapted forest as in step3.\n");
c"t8_global_productionf"(" [step5] \n");

forest = t8_step5_build_forest(comm, level);
c"t8_forest_write_vtk"(forest, prefix_forest);
c"t8_global_productionf"(" [step5] Wrote forest to vtu files: %s*\n", prefix_forest);

# /*
# * Build data array and gather data for the local elements.
# */
data = t8_step5_create_element_data(forest);

c"t8_global_productionf"(" [step5] Computed level and volume data for local elements.\n");
if c"t8_forest_get_local_num_elements"(forest) > 0
    # /* Output the stored data of the first local element (if it exists). */
    c"t8_global_productionf"(" [step5] Element 0 has level %i and volume %e.\n",data[1].level, data[1].volume);
end

# /*
# * Exchange the data values of the ghost elements
# */
t8_step5_exchange_ghost_data(forest, data);
c"t8_global_productionf"(" [step5] Exchanged ghost data.\n");

if c"t8_forest_get_num_ghosts"(forest) > 0
  # /* output the data of the first ghost element (if it exists) */
  first_ghost_index = c"t8_forest_get_local_num_elements"(forest);
  c"t8_global_productionf"(" [step5] Ghost 0 has level %i and volume %e.\n",
                       data[first_ghost_index+1].level,
                       data[first_ghost_index+1].volume);
end

# /*
# * Output the volume data to vtu.
# */
t8_step5_output_data_to_vtu(forest, data, prefix_forest_with_data);
c"t8_global_productionf"(" [step5] Wrote forest and volume data to %s*.\n", prefix_forest_with_data);

# /*
# * clean-up
# */

# /* Destroy the forest. */
c"t8_forest_unref"(Ref(forest));
c"t8_global_productionf"(" [step5] Destroyed forest.\n");

c"sc_finalize"();
c"sc_MPI_Finalize"();
