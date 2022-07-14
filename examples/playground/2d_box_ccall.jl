# vim: tabstop=2:softtabstop=2

using MPI
using Printf
# using T8code
# using CBinding

struct t8_forest end

Cptr = Ptr{Cvoid}
CChar = UInt8
MPI_Comm_t = Cptr
t8_locidx_t = Int32

T8DIR = ENV["JULIA_T8CODE_PATH"]

libt8 = "$(T8DIR)/lib/libt8.so"
libsc = "$(T8DIR)/lib/libsc.so"
# libp4 = "/home/jmark/install/t8code/lib/libt8.so"

const SC_LP_DEFAULT    =  -1     # /**< this selects the SC default threshold */
const SC_LP_ALWAYS     =   0     # /**< this will log everything */
const SC_LP_TRACE      =   1     # /**< this will prefix file and line number */
const SC_LP_DEBUG      =   2     # /**< any information on the internal state */
const SC_LP_VERBOSE    =   3     # /**< information on conditions, decisions */
const SC_LP_INFO       =   4     # /**< the main things a function is doing */
const SC_LP_STATISTICS =   5     # /**< important for consistency/performance */
const SC_LP_PRODUCTION =   6     # /**< a few lines for a major api function */
const SC_LP_ESSENTIAL  =   7     # /**< this logs a few lines max per program */
const SC_LP_ERROR      =   8     # /**< this logs errors only */
const SC_LP_SILENT     =   9     # /**< this never logs anything */

@enum t8_vtk_data_type_t begin
  T8_VTK_SCALAR                # /* One double value per element */
  T8_VTK_VECTOR                # /* 3 double values per element */
end

# P4EST_QUADRANT_LEN(l) = 1 << (P4EST_MAXLEVEL - l)

macro SC_ASSERT(q)
  :( $(esc(q)) ? nothing : throw(AssertionError($(string(q)))) )
end

macro P4EST_ASSERT(q)
  :( @SC_ASSERT($(esc(q)) ) )
end

macro T8_ASSERT(q)
  :( @SC_ASSERT($(esc(q)) ) )
end

macro T8_CCALL(fname, rtype, args...)
  # This macro is not supposed to be understood by weaklings ...
  quote
    function $(esc(fname))($((typeof(arg.args[1]) == Symbol ? esc(arg.args[1]) : esc(Expr(:kw, arg.args[1].args[1], arg.args[2])) for arg in args)...))
      ccall(($(string(fname)), libt8), 
        $(esc(rtype)), ($((esc(typeof(arg.args[1]) == Symbol ? arg.args[2] : arg.args[1].args[2]) for arg in args)...),), 
        $((esc(typeof(arg.args[1]) == Symbol ? arg.args[1] : arg.args[1].args[1]) for arg in args)...))
    end
  end
end

@T8_CCALL(t8_init, Cvoid, log_threshold :: Cint = SC_LP_PRODUCTION)

@T8_CCALL(sc_init, Cvoid, comm :: MPI_Comm_t, catch_signals :: Cint, print_backtrace :: Cint, log_handler :: Cptr, log_threshold :: Cint)

function sc_finalize()
  return ccall(("sc_finalize", libsc), Cvoid, ())
end

function sc_MPI_Finalize()
  return ccall(("sc_MPI_Finalize", libsc), Cvoid, ())
end

# t8_locidx_t t8_forest_get_num_ghosts (t8_forest_t forest);
function t8_forest_get_num_ghosts(forest)
  return ccall(("t8_forest_get_num_ghosts", libt8), t8_locidx_t, (Cptr,), forest)
end

@T8_CCALL(t8_forest_get_num_local_trees, t8_locidx_t, forest :: Cptr)

@T8_CCALL(t8_forest_get_local_num_elements, t8_locidx_t, forest :: Cptr)

# t8_locidx_t t8_forest_get_tree_num_elements (t8_forest_t forest, t8_locidx_t ltreeid);
function t8_forest_get_tree_num_elements(forest, ltreeid)
  return ccall(("t8_forest_get_tree_num_elements", libt8), t8_locidx_t, (Cptr,t8_locidx_t), forest, ltreeid)
end

# t8_eclass_t t8_forest_get_tree_class (t8_forest_t forest, t8_locidx_t ltreeid);
function t8_forest_get_tree_class(forest, ltreeid)
  return ccall(("t8_forest_get_tree_class", libt8), Cptr, (Cptr, t8_locidx_t), forest, ltreeid)
end

# t8_eclass_scheme_c *t8_forest_get_eclass_scheme (t8_forest_t forest, t8_eclass_t eclass);
function t8_forest_get_eclass_scheme(forest, eclass)
  return ccall(("t8_forest_get_eclass_scheme", libt8), Cptr, (Cptr, Cptr), forest, eclass)
end

function t8_cmesh_new_periodic(comm, ndim) :: Cptr
  return ccall(("t8_cmesh_new_periodic", libt8), Cptr, (MPI_Comm_t,Cint), comm, ndim)
end

@T8_CCALL(t8_cmesh_new_periodic_hybrid, Cptr, comm :: MPI_Comm_t)

function t8_cmesh_unref(cmesh)
  ccall(("t8_cmesh_unref",libt8), Cvoid, (Cptr,), pcmesh)
end

# sc_MPI_Comm t8_forest_get_mpicomm (t8_forest_t forest);
function t8_forest_get_mpicomm(forest)
  return ccall(("t8_forest_get_mpicomm",libt8), Cptr, (Cptr,), forest)
end

# void t8_forest_unref (t8_forest_t *pforest);
function t8_forest_unref(pforest)
  ccall(("t8_forest_unref",libt8), Cvoid, (Cptr,), pforest)
end

function t8_scheme_cxx_unref(pscheme)
  # void t8_scheme_cxx_unref (t8_scheme_cxx_t **pscheme);
  ccall(("t8_scheme_cxx_unref",libt8), Cvoid, (Cptr,), pscheme)
end

# function t8_cmesh_new_from_p4est(conn, do_partition = 0) :: Cptr
#   # cmesh = c"t8_cmesh_new_from_p4est"(conn, comm, do_partition)
#   # cmesh = C_NULL
#   # return cmesh
#   return C_NULL
# 
#   return ccall(("t8_cmesh_new_from_p4est", libt8), Cptr, (Cptr, MPI_Comm_t, Cint), conn, comm, do_partition)
# end

# function t8_scheme_new_default() :: Cptr
#   # t8_scheme_cxx_t    *t8_scheme_new_default_cxx (void);
#   return ccall(("t8_scheme_new_default_cxx", libt8), Cptr, ())
# end

@T8_CCALL(t8_scheme_new_default_cxx, Cptr)

# t8_forest_t t8_forest_new_uniform (t8_cmesh_t cmesh,
#                                    t8_scheme_cxx_t *scheme,
#                                    int level, int do_face_ghost,
#                                    sc_MPI_Comm comm);
# function t8_forest_new_uniform(cmesh,scheme,level,do_face_ghost,comm) :: Cptr
#   return ccall(("t8_forest_new_uniform",libt8), Cptr, 
#     (Cptr, Cptr, Cint, Cint, MPI_Comm_t), cmesh, scheme, level, do_face_ghost, comm)
# end
@T8_CCALL(t8_forest_new_uniform, Cptr, cmesh :: Cptr, scheme :: Cptr, level :: Cint, do_face_ghost :: Cint, comm :: MPI_Comm_t)

function string2ASCII(s,BUFSIZ=BUFSIZ)
  n = min(length(s),BUFSIZ)
  return NTuple{BUFSIZ,UInt8}(Vector{UInt8}(s[1:n] * "\0"^(BUFSIZ-n)))
end

@T8_CCALL(t8_forest_init, Cvoid, pforest :: Cptr)

function t8_forest_init()
  forest = Ptr{t8_forest}(0);
  forest_ref = Ref(forest);
  t8_forest_init(forest_ref);
  return forest_ref[];
end

# void t8_forest_set_user_data (t8_forest_t forest, void *data);
function t8_forest_set_user_data(forest, data)
  ccall(("t8_forest_set_user_data",libt8), Cvoid, (Cptr, Cptr), forest, data)
end

# void *t8_forest_get_user_data (t8_forest_t forest);
function t8_forest_get_user_data(forest)
  return ccall(("t8_forest_get_user_data",libt8), Cptr, (Cptr,), forest)
end

# int t8_forest_is_committed (t8_forest_t forest);
function t8_forest_is_committed(forest)
  return ccall(("t8_forest_is_committed",libt8), Cint, (Cptr,), forest)
end

# void t8_forest_set_adapt (t8_forest_t forest,
#                           const t8_forest_t set_from,
#                           t8_forest_adapt_t adapt_fn,
#                           int recursive);
function t8_forest_set_adapt(forest, set_from, adapt_fn, recursive=0)
  ccall(("t8_forest_set_adapt",libt8), Cvoid, (Cptr, Cptr, Cptr, Cint), forest, set_from, adapt_fn, recursive)
end

# void t8_forest_set_partition (t8_forest_t forest,
#                               const t8_forest_t set_from,
#                               int set_for_coarsening);
function t8_forest_set_partition(forest, set_from, set_for_coarsening=0)
  ccall(("t8_forest_set_partition",libt8), Cvoid, (Cptr, Cptr, Cint), forest, set_from, set_for_coarsening)
end

# void t8_forest_set_balance (t8_forest_t forest,
#                             const t8_forest_t set_from,
#                             int no_repartition);
function t8_forest_set_balance(forest, set_from, no_repartion=0)
  ccall(("t8_forest_set_balance",libt8), Cvoid, (Cptr, Cptr, Cint), forest, set_from, no_repartion)
end

# void t8_forest_commit (t8_forest_t forest);
function t8_forest_commit(forest)
  return ccall(("t8_forest_commit",libt8), Cvoid, (Cptr,), forest)
end

# int t8_element_level (t8_eclass_scheme_c *ts, const t8_element_t *elem);
function t8_element_level(ts, elem)
  return ccall(("t8_element_level",libt8), Cint, (Cptr,Cptr), ts, elem)
end

# double              t8_forest_element_volume (t8_forest_t forest,
#                                               t8_locidx_t ltreeid,
#                                               const t8_element_t *element);
function t8_forest_element_volume(forest, ltreeid, element)
  return ccall(("t8_forest_element_volume",libt8), Cdouble, (Cptr,t8_locidx_t, Cptr), forest, ltreeid, element)
end

# double *t8_forest_get_tree_vertices (t8_forest_t forest, t8_locidx_t ltreeid);
function t8_forest_get_tree_vertices(forest, ltreeid)
  return ccall(("t8_forest_get_tree_vertices",libt8), Ptr{Cdouble}, (Cptr,t8_locidx_t), forest, ltreeid)
end

# void t8_forest_element_centroid (t8_forest_t forest,
#                                  t8_locidx_t ltreeid,
#                                  const t8_element_t *element,
#                                  double *coordinates);
function t8_forest_element_centroid(forest, ltreeid, element, coordinates)
  ccall(("t8_forest_element_centroid",libt8), Cvoid, (Cptr,t8_locidx_t,Cptr,Ptr{Cdouble}), forest, ltreeid, element, coordinates)
end

# t8_element_t       *t8_forest_get_element_in_tree (t8_forest_t forest,
#                                                    t8_locidx_t ltreeid,
#                                                    t8_locidx_t leid_in_tree);
function t8_forest_get_element_in_tree(forest, ltreeid, leid_in_tree)
  return ccall(("t8_forest_get_element_in_tree",libt8), Cptr, (Cptr,t8_locidx_t,t8_locidx_t), forest, ltreeid, leid_in_tree)
end

# double t8_vec_dist (const double vec_x[3], const double vec_y[3]);
function t8_vec_dist(vec_x, vec_y)
  return ccall(("t8_vec_dist",libt8), Cdouble, (Ptr{Cdouble},Ptr{Cdouble}), vec_x, vec_y)
end

# int                 t8_forest_write_vtk_ext (t8_forest_t forest,
#                                              const char *fileprefix,
#                                              int write_treeid,
#                                              int write_mpirank,
#                                              int write_level,
#                                              int write_element_id,
#                                              int write_ghosts,
#                                              int write_curved,
#                                              int do_not_use_API,
#                                              int num_data,
#                                              t8_vtk_data_field_t *data);
@T8_CCALL(t8_forest_write_vtk_ext, Cint, 
  forest            :: Cptr, 
  fileprefix        :: Cstring, 
  write_treeid      :: Cint,
  write_mpirank     :: Cint,
  write_level       :: Cint,
  write_element_id  :: Cint,
  write_ghosts      :: Cint,
  write_curved      :: Cint,
  do_not_use_API    :: Cint,
  num_data          :: Cint,
  data              :: Cptr)

# function t8_global_productionf(fmtstr :: String, args...)
#   Printf.@printf( "[t8] "*fmtstr, args...)
# 
#   # j2c = Dict(String => Cstring, Int => Clong, Int32 => Cint)
#   # # @ccall libt8.t8_global_productionf(fmtstr :: Cstring, args :: Cstring) :: Cvoid
#   # 
#   # # Targs = (Cstring)
#   # # println((Cstring, Targs...))
#   # #
#   # Targs = :(Cstring,Tuple(j2c[typeof(arg)] for arg in  $args)...)
# 
#   # println(eval(Targs))
# 
#   # # Targs = :(Cstring, Cstring)
#   # # println(typeof(Targs))
# 
#   # # ccall(("t8_global_productionf",libt8), Cvoid, (Cstring, Targs...), fmtstr, args...)
#   # # ccall(("t8_global_productionf",libt8), Cvoid, (Cstring, Cstring), fmtstr, args...)
#   # # @eval ccall(("t8_global_productionf",libt8), Cvoid, $Targs, $fmtstr, $args...)
#   # println(:(ccall(("t8_global_productionf",$libt8), Cvoid, $Targs, $fmtstr, $(args...))))
#   # # eval(:(ccall(("t8_global_productionf",$libt8), Cvoid, ($(Targs...)), $fmtstr, $(args...))))
# end

# /* This is our own defined data that we will pass on to the
#  * adaptation callback. */
struct t8_step3_adapt_data_t
  midpoint                    :: Vector{Float64}
  refine_if_inside_radius     :: Float64
  coarsen_if_outside_radius   :: Float64
end

# /* The data that we want to store for each element.
#  * In this example we want to store the element's level and volume. */
struct t8_step5_data_per_element_t
  level   :: Cint
  volume  :: Cdouble
  state   :: Cdouble
end

const BUFSIZ = 8192

struct t8_vtk_data_field_t
  type        :: Int32
  description :: NTuple{BUFSIZ,UInt8}
  data        :: Ptr{Cvoid}
end

function _adapt_callback(forest,
                         forest_from,
                         which_tree,
                         lelement_id,
                         ts,
                         is_family, 
                         num_elements,
                         elements) :: Cint

  # tree_class = c"t8_forest_get_tree_class"(forest_from, which_tree);
  # eclass_scheme = c"t8_forest_get_eclass_scheme"(forest_from, tree_class);
  
  element = unsafe_load(elements,1)

  level = t8_element_level(ts, element)

  # /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
  #  * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. */
  centroid = Array{Float64}(undef,3);      # /* Will hold the element midpoint. */
  
  # /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
  #  * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
  #  * and we can now access it with t8_forest_get_user_data (forest). */
  adapt_data_ptr = Ptr{t8_step3_adapt_data_t}(t8_forest_get_user_data(forest))
  
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
  
  tree_vertices = t8_forest_get_tree_vertices(forest_from, which_tree);
  
  # /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid(forest_from, which_tree, element, pointer(centroid));

  # /* Compute the distance to our sphere midpoint. */
  dist = t8_vec_dist(pointer(centroid), pointer(adapt_data.midpoint));

  if level < level_max && dist < adapt_data.refine_if_inside_radius
    # /* Refine this element. */
    return 1;
  elseif level > level_min && num_elements > 1 && dist > adapt_data.coarsen_if_outside_radius
    # /* Coarsen this family. Note that we check for num_elements > 1 before, since returning < 0
    # * if we do not have a family as input is illegal. */
    return -1;
  end
  
  # /* Do not change this element. */
  return 0;
end

function _create_element_data(forest)
  # t8_locidx_t         num_local_elements;
  # t8_locidx_t         num_ghost_elements;
  # struct t8_step5_data_per_element *element_data;

  # /* Check that forest is a committed, that is valid and usable, forest. */
  @T8_ASSERT (t8_forest_is_committed(forest) != 0);

  # /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements(forest);
  # /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts(forest);

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
  num_local_trees = t8_forest_get_num_local_trees(forest);

  current_index = 0

  for itree = 0:num_local_trees-1
    # /* This loop iterates through all local trees in the forest. */
    # /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
    #  * also a different way to interpret its elements. In order to be able to handle elements
    #  * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */

    tree_class = t8_forest_get_tree_class(forest, itree);
    eclass_scheme = t8_forest_get_eclass_scheme(forest, tree_class);

    tree_vertices = t8_forest_get_tree_vertices(forest, itree);

    # /* Get the number of elements of this tree. */
    num_elements_in_tree = t8_forest_get_tree_num_elements(forest, itree);

    for ielement = 0:num_elements_in_tree-1
      current_index += 1
      # /* This loop iterates through all the local elements of the forest in the current tree. */
      # /* We can now write to the position current_index into our array in order to store
      #  * data for this element. */
      # /* Since in this example we want to compute the data based on the element in question,
      #  * we need to get a pointer to this element. */

      element = t8_forest_get_element_in_tree(forest, itree, ielement);

      # /* We want to store the elements level and its volume as data. We compute these
      #  * via the eclass_scheme and the forest_element interface. */
      level = t8_element_level(eclass_scheme, element)
      volume = t8_forest_element_volume(forest, itree, element);

      centroid = Array{Float64}(undef,3);      # /* Will hold the element midpoint. */
      t8_forest_element_centroid(forest, itree, element, pointer(centroid));

      # println("centroid = ", centroid)
      # println("volume = ", volume)

      # state = centroid[1]^2 - centroid[2]^2
      state = sin(7*centroid[1]) + cos(2*centroid[2])

      element_data[current_index] = t8_step5_data_per_element_t(level,volume,state)
    end # for
  end # for

  # println("centroid = ", centroid)

  return element_data;
end

# # /* Each process has computed the data entries for its local elements.
# #  * In order to get the values for the ghost elements, we use t8_forest_ghost_exchange_data.
# #  * Calling this function will fill all the ghost entries of our element data array with the
# #  * value on the process that owns the corresponding element. */
# function t8_step5_exchange_ghost_data(forest, data)
#   num_elements = c"t8_forest_get_local_num_elements"(forest);
#   num_ghosts = c"t8_forest_get_num_ghosts"(forest);
# 
#   # /* t8_forest_ghost_exchange_data expects an sc_array (of length num_local_elements + num_ghosts).
#   #  * We wrap our data array to an sc_array. */
#   sc_array_wrapper = c"sc_array_new_data"(Ref(data), sizeof(t8_step5_data_per_element_t), num_elements + num_ghosts);
# 
#   # /* Carry out the data exchange. The entries with indices > num_local_elements will get overwritten. */
#   c"t8_forest_ghost_exchange_data"(forest, sc_array_wrapper);
# 
#   # /* Destroy the wrapper array. This will not free the data memory since we used sc_array_new_data. */
#   c"sc_array_destroy"(sc_array_wrapper);
# end

# /* Write the forest as vtu and also write the element's volumes in the file.
#  * 
#  * t8code supports writing element based data to vtu as long as its stored
#  * as doubles. Each of the data fields to write has to be provided in its own
#  * array of length num_local_elements.
#  * We support two types: T8_VTK_SCALAR - One double per element
#  *                  and  T8_VTK_VECTOR - 3 doubles per element
#  */
function _output_data_to_vtu(forest, data, prefix)
  num_elements = t8_forest_get_local_num_elements(forest);

  # /* Copy the elment's volumes from our data array to the output array. */
  element_volumes = Vector{Cdouble}([data[i].volume for i in 1:num_elements])
  element_states  = Vector{Cdouble}([data[i].state for i in 1:num_elements])

  # # /* The number of user defined data fields to write. */
  num_data = 2;

  vtk_data = t8_vtk_data_field_t[
    t8_vtk_data_field_t(Cint(T8_VTK_SCALAR),string2ASCII("Element Volume"),pointer(element_volumes)),
    t8_vtk_data_field_t(Cint(T8_VTK_SCALAR),string2ASCII("Element State"),pointer(element_states))
  ]

  # /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
  #  * from t8_forest_vtk.h. Despite writing user data, it also offers more control over which 
  #  * properties of the forest to write. */
  write_treeid      = 1;
  write_mpirank     = 1;
  write_level       = 1;
  write_element_id  = 1;
  write_ghosts      = 0;

  t8_forest_write_vtk_ext(forest, prefix, write_treeid, write_mpirank,
                          write_level, write_element_id, write_ghosts,
                          0, 0, num_data, pointer(vtk_data));
end

# /* The prefix for our output files. */
const prefix_forest           = "2d_forest";
const prefix_forest_with_data = "2d_forest_with_volume_data";

# /* The uniform refinement level of the forest. */
const level_min = 4;
const level_max = 7;
const ndim = 2;

# /* Initialize MPI. This has to happen before we initialize sc or t8code. */
MPI.Init()
comm = MPI.COMM_WORLD.val

# /* Initialize the sc library, has to happen before we initialize t8code. */
sc_init(comm, 1, 1, C_NULL, SC_LP_VERBOSE);
# /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
# t8_init"(SC_LP_PRODUCTION);
t8_init(SC_LP_VERBOSE);

# cmesh = t8_cmesh_new_periodic(comm, ndim);
# cmesh = c"t8_cmesh_new_periodic_tri"(comm);
cmesh = t8_cmesh_new_periodic_hybrid(comm);
scheme = t8_scheme_new_default_cxx();

# /* Start with a uniform forest. */
root_forest = t8_forest_new_uniform(cmesh, scheme, level_min, 0, comm);

# /* Adapt, partition, balance and create ghost elements all in the same step. */
forest = t8_forest_init();

adapt_data = t8_step3_adapt_data_t(
  [0.5, 0.5, 0.0],            # /* Midpoints of the sphere. */
  0.2,                        # /* Refine if inside this radius. */
  0.4                         # /* Coarsen if outside this radius. */
)

# typedef int         (*t8_forest_adapt_t) (t8_forest_t forest,
#                                           t8_forest_t forest_from,
#                                           t8_locidx_t which_tree,
#                                           t8_locidx_t lelement_id,
#                                           t8_eclass_scheme_c *ts,
#                                           const int is_family,
#                                           const int num_elements,
#                                           t8_element_t *elements[]);
c_adapt_callback = @cfunction(_adapt_callback, Cint, 
  (Ptr{t8_forest}, Ptr{t8_forest}, t8_locidx_t, t8_locidx_t, Cptr, Cint, Cint, Ptr{Cptr}))

t8_forest_set_user_data(forest, Ref(adapt_data));

t8_forest_set_adapt(forest, root_forest, c_adapt_callback, 1);

t8_forest_set_partition(forest, C_NULL, 0);
t8_forest_set_balance(forest, C_NULL, 0);
# # c"t8_forest_set_ghost(forest, 1, c"T8_GHOST_FACES");
t8_forest_commit(forest);

# /*
# * Build data array and gather data for the local elements.
# */
data = _create_element_data(forest);

# t8_global_productionf(" [step5] Computed level and volume data for local elements.\n");
# if t8_forest_get_local_num_elements(forest) > 0
#     # /* Output the stored data of the first local element (if it exists). */
#     c"t8_global_productionf"(" Element 0 has level %i and volume %e.\n",data[1].level, data[1].volume);
# end

# # /*
# # * Exchange the data values of the ghost elements
# # */
# t8_step5_exchange_ghost_data(forest, data);
# c"t8_global_productionf"(" [step5] Exchanged ghost data.\n");
# 
# if c"t8_forest_get_num_ghosts"(forest) > 0
#   # /* output the data of the first ghost element (if it exists) */
#   first_ghost_index = c"t8_forest_get_local_num_elements"(forest);
#   c"t8_global_productionf"(" [step5] Ghost 0 has level %i and volume %e.\n",
#                        data[first_ghost_index+1].level,
#                        data[first_ghost_index+1].volume);
# end

# /*
# * Output the volume data to vtu.
# */
_output_data_to_vtu(forest, data, prefix_forest);
# c"t8_global_productionf"(" Wrote forest and volume data to %s*.\n", prefix_forest);

# /*
# * clean-up
# */

# /* Destroy the forest. */
t8_forest_unref(Ref(forest));
# t8_forest_unref(Ref(root_forest));

# t8_cmesh_unref(Ref(cmesh));
# t8_scheme_cxx_unref(Ref(scheme));

sc_finalize();
# sc_MPI_Finalize();
