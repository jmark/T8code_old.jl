module T8code

using CBinding;

export @P4EST_ASSERT
export @T8_ASSERT

# Old hack: deprecated
# using Base.Libc.Libdl
# p4 = dlopen("$(LIBDIR)/libp4est.so",RTLD_GLOBAL)

T8DIR = ENV["JULIA_T8CODE_PATH"]

INCDIR = "$(T8DIR)/include"
LIBDIR = "$(T8DIR)/lib"

c`
    -Wall
    -DSC_ENABLE_MPI=1
    -I$(INCDIR)
    -I$(INCDIR)/src
    -L$(LIBDIR) -lt8`

    # -L$(LIBDIR) -lsc -lp4est -lt8`

const c"int8_t" = Int8
const c"int16_t" = Int16
const c"int32_t" = Int32
const c"int64_t" = Int64
const c"uint8_t" = UInt8
const c"uint16_t" = UInt16
const c"uint32_t" = UInt32
const c"uint64_t" = UInt64
const c"size_t" = Csize_t
const c"wchar_t" = Cwchar_t
const c"wchar_t" = Cwchar_t

# define them as Cvoid since they are usually used as opaque
# const c"FILE"       = Cvoid
# const c"va_list"    = Cvoid
# const c"MPI_Comm"   = Cvoid

c"""
# include <stdint.h>
# include <ctype.h>
# include <stdarg.h>
# include <stddef.h>
# include <stdio.h>
# include <mpi.h>

# include <sc.h>
# include <sc_shmem.h>
# include <sc_refcount.h>
# include <p4est.h>
# include <p4est_vtk.h>
# include <p4est_extended.h>
# include <p8est.h>
# include <p8est_vtk.h>
# include <p8est_extended.h>
# include <t8.h>
# include <t8_element.h>
# include <t8_cmesh.h>
# include <t8_cmesh_vtk.h>
# include <t8_forest.h>
# include <t8_cmesh/t8_cmesh_examples.h>

# include <t8_forest_vtk.h>

# include <t8_element_c_interface.h>
# include <t8_vec.h>
"""i

c"t8_scheme_cxx_t    *t8_scheme_new_default_cxx (void);"

# c"""
# typedef struct {
#     int dummy;
# } t8_element_t;"""

# c"""
#     typedef struct t8_element t8_element_t;
# """

# c"""
# int t8_forest_adapt_t (t8_forest_t forest,
#                           t8_forest_t forest_from,
#                           t8_locidx_t which_tree,
#                           t8_locidx_t lelement_id,
#                           t8_eclass_scheme_c * ts,
#                           int num_elements,
#                           t8_element_t ** elements);
# """

macro SC_ASSERT(q)
    :( $(esc(q)) ? nothing : throw(AssertionError($(string(q)))) )
end

macro P4EST_ASSERT(q)
    :( @SC_ASSERT($(esc(q)) ) )
end

P4EST_QUADRANT_LEN(l) = 1 << (c"P4EST_MAXLEVEL" - l)

macro T8_ASSERT(q)
    :( @SC_ASSERT($(esc(q)) ) )
end

end # module
