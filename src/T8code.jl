module T8code

using CBinding;

# Old hack: deprecated
# using Base.Libc.Libdl
# p4 = dlopen("$(LIBDIR)/libp4est.so",RTLD_GLOBAL)

T8DIR = ENV["JULIA_T8CODE_PATH"]

INCDIR = "$(T8DIR)/include"
LIBDIR = "$(T8DIR)/lib"

c`
    -Wall
    -fparse-all-comments
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
# include <p4est_extended.h>
# include <p8est.h>
# include <p8est_extended.h>
# include <t8.h>
# include <t8_cmesh.h>
# include <t8_cmesh_vtk.h>
# include <t8_forest.h>
# include <t8_forest_vtk.h>
"""i

c"t8_scheme_cxx_t    *t8_scheme_new_default_cxx (void);"

end # module
