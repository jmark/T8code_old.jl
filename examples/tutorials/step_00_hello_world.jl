
using MPI
using T8code
using CBinding
using Base.Libc.Libdl

HOME = ENV["HOME"];

T8DIR = "$(HOME)/opt/t8code"

INCDIR = "$(T8DIR)/include"
LIBDIR = "$(T8DIR)/lib"

p4 = dlopen("$(LIBDIR)/libp4est.so",RTLD_GLOBAL)

# /* Initialize MPI. This has to happen before we initialize sc or T8. */
MPI.Init()
comm = MPI.COMM_WORLD.val

# /* Initialize the sc library, has to happen before we initialize T8. */
c"sc_init"(comm, 1, 1, C_NULL, c"SC_LP_VERBOSE")

# /* Initialize T8 with log level SC_LP_PRODUCTION. See sc.h for more info on the leg levels. */
c"t8_init"(c"SC_LP_VERBOSE")

# # /* Print a message on the root process. */
c"t8_global_productionf"(" [step0] \n")
c"t8_global_productionf"(" [step0] Hello, this is T8 :)\n")
c"t8_global_productionf"(" [step0] \n")

c"sc_finalize"()
