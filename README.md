# T8code.jl

**T8code.jl** is a Julia package that wraps
[t8code](https://github.com/holke/t8code), a C library to manage multiple
connected adaptive binary trees, quadtrees, and octrees in parallel.


## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). T8code.jl works
with Julia v1.6 or newer.

T8code.jl depends on [t8code](https://github.com/holke/t8code) and a MPI
distribution provided by your system.

### Setup
For example, if your t8code library is installed to `/opt/t8code`,
use it from T8code.jl by executing
```bash
julia --project -e 'ENV["JULIA_T8CODE_PATH"] = "/opt/t8code";'
```

### Usage
Check out the `examples/tutorials` directory on how to use this package.

## Authors
T8code.jl is maintained by
[Johannes Markert](https://www.hlrs.de/people/schlottke-lakemper)
(German Aerospace Center (DLR), Germany),
The [t8code](https://github.com/holke/t8code) library itself is written by
Johannes Holke, Carsten Burstedde, and many others.


## License and contributing
T8code.jl is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
[t8code](https://github.com/holke/t8code) itself is licensed under the GNU
General Public License, version 2.
