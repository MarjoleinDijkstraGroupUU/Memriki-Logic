# memriki-logic

This projects is a collection of simulation scripts for electrical circuits based on memristive[^1] microfluidic ionic channels[^2]. Shinriki inspired osciallators [^3] are build and then connected together to form logic gates. These scripts are the basis for the paper [Neuromorphic Computing with Microfluidic Memristors]().

To reproduce the plots of the paper follow these instructions:

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> memriki-logic

To (locally) reproduce this project, do the following:

0. Download this code base.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "memriki-logic"
```
which auto-activate the project and enable local path handling from DrWatson.

The conductance-voltage and current-voltage plots of a single memristor can be reproduced by running `scripts/memristor.jl`.
Plots showing the dynamics of the oscillator can be generated with `scripts/dynamics.jl`
The `NAND` and `XOR` gates can be reproduced by `scripts/nand-xor.jl`.
The code for `AND` and `OR` gates assembled from `NAND` gates can be found in `scripts/combining-nands.jl`.
To test the repeatability of the gates run `scripts/repeatability.jl`


[^1]: L. Chua, Memristor-the missing circuit element, IEEE Transactions on Circuit Theory 18, 507 (1971) [https://doi.org/10.1109/TCT.1971.1083337](https://doi.org/10.1109/TCT.1971.1083337)
[^2]: T. Kamsma, W. Boon, T. ter Rele, C. Spitoni, and R. van Roij, Iontronic neuromorphic signaling with conical microfluidic memristors, Physical Review Letters 130, 268401 (2023) [https://doi.org/10.1103/PhysRevLett.130.268401](https://doi.org/10.1103/PhysRevLett.130.268401)
[^3]: M. Shinriki, M. Yamamoto, and S. Mori, Multimode oscillations in a modified van der pol oscillator containing a positive nonlinear conductance, Proceedings of the IEEE 69, 394 (1981) [https://doi.org/10.1109/proc.1981.11973](https://doi.org/10.1109/proc.1981.11973)
