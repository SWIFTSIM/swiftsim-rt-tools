GEAR-RT: Use of boundary particles
----------------------------------------------------

In some tests, an additional treatment of boundary particles is necessary. E.g.
for the Iliev Test 2, the gas density is too low to absorb a sufficient amount
of radiation before it reaches the box edge, and the additional radiation
influences the outer regions of the Stromgren sphere.

In order to alleviate the issues, some special treatment of boundary particles
is necessary. Three things are required:

1)  Define which particles are the boundary particles through their ID. Particles
    with an ID above the specified threshold ID are treated as boundary
    particles. You need to make sure to assign the correct IDs while generating
    ICs. Note that if you're generating ICs with gas and star particles using
    swiftsimio, you need to assign particle IDs manually to *all* particle
    types, otherwise swiftsimio will overwrite them.
    
2)  Keep the boundary particles static. You can either declare all particles in 
    the simulation boundary particles at compile time, or you can activate the 
    `#GIZMO_FIXED_PARTICLES` macro in `src/const.h`. Alternatively, you can
    configure swift to use boundary particles above a provided particle ID,
    while the other particles are treated normally. Consult the swift 
    documentation on how to do that.

3)  Add the special treatment of boundary particles to the code. It should
    suffice to call the `rt_tchem_set_boundary_particles_for_test(p)` funciton
    during the thermochemistry part. You can then set manually what to do with
    boundary particles.
