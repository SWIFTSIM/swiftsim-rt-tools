Some grackle standalone programs.

Note: We generally use the 64 bit float version of grackle. This should be the
default setup when compiling it, and you'd need to manually change that during
installation. 

Note: The grackle library is experiencing active development (state 2023). The API
might change in the future. We keep a frozen version forked on https://github.com/mladenivkovic/grackle-swift .
This version is guaranteed to work with swift (and with this repository).

- CoolingTest: (Re)produce the results from `examples/RadiativeTransferTests/CoolingTest`

- CosmoCoolingTest: Same as CoolingTest, but including cosmological expansion.

- Iliev06Test0Part3: (Re)produce the results from Iliev et al 2006, Test 0, part 3.
  (In which a gas parcel is first heated, and then cooled.)

- include: some global grackle related includes that are used in each example.
