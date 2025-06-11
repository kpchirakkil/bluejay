# Verification of Jacobian and Rate Functions

This note documents that the multicolumn model keeps a consistent ordering of species, altitude, and horizontal columns when building transport matrices.

## Indexing order
- `flatten_atm` builds vectors by species → altitude → horizontal.
- `unflatten_atm` reshapes back with `(length(species_list), GV.num_layers, n_horiz)`.

## Horizontal offsets
- In `chemJmat`, the base index for a given column/altitude uses `(ihoriz-1)*(length(GV.active_longlived)*GV.num_layers)`.
- When connecting neighboring columns, the Jacobian indices shift by `GV.num_layers*length(GV.active_longlived)`.

## Reshaping in rate functions
- `ratefn` reshapes inputs using `(species, altitude, horizontal)` ordering before looping over columns and altitudes.

## Temperature arrays
- `MODEL_SETUP.jl` initializes `Tn_temp`, `Ti_temp`, and `Te_temp` with dimensions `(n_horiz, num_layers+2)` and accesses them as `Tn_temp[ihoriz, :]` etc.

These observations verify that the original Julia code consistently follows the species → altitude → horizontal convention across the major files.