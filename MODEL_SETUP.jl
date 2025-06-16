################################################################################
# MODEL_SETUP.jl
# DESCRIPTION: Sets up variables, constants, etc. that typically don't require
# checking by the user before each run but change depending on the type of 
# model being run.
# 
# Eryn Cangi
# Created April 2024
# Last edited: May 2024
# Currently tested for Julia: 1.8.5
################################################################################

using DataFrames
using DoubleFloats

# First do some error checking
if (special_seasonal_case!=nothing) & (exp_type=="all")
    throw("Only use exp_type='all' when special_seasonal_case != nothing")
end

# ***************************************************************************************************** #
#                                                                                                       #
#                                       Planet-dependent constants                                      #
#                                                                                                       #
# ***************************************************************************************************** #
# Yeah ok they're not really constants but this is the easiest way to manage it because they
# are used in the photochemistry module, which gets loaded before the parameter/constant files, and which
# cannot load something based on a variable defined in the enclosing scope. 
const M_P = Dict( # Planetary mass in g 
                 "Mars"=>0.1075*5.972e27, 
                 "Venus"=>4.867e27
                )[planet]
const R_P = Dict( # Planetary radius in cm
                 "Mars"=>3396e5, 
                 "Venus"=>6050e5
                )[planet] 
const DH = Dict( # Atmospheric D/H ratio 
                "Mars"=>5.5 * 1.6e-4, # Yung 1988
                "Venus"=>240 * 1.6e-4,
               )[planet]
const sol_in_sec = Dict(
                        "Mars"=>88775.2438,   # One Mars sol in seconds
                        "Venus"=>2.09968e7
                       )[planet]
const season_in_sec = Dict(
                           "Mars"=>1.4838759e7,
                           "Venus"=>4.8535373e6
                          )[planet]
const g = bigG * M_P / (R_P^2);
const SA = 4*pi*(R_P)^2 # cm^2
 

# ***************************************************************************************************** #
#                                                                                                       #
#              Model species, Jrate lists, and lists of chemistry/transport species                     #
#                                                                                                       #
# ***************************************************************************************************** #

remove_ignored_species = true # Whether to use a slightly smaller list of species and reactions (removing minor species that Roger had in his model)
ignored_species = [:CNpl,:HCNpl,:HCNHpl,:HN2Opl,:NH2pl,:NH3pl,:N2Opl,:NO2pl,:CH,:CN,:HCN,:HNO,:NH,:NH2,:HD2pl]#:N2O,:NO2

#                                        Neutrals
# =======================================================================================================
const orig_neutrals = [:O,
                    #   :HO2, :HOCO, 
                        :O2, :O3,
                    #   :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD,

                       # Turn these off for minimal ionosphere:
                     #  :C, :HCN, :HCO, :N, :NO, 
                       ]; 
const conv_neutrals = remove_ignored_species==true ? setdiff(orig_neutrals, ignored_species) : orig_neutrals
const new_neutrals = [];
const neutral_species = [conv_neutrals..., new_neutrals...];

#                                          Ions
# =======================================================================================================
const orig_ions = [:Opl] #, :Opl, :O2pl, :COpl];# Nair minimal ionosphere 
               #    :Arpl, :ArHpl, :ArDpl, 
               #    :Cpl, :CHpl,  :COpl, 
               #    :Hpl, :Dpl, :H2pl, :HDpl, :H3pl, :H2Dpl, :HD2pl, 
               #    :H2Opl,  :HDOpl, :H3Opl, :H2DOpl, 
               #    :HO2pl, :HCOpl, :DCOpl, :HOCpl, :DOCpl, :DCO2pl, 
               #    :HNOpl,   
               #    :Npl, :NHpl, :N2pl, :N2Hpl, :N2Dpl, :NOpl,
               #    :OHpl, :ODpl];
const new_ions = [];
const ion_species = remove_ignored_species==true ? setdiff([orig_ions..., new_ions...], ignored_species) : [orig_ions..., new_ions...]
const nontherm = ions_included==true ? true : false   # whether to do non-thermal escape; this has to be here because it's needed in short order to do Jrates.
 
#                                     Full species list
# =======================================================================================================
const all_species = [neutral_species..., ion_species...];


#                        Photolysis and Photoionization rate symbol lists 
# =======================================================================================================

const conv_Jrates, newJrates = format_Jrates(reaction_network_spreadsheet, all_species, "Jratelist"; ions_on=ions_included, hot_atoms=nontherm)
const Jratelist = [conv_Jrates..., newJrates...];

# These dictionaries specify the species absorbing a photon for each J rate, and the products of the reaction.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])


#                               Miscellaneous logical groupings
# =======================================================================================================
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = get_deuterated(all_species)
const D_ions = get_deuterated(ion_species) #[s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];


#                             Define short- and long-lived species
# =======================================================================================================

# Short lived species, whose chemical lifetime is << diffusion timescale ------- #
const short_lived_species = [];# technically shortlived but count as longlived: :CH, :HCO, :HO2, :O3, :OH, :O1D, :DO2, :OD...
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2, :C, :CH])
    append!(short_lived_species, ion_species)
end

# Long lived species 
const long_lived_species = setdiff(all_species, short_lived_species)


#                         Species participating in chemistry and transport
# =======================================================================================================

no_chem_species = []; 
no_transport_species = [];

# Fixed species (Densities don't update)
# -------------------------------------------------------------------
for s in dont_compute_chemistry
    if ~(s in no_chem_species)
        push!(no_chem_species, s)
    end
end 
for s in dont_compute_transport
    if ~(s in no_transport_species)
        push!(no_transport_species, s)
    end
end

for s in dont_compute_either_chem_or_transport
    push!(no_chem_species, s)
    push!(no_transport_species, s)
end

# Chemistry and transport participants
# -------------------------------------------------------------------
if converge_which == "neutrals"
    append!(no_chem_species, union(conv_ions, N_neutrals)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_neutrals))
    append!(no_transport_species, setdiff(conv_neutrals, N_neutrals))
elseif converge_which == "both"
    append!(no_transport_species, short_lived_species)
end

# Disallow transport and/or chemistry if the appropriate setting is toggled
# TODO: Currently not working :(
if do_trans==false
    append!(no_transport_species, all_species)
end

if do_chem==false
    append!(no_chem_species, all_species)
end

const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species 
# -------------------------------------------------------------------
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# Sort name lists created here
# -------------------------------------------------------------------
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)
sort!(active_species)
sort!(inactive_species)
sort!(short_lived_species)
sort!(long_lived_species)
sort!(active_longlived)
sort!(active_shortlived)
sort!(chem_species)
sort!(transport_species)
sort!(no_chem_species)
sort!(no_transport_species)
sort!(D_bearing_species)
sort!(D_ions)
sort!(N_neutrals)


# ***************************************************************************************************** #
#                                                                                                       #
#                                        Atmospheric setup                                              #
#                                                                                                       #
# ***************************************************************************************************** #

#                             Ion chemistry and non-thermal escape
# =======================================================================================================

const e_profile_type = ions_included==true ? "quasineutral" : "none" 
    # OPTIONS: 
    # "quasineutral" - n_e = sum of all the ion densities; PREFERRED
    # "O2+" - n_e = sum(n_O2pl)
    # "constant" - n_e set to some constant value which I believe is 1e5.
    # "none" - no electrons, for use in neutral-only models

#                                       Altitude grid                  
# =======================================================================================================
const zmin = Dict("Venus"=>90e5, "Mars"=>0.)[planet]
const dz = 2e5  # Discretized layer thickness
const dx = 2e5  # Width of one vertical column -- could be used to calculate horizontal transport flux boundary condition (if it wasn't 0)
const zmax = 106e5  # Top altitude (cm)
const alt = convert(Array, (zmin:dz:zmax)) # These are the layer centers.
const n_all_layers = length(alt)
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(non_bdy_layers) # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])
const hygropause_alt = 40e5  # Location of the hygropause

#                              Temperature profile construction                      
# =======================================================================================================

# Establish the options for controltemps[3]
const Texo_opts = Dict("Mars"=>Dict("min-P2"=>190., "mean-P2"=>210., "max-P2"=>280.,   # These are based on solar min, mean, max.
                                    "min"=>175., "mean"=>225., "max"=>275.,   # These are based on solar min, mean, max.
                                    "meansundist"=>225., "aphelion"=>225., "perihelion"=>225.),
                       "Venus"=>Dict("min"=>260., "mean"=>290., "max"=>320.))

const Texo_inclusive_opts = Dict("inclusive-ap"=>175., 
                                 "inclusive-mean"=>225., 
                                 "inclusive-peri"=>275.)

const Tsurf = Dict("Mars"=>230., "Venus"=>735.)
const Tmeso = Dict("Mars"=>130., "Venus"=>170.)

# Create the temperature profile control array
const controltemps = [Tsurf[planet], Tmeso[planet], Texo_opts[planet]["mean"]]
if planet=="Venus"
    const meantemps = [Tsurf[planet], Tmeso[planet], Texo_opts[planet]["min"]] # Used for saturation vapor pressure. DON'T CHANGE!
elseif planet=="Mars"
    const meantemps = [Tsurf[planet], Tmeso[planet], Texo_opts[planet]["mean"]] # Used for saturation vapor pressure. DON'T CHANGE!
end


# Modify the settings if doing a special isothermal atmosphere.
if temp_scenario=="isothermal"
    controltemps .= [225., 225., 225.]
    meantemps .= [225., 225., 225.] # Used for saturation vapor pressure. DON'T CHANGE!
else # Set the exobase temp according to the temp scenario.
    controltemps[3] =  Texo_opts[planet][temp_scenario]
end

# Modify the array for the special case where multiple parameters are changed for the seasonal model
if special_seasonal_case!=nothing
    controltemps .= [Tsurf[planet], Tmeso[planet], Texo_inclusive_opts[special_seasonal_case]]
end

# MULTICOL ADDED: allocate 2D arrays for Tn, Ti, Te

# Initialize the mean temperature profile for SVP, remains 1-D
local Tn_meanSVP_temp
if planet == "Mars"
    Tn_meanSVP_temp = T_Mars(meantemps...; alt)["neutrals"]
elseif planet=="Venus"
    Tn_meanSVP_temp = T_Venus(meantemps..., "Venus-Inputs/FoxandSung2001_temps_mike.txt"; alt)["neutrals"]
end
const Tn_meanSVP = Tn_meanSVP_temp  # This stays as 1-D

# Initialize the 2-D temperature arrays [n_horiz, num_layers+2]
Tn_temp = zeros(n_horiz, num_layers+2)
Ti_temp = zeros(n_horiz, num_layers+2)
Te_temp = zeros(n_horiz, num_layers+2)

# Loop over horizontal columns to create temperature profiles for each column

# IDENTICAL temperatures
for ihoriz in 1:n_horiz
    # For each column, call T_Mars or T_Venus as if we had a single-column scenario:
    local T_array_dict
    if planet == "Mars"
        T_array_dict = T_Mars(controltemps[1], controltemps[2], controltemps[3]; alt=alt)
    elseif planet == "Venus"
        T_array_dict = T_Venus(controltemps[1], controltemps[2], controltemps[3], 
                               "Venus-Inputs/FoxandSung2001_temps_mike.txt"; alt=alt)
    end

    # Assign temperature profiles to the respective horizontal column (ihoriz, :)
    Tn_temp[ihoriz, :] = T_array_dict["neutrals"]
    Ti_temp[ihoriz, :] = T_array_dict["ions"]
    Te_temp[ihoriz, :] = T_array_dict["electrons"]
end

# DIFFERENT temperatures
# for ihoriz in 1:n_horiz
#     local T_array_dict
#     if planet == "Mars"
#         T_array_dict = T_Mars(controltemps[1], controltemps[2], controltemps[3]; alt=alt)
#     elseif planet == "Venus"
#         T_array_dict = T_Venus(controltemps[1], controltemps[2], controltemps[3], 
#                                "Venus-Inputs/FoxandSung2001_temps_mike.txt"; alt=alt)
#     end

#     # Assign ion and electron temperature profiles (same for all columns)
#     Ti_temp[ihoriz, :] = T_array_dict["ions"]
#     Te_temp[ihoriz, :] = T_array_dict["electrons"]

#     # Assign neutral temperature profiles with column-specific increments
#     if ihoriz == 1
#         # First column unchanged
#         Tn_temp[ihoriz, :] = T_array_dict["neutrals"]
#     elseif ihoriz == 2
#         # Second column increased by 20 K
#         Tn_temp[ihoriz, :] = T_array_dict["neutrals"] .+ 20.0
#     elseif ihoriz == 3
#         # Third column increased by 40 K
#         Tn_temp[ihoriz, :] = T_array_dict["neutrals"] .+ 40.0
#     end
# end

# Final assignment to global constants (2-D arrays)
const Tn_arr = Tn_temp
const Ti_arr = Ti_temp
const Te_arr = Te_temp

# @show size(Tn_arr), size(Ti_arr), size(Te_arr)

const Tplasma_arr = Ti_arr .+ Te_arr;
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)

#                              Horizontal winds construction                      
# =======================================================================================================
# const horiz_wind_v = [zeros(9) for ihoriz in 1:n_horiz]   # horizontal wind profiles; one altitude profile for each vertical column; units cm/s # MULTICOL WARNING hardcoded values

# Provide a simple default horizontal wind profile for each column.  These
# values may be overwritten by user supplied profiles in `INPUT_PARAMETERS.jl`.
# A small constant wind of 10 cm/s is used at all altitudes so that
# horizontal transport terms are non-zero by default.

# const horiz_wind_v = [fill(0.0, length(alt)) for ihoriz in 1:n_horiz]
const horiz_wind_v = [fill(10.0, length(alt)) for ihoriz in 1:n_horiz]

#                                      Water profile settings
# =======================================================================================================

if dust_storm_on==true
    const excess_peak_alt = 60
end

# Water mixing ratios to use
if planet=="Mars"
    const water_MRs = Dict("loweratmo"=>Dict("standard"=>1.3e-4, "low"=>0.65e-4, "high"=>2.6e-4), 
                       "mesosphere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4), 
                       "everywhere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4))
    const water_mixing_ratio = water_MRs[water_loc][water_case]
elseif planet=="Venus"
    const water_mixing_ratio = Dict("standard"=>1e-6)[water_case]  # parse(Float64, water_case) 
end

# Whether to install a whole new water profile or just use the initial guess with modifications (for seasonal model)
if planet=="Venus"
    const reinitialize_water_profile = false 
elseif planet=="Mars"
    const reinitialize_water_profile = seasonal_cycle==true ? false : true # should be off if trying to run simulations for seasons
end

const update_water_profile = seasonal_cycle==true ? true : false # this is for modifying the profile during cycling, MAY be fixed?
const modified_water_alts = "below fixed point"

# altitude at which to add the extra water -- applies to both dust storm parcels and the tanh profile
const add_water_alt_opts = Dict("low"=>45, "standard"=>60, "high"=>65)
const f_fac_opts = Dict("low"=>0.005, "standard"=>10, "high"=>100) # a parameter named f which helps manipulate the water profile. Not related to any other f.

# Set the saturation vapor pressure curves, used in boundary conditions
const H2Osat = map(x->Psat(x), Tn_meanSVP) # Using this function keeps SVP fixed 
const HDOsat = map(x->Psat_HDO(x), Tn_meanSVP)

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position within the active species vector - these are used later in chemJ_mat.
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

# Altitude at which water transitions from fixed to freely solved for
# Calculate the saturation vapor pressure fraction for each horizontal column
initial_atm = get_ncurrent(initial_atm_file, n_horiz)
H2Osatfrac = hcat([H2Osat ./ map(z -> n_tot(initial_atm, z, ihoriz; all_species, n_alt_index), alt) for ihoriz in 1:n_horiz]...)
const upper_lower_bdy = alt[something(findfirst(isequal(minimum(H2Osatfrac[:, 1])), H2Osatfrac[:, 1]), 0)] # in cm
const upper_lower_bdy_i = n_alt_index[upper_lower_bdy]  # the uppermost layer at which water will be fixed, in cm
# Control whether the removal of rates etc at "Fixed altitudes" runs. If the boundary is 
# the bottom of the atmosphere, we shouldn't do it at all.
const remove_rates_flag = true
if upper_lower_bdy == zmin
    const remove_rates_flag = false 
end

#                              Species-specific scale heights
# =======================================================================================================
const Hs_dict = Dict{Symbol, Vector{Vector{Float64}}}([sp => [scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)][ihoriz, :]; molmass, M_P, R_P) for ihoriz in 1:n_horiz] for sp in all_species])
# @show typeof(Hs_dict[:O][1]), size(Hs_dict[:O][1])

#                                     Boundary conditions (lower and upper)
# =======================================================================================================
# "n": density boundary condition; "f": flux bc; "v": velocity bc; 
# "see boundaryconditions()" -- nonthermal escape depends on the dynamic density of the
# atmosphere, so it can't be imposed as a constant here and is calculated on the fly.
if planet=="Mars"
    co2_lower = fill(2.1e17, n_horiz)
    ar_lower  = fill(2.0e-2*2.1e17, n_horiz)
    n2_lower  = fill(1.9e-2*2.1e17, n_horiz)
    h2o_lower = fill(H2Osat[1], n_horiz)
    hdo_lower = fill(HDOsat[1], n_horiz)
    const speciesbclist = Dict(
                        :CO2=>Dict("n"=>[[co2_lower[i], NaN] for i in 1:n_horiz], "f"=>[[NaN, 0.] for _ in 1:n_horiz]),
                        :Ar=>Dict("n"=>[[ar_lower[i], NaN] for i in 1:n_horiz], "f"=>[[NaN, 0.] for _ in 1:n_horiz]),
                        :N2=>Dict("n"=>[[n2_lower[i], NaN] for i in 1:n_horiz], "f"=>[[NaN, 0.] for _ in 1:n_horiz]),
                        #:C=>Dict("f"=>[NaN, 4e5]), # NEW: Based on Lo 2021
                        :H2O=>Dict("n"=>[[h2o_lower[i], NaN] for i in 1:n_horiz], "f"=>[[NaN, 0.] for _ in 1:n_horiz]),
                        :HDO=>Dict("n"=>[[hdo_lower[i], NaN] for i in 1:n_horiz], "f"=>[[NaN, 0.] for _ in 1:n_horiz]),
                        :O=> Dict("f"=>[[0., 1.2e8] for ihoriz in 1:n_horiz]),
                        :H2=>Dict("f"=>[[0., NaN] for ihoriz in 1:n_horiz], "v"=>[[NaN, effusion_velocity(Tn_arr[ihoriz, end], 2.0; M_P, R_P, zmax)] for ihoriz in 1:n_horiz], "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),
                        :HD=>Dict("f"=>[[0., NaN] for ihoriz in 1:n_horiz], "v"=>[[NaN, effusion_velocity(Tn_arr[ihoriz, end], 3.0; M_P, R_P, zmax)] for ihoriz in 1:n_horiz], "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),
                        :H=> Dict("f"=>[[0., NaN] for ihoriz in 1:n_horiz], "v"=>[[NaN, effusion_velocity(Tn_arr[ihoriz, end], 1.0; M_P, R_P, zmax)] for ihoriz in 1:n_horiz], "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),
                        :D=> Dict("f"=>[[0., NaN] for ihoriz in 1:n_horiz], "v"=>[[NaN, effusion_velocity(Tn_arr[ihoriz, end], 2.0; M_P, R_P, zmax)] for ihoriz in 1:n_horiz], "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),
                       );
elseif planet=="Venus"
    const ntot_at_lowerbdy = 9.5e15 # at 90 km
    # const KoverH_lowerbdy = Keddy([zmin], [ntot_at_lowerbdy]; planet)[1]/scaleH_lowerboundary(zmin, Tn_arr[1]; molmass, M_P, R_P, zmin)
    const KoverH_lowerbdy = Keddy([zmin], [ntot_at_lowerbdy]; planet=planet)[1] / scaleH_lowerboundary(zmin, Tn_arr[1]; molmass, M_P, R_P, zmin)
    const manual_speciesbclist=Dict(# major species neutrals at lower boundary (estimated from Fox&Sung 2001, Hedin+1985, agrees pretty well with VIRA)
                                    :CO2=>Dict("n"=>[[0.965*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    :Ar=>Dict("n"=>[[5e11, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    :CO=>Dict("n"=>[[4.5e-6*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    :O2=>Dict("n"=>[[3e-3*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    # :O2 => Dict("n" => [[3e-3 * ntot_at_lowerbdy, NaN], [2.9e-3 * ntot_at_lowerbdy, NaN], [3.1e-3 * ntot_at_lowerbdy, NaN]], "f" => [[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    :N2=>Dict("n"=>[[0.032*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz]),

                                    # water mixing ratio is fixed at lower boundary
                                    :H2O=>Dict("n"=>[[water_mixing_ratio*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),
                                    # we assume HDO has the bulk atmosphere ratio with H2O at the lower boundary, ~consistent with Bertaux+2007 observations
                                    :HDO=>Dict("n"=>[[2*DH*water_mixing_ratio*ntot_at_lowerbdy, NaN] for ihoriz in 1:n_horiz], "f"=>[[NaN, 0.] for ihoriz in 1:n_horiz]),

                                    # atomic H and D escape solely by photochemical loss to space, can also be mixed downward
                                    :H=> Dict("v"=>[[-KoverH_lowerbdy, effusion_velocity(Tn_arr[ihoriz, end], 1.0; zmax, M_P, R_P)] for ihoriz in 1:n_horiz],
                                                    #                 ^^^ other options here:
                                                    #                 effusion_velocity(Tn_arr[ihoriz, end], 1.0; zmax) # thermal escape, negligible
                                                    #                 100 # representing D transport to nightside, NOT escape
                                                    #                 NaN # No thermal escape to space, appropriate for global average model
                                              "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),
                                    :D=> Dict("v"=>[[-KoverH_lowerbdy, effusion_velocity(Tn_arr[ihoriz, end], 2.0; zmax, M_P, R_P)] for ihoriz in 1:n_horiz],
                                                    #                 ^^^ other options here:
                                                    #                  effusion_velocity(Tn_arr[ihoriz, end], 2.0; zmax) # thermal escape, negligible
                                                    #                 100 # representing D transport to nightside, NOT escape
                                                    #                 NaN # No thermal escape to space, appropriate for global average model
                                              "ntf"=>[[NaN, "see boundaryconditions()"] for ihoriz in 1:n_horiz]),

                                    # # H2 mixing ratio at lower boundary adopted from Yung&DeMore1982 as in Fox&Sung2001
                                    # :H2=>Dict("n"=>[1e-7*ntot_at_lowerbdy, NaN],
                                    #           "v"=>[NaN, effusion_velocity(Tn_arr[ihoriz, end], 2.0; zmax)],
                                    #           "ntf"=>[NaN, "see boundaryconditions()"]),
                                    # :HD=>Dict("n"=>[DH*1e-7*ntot_at_lowerbdy, NaN],
                                    #           "v"=>[NaN, effusion_velocity(Tn_arr[ihoriz, end], 3.0; zmax)],
                                    #           "ntf"=>[NaN, "see boundaryconditions()"]),

                                    # unusued neutral boundary conditions
                                    #:O=> Dict("v"=>[-KoverH_lowerbdy, NaN], "f"=>[NaN, 0.#=1.2e6=#]), # no effect on O profile
                                    #:N=>Dict("v"=>[-KoverH_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    #:NO=>Dict("v"=>[-KoverH_lowerbdy, NaN], #="n"=>[3e8, NaN],=# #="n"=>[5.5e-9*ntot_at_lowerbdy, NaN], =# "f"=>[NaN, 0.]),

                                    # assume no ion loss, appropriate for global average and small observed rates
                                    #:Hpl=>Dict("v"=>[-KoverH_lowerbdy, 0.0 #=effusion_velocity(Ti_arr[ihoriz, end], 1.0; zmax)=#]),#, "f"=>[NaN, 1.6e7]),
                                    #:H2pl=>Dict("v"=>[-KoverH_lowerbdy, 0.0 #=effusion_velocity(Ti_arr[ihoriz, end], 2.0; zmax)=#]),#, "f"=>[NaN, 2e5]),
                                    #:Opl=>Dict("v"=>[-KoverH_lowerbdy, 2e5], ), # "f"=>[NaN, 2.1e8] # tends to cause hollowing out of atmosphere
                                    );

    # add in downward mixing velocity boundary condition for all other species
    auto_speciesbclist = Dict()
    for sp in all_species
        if sp in keys(manual_speciesbclist)
            auto_speciesbclist[sp] = manual_speciesbclist[sp]
        else
            auto_speciesbclist[sp] = Dict("v"=>[[-KoverH_lowerbdy, 0.0] for ihoriz in 1:n_horiz])
        end
    end

    const speciesbclist = deepcopy(auto_speciesbclist)
end

#                                     Boundary conditions (back edge and front edge)
# =======================================================================================================
# "n": density boundary condition; "f": flux bc; "v": velocity bc; 
# The default boundary conditions are zero flux boundary conditions at the back edge and the front edge. If different boundary conditions are required, they will need to be implemented here and in the boundaryconditions_horiz function.
# The zero flux boundary conditions will be input as vectors with altitude.
# For each species, there are two vectors of length num_layers. The first directs the BC values at the first vertical column and the second directs the BC values at the last vertical column

# The dictionary `speciesbclist_horiz` sets horizontal fluxes at the back and
# front edges.  By default both profiles are zero, representing closed
# boundaries, but you can modify the values below (or in a separate script)
# to impose an influx or outflux for any species.  Each entry contains two
# vectors of length `num_layers` giving the flux [#/cm²/s] at the back and front
# edges respectively.

# add in zero flux edge boundary conditions on both edges for all species
auto_speciesbclist_horiz = Dict()
for sp in all_species
    auto_speciesbclist_horiz[sp] = Dict("f"=>[[0.0 for ialt in 1:num_layers] for c in 1:2])
end

const speciesbclist_horiz = deepcopy(auto_speciesbclist_horiz)

# Example modification: Set non-zero flux boundary conditions for O
# speciesbclist_horiz[:O] = Dict(
#     "f" => [
#         fill(1e7, num_layers),   # Influx at the back edge (cm⁻² s⁻¹)
#         fill(-1e7, num_layers)   # Outflux at the front edge (cm⁻² s⁻¹)
#     ]
# )

# ***************************************************************************************************** #
#                                                                                                       #
#                         Set up simulation filenames and define input files                            #
#                                                                                                       #
# ***************************************************************************************************** #

# Crosssection file sources
# -------------------------------------------------------------------
const photochem_data_files = Dict(:CO2=>Dict("main"=>"CO2.dat"), 
                                   :H2O=>Dict("main"=>"h2oavgtbl.dat"), 
                                   :HDO=>Dict("main"=>"HDO.dat"), 
                                   :H2O2=>Dict("main"=>"H2O2.dat"), 
                                   :HDO2=>Dict("main"=>"H2O2.dat"), 
                                   :O3=>Dict("main"=>"O3.dat", "chapman"=>"O3Chap.dat"), 
                                   :O2=>Dict("main"=>"O2.dat", "schr_short"=>"130-190.cf4", "schr_mid"=>"190-280.cf4", "schr_long"=>"280-500.cf4"), 
                                   :H2=>Dict("main"=>"binnedH2.csv"), 
                                   :HD=>Dict("main"=>"binnedH2.csv"), 
                                   :OH=>Dict("main"=>"binnedOH.csv", "O1D+H"=>"binnedOHo1D.csv"), 
                                   :OD=>Dict("main"=>"OD.csv"))

# Filename tags and codes
# -------------------------------------------------------------------

extra_str = seasonal_cycle==true ? "cycle" : "eq"

if seasonal_cycle == true
    # folder naming scheme
    filetag = Dict("temperature"=>"seasons_temp$(extra_str)_Texo=$(Int64(controltemps[3]))_$(results_version)",
                   "water"=>"seasons_water$(extra_str)_$(water_case)_$(water_loc)_$(results_version)",
                   "insolation"=>"seasons_insolation_$(solar_scenario)_$(results_version)",)

    if special_seasonal_case != nothing
        const tag = "seasons_$(special_seasonal_case)_$(results_version)"
    else 
        const tag = filetag[exp_type]
    end
else # not a seasonal cycle experiment
    const tag = "eqrun_$(exp_type)_$(results_version)"
end

# Tags, shortcodes, and filenames
# -------------------------------------------------------------------
# The shortcodes provide unique identifiers for a simulation. Necessary because you end up running the model many times...
const hrshortcode, rshortcode = generate_code(ions_included, controltemps[1], controltemps[2], controltemps[3], water_case, solar_scenario)
const sim_folder_name = "$(hrshortcode)_$(rshortcode)_$(tag)"


# ***************************************************************************************************** #
#                                                                                                       #
#                           Algorithm, solver, and float type settings                                  #
#                                                                                                       #
# ***************************************************************************************************** #

# Simulation run time and timestep size  
const season_length_in_sec = seasonal_cycle==true ? season_in_sec : 1e16
const maxlogdt = seasonal_cycle==true ? 5 : 16 # simulation will run until dt = 10^maxlogdt seconds
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-3, maxlogdt])
const timestep_type = seasonal_cycle==true ? "log-linear" : "dynamic-log" 
    # OPTIONS: "static-log": Logarithmically spaced timesteps that are imposed and don't adjust.
    #                        Should basically never be used, but can be used for testing.
    # "dynamic-log": Logarithmically spaced timesteps which can be dynamically adjusted
    #                as the model runs to improve stability.
    # "log-linear": A hybrid timestep that uses logarithmically-spaced timesteps 
    #               at small elapsed t to get the model going, but then switches
    #               to linearly spaced timesteps of ~1 week. This is used exclusively
    #               with the seasonal model, so that the output can be saved 
    #               at every week and end specifically at 1 season.
    

#                                        Solver algorithm type 
# =======================================================================================================
const problem_type = "Gear" 
    # OPTIONS: 
    # "SS": Julia solver SteadyState.
    # "ODE": Julia ODE solver.
    # "Gear": Preferred solver; Mike's hand-coded Gear method.
# for static timesteps:
const n_steps = 800 # Used with timestep_type="static-log"
# for dynamic timesteps:
const dt_incr_factor = 1.5
const dt_decr_factor = 10

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be
if problem_type == "Gear"
    const assume_photochem_eq = false
else # In using the Julia-provided solvers, it was necessary to assume photochemical equilibrium for short-lived species.
    const assume_photochem_eq = converge_which == "both" ? true : false
end

#                Nitty gritty additions that were included to improve model stability
# =======================================================================================================
# whether to include differentiation terms in Jacobian with respect to electron density 
# or generic thirdbody M. After much testing, these were determined to not be necessary.
const ediff = false # true 
const mdiff = false # true 
const error_checking_scheme = "new"
    # OPTIONS: "new", "old" 
    # We had to work on the method of checking error a lot when dealing with stability issues.
    # The old version has been left in for testing purposes just in case problems are ever 
    # encountered again.

#                                            Float types
# =======================================================================================================
# See CONSTANTS.jl for the setting. no, it's not ideal, but it's the best option 
# right now.
# this means this file must be loaded after CONSTANTS.

# Logic to require Double64 when using the Gear solver. Currently off as Doubles 
# are not being used.
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in MODEL_SETUP.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end

# ***************************************************************************************************** #
#                                                                                                       #
#                         Misc. things that depend on things defined above                              #
#                                                                                                       #
# ***************************************************************************************************** #

# ***************************************************************************************************** #
#                                                                                                       #
#                          Create a parameter dataframe for logging ease                                #
#                                                                                                       #
# ***************************************************************************************************** #

PARAMETERS_GEN = DataFrame(Field=[], Value=[])

push!(PARAMETERS_GEN, ("PLANET", planet));
push!(PARAMETERS_GEN, ("M_P", M_P));
push!(PARAMETERS_GEN, ("R_P", R_P));
push!(PARAMETERS_GEN, ("HRSHORTCODE", hrshortcode));
push!(PARAMETERS_GEN, ("RSHORTCODE", rshortcode));
push!(PARAMETERS_GEN, ("VARIED_PARAM", exp_type))
push!(PARAMETERS_GEN, ("INITIAL_ATM", initial_atm_file));
push!(PARAMETERS_GEN, ("RXN_SOURCE", reaction_network_spreadsheet));
push!(PARAMETERS_GEN, ("IONS", ions_included ));
push!(PARAMETERS_GEN, ("CONVERGE", converge_which));
push!(PARAMETERS_GEN, ("NONTHERMAL_ESC", nontherm));
push!(PARAMETERS_GEN, ("SOLAR_SCENARIO", solar_scenario));
push!(PARAMETERS_GEN, ("SOLARFILE", solarfile));
push!(PARAMETERS_GEN, ("ELECTRON_PROF", e_profile_type));
push!(PARAMETERS_GEN, ("EDIFF", ediff));
push!(PARAMETERS_GEN, ("MDIFF", mdiff));
push!(PARAMETERS_GEN, ("DH", DH));
push!(PARAMETERS_GEN, ("AMBIPOLAR_DIFFUSION_ON", use_ambipolar));
push!(PARAMETERS_GEN, ("MOLEC_DIFFUSION_ON", use_molec_diff));

# Log altitude grid so we can avoid loading this very file when doing later analysis.
PARAMETERS_ALTGRID = DataFrame(Alt=alt)

PARAMETERS_CONDITIONS = DataFrame(Field=[], Value=[], Unit=[]);

push!(PARAMETERS_CONDITIONS, ("SZA", SZA, "deg"));
push!(PARAMETERS_CONDITIONS, ("TSURF", controltemps[1], "K"));
push!(PARAMETERS_CONDITIONS, ("TMESO", controltemps[2], "K"));
push!(PARAMETERS_CONDITIONS, ("TEXO", controltemps[3], "K"));
push!(PARAMETERS_CONDITIONS, ("MEAN_TEMPS", join(meantemps, " "), "K"));
push!(PARAMETERS_CONDITIONS, ("WATER_MR", water_mixing_ratio, "mixing ratio"));
push!(PARAMETERS_CONDITIONS, ("WATER_CASE", water_case, "whether running with 10x, 1/10th, or standard water in middle/upper atmo"));

waterbdy = :H2O in inactive_species ? zmax : upper_lower_bdy/1e5
push!(PARAMETERS_CONDITIONS, ("WATER_BDY", waterbdy, "km"))

# This is so ugly because the XLSX package won't write columns of different lengths, so I have to pad all the shorter lists
# with blanks up to the length of the longest list and also transform all the symbols into strings. 
L = max(length(all_species), length(neutral_species), length(ion_species), length(no_chem_species), length(no_transport_species), length(Jratelist))
PARAMETERS_SPLISTS = DataFrame(AllSpecies=[[string(a) for a in all_species]..., ["" for i in 1:L-length(all_species)]...], 
                               Neutrals=[[string(n) for n in neutral_species]..., ["" for i in 1:L-length(neutral_species)]...], 
                               Ions=[[string(i) for i in ion_species]..., ["" for i in 1:L-length(ion_species)]...],
                               NoChem=[[string(nc) for nc in no_chem_species]..., ["" for i in 1:L-length(no_chem_species)]...],
                               NoTransport=[[string(nt) for nt in no_transport_species]..., ["" for i in 1:L-length(no_transport_species)]...],
                               Jratelist=[[string(j) for j in Jratelist]..., ["" for i in 1:L-length(Jratelist)]...]);
PARAMETERS_SOLVER = DataFrame(Field=[], Value=[]);
PARAMETERS_XSECTS = DataFrame(Species=[], Description=[], Filename=[]);
# PARAMETERS_BCS = DataFrame(Species=[], Type=[], Lower=[], Upper=[]);
# PARAMETERS_BCS_HORIZ = DataFrame(Species=[], Type=[], BackEdge=[], FrontEdge=[]);
PARAMETERS_BCS = DataFrame(Species=[], Type=[], Column=Int[], Lower=[], Upper=[]);
PARAMETERS_BCS_HORIZ = DataFrame(Species=[], Type=[], Altitude=Int[], BackEdge=[], FrontEdge=[]);

# LOG THE TEMPERATURES
PARAMETERS_TEMPERATURE_ARRAYS = DataFrame(Neutrals = vec(Tn_arr), Ions = vec(Ti_arr), Electrons = vec(Te_arr))
