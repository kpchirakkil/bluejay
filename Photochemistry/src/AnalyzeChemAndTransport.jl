# **************************************************************************** #
#                                                                              #
#          Non-core functions for analyzing the chemistry and transport        #
#                                                                              #
# These functions are generally used after the model has run and you want to   #
# Do some analysis on the output. I don't believe any of them are called       #
# inside of the core of the model (Core.jl). If so, they should be in there.   #
# If such functions are found, please open an issue on Github.                 #
# **************************************************************************** #

#===============================================================================#
#                           Chemistry functions                                 #
#===============================================================================#

function chemical_lifetime(s::Symbol, atmdict; globvars...)
    #=
    Calculates chemical lifetime of a molecule s in the atmosphere atmdict. 
    Good for comparing with the results of diffusion_timescale.
    =#
    GV = values(globvars)
    required = [:all_species, :Jratelist, :n_alt_index, :ion_species, :num_layers, :reaction_network, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    # We'll accumulate chemical lifetimes in a matrix, one column per horizontal slice
    chem_lt = zeros(Float64, n_horiz, GV.num_layers)  # (n_horiz, num_layers)

    for ihoriz in 1:n_horiz
        # For the interior (bulk) layers, we slice Tn, Ti, Te
        # Tn_col, Ti_col, Te_col each has length = num_layers (since we skip boundary)
        Tn_col = GV.Tn[ihoriz, :]
        Ti_col = GV.Ti[ihoriz, :]
        Te_col = GV.Te[ihoriz, :]

        loss_all_rxns, ratecoefs = get_volume_rates(
            s, atmdict, n_horiz;
            species_role="reactant", 
            which="all",
            remove_sp_density=true,
            # pass the sliced arrays so code sees only the interior layers:
            Tn=Tn_col[2:end-1], Ti=Ti_col[2:end-1], Te=Te_col[2:end-1],
            globvars...
        )
    
        total_loss_by_alt = zeros(size(Tn_col))  # length = num_layers
    
        for k in keys(loss_all_rxns)
            total_loss_by_alt .+= loss_all_rxns[k][ihoriz]
        end
    
        chem_lt[ihoriz, :] = 1.0 ./ total_loss_by_alt
    end

    return chem_lt
end


function get_column_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, 
                          ihoriz::Int64; which="all", sp2=nothing, role="product", 
                          startalt_i=1, returntype="df", globvars...)
    #=
    Input:
        sp: species for which to search for reactions
        atmdict: the present atmospheric state to calculate on
        ihoriz: horizontal column index to extract correct column data
        Tn, Ti, Te: Arrays of the temperature profiles including boundary layers
        which: whether to do photochemistry, or just bimolecular reactions. "all", "Jrates" or "krates"
        sp2: optional second species to include, i.e. usually sp's ion.
        role: "product" or "reactant" only.
        startalt_i: Index of the altitude at which to start. 
        returntype: whether to return a dataframe ("df") or a dictionary "dict"
    Output:
        sorted: Total column rates for all reactions of species sp. 
                Sorted, in order of largest rate to smallest.
    =#
    GV = values(globvars)
    required = [:Tn, :Ti, :Te, :all_species, :ion_species, :reaction_network, :num_layers, :dz]
    check_requirements(keys(GV), required)

    # Extract FULL column temperatures explicitly
    Tn_col = GV.Tn[ihoriz, :]
    Ti_col = GV.Ti[ihoriz, :]
    Te_col = GV.Te[ihoriz, :]

    # Internal slicing for bulk altitudes
    # Tn_col = Tn_col_full[2:end-1]
    # Ti_col = Ti_col_full[2:end-1]
    # Te_col = Te_col_full[2:end-1]

    # Debug: confirm shapes
    # @assert length(Tn_col) == GV.num_layers + 2 "Tn_col length should match num_layers+2"
    # @assert length(Ti_col) == GV.num_layers + 2 "Ti_col length should match num_layers+2"
    # @assert length(Te_col) == GV.num_layers + 2 "Te_col length should match num_layers+2"

    # Now we compute reaction rates for just this column
    rxd, coefs = get_volume_rates(sp, atmdict, n_horiz;
                                  species_role=role,
                                  which=which,
                                  globvars...,
                                  # pass the column slices for Tn, Ti, Te
                                  Tn=Tn_col[2:end-1], Ti=Ti_col[2:end-1], Te=Te_col[2:end-1])

    # Sum up the column rates for each reaction
    columnrate = Dict{String, Float64}()

    for k in keys(rxd)
        # skip altitudes below startalt_i if requested
        columnrate[k] = sum(rxd[k][startalt_i:end] .* GV.dz)
    end

    # If second species is requested
    if sp2 != nothing
        rxd2, _ = get_volume_rates(sp2, atmdict, n_horiz;
                                   species_role=role,
                                   which=which,
                                   globvars...,
                                   Tn=Tn_col[2:end-1], Ti=Ti_col[2:end-1], Te=Te_col[2:end-1])
        columnrate2 = Dict{String, Float64}()

        for k in keys(rxd2)
            columnrate2[k] = sum(rxd2[k][startalt_i:end] .* GV.dz)
        end
        colrate_dict = merge(columnrate, columnrate2)
    else
        colrate_dict = columnrate
    end

    # sort them from largest to smallest
    sorted = sort(collect(colrate_dict), by=x->x[2], rev=true)

    if returntype=="df"
        return DataFrame([[names(DataFrame(sorted))]; collect.(eachrow(DataFrame(sorted)))],
                         [:Reaction; :ColumnRate])
    else
        return sorted
    end
end


function get_volume_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, 
                          n_horiz::Int64; species_role="both", which="all", 
                          remove_sp_density=false, globvars...)

    #=
    Input:
        sp: Species name
        atmdict: Present atmospheric state dictionary
        n_horiz: Number of vertical columns in the simulation
        Tn, Ti, Te: 2D temperature arrays each of shape (n_horiz, num_layers+2)
        species_role: whether to look for sp as a "reactant", "product", or "both"
        which: "all", "Jrates", "krates"
        remove_sp_density: if set to true, remove sp from the reactants
    Output:
        rxn_dat: dictionary reaction => Vector{Array{ftype_ncur}} 
                 where each entry is the alt-dependent rate for each column
        rate_coefs: dictionary of reaction => rate coefficient
    =#
    GV = values(globvars)
    required = [:all_species, :ion_species, :num_layers, :reaction_network, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    rxn_dat = Dict{String, Vector{Array{ftype_ncur}}}()
    rate_coefs = Dict{String, Vector{Array{ftype_ncur}}}()

    # filter the reaction network to only those relevant
    filtered_rxn_list = filter_network(sp, which, species_role; GV.reaction_network)

    # Now we loop over each horizontal column
    for ihoriz in 1:n_horiz
        # Extract FULL temperature arrays for this column (length=num_layers+2)
        Tn_col = GV.Tn[ihoriz, :]
        Ti_col = GV.Ti[ihoriz, :]
        Te_col = GV.Te[ihoriz, :]

        # Internally slice for bulk altitudes (length=num_layers)
        # Tn_col = Tn_col_full[2:end-1]
        # Ti_col = Ti_col_full[2:end-1]
        # Te_col = Te_col_full[2:end-1]

        # println("[DEBUG] Inside loop for ihoriz = ", ihoriz)
        # println("length(Tn_col) = ", length(Tn_col))
        # println("length(Ti_col) = ", length(Ti_col))
        # println("length(Te_col) = ", length(Te_col))
        # @assert length(Tn_col) == GV.num_layers + 2 "Tn_col length mismatch"

        # For each reaction in the filtered list:
        for rxn in filtered_rxn_list
            rxn_str = format_chemistry_string(rxn[1], rxn[2])

            if typeof(rxn[3]) == Symbol
                # photodissociation
                if remove_sp_density == false
                    # multiply [Molecule] * Jrate
                    # store alt profiles for each column
                    rxn_dat[rxn_str] = [atmdict[rxn[1][1]][ihoriz] .* vec(atmdict[rxn[3]][ihoriz]) for ihoriz in 1:n_horiz]
                else
                    rxn_dat[rxn_str] = 1 .* atmdict[rxn[3]]
                end
                rate_coefs[rxn_str] = vec(atmdict[rxn[3]])
            else
                # bi- or termolecular
                remove_me = (remove_sp_density == true) ? sp : nothing
                density_prod = [reactant_density_product(atmdict, rxn[1], ihoriz; removed_sp=remove_me, globvars...)
                                for ihoriz in 1:n_horiz]

                thisrate = (typeof(rxn[3]) != Expr) ? :($rxn[3] + 0) : rxn[3]
                rate_coef = [eval_rate_coef(atmdict, thisrate, ihoriz; globvars...)
                             for ihoriz in 1:n_horiz]

                rxn_dat[rxn_str] = [density_prod[ihoriz] .* rate_coef[ihoriz] for ihoriz in 1:n_horiz]

                if typeof(rate_coef) == Vector{Float64} && typeof(rate_coef[1]) == Float64
                    rate_coef = [rate_coef[ihoriz] .* ones(GV.num_layers) for ihoriz in 1:n_horiz]
                end
                rate_coefs[rxn_str] = rate_coef
            end
        end
    end

    return rxn_dat, rate_coefs
end


function get_volume_rates(sp::Symbol, source_rxn::Vector{Any}, source_rxn_rc_func, 
                          atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, 
                          Mtot, ihoriz::Int64; globvars...)
    #=
    Override to call for a single reaction. Useful for doing non-thermal flux boundary conditions.
    Input:
        sp: species name
        source_rxn: chemical reaction for which to get the volume rate 
        atmdict: Present atmospheric state dictionary
        Mtot: total atmospheric density at each altitude
        ihoriz: vertical column index
    Output: 
        vol_rates: alt-dependent rates (#/cm^3/s) for the reaction
    =#
    GV = values(globvars)
    required = [:all_species, :ion_species, :Jratedict, :num_layers, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    # println("\n[DEBUG - second (override) get_volume_rates called]")
    # println("[DEBUG] Species = ", sp, ", ihoriz = ", ihoriz)
    # @show size(GV.Tn), GV.num_layers, ihoriz

    # Now Tn_col, Ti_col, Te_col each is length num_layers+2
    Tn_col = GV.Tn[ihoriz, :]
    Ti_col = GV.Ti[ihoriz, :]
    Te_col = GV.Te[ihoriz, :]

    # println("length(Tn_col) at column $(ihoriz) = ", length(Tn_col))
    # @assert length(Tn_col) == GV.num_layers + 2 "Tn_col length mismatch"

    # Fill in the rate * density
    if typeof(source_rxn[3]) == Symbol
        # photodissociation
        # multiply [molecule] * Jrate
        vol_rates = atmdict[source_rxn[1][1]][ihoriz] .* GV.Jratedict[source_rxn[3]][ihoriz]
    else
        # bi- or termolecular reaction
        # evaluate the rate coef for just the interior alt range
        rate_coef = source_rxn_rc_func(Tn_col[2:end-1], Ti_col[2:end-1], Te_col[2:end-1], Mtot)
        vol_rates = reactant_density_product(atmdict, source_rxn[1], ihoriz; globvars...) .* rate_coef
    end

    return vol_rates
end

function make_chemjac_key(fn, fpath, list1, list2) 
    #=
    This somewhat superfluous function makes a key to the chemical jacobian,
    telling which index corresponds to which species. But really it just gives the 
    indices of the entries in all_species, because that's how the jacobian is ordered,
    but this function is written agnostically so that could technically change and this
    function would still work.

    fn: filename to save the key to
    fpath: where to save fn
    list1: jacobian row indices
    list2: jacobian col indices
    =#
    dircontents = readdir(fpath)
    if !(fn in dircontents)
        println("Creating the chemical jacobian row/column key")
        f = open(fpath*"/"*fn, "w")
        write(f, "Chemical jacobian rows (i):\n")
        write(f, "$([i for i in 1:length(list1)])\n")
        write(f, "$(list1)\n\n")
        write(f, "Chemical jacobian cols (j):\n")
        write(f, "$([j for j in 1:length(list2)])\n")
        write(f, "$(list2)\n\n")
        close(f)
    end
end

function reactant_density_product(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, reactants, ihoriz::Int64; removed_sp=nothing, globvars...)
    #=
    Calculates the product of all reactant densities for a chemical reaction for the whole atmosphere, 
    i.e. for A + B --> C + D, return n_A * n_B.

    Input:
        atmdict: the atmospheric state dictionary
        reactants: a list of reactant symbols.
	ihoriz: Vertical column index
    Output: 
        density_product: returns n_A * n_B for all altitudes for the reaction A + B for column ihoriz--> ...
    =#

    GV = values(globvars)
    required =  [:all_species, :ion_species, :num_layers]
    check_requirements(keys(GV), required)

    if removed_sp != nothing # remove reactant if requested - useful for calculating chemical lifetimes
        deleteat!(reactants, findfirst(x->x==removed_sp, reactants))
    end

    density_product = ones(GV.num_layers)
    for r in reactants
        if r != :M && r != :E
           # species densities by altitude
           density_product .*= atmdict[r][ihoriz]  # multiply by each reactant density
        elseif r == :M
           density_product .*= sum([atmdict[sp][ihoriz] for sp in GV.all_species])
        elseif r == :E
            density_product .*= sum([atmdict[sp][ihoriz] for sp in GV.ion_species])
        else
            throw("Got an unknown symbol in a reaction rate: $(r)")
        end
    end

    return density_product 
end 

function volume_rate_wrapper(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot, ihoriz::Int64; returntype="array", globvars...)
    #=
    Gets altitude-dependent volume production or loss of species sp due to reactions in source_rxns.
    Does NOT care if it is production or loss. 
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :collision_xsect, :ion_species, :Jratedict, :molmass, :non_bdy_layers, :num_layers,  
                :n_alt_index, :Tn, :Ti, :Te, :dz, :zmax]
    check_requirements(keys(GV), required)

    # each column corresponds to a reaction and each row corresponds to an
    # altitude. This orientation matches how we populate the array below.
    rates = Array{ftype_ncur}(undef, GV.num_layers, length(source_rxns))

    # Extract the temperature for this horizontal column explicitly
    Tn_col = GV.Tn[ihoriz, :]
    Ti_col = GV.Ti[ihoriz, :]
    Te_col = GV.Te[ihoriz, :]

    # Loop through reactions
    for (i, source_rxn) in enumerate(source_rxns)
        rate_coef = source_rxn_rc_funcs[source_rxn](Tn_col[2:end-1], Ti_col[2:end-1], Te_col[2:end-1], Mtot)
        rates[:, i] = reactant_density_product(atmdict, source_rxn[1], ihoriz; globvars...) .* rate_coef
    end

    # Returns an array where rows represent altitudes and columns are reactions.
    if returntype=="by rxn"
        return sum(rates, dims=1)
    elseif returntype=="by alt"
        return sum(rates, dims=2)
    elseif returntype=="array"
        return rates
    elseif returntype=="df" # Useful if you want to look at the arrays yourself.
        ratesdf = DataFrame(rates, vec([format_chemistry_string(r[1], r[2]) for r in source_rxns]))
        return ratesdf
    end
end

#===============================================================================#
#                   Transport and escape functions                              #
#===============================================================================#

function diffusion_timescale(s::Symbol, T_arr::Array, atmdict; globvars...)
    #=
    Inputs:
        s: species symbol
        T_arr: temperature array for the given species 
        atmdict: atmospheric state dict
        Dcoef_template: array of 0s to use in Dcoef!
    Output: Molecular and eddy diffusion timescale (s) by alt
    =#
    
    GV = values(globvars)
    required = [:all_species, :alt, :molmass, :n_alt_index, :neutral_species, :polarizability, :planet, :q, :speciesbclist, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    # Get diffusion coefficient array template
    Dcoef_template = zeros(size(T_arr)) 

    # Other stuff
    Hs = scaleH(GV.alt, s, T_arr; GV.molmass)
    ncur_with_bdys =  ncur_with_boundary_layers(atmdict, n_horiz; GV.all_species, GV.n_alt_index)
    
    # Molecular diffusion timescale: H_s^2 / D, scale height over diffusion constant
    D = Dcoef!(Dcoef_template, T_arr, s, ncur_with_bdys; globvars...)
    molec_or_ambi_timescale = (Hs .^ 2) ./ D
   
    # Eddy timescale... this was in here only as scale H... 
    K = [Keddy(alt, n_tot(ncur_with_bdys, ihoriz; GV.all_species, GV.molmass); GV.planet) for ihoriz in 1:n_horiz]
    #eddy_timescale = (Hs .^ 2) ./ K # original
    eddy_timescale = ([Hs for ihoriz in 1:n_horiz] .^ 2) ./ K

    # Combined timescale?!??
    combined_timescale = ([Hs for ihoriz in 1:n_horiz] .^ 2) ./ (K .+ D)

    return molec_or_ambi_timescale, eddy_timescale, combined_timescale
end

function final_escape(thefolder, thefile, n_horiz::Int64; globvars...)
    #=
    thefolder: Folder in which an atmosphere file lives
    thefile: the file containing an atmosphere for which you'd like to calculate the final escape fluxes of H and D.
    n_horiz: Number of vertical columns in the simulation
    =#
    
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
                                   :q, :molmass, :polarizability, :collision_xsect,
                                   # From CUSTOMIZATIONS.jl
                                   :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, 
                                   # Simulation-unique stuff 
                                    :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)
    
    # First load the atmosphere and associated variables.
    atmdict = get_ncurrent(thefolder*thefile, n_horiz);

    vardict = load_from_paramlog(thefolder; globvars...);
    
    # Get Jrate list 
    Jratelist = format_Jrates(thefolder*"active_rxns.xlsx", GV.all_species, "Jratelist"; hot_atoms=true, ions_on=true)[1];
    Jratedict = Dict([j=>atmdict[j] for j in Jratelist])
    
    # Make a dataframe to store things
    escdf = DataFrame("EscapeType"=>["Thermal", "Nonthermal", "Total"], 
                      "H"=>[0, 0, 0], "D"=>[0, 0, 0], "H2"=>[0, 0, 0], "HD"=>[0, 0, 0])

    # Now collect non-thermal and thermal fluxes for each species. 
    for s in ["H", "D", "H2", "HD"]
        nonthermal_esc, thermal_esc = get_transport_PandL_rate(Symbol(s), atmdict, n_horiz; returnfluxes=true, Jratedict, zmax=GV.alt[end],
                                                               hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, hot_HD_network=GV.hHDnet,
                                                               hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, hot_HD_rc_funcs=GV.hHDrc, 
                                                               Hs_dict=vardict["Hs_dict"], ion_species=vardict["ion_species"], neutral_species=vardict["neutral_species"],
                                                               speciesbclist=vardict["speciesbclist"],
                                                               Tprof_for_Hs=vardict["Tprof_for_Hs"], Tprof_for_diffusion=vardict["Tprof_for_diffusion"], 
                                                               transport_species=vardict["transport_species"], 
                                                               Tn=vardict["Tn_arr"], Ti=vardict["Ti_arr"], Te=vardict["Te_arr"], Tp=vardict["Tplasma_arr"],
                                                               globvars...)

        escdf.:($s) = [thermal_esc, nonthermal_esc, thermal_esc+nonthermal_esc]
    end
    
    # Calculate total H atoms lost
    escdf."TotalHAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :H2, :D])])) .+ 2 .* escdf[:, :H2]
    # Calculate total H atoms lost
    escdf."TotalDAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :H2, :H, :TotalHAtomsLost])])) # adds up 1 * D and 1 * HD
    # Calculate total atoms lost
    escdf."TotalHnDAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :D, :H, :H2, :HD])]))
    
    return escdf
end

function fractionation_factor(esc_df, h2o_0, hdo_0; ftype="total")
    #=
    Calculates fractionation factor if given the H2O and HDO at the bottom of the atmosphere, as well as a dataframe
    of escape rates as generated by final_escape.
    =#
    flux_t_D = df_lookup(esc_df, "EscapeType", "Thermal", "D")[1] + df_lookup(esc_df, "EscapeType", "Thermal", "HD")[1]
    flux_t_H = df_lookup(esc_df, "EscapeType", "Thermal", "H")[1] + 2*df_lookup(esc_df, "EscapeType", "Thermal", "H2")[1] + df_lookup(esc_df, "EscapeType", "Thermal", "HD")[1]

    flux_nt_D = df_lookup(esc_df, "EscapeType", "Nonthermal", "D")[1] + df_lookup(esc_df, "EscapeType", "Nonthermal", "HD")[1]
    flux_nt_H = df_lookup(esc_df, "EscapeType", "Nonthermal", "H")[1] + 2*df_lookup(esc_df, "EscapeType", "Nonthermal", "H2")[1] + df_lookup(esc_df, "EscapeType", "Nonthermal", "HD")[1]

    if ftype=="thermal"
        flux_nt_D = 0
        flux_nt_H = 0
    elseif ftype=="nonthermal"
        flux_t_D = 0
        flux_t_H = 0
    end
    
    return f = ((flux_t_D + flux_nt_D) / (flux_t_H + flux_nt_H)) / (hdo_0 / (2 * h2o_0))
end

function get_transport_PandL_rate(sp::Symbol, atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, n_horiz::Int64; returnfluxes=false, nonthermal=true, globvars...)
    #=
    Input:
        sp: species for which to return the transport production and loss
        atmdict: species number density by altitude for each vertical column
        returnfluxes: whether to return fluxes (thermal and nonthermal) instead of production/loss
        nonthermal: whether to consider nonthermal escape
	n_horiz: Number of vertical columns in simulation
    Output
        Array of production and loss (#/cm³/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :dz, :Hs_dict, :molmass,  :n_alt_index,
                                   :neutral_species, :num_layers, :polarizability, :q, :speciesbclist, :Te, :Ti, :Tn, :Tp, 
                                   :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    if nonthermal
        required = [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                                       :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict]
        check_requirements(keys(GV), required)
    end

    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, n_horiz; globvars...)
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict, n_horiz; globvars...)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>[fluxcoefs_all[s][ihoriz][2:end-1, :] for ihoriz in 1:n_horiz] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, sum([atmdict[sp] for sp in GV.all_species]), n_horiz; nonthermal=nonthermal, globvars...)

    # each element in thesebcs has the format [downward, upward]
    thesebcs = bc_dict[sp]

    # Fill array 
    transport_PL = [fill(convert(ftype_ncur, NaN), GV.num_layers) for ihoriz in 1:n_horiz]

    # These are the derivatives, which should be what we want (check math)
    for ihoriz in [1:n_horiz;]
    	transport_PL[ihoriz][1] = ((atmdict[sp][ihoriz][2]*fluxcoefs_bulk_layers[sp][ihoriz][2, 1]  # in from layer above
                        -atmdict[sp][ihoriz][1]*fluxcoefs_bulk_layers[sp][ihoriz][1, 2]) # out to layer above
                    +(-atmdict[sp][ihoriz][1]*thesebcs[ihoriz][1, 1] # out to boundary layer
                      +thesebcs[ihoriz][1, 2])) # in from the boundary layer
        for ialt in 2:length(transport_PL) - 1
            transport_PL[ihoriz][ialt] = ((atmdict[sp][ihoriz][ialt+1]*fluxcoefs_bulk_layers[sp][ihoriz][ialt+1, 1]  # coming in from above    
                               -atmdict[sp][ihoriz][ialt]*fluxcoefs_bulk_layers[sp][ihoriz][ialt, 2])    # leaving out to above layer
                             +(-atmdict[sp][ihoriz][ialt]*fluxcoefs_bulk_layers[sp][ihoriz][ialt, 1]     # leaving to the layer below
                               +atmdict[sp][ihoriz][ialt-1]*fluxcoefs_bulk_layers[sp][ihoriz][ialt-1, 2]))  # coming in from below
        end
    	transport_PL[ihoriz][end] = ((thesebcs[ihoriz][2, 2] # in from upper boundary layer - (non-thermal loss from flux bc)
                          - atmdict[sp][ihoriz][end]*thesebcs[ihoriz][2, 1]) # (#/cm³) * (#/s) out to space from upper bdy (thermal loss from velocity bc)
                        + (-atmdict[sp][ihoriz][end]*fluxcoefs_bulk_layers[sp][ihoriz][end, 1] # leaving out to layer below
                           +atmdict[sp][ihoriz][end-1]*fluxcoefs_bulk_layers[sp][ihoriz][end-1, 2])) # coming in to top layer from layer below
    end

    # Use these for a sanity check if you like. 
    # println("Activity in the top layer for sp $(sp) AS FLUX:")
    # println("Flux calculated from flux bc. for H and D, this should be the nonthermal flux: $([thesebcs[ihoriz][2, 2]*GV.dz for ihoriz in 1:n_horiz])")
    # println("Calculated flux from velocity bc. For H and D this should be thermal escape: $([atmdict[sp][ihoriz][end]*thesebcs[ihoriz][2, 1]*GV.dz for ihoriz in 1:n_horiz])")
    # println("Down to layer below: $(-atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")
    # println("In from layer below: $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz)")
    if returnfluxes
        tflux = zeros(0)
        if nonthermal
            ntflux = zeros(0)
        end
        for ihoriz in [1:n_horiz;]
            append!(tflux, atmdict[sp][ihoriz][end]*thesebcs[ihoriz][2, 1]*GV.dz)
            if nonthermal
                append!(ntflux, thesebcs[ihoriz][2, 2]*GV.dz)
                if sp in [:H, :D, :H2, :HD]
                    ntflux[ihoriz] = ntflux[ihoriz] < 0 ? abs(ntflux[ihoriz]) : throw("I somehow got a positive nonthermal flux, meaning it's going INTO the atmosphere? for $(sp)")
                else 
                    ntflux[ihoriz] = 0 
                end
            end
        end
        if nonthermal
            return ntflux, tflux
        else 
            return tflux
        end
    else 
        return transport_PL
    end
end

function get_directional_fluxes(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; nonthermal=true, return_up_n_down=false, globvars...)
    #=
    Returns the flux up and down from each atmospheric cell. 

    Input:
        sp: species for which to return the transport production and loss
        atmdict: species number density by altitude
        returnfluxes: whether to return fluxes (thermal and nonthermal) instead of production/loss
        nonthermal: whether to consider nonthermal escape 
    Output
        Array of production and loss (#/cm³/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :dz, :Hs_dict, :molmass, :n_alt_index,
               :neutral_species, :num_layers, :polarizability, :q, :speciesbclist, :Te, :Ti, :Tn, :Tp, 
               :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    if nonthermal
        required = [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                    :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict]
        check_requirements(keys(GV), required)
    end

    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, n_horiz; globvars...)
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict, n_horiz; globvars...)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>fluxcoefs_all[s][2:end-1, :] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, sum([atmdict[sp] for sp in GV.all_species]), n_horiz; nonthermal=nonthermal, globvars...)

    # each element in each vector within thesebcs has the format [downward, upward]
    thesebcs = bc_dict[sp]

    # Fill array 
    flux = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1)) # will store positive and negative values at each alt, with + meaning up, - meaning down.
    up = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1)) # Bethan asked for just up fluxes on time so now I also track them separately. 
    down = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1))
    flux[1] = NaN # 0 alt
    up[1] = NaN
    down[1] = NaN

    # Lower boundary 
    up2 = thesebcs[1, 2] # in from the boundary layer (upwards)
    down2 = atmdict[sp][1]*thesebcs[1, 1] # out to boundary layer (downwards)
    flux[2] = up2 - down2
    up[2] = up2
    down[2] = down2

    # println("Alt $(alt[2]): up $(up2*dz), down $(down2*dz), net $(up2*dz-down2*dz)")

    for ialt in 3:length(flux) - 2
        up_i = atmdict[sp][ialt-1]*fluxcoefs_bulk_layers[sp][ialt-1, 2] # in from layer below (upwards)
        down_i = atmdict[sp][ialt]*fluxcoefs_bulk_layers[sp][ialt, 1]     # out to the layer below [downwards]
        # println("Alt $(GV.alt[ialt]): up $(up*dz), down $(down*dz), net $(up*dz-down*dz)")
        flux[ialt] = up_i - down_i
        up[ialt] = up_i
        down[ialt] = down_i
    end

    up_penult = atmdict[sp][end]*thesebcs[2, 1] # (#/cm³) * (#/s) out to space from upper bdy (thermal loss from velocity bc) (upwards)
    down_penult = thesebcs[2, 2] # in from upper boundary layer (non-thermal loss from flux bc) (downwards) (has to be negative to come out positive)
    flux[end-1] = up_penult - down_penult
    up[end-1] = up_penult
    down[end-1] = down_penult
    # println("Alt $(alt[end-1]): up $(up_penult*dz), down $(down_penult*dz), net $(up_penult*dz-down_penult*dz)")

    flux[end] = NaN
    up[end] = NaN 
    down[end] = NaN

    # mult by dz so the array is truly a flux
    flux = flux .* GV.dz
    up = up .* GV.dz 
    down = down .* GV.dz

    # Use these for a sanity check if you like. 
    # println("Activity in the top layer for sp $(sp) AS FLUX:")
    # println("Flux calculated from flux bc. for H and D, this should be the nonthermal flux: $(thesebcs[2, 2]*GV.dz)")
    # println("Calculated flux from velocity bc. For H and D this should be thermal escape: $(atmdict[sp][end]*thesebcs[2, 1]*GV.dz)")
    # println("Down to layer below: $(-atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")
    # println("In from layer below: $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz)")
    # println("net in the second to last layer $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz -atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")

    if !return_up_n_down
        return flux
    else 
        return flux, up, down
    end
end

function flux_pos_and_neg(fluxarr) 
    #=
    Input:
        fluxarr: the output of function get_flux. 
    Outputs: 
        This generates two arrays, one with the positive flux
        and one with the negative flux, but all values are positive. This is just so 
        you can easily plot flux on a log axis with different markers for positive and negative.
    =#
    pos = []
    abs_val_neg = []

    for f in fluxarr
        if f > 0
            append!(pos, f)
            append!(abs_val_neg, NaN)
        else
            append!(abs_val_neg, abs(f))
            append!(pos, NaN)
        end
    end
    return pos, abs_val_neg
end

# Note: These functions are probably misleading. They were used to create plots that never made it 
# to publication.
function limiting_flux(sp, atmdict, T_arr; treat_H_as_rare=false, full_equation=true, globvars...)
    #=
    Calculate the limiting upward flux (Hunten, 1973; Zahnle, 2008). 
    Inputs:
        sp: A species that is traveling upwards
        atmdict: present atmospheric state
        T_arr: Array of neutral temperatures
    Output:
        Φ, limiting flux for a hydrostatic atmosphere
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :dz, :non_bdy_layers, :molmass, :n_alt_index, :n_all_layers]
    check_requirements(keys(GV), required)
    
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s
    if treat_H_as_rare==true
        if sp==:H
            thedensity = atmdict[:D]
        elseif sp==:H2 
            thedensity = atmdict[:HD]
        end
    elseif treat_H_as_rare==false
        thedensity = atmdict[sp]
    end

    if length(T_arr)==length(GV.alt)
        T_arr = T_arr[2:end-1]
    end

    Ha = scaleH(atmdict, T_arr, n_horiz; ignore=[sp], globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr) # AT^s

    if full_equation
        dTdz = zeros(GV.n_all_layers)
        dTdz = dTdz[2:end] = @. (T_arr[2:end] - T_arr[1:end-1]) / GV.dz # make the temp gradient
        print(dTdz)
        fi = thedensity ./ n_tot(atmdict, 1; ignore=[sp], globvars...)
        ma = meanmass(atmdict, n_horiz; ignore=[sp], globvars...) 

        return @. ((bi*fi)/(1+fi)) * ( mH*(ma - GV.molmass[sp]) * (g/(kB*T_arr)) - (thermaldiff(sp)/T_arr) * dTdz[1:end-1])
    else
        D = Dcoef_neutrals(non_bdy_layers, sp, bi, atmdict; globvars...)    
        return (D .* atmdict[sp] ./ Ha) .* (1 .- GV.molmass[sp] ./ meanmass(atmdict, n_horiz; ignore=[sp], globvars...))
    end
end

# function limiting_flow_velocity(sp, atmdict, T_arr; globvars...)
function limiting_flow_velocity(sp, atmdict, T_arr; ihoriz::Int=1, globvars...)
    #=
    Calculate the limiting upward flux (Hunten, 1973; Zahnle, 2008). 
    Inputs:
        sp: A species that is traveling upwards
        atmdict: present atmospheric state
        T_arr: Array of neutral temperatures
    Output:
        Φ, limiting flux for a hydrostatic atmosphere
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :n_alt_index, :non_bdy_layers, :molmass]
    check_requirements(keys(GV), required)
    
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s
    Ha = scaleH(atmdict, T_arr[2:end-1]; ignore=[sp], globvars..., alt=GV.non_bdy_layers)
    Hi = scaleH(GV.non_bdy_layers, sp, T_arr[2:end-1]; GV.molmass)
    bi = binary_dcoeff_inCO2(sp, T_arr[2:end-1]) # AT^s
    na = n_tot(atmdict, ihoriz; ignore=[sp], GV.all_species)

    return @. (bi / na) * (1/Ha - 1/Hi)
end

# function limiting_flux_molef(sp, atmdict, T_arr; globvars...)
function limiting_flux_molef(sp, atmdict, T_arr; ihoriz::Int=1, globvars...)
    #=
    Roger requested the limiting flux in in mole fraction. This is actually the same result as above. But this way we're sure
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :non_bdy_layers, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    avogadro = 6.022e23

    X = (atmdict[sp] ./ avogadro) ./ (n_tot(atmdict, ihoriz; globvars...) ./ avogadro)
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s

    Ha = scaleH(atmdict, T_arr, n_horiz; globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr)
    Hi = scaleH(non_bdy_layers, sp, T_arr, n_horiz; globvars...)

    return bi .* X .* (1 ./ Ha - 1 ./ Hi), X
end
