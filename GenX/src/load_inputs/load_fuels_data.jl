"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
    load_fuels_data!(setup::Dict, path::AbstractString, inputs::Dict)

Read input parameters related to fuel costs and CO$_2$ content of fuels
"""
function load_fuels_data!(setup::Dict, path::AbstractString, inputs::Dict)

    # Fuel related inputs - read in different files depending on if time domain reduction is activated or not
    data_directory = joinpath(path, setup["TimeDomainReductionFolder"])
    if setup["TimeDomainReduction"] == 1  && time_domain_reduced_files_exist(data_directory)
        my_dir = data_directory
    else
        my_dir = path
    end
    filename = "Fuels_data.csv"
    file_path = joinpath(my_dir, filename)
    fuels_in = DataFrame(CSV.File(file_path, header=true), copycols=true)

	use_minmax_supply = false
	start_ind = 2
	if (fuels_in[1, :Time_Index] == -2) && (fuels_in[2, :Time_Index] == -1)
		println("\tMinimum and Maximum supply constraints found")
		start_ind += 2
		use_minmax_supply = true
	end

    # Fuel costs & CO2 emissions rate for each fuel type
	fuels = names(fuels_in)[2:end] # fuel type indexes
	costs = Matrix(fuels_in[start_ind:end,2:end])
	CO2_content = fuels_in[start_ind-1,2:end] # tons CO2/MMBtu
	fuel_costs = Dict{AbstractString,Array{Float64}}()
	fuel_CO2 = Dict{AbstractString,Float64}()
	Minimum_Supply_MMBTU = Vector{Float64}(undef, length(fuels))
	Maximum_Supply_MMBTU = Vector{Float64}(undef, length(fuels))
	if use_minmax_supply
		Minimum_Supply_MMBTU_content = fuels_in[1,2:end] # New addition for minimum supply constraint
		Maximum_Supply_MMBTU_content = fuels_in[2,2:end] # New addition for maximum supply constraint
	end

    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    for i = 1:length(fuels)
            fuel_costs[fuels[i]] = costs[:,i] / scale_factor
            # fuel_CO2 is kton/MMBTU with scaling, or ton/MMBTU without scaling.
            fuel_CO2[fuels[i]] = CO2_content[i] / scale_factor
			if use_minmax_supply
				Minimum_Supply_MMBTU[i] = Minimum_Supply_MMBTU_content[i] / scale_factor 
				Maximum_Supply_MMBTU[i] = Maximum_Supply_MMBTU_content[i] / scale_factor 
            end
    end
    
    inputs["fuels"] = fuels
    inputs["fuel_costs"] = fuel_costs
    inputs["fuel_CO2"] = fuel_CO2
	inputs["Minimum_Supply_MMBTU"] = Minimum_Supply_MMBTU
	inputs["Maximum_Supply_MMBTU"] = Maximum_Supply_MMBTU

    println(filename * " Successfully Read!")

    return fuel_costs, fuel_CO2
end
