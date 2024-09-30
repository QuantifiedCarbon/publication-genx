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
	load_transmission_variability(setup::Dict, path::AbstractString, inputs_genvar::Dict)

Function for reading input parameters related to hourly transmission capacity availabile for all network lines.
"""
function load_transmission_variability(setup::Dict, path::AbstractString, inputs_transvar::Dict)

	# Hourly capacity factors
	#data_directory = chop(replace(path, pwd() => ""), head = 1, tail = 0)
	data_directory = joinpath(path, setup["TimeDomainReductionFolder"])
	# TimeDomainReduction for this Transmission_variability has not been implemented!!!
	if setup["TimeDomainReduction"] == 1  && isfile(joinpath(data_directory,"Load_data.csv")) && isfile(joinpath(data_directory,"Generators_variability.csv")) && isfile(joinpath(data_directory,"Transmission_variability.csv")) && isfile(joinpath(data_directory,"Fuels_data.csv")) # Use Time Domain Reduced data for GenX
		trans_var = DataFrame(CSV.File(joinpath(data_directory,"Transmission_variability.csv"), header=true), copycols=true)
	else # Run without Time Domain Reduction OR Getting original input data for Time Domain Reduction
		trans_var = DataFrame(CSV.File(joinpath(path,"Transmission_variability.csv"), header=true), copycols=true)
	end

	# DOES THIS NEED DOING?
	# Reorder DataFrame to R_ID order (order provided in Generators_data.csv)
	#select!(trans_var, [:Time_Index; Symbol.(inputs_genvar["RESOURCES"]) ])

	# Maximum power output and variability of each energy resource
	inputs_transvar["pTransVar_Max"] = transpose(Matrix{Float64}(trans_var[1:inputs_transvar["T"],2:(2*inputs_transvar["L"]+1)]))

	println("Transmission_variability.csv Successfully Read!")

	return inputs_transvar
end
