## Importing packages
using Pkg
Pkg.activate("EXPRES")
Pkg.instantiate()

using JLD2
import StellarSpectraObservationFitting as SSOF
using Dates
## Setting up necessary variables

stars = ["10700", "26965", "34411"]
# orders_list = [1:85, 1:85, 1:85]
orders_list = [1:85, 1:85, 1:85]
include("data_locs.jl")  # defines expres_data_path and expres_save_path
# prep_str = "noreg_"
prep_str = ""
cutoff = now() - Week(1)
input_ind = SSOF.parse_args(1, Int, 0)
delete = SSOF.parse_args(2, Bool, false)

function clean(order::Int, star::String)
    dir = expres_save_path*star*"/$(order)/"
    if isdir(dir)
        ls = readdir(dir)
        println(order)
        for file in ls
            if file != "data.jld2" && !isdir(dir * file) && (mtime(dir * file) < datetime2unix(cutoff))
                println(file)
                if delete; rm(dir * file) end
            end
        end
    else
        println("couldn't find " * dir)
    end
end

input_ind == 0 ? 1:length(stars) : star_inds = input_ind
for star_ind in star_inds
    star = stars[star_ind]
    orders = orders_list[star_ind]
    n_ord = length(orders)
    for i in 1:n_ord
        clean(orders[i], star)
    end
end
