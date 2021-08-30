## Importing packages
using Pkg
Pkg.activate("EXPRES")
Pkg.instantiate()

using JLD2
import StellarSpectraObservationFitting; SSOF = StellarSpectraObservationFitting

## Setting up necessary variables

stars = ["10700", "26965", "34411"]
orders_list = [1:85, 1:85, 1:85]
include("data_locs.jl")  # defines expres_data_path and expres_save_path
# prep_str = "noreg_"
prep_str = ""

function retrieve(order::Int, star::String)
    @load expres_save_path*star*"/$(order)/$(prep_str)results.jld2" model rv_errors
    rvs_notel_opt = SSOF.rvs(model)
    return rvs_notel_opt, rv_errors
end

function retrieve_md(order::Int, star::String)
    @load expres_save_path*star*"/$(order)/$(prep_str)model_decision.jld2" comp_ls ℓ aic bic ks test_n_comp_tel test_n_comp_star
    ans_aic = argmin(aic)
    ans_bic = argmin(bic)
    n_comps = [test_n_comp_tel[ans_aic[1]], test_n_comp_star[ans_aic[2]]]
    n_comps_bic = [test_n_comp_tel[ans_bic[1]], test_n_comp_star[ans_bic[2]]]
    return n_comps, n_comps_bic, ans_aic==ans_bic
end

input_ind = SSOF.parse_args(1, Int, 0)
input_ind == 0 ? star_inds = (1:3) : star_inds = input_ind
for star_ind in star_inds
    star = stars[star_ind]
    orders = orders_list[star_ind]

    @load expres_save_path*star*"/$(orders[1])/data.jld2" n_obs times_nu airmasses
    n_ord = length(orders)
    rvs = zeros(n_ord,  n_obs)
    rvs_σ = zeros(n_ord, n_obs)
    n_comps = zeros(Int, n_ord, 2)
    n_comps_bic = zeros(Int, n_ord, 2)
    robust = zeros(Bool, n_ord)
    for i in 1:n_ord
        try
            rvs[i, :], rvs_σ[i, :] = retrieve(orders[i], star)
            n_comps[i, :], n_comps_bic[i, :], robust[i] = retrieve_md(orders[i], star)
        catch
            rvs_σ[i, :] .= Inf
            n_comps[i, :] .= -1
            println("order $(orders[i]) is missing")
        end
    end

    @save "$(prep_str)$(star)_rvs.jld2" rvs rvs_σ n_obs times_nu airmasses n_ord
    @save "$(prep_str)$(star)_md.jld2" n_comps n_comps_bic robust
end
