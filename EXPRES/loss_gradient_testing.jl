## Importing packages
using Pkg
Pkg.activate("EXPRES")

import StellarSpectraObservationFitting; SSOF = StellarSpectraObservationFitting
using JLD2
using Statistics
import StatsBase

## Setting up necessary variables

stars = ["10700", "26965", "34411"]
star = stars[SSOF.parse_args(1, Int, 2)]
interactive = length(ARGS) == 0
save_plots = true
include("data_locs.jl")  # defines expres_data_path and expres_save_path
desired_order = SSOF.parse_args(2, Int, 68)  # 68 has a bunch of tels, 47 has very few
use_reg = SSOF.parse_args(3, Bool, true)

## Loading in data and initializing model
save_path = expres_save_path * star * "/$(desired_order)/"
@load save_path * "data.jld2" n_obs data times_nu airmasses
if !use_reg
    save_path *= "noreg_"
end

## it appears that Array = BandedMatrix << SparseArray for basic gradient timing

using Zygote
using ParameterHandling
using BenchmarkTools
using SparseArrays
using MKLSparse
using BandedMatrices

inds = 1:1000
y = SSOF.FullLinearModel(ones(length(inds),2), ones(2,114), ones(length(inds)))
include("lsf.jl")
xs = EXPRES_lsf(exp.(data.log_λ_obs[inds,:]))
xb = BandedMatrix(xs)
xa = Array(xb)
# y = data.flux[inds,:]
function g_test(inps, l)
	p0, unflatten = flatten(inps)  # unflatten returns NamedTuple of untransformed params
	f = l ∘ unflatten
	return unflatten(only(Zygote.gradient(f, p0))), f, p0, unflatten
end

fb(y) = sum(xb * y())
_, fx, p0 = g_test(y, fb)
@btime Zygote.gradient(fx, p0)
fs(y) = sum(xs * y())
_, fx, p0 = g_test(y, fs)
Zygote.gradient(fx, p0)
fa(y) = sum(xa * y())
_, fx2, p02 = g_test(y, fa)
@btime Zygote.gradient(fx, p0)

using Juno
@profiler Zygote.gradient(fx, p0)
@profiler Zygote.gradient(fx2, p02)



# Diagonal
ProjectTo(x::Diagonal) = ProjectTo{Diagonal}(; diag=ProjectTo(x.diag))
(project::ProjectTo{Diagonal})(dx::AbstractMatrix) = Diagonal(project.diag(diag(dx)))
(project::ProjectTo{Diagonal})(dx::Diagonal) = Diagonal(project.diag(dx.diag))

# BandedMatrices
ProjectTo(x::Diagonal) = ProjectTo{BandedMatrix}(; diag=ProjectTo(x.diag))
(project::ProjectTo{BandedMatrix})(dx::AbstractMatrix) = Diagonal(project.diag(diag(dx)))
(project::ProjectTo{BandedMatrix})(dx::BandedMatrix) = Diagonal(project.diag(dx.diag))


## Nabla testing

using Nabla
using BenchmarkTools

yy = [y.M, y.s, y.μ]
ke
dump(y)

x = ones(3000, 3000)
y = ones(3000)
z = zeros(3000,3000)
zz = zeros(3000)
@time x, y[:] = zeros(3000,3000), zz
z[1] = 1
x[1]

fb(yy) = sum(xb * ((yy[1] * yy[2]) .+ yy[3]))
_, fx, p0, unf = g_test(yy, fb)
@btime flatten(∇(fb)(yy))
@btime only(Zygote.gradient(fx, p0))

fs(yy) = sum(xs * ((yy[1] * yy[2]) .+ yy[3]))
_, fx, p0, unf = g_test(yy, fs)
@btime flatten(∇(fs)(yy))
@btime only(Zygote.gradient(fx, p0))

83.3/4.465
494.9/1.364



##
inds = 1:100
x = data.log_λ_obs_bounds[1:inds[end]+1, :]
log_λ, _ = SSOF.create_λ_template(x; upscale=sqrt(2))
y = SSOF.FullLinearModel(ones(length(log_λ),2), ones(2,114), ones(length(log_λ)))
zs = SSOF.oversamp_interp_helper(x, log_λ)
za = [Array(i) for i in zs]


heatmap(SSOF.spectra_interp(y(), zs))
heatmap(SSOF.spectra_interp(y(), za))

ss = SSOF.spectra_interp(y(), zs)
ss[1,1]
zs[1] * y()


fs(y) = sum(SSOF.spectra_interp(y(), zs))
_, fx, p0 = g_test(y, fs)
@btime Zygote.gradient(fx, p0)

spectra_interp(y(), zs)
sum(xb * y())








# 7020, 114
ind_λ = 1:1200; ind_t = 1:114
data_small = SSOF.LSFData(data.flux[ind_λ, ind_t], data.var[ind_λ, ind_t], data.log_λ_obs[ind_λ, ind_t], data.log_λ_star[ind_λ, ind_t], data.lsf_broadener[ind_λ,ind_λ])

function lsf_broadener(λ::AbstractVector; safe::Bool=true)
    max_w = 0
    for i in 1:length(nwn)
        lo, mid, hi = SSOF.searchsortednearest(nwn, [nwn[i] - 3 * σs[i], nwn[i], nwn[i] + 3 * σs[i]])
        lsf = Normal(wn[i], σs[i])
        holder[i, lo:hi] = pdf.(lsf, wn[lo:hi])
        holder[i, lo:hi] ./= sum(view(holder, i, lo:hi))
        max_w = max(max_w, max(hi-mid, mid-lo))
    end
    return BandedMatrix(holder, (max_w, max_w))
end
if false#isfile(save_path*"results.jld2")
    @load save_path*"results.jld2" model rvs_naive rvs_notel
    if model.metadata[:todo][:err_estimated]
        @load save_path*"results.jld2" rv_errors
    end
    if model.metadata[:todo][:downsized]
        @load save_path*"model_decision.jld2" comp_ls ℓ aic bic ks test_n_comp_tel test_n_comp_star
    end
else
    model_upscale = sqrt(2)
    # model_upscale = 2 * sqrt(2)
    @time model = SSOF.OrderModel(data_small, "EXPRES", desired_order, star; n_comp_tel=3, n_comp_star=3, upscale=model_upscale)
    @time rvs_notel, rvs_naive, _, _ = SSOF.initialize!(model, data_small; use_gp=true)
    if !use_reg
        SSOF.zero_regularization(model)
        model.metadata[:todo][:reg_improved] = true
    end
    # @save save_path*"results.jld2" model rvs_naive rvs_notel
end


## Creating optimization workspace
workspace, loss = SSOF.OptimWorkspace(model, data_small; return_loss_f=true)

## loss testing

using Zygote
using ParameterHandling
using BenchmarkTools
ts = workspace.telstar
@btime ts.obj.f(ts.p0)
@btime Zygote.gradient(ts.obj.f, ts.p0)

nλ = size(data_small.flux, 1)

# No LSF and χ²
# _loss(tel::AbstractMatrix, star::AbstractMatrix, rv::AbstractMatrix) =
#     sum(((SSOF.total_model(tel, star, rv) .- workspace.d.flux) .^ 2) ./ workspace.d.var)
# LSF and χ²
# _loss(tel::AbstractMatrix, star::AbstractMatrix, rv::AbstractMatrix) =
#     mapreduce(i -> sum((((workspace.d.lsf_broadener[i] * (view(tel, :, i) .* (view(star, :, i) .+ view(rv, :, i)))) .- view(workspace.d.flux, :, i)) .^ 2) ./ view(workspace.d.var, :, i)), +, 1:size(tel, 2))
# _loss(tel::AbstractMatrix, star::AbstractMatrix, rv::AbstractMatrix) =
#     mapreduce(i -> sum((((workspace.d.lsf_broadener * (view(tel, :, i) .* (view(star, :, i) .+ view(rv, :, i)))) .- view(workspace.d.flux, :, i)) .^ 2) ./ view(workspace.d.var, :, i)), +, 1:size(tel, 2))
_loss(tel::AbstractMatrix, star::AbstractMatrix, rv::AbstractMatrix) =
    sum((((workspace.d.lsf_broadener * (tel .* (star .+ rv))) .- workspace.d.flux) .^ 2) ./ workspace.d.var)

# Interpolation
loss2(om::SSOF.OrderModel; tel::SSOF.LinearModel=om.tel.lm, star::SSOF.LinearModel=om.star.lm) =
	_loss(SSOF.tel_model(om; lm=tel), SSOF.star_model(om; lm=star), workspace.o.rv)
# no Interpolation
# loss2(om::SSOF.OrderModel; tel::SSOF.LinearModel=om.tel.lm, star::SSOF.LinearModel=om.star.lm) =
# 	_loss(tel()[1:nλ, :], star()[1:nλ, :], workspace.o.rv)

l_telstar(nt::NamedTuple{(:tel, :star,),<:Tuple{SSOF.LinearModel, SSOF.LinearModel}}) =
	loss2(workspace.om; tel=nt.tel, star=nt.star) + SSOF.model_prior(nt.tel, workspace.om.reg_tel) + SSOF.model_prior(nt.star, workspace.om.reg_star)

p0, unflatten = ParameterHandling.flatten(ts.θ)  # unflatten returns NamedTuple of untransformed params
f = l_telstar ∘ unflatten
@btime f(p0)
@btime Zygote.gradient(f, p0)


########## ONLY USELESS GARBAGE BELOW
## one parameter
function g_test(inps, l)
	p0, unflatten = flatten(inps)  # unflatten returns NamedTuple of untransformed params
	f = l ∘ unflatten
	return unflatten(only(Zygote.gradient(f, p0)))#, f, p0
end

f2(lm) = sum(lm())
g_test(model.tel.lm, f2)

x = copy(model.tel.lm)
x.M .= 0; x.μ .= 0; x.M[1000] = 1;
# x.M .= 1; x.μ .= 0;

model.tel.lm

f2(lm) = sum(lm())
g_test(model.tel.lm, f2).M[1000]
sum(x.M * x.s)
x.s[]
x.s

f2(lm) = sum(model.t2o[1] * lm())
g_test(model.tel.lm, f2).M[1000]
sum(model.t2o[1] * (x.M * x.s))

sm = SSOF.star_model(model) + SSOF.rv_model(model)
f2(lm) = sum((model.t2o[1] * lm() .* sm))
g_test(model.tel.lm, f2).M[1000]
sum((model.t2o[1] * (x.M * x.s)) .* sm)

f2(lm) = sum(data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm))
g_test(model.tel.lm, f2).M[1000]
sum(data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm))

f2(lm) = sum(data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux)
g_test(model.tel.lm, f2).M[1000]
sum(data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm))

f2(lm) = sum((data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux) .^ 2)
g_test(model.tel.lm, f2).M[1000]
β = (data_small.lsf_broadener[1] * (model.t2o[1] * model.tel.lm() .* sm) - data_small.flux)
dβ = data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm)
sum(2 .* β .* dβ)

f2(lm) = sum(((data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux) .^ 2) ./ data_small.var)
g_test(model.tel.lm, f2).M[1000]
sum(2 .* β .* dβ ./ data_small.var)

## all M

x = copy(model.tel.lm)
x.μ .= 0; x.M .= 0; x.M[116,3] = 1;

(x.M * x.s)
(x.M * x.s)[116,:]
x.s[3, :]
f2(lm) = sum(lm())
@time g_test(model.tel.lm, f2).M
@time repeat(sum(x.s, dims=2)', size(x.M, 1))


f2(lm) = sum(model.t2o[1] * lm())
@btime g_test(model.tel.lm, f2).M
# sum(model.t2o[1] * (x.M * x.s))
dM = zeros(size(x.M))
mapreduce((i,j) -> i*j, +, 1:size(dM, 1), 1:size(dM, 2))

repeat(1:size(dM, 1), size(dM, 2))

1:size(dM, 2)
1:size(dM, 1)
sum([i^2 for i in 1:8])

hcat([interp_helper[i] * view(model, :, i) for i in 1:size(model, 2)]...)
@btime for i in 1:size(dM, 1)
	for j in 1:size(dM, 2)
		dM[i,j] = sum(kron(view(model.t2o[1], :, i),view(x.s, j, :)))
	end
end
dM

i=100;j=3;
@time sum(view(model.t2o[1], :, i)*transpose(view(x.s, j, :)))
@time sum(view(model.t2o[1], :, i)*view(x.s, j, :)')
@time sum(view(model.t2o[1], :, i).*view(x.s, j, :)')
@time sum(kron(view(model.t2o[1], :, i),view(x.s, j, :)))


sm = SSOF.star_model(model) + SSOF.rv_model(model)
f2(lm) = sum((model.t2o[1] * lm() .* sm))
@time g_test(model.tel.lm, f2).M[1000]
sum((model.t2o[1] * (x.M * x.s)) .* sm)

f2(lm) = sum(data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm))
@time g_test(model.tel.lm, f2).M[1000]
sum(data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm))

f2(lm) = sum(data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux)
@time g_test(model.tel.lm, f2).M[1000]
sum(data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm))

f2(lm) = sum(abs2, data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux)
@time g_test(model.tel.lm, f2).M[1000]
β = (data_small.lsf_broadener[1] * (model.t2o[1] * model.tel.lm() .* sm) - data_small.flux)
dβ = data_small.lsf_broadener[1] * ((model.t2o[1] * (x.M * x.s)) .* sm)
sum(2 .* β .* dβ)

f2(lm) = sum(((data_small.lsf_broadener[1] * (model.t2o[1] * lm() .* sm) - data_small.flux) .^ 2) ./ data_small.var)
@time _, f3, p0 = g_test(model.tel.lm, f2)
sum(2 .* β .* dβ ./ data_small.var)

@time Zygote.gradient(f3, p0)
