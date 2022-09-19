import StellarSpectraObservationFitting as SSOF
using JLD2
using CSV, DataFrames, Query, StatsBase, Statistics, Dates
using FITSIO
using Statistics
using Plots

function n2s(i)
	@assert -1 < i < 1000
	ans = string(i)
	return "0"^(3-length(ans))*ans
end

function neid_extras(df_files::DataFrame, save_path_base::String)
	n_obs = nrow(df_files)
	f = FITS(df_files.Filename[1])
	neid_tel = zeros(n_obs, size(f[11], 1), size(f[11], 2))
	neid_tel_wvapor = zeros(n_obs)
	neid_tel_zenith = zeros(n_obs)
	_df = DataFrame(f[14])
	inds = _df.INDEX
	lcs = _df.LINE_CENTER
	d_lcs = Dict()
	for i in 1:length(lcs)
		_lcs = split(lcs[i], ", ")
		for j in _lcs
			if haskey(d_lcs, j)
				append!(d_lcs[j], [inds[i]])
			else
				d_lcs[j] = [inds[i]]
			end
		end
	end
	Hα = "6562.808"
	if Hα in keys(d_lcs) && "Ha06_2" in d_lcs[Hα]; d_lcs[Hα] = ["Ha06_2"] end
	df_cols = String[]
	for i in inds
		append!(df_cols, [i, i*"_σ"])
	end
	df_act = zeros(n_obs, length(df_cols))
	neid_time = zeros(n_obs)
	neid_rv = zeros(n_obs)
	neid_rv_σ = zeros(n_obs)
	neid_order_rv = zeros(n_obs, 118)
	for i in 1:n_obs # at every time
		f = FITS(df_files.Filename[i])
		driftfun = read_header(f[1])["DRIFTFUN"]
		if driftfun != "dailymodel0"
			println("spectrum $i ($(df_files.Filename[i])) has a wavelength calib. drift func \"$driftfun\" instead of \"dailymodel0\", consider removing it from your analysis")
		end
		neid_tel[i, :, :] .= read(f[11])
		_df_h = read_header(f[11])
		neid_tel_wvapor[i] = _df_h["WVAPOR"]
		neid_tel_zenith[i] = _df_h["ZENITH"]
		_df = DataFrame(f[14])
		df_act[i, 1:2:end] = _df.VALUE
		df_act[i, 2:2:end] = _df.UNCERTAINTY
		ccf_header = read_header(f[13])
		neid_time[i] = ccf_header["CCFJDMOD"]
		neid_rv[i] = ccf_header["CCFRVMOD"] * 1000  # m/s
		neid_rv_σ[i] = ccf_header["DVRMSMOD"] * 1000  # m/s
		ccfs_exist = vec(.!all(iszero.(read(f[13])); dims=1))
		# neid_rv_ords = [j for j in 1:length(ccfs_exist) if ccfs_exist[j]]
		# neid_rv_ords = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 72, 73, 74, 75, 76, 77, 78, 79, 81, 82, 83, 84, 85, 91, 92, 93, 95, 96]
		# println(ccf_header)
		for j in 1:length(ccfs_exist)
		    if ccfs_exist[j]
				# neid_order_rv[i, j] = ccf_header["CCFRV"*n2s(j+51)] * 1000  # m/s
				neid_order_rv[i, j] = ccf_header["CCFRV"*n2s(174-j)] * 1000  # m/s
		    end
		end
	end
	d_act_tot = Dict()
	for i in 1:2:length(df_cols)
		d_act_tot[df_cols[i]] = df_act[:, i]
		d_act_tot[df_cols[i+1]] = df_act[:, i+1]
	end
	@save save_path_base * "/neid_pipeline.jld2" neid_time neid_rv neid_rv_σ neid_order_rv d_act_tot neid_tel d_lcs neid_tel_wvapor neid_tel_zenith
end

function neid_activity_indicators(pipeline_path::String, data::SSOF.Data)
	@load pipeline_path d_act_tot d_lcs
	lo, hi = exp.(quantile(vec(data.log_λ_star), [0.05, 0.95]))
	df_act = Dict()
	for wv in keys(d_lcs)
		if lo < parse(Float64, wv) < hi
			for key in d_lcs[wv]
				df_act[key] = d_act_tot[key]
				df_act[key*"_σ"] = d_act_tot[key*"_σ"]
			end
		end
	end
	return df_act
end

function neid_plots(mws::SSOF.ModelWorkspace,
	airmasses::AbstractVector,
	times_nu::AbstractVector,
	rvs::AbstractVector,
	rv_errors::AbstractVector,
	star::String,
	base_path::String,
	pipeline_path::String,
	desired_order::Int;
	mask = :,
	display_plt::Bool=false,
	tel_errors::Union{AbstractMatrix, Nothing}=nothing,
	star_errors::Union{AbstractMatrix, Nothing}=nothing,
	df_act::Dict=Dict(),
	kwargs...)

	@load pipeline_path neid_time neid_rv neid_rv_σ neid_order_rv
	neid_time .-= 2400000.5

	# Compare RV differences to actual RVs from activity
	plt = plot_model_rvs(view(times_nu, mask), view(rvs, mask), view(rv_errors, mask), view(neid_time, mask), view(neid_rv, mask), view(neid_rv_σ, mask); display_plt=display_plt, title="$star (median σ: $(round(median(vec(view(rv_errors, mask))), digits=3)))");
	png(plt, base_path * "model_rvs.png")

	save_model_plots(mws, airmasses, times_nu, base_path; display_plt=display_plt, tel_errors=tel_errors, star_errors=star_errors, df_act=df_act, kwargs...);

	if all(.!iszero.(view(neid_order_rv, :, desired_order)))
	    plt = plot_model_rvs(view(times_nu, mask), view(rvs, mask), view(rv_errors, mask), view(neid_time, mask), view(neid_order_rv, mask, desired_order), zeros(length(view(neid_time, mask))); display_plt=display_plt, title="$star (median σ: $(round(median(vec(rv_errors)), digits=3)))");
	    png(plt, base_path * "model_rvs_order.png")
	end
end


# used to easily identify how far to mask https://apps.automeris.io/wpd/
function neid_order_masks!(data::SSOF.Data, order::Int, star::String)
	if star=="26965"
		if order==31
		    SSOF.mask_stellar_features!(data, log(4326.7), 100)
		elseif order==38
		    SSOF.mask_stellar_features!(data, log(4549.5), 100)
		elseif order==41
		    SSOF.mask_stellar_features!(data, log(4651.9), 100)
		elseif order==47
		    SSOF.mask_stellar_features!(data, log(4871.6), 100)
		elseif order==48
		    SSOF.mask_stellar_features!(data, log(4909.8), 100)
		elseif order==60
		    SSOF.mask_tellurics!(data, 0, log(5326.2))
		elseif order==61
			# not sure which it is
		    SSOF.mask_stellar_features!(data, 0, log(5374))
			SSOF.mask_tellurics!(data, 0, log(5374.5))
		elseif order==95
		    SSOF.mask_stellar_features!(data, log(7832.6), 100)
		end
	end
end
