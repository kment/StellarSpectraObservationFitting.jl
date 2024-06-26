# Heavily inspired (if not blatantly ripped of from) the functions at
# https://github.com/megbedell/wobble/blob/master/wobble/data.py
using LinearAlgebra
using Statistics
import StatsBase: winsor

_high_quantile_default = 0.9

"""
	fit_continuum(x, y, σ²; ignore_weights=false, order=6, nsigma=[0.5,3.0], maxniter=50, edge_mask::Int=0)

Fit the continuum of `y` with a polynomial of order `order` using assymetric sigma-clipping, similar to `wobble`
"""
function fit_continuum(x::AbstractVector, y::AbstractVector, σ²::AbstractVector; ignore_weights::Bool=false, order::Int=6, nsigma::Vector{<:Real}=[0.5,3.0], maxniter::Int=50, edge_mask::Int=0)
    @assert 0 <= order < length(x)
    @assert length(x) == length(y) == length(σ²)
    @assert length(nsigma) == 2

    A = _vander(x .- mean(x), order)  # the constant design matrix
    m = fill(true, length(x))  # the sigma-clipped mask
	m[y .< 0.5] .= false  # mask out the low pixels
	m[σ² .== Inf] .= false  # mask out the bad pixels

	# masking out edge if desired
	if edge_mask > 0
		# m[edge_pad+1:edge_mask+edge_pad] .= false
		# m[end-edge_mask-edge_pad+1:end-edge_pad] .= false
		hold_left = y[1:edge_mask]
		hold_right = y[end-edge_mask+1:end]
		y[1:edge_mask] .= 1
		y[end-edge_mask+1:end] .= 1
	end

	# preallocating arrays
    μ = ones(length(x))
	w = Array{Float64}(undef, order+1)

	# keep fitting polynomial and sigma clipping until the sigma-clipped mask stops changing
    for i in 1:maxniter
		if ignore_weights
			w[:] = general_lst_sq(view(A, m, :), view(y, m))
        else
			w[:] = general_lst_sq(view(A, m, :), view(y, m), view(σ², m))
		end
        μ[:] = A * w
        resid = y - μ
        # sigma = median(abs.(resid))
		sigma = std(resid)
        m_new = (-nsigma[1]*sigma) .< resid .< (nsigma[2]*sigma)
		m_new[y .< 0.5] .= false
		m_new[σ² .== Inf] .= false  # mask out the bad pixels
        if sum(m) == sum(m_new); break end
        m = m_new
    end

	# restoring edge
	if edge_mask > 0
		y[1:edge_mask] .= hold_left
		y[end-edge_mask+1:end] .= hold_right
	end

    return μ, w
end


"""
	continuum_normalize!(d; order=6, kwargs...)

Divide `d.flux` (and `d.var`) by its sigma-clipped polynomial continuum (or its square)
"""
function continuum_normalize!(d::Data; order::Int=6, kwargs...)
	continuum = ones(size(d.log_λ_obs, 1))
	w = Array{Float64}(undef, order + 1, size(d.flux, 2))

	# fit the continuum at each time
	for i in axes(d.log_λ_obs, 2)
		continuum[:], w[:, i] = fit_continuum(view(d.log_λ_obs, :, i), view(d.flux, :, i), view(d.var, :, i); order=order, kwargs...)
		d.flux[:, i] ./= continuum
		d.var[:, i] ./= continuum .* continuum
		d.var_s[:, i] ./= continuum .* continuum
	end
	return w
end

_min_flux_default = 0.
_max_flux_default = 2.


or(a::Bool, b::Bool) = a || b

"""
	mask_infinite_pixels!(y, σ², log_λ_star, log_λ_star_bounds; padding=0, include_bary_shifts=false, verbose=true, kwargs...)

Mask out pixels whose `y` or `σ²` values are infinite
"""
function mask_infinite_pixels!(y::AbstractMatrix, σ²::AbstractMatrix, log_λ_star::AbstractMatrix, log_λ_star_bounds::AbstractMatrix; padding::Int=0, include_bary_shifts::Bool=false, verbose::Bool=true, kwargs...)
	i = findall(vec(all(or.(.!isfinite.(y), .!isfinite.(σ²)); dims=2)))
	if length(i) > 0
		if verbose; println("Instrumental pipeline already masked out pixels $i at all times") end
		if include_bary_shifts
			return mask_stellar_pixel!(σ², log_λ_star, log_λ_star_bounds, i; padding=padding, verbose=false)
		else
			return mask!(σ², i; padding=padding)
		end
	end
	return Int[]
end
function mask_infinite_pixels!(d::Data; kwargs...)
	mask_infinite_pixels!(d.flux, d.var_s, d.log_λ_star, d.log_λ_star_bounds; include_bary_shifts=true, kwargs...)
	return mask_infinite_pixels!(d.flux, d.var, d.log_λ_star, d.log_λ_star_bounds; include_bary_shifts=false, verbose=false, kwargs...)
end

"""
	mask_low_pixels!(bad, y, σ²; padding=0, min_flux=0., padding=2, using_weights=false)

Mask out pixels whose `y` values are below `min_flux`
"""
function mask_low_pixels!(bad::AbstractVector, y::AbstractVector, σ²::AbstractVector; min_flux::Real=_min_flux_default, padding::Int=2, using_weights::Bool=false)
	bad[:] = y .< min_flux
	# y[bad] .= min_flux
	l = length(bad)
	for i in findall(bad)
		bad[max(1, i - padding):min(i + padding, l)] .= true
	end
	return mask!(σ², bad; using_weights=using_weights)
end
function mask_low_pixels!(y::AbstractVector, σ²::AbstractVector; kwargs...)
	bad = Array{Bool}(undef, length(y))
	mask_low_pixels!(bad, y, σ²; kwargs...)
end
function mask_low_pixels!(y::AbstractMatrix, σ²::AbstractMatrix; kwargs...)
	affected = Int64[]
	bad = Array{Bool}(undef, size(y, 1))
	for i in axes(y, 2)
		affected2 = mask_low_pixels!(view(y, :, i), view(σ², :, i); kwargs...)
		affected_pixels!(affected, affected2)
	end
	return affected
end
function mask_low_pixels!(d::Data; kwargs...)
	mask_low_pixels!(d.flux, d.var_s; kwargs...)
	return mask_low_pixels!(d.flux, d.var; kwargs...)
end


"""
	mask_high_pixels!(bad, y, σ²; padding=0, max_flux=2., padding=2, using_weights=false)

Mask out pixels whose `y` values are above `max_flux`
"""
function mask_high_pixels!(bad::AbstractVector, y::AbstractVector, σ²::AbstractVector; max_flux::Real=_max_flux_default, padding::Int=2, kwargs...)
	bad[:] = y .> max_flux
	# y[bad] .= max_flux
	l = length(bad)
	for i in findall(bad)
		bad[max(1, i - padding):min(i + padding, l)] .= true
	end
	return mask!(σ², bad; kwargs...)
end
function mask_high_pixels!(y::AbstractVector, σ²::AbstractVector; kwargs...)
	bad = Array{Bool}(undef, length(y))
	mask_high_pixels!(bad, y, σ²; kwargs...)
end
function mask_high_pixels!(y::AbstractMatrix, σ²::AbstractMatrix; kwargs...)
	affected = Int64[]
	bad = Array{Bool}(undef, size(y, 1))
	for i in axes(y, 2)
		affected2 = mask_high_pixels!(bad, view(y, :, i), view(σ², :, i); kwargs...)
		affected_pixels!(affected, affected2)
	end
	return affected
end
function mask_high_pixels!(d::Data; kwargs...)
	mask_high_pixels!(d.flux, d.var_s; kwargs...)
	return mask_high_pixels!(d.flux, d.var; kwargs...)
end


"""
	mask_bad_edges!(y, σ², log_λ_star, log_λ_star_bounds; window_width=128, min_snr=8., verbose=true, always_mask_something=false, edges=nothing, kwargs...)

Mask out edges of `y` (either given by `edges` or where `y`/`σ` is on average below `min_snr`)
"""
function mask_bad_edges!(y::AbstractMatrix, σ²::AbstractMatrix, log_λ_star::AbstractMatrix, log_λ_star_bounds::AbstractMatrix; window_width::Int=128, min_snr::Real=8., verbose::Bool=true, always_mask_something::Bool=false, edges=nothing, kwargs...)
	
	# if `edges` isn't given, find how much to mask based on SNR in a sliding window
	if isnothing(edges)
		n_pix = size(y, 1)
		window_start_tot = 0
		window_end_tot = n_pix + 1
		for i in axes(y, 2)
			for window_start in 1:Int(floor(window_width/10)):(n_pix - window_width)
				window_end = window_start + window_width
				mean_snr = sqrt(mean((y[window_start:window_end, i] .^2) ./ abs.(σ²[window_start:window_end, i])))
				if mean_snr > min_snr
					window_start_tot = max(window_start_tot, window_start)
					break
				end
			end
		end
		for i in axes(y,2)
			for window_end in n_pix:-Int(floor(window_width/10)):(window_width + 1)
				window_start = window_end - window_width
				mean_snr = sqrt(mean((y[window_start:window_end] .^2) ./ abs.(σ²[window_start:window_end])))
				if mean_snr > min_snr
					window_end_tot = min(window_end_tot, window_end)
					break
				end
			end
		end
	else
		window_start_tot, window_end_tot = edges
	end

	# do the masking
	if always_mask_something || window_start_tot > 1
		# σ²[1:window_start_tot, :] .= Inf # trim everything to left of window
		if verbose; println("masking out low SNR edge from 1:$window_start_tot") end
		log_λ_low = maximum(view(log_λ_star_bounds, max(1, window_start_tot), :))
	else
		log_λ_low = -Inf
	end
	if always_mask_something || window_end_tot < n_pix
		# σ²[window_end_tot:end, :] .= Inf # trim everything to right of window
		log_λ_high = minimum(view(log_λ_star_bounds, min(window_end_tot, size(log_λ_star_bounds, 1)), :))
		if verbose; println("masking out low SNR edge from $window_end_tot:end") end
	else
		log_λ_high = Inf
	end
	if always_mask_something || window_start_tot > 1 || window_end_tot < n_pix
		return mask_stellar_feature!(σ², log_λ_star, log_λ_low, log_λ_high; inverse=true, verbose=false, kwargs...), [window_start_tot, window_end_tot]
	else
		return Int[], [window_start_tot, window_end_tot]
	end
end
function mask_bad_edges!(d::Data; kwargs...)
	affected, edges = mask_bad_edges!(d.flux, d.var, d.log_λ_star, d.log_λ_star_bounds; kwargs...)
	mask_bad_edges!(d.flux, d.var_s, d.log_λ_star, d.log_λ_star_bounds; always_mask_something=true, verbose=false, edges=edges, kwargs...)
	return affected
end

"""
	flat_normalize!(d; kwargs...)

Divide `d.flux` (and `d.var`) by its `0.9` quantile value (or its square)
"""
function flat_normalize!(d::Data; kwargs...)
	for i in axes(d.log_λ_obs, 2)
		continuum = quantile(view(d.flux, .!(isnan.(view(d.flux, :, i))), i), _high_quantile_default)
		d.flux[:, i] ./= continuum
		d.var[:, i] ./= continuum * continuum
		d.var_s[:, i] ./= continuum * continuum
	end
end


and(a::Bool, b::Bool) = a && b


"""
	outlier_mask(v; thres=10, prop=0.2, return_stats=false, only_low::Bool=false)

Mask outliers in `v` that are outside `thres`-σ after winsorizing out the top and bottom `prop` quantiles
"""
function outlier_mask(v::AbstractVecOrMat; thres::Real=10, prop::Real=0.2, return_stats::Bool=false, only_low::Bool=false)
	wv = winsor(v; prop=prop)
	μ = mean(wv)
	σ = stdm(wv, μ)

	if only_low
		mask = v .> (μ - thres * σ)
		if return_stats
			return mask, (v .- μ) ./ σ
		end
	else
		mask = and.(v .< (μ + thres * σ), v .> (μ - thres * σ))
	end

	if return_stats
		return mask, abs.((v .- μ) ./ σ)
	end
	return mask
end

# function recognize_bad_normalization!(d::Data; kwargs...)
# 	mask = outlier_mask([mean(view(d.var, isfinite.(view(d.var, :, i)), i)) for i in axes(d.var, 2)]; kwargs...) .|| outlier_mask(vec(std(d.flux; dims=1)); kwargs...)
# 	for i in axes(d.log_λ_obs, 2)
# 		if !mask[i]
# 			# d.var[:, i] .= Inf
# 			println("spectrum $i has a weird continuum normalization, consider removing it from your analysis")
# 		end
# 	end
# end

# function recognize_bad_drift!(d::Data; kwargs...)
# 	mask = outlier_mask(d.log_λ_obs[1, :]; kwargs...)
# 	for i in axes(d.log_λ_obs, 2)
# 		if !mask[i]
# 			# d.var[:, i] .= Inf
# 			println("spectrum $i has a weird drift, consider removing it from your analysis")
# 		end
# 	end
# end


"""
	snap(y)

Estimate δ⁴`y`/δ⁴x (a.k.a. the snap) with finite differences. See https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences
"""
function snap(y::AbstractMatrix)
	snp = Array{Float64}(undef, size(y, 1), size(y, 2))
	l = size(snp, 1)
	snp[1:2, :] .= 0
	snp[end-1:end, :] .= 0
	snp[3:end-2, :] .= abs.(view(y, 5:l, :) - 4view(y, 4:(l-1), :) + 6view(y, 3:(l-2), :) - 4view(y, 2:(l-3), :) + view(y, 1:(l-4), :))
	return snp
end
function snap(y::AbstractMatrix, σ²::AbstractMatrix)
	@assert size(y) == size(σ²)
	snp = Array{Float64}(undef, size(y, 1), size(y, 2))
	m = isfinite.(σ²)
	l = size(m, 1)
	# snp[1, :] = abs.(view(y, 3, :) - 4view(y, 2, :) + 6view(y, 1, :))
	# snp[2, :] = abs.(view(y, 4, :) - 4view(y, 3, :) + 6view(y, 2, :) - 4view(y, 1, :))
	# snp[end-1, :] = abs.(-4view(y, l, :) + 6view(y, l-1, :) - 4view(y, l-2, :) + view(y, l-3, :))
	# snp[end, :] = abs.(6view(y, l, :) - 4view(y, l-1, :) + view(y, l-2, :))
	snp[1:2, :] .= 0
	snp[end-1:end, :] .= 0
	snp[3:end-2, :] .= abs.(view(y, 5:l, :) - 4view(y, 4:(l-1), :) + 6view(y, 3:(l-2), :) - 4view(y, 2:(l-3), :) + view(y, 1:(l-4), :)) .* view(m, 5:l, :) .* view(m, 4:(l-1), :) .* view(m, 2:(l-3), :) .* view(m, 1:(l-4), :)
	# snp[.!m] .= 0
	snp[all(.!m; dims=2), :] .= 0
	return snp
end
function _snap!(snp::AbstractVector, y::AbstractVector; def_val::Real=1.)
	l = length(snp)
	@assert length(y) == l
	snp[1] = abs.(y[3] - 4y[2] + 6y[1] - 4def_val + def_val)
	snp[2] = abs.(y[4] - 4y[3] + 6y[2] - 4y[1] + def_val)
	snp[end-1] = abs.(def_val - 4y[l] + 6y[l-1] - 4y[l-2] + y[l-3])
	snp[end] = abs.(def_val - 4def_val + 6y[l] - 4y[l-1] + y[l-2])
	snp[3:end-2] .= abs.(view(y, 5:l) - 4view(y, 4:(l-1)) + 6view(y, 3:(l-2)) - 4view(y, 2:(l-3)) + view(y, 1:(l-4)))
	return snp
end
function _snap(y::AbstractVector; kwargs...)
	snp = Array{Float64}(undef, length(y))
	_snap!(snp, y; kwargs...)
	return snp
end


"""
	bad_pixel_flagger(y, σ²; prop=.005, thres=8)

Find anomalous pixels where δ⁴`y`/δ⁴x is high, passing `prop` and `thres` to `outlier_mask()`
"""
function bad_pixel_flagger(y::AbstractMatrix, σ²::AbstractMatrix; prop::Real=.005, thres::Real=8)
	snp = snap(y, σ²)
	snp = vec(mean(snp; dims=2))
	high_snap_pixels = find_modes(snp)
	return high_snap_pixels[.!outlier_mask(snp[high_snap_pixels]; prop=prop, thres=thres)]
end


"""
	mask_bad_pixel!(y, σ², log_λ_star, log_λ_star_bounds; padding=2, include_bary_shifts=false, verbose=true, bad_pixels=nothing, kwargs...)

Mask anomalous pixels where δ⁴`y`/δ⁴x is high
"""
function mask_bad_pixel!(y::AbstractMatrix, σ²::AbstractMatrix, log_λ_star::AbstractMatrix, log_λ_star_bounds::AbstractMatrix; padding::Int=2, include_bary_shifts::Bool=false, verbose::Bool=true, bad_pixels=nothing, kwargs...)
	if isnothing(bad_pixels); bad_pixels = bad_pixel_flagger(y, σ²; kwargs...) end
	if length(bad_pixels) > 0
		if length(bad_pixels) > 15
			if verbose; println("lots of snappy pixels, investigate?") end
		else
			if verbose; println("masked out high snap pixels $bad_pixels at all times") end
		end
		if include_bary_shifts
			return mask_stellar_pixel!(σ², log_λ_star, log_λ_star_bounds, bad_pixels; padding=padding, verbose=false), bad_pixels
		else
			return mask!(σ², bad_pixels; padding=padding), bad_pixels
		end
	end
	return Int[], bad_pixels
end
function mask_bad_pixel!(d::Data; kwargs...)
	affected, bad_pixels = mask_bad_pixel!(d.flux, d.var, d.log_λ_star, d.log_λ_star_bounds; include_bary_shifts=false, verbose=false, kwargs...)
	mask_bad_pixel!(d.flux, d.var_s, d.log_λ_star, d.log_λ_star_bounds; include_bary_shifts=true, bad_pixels=bad_pixels, kwargs...)
	return affected
end


"""
	mask_isolated_pixels!(σ²; neighbors_required::Int=29, verbose::Bool=true)

Mask between two flagged pixels if there are fewer than `neighbors_required` pixels between them
"""
function mask_isolated_pixels!(σ²::AbstractMatrix; neighbors_required::Int=29, verbose::Bool=true)
	affected = Int[]
	lo = 1
	hi = 1
	m = vec(all(isinf.(σ²); dims=2))
	l = length(m)
	while lo <= l

		# if pixel is already flagged, keep going
		if m[lo]
			lo += 1

		# if pixel is not flagged, go until the next flagged pixel is found
		# and make sure the distance to it is <`neighbors_required`
		else
			hi = lo + 1
			while hi <= l && !m[hi]
				hi += 1
			end
			hi -= 1
			if hi-lo < neighbors_required
				σ²[lo:hi, :] .= Inf
				if verbose; println("masked isolated pixels $lo:$hi") end
				affected_pixels!(affected, lo:hi)
			end
			lo = hi + 1
		end
	end
	return affected
end
function mask_isolated_pixels!(d::Data; kwargs...)
	mask_isolated_pixels!(d.var_s; kwargs...)
	return mask_isolated_pixels!(d.var; verbose=false, kwargs...)
end


"""
	process!(d; λ_thres=4000., min_snr=8, kwargs...)

Doing all of the data preprocessing, including continuum-fitting and masking of infinite points, anomalously snappy or high pixels, and low SNR edges.
"""
function process!(d; λ_thres::Real=4000., min_snr::Real=8, kwargs...)
	
	flat_normalize!(d)

	## masking
	# mask_low_pixels_all_times!(d; padding=2)
	bad_inst = mask_infinite_pixels!(d; padding=1)
	bad_snap = mask_bad_pixel!(d; padding=1)  # thres from 4-11 seems good
	bad_edge = mask_bad_edges!(d; min_snr=min_snr)
	bad_high = mask_high_pixels!(d; padding=1)
	bad_isol = mask_isolated_pixels!(d)

	## why was each pixel masked
	filter_bads(bad) = [i for i in bad if !(i in bad_edge)]
	bad_inst = filter_bads(bad_inst)
	bad_high = filter_bads(bad_high)
	bad_snap = filter_bads(bad_snap)
	bad_isol = filter_bads(bad_isol)

	## continuum normalizing
	# λ_thres = 4000 # is there likely to even be a continuum (neid index order 23+)
	# λ_thres = 6200 # where neid blaze correction starts to break down (neid index order 77+)
	red_enough = minimum(d.log_λ_obs) > log(λ_thres)
	# enough_points = (sum(isinf.(d.var)) / length(d.var)) < 0.5
	enough_points = true
	if (red_enough && enough_points)
		println("normalizing")
		w = continuum_normalize!(d; kwargs...)
	else
		w = nothing
	end

	return bad_inst, bad_high, bad_snap, bad_edge, bad_isol
end
