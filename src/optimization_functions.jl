_loss(tel::AbstractMatrix, star::AbstractMatrix, rv::AbstractMatrix, tfd::TFData) =
    sum((((tel .* (star + rv)) - tfd.flux) .^ 2) ./ tfd.var)
loss(tfo::TFOutput, tfd) = _loss(tfo.tel, tfo.star, tfo.rv, tfd)
loss(tfm::TFModel, tfd) = _loss(tel_model(tfm), star_model(tfm), rv_model(tfm), tfd)
loss_tel(tfo::TFOutput, tfm::TFModel, tfd) = _loss(tel_model(tfm), tfo.star, tfo.rv, tfd) + tel_prior(tfm)
loss_star(tfo::TFOutput, tfm::TFModel, tfd) = _loss(tfo.tel, star_model(tfm), tfo.rv, tfd) + star_prior(tfm)
loss_telstar(tfo::TFOutput, tfm::TFModel, tfd) = _loss(tel_model(tfm), star_model(tfm), tfo.rv, tfd) + star_prior(tfm) + tel_prior(tfm)
loss_rv(tfo::TFOutput, tfm::TFModel, tfd) = _loss(tfo.tel, tfo.star, rv_model(tfm), tfd)
function loss_funcs(tfo::TFOutput, tfm::TFModel, tfd::TFData)
    l() = loss(tfo, tfd)
    l_tel() = loss_tel(tfo, tfm, tfd)
    l_star() = loss_star(tfo, tfm, tfd)
    l_telstar() = loss_telstar(tfo, tfm, tfd)
    l_rv() = loss_rv(tfo, tfm, tfd)
    return l, l_tel, l_star, l_telstar, l_rv
end


struct TFOptimSubWorkspace
    θ::Flux.Params
    obj::OnceDifferentiable
    opt::Optim.FirstOrderOptimizer
    optstate::Optim.AbstractOptimizerState
    p0::Vector
    function TFOptimSubWorkspace(θ::Flux.Params, loss::Function)
        _, _, _, p0, obj = opt_funcs(loss, θ)
        opt = LBFGS()
        # initial_state(method::LBFGS, ...) doesn't use the options for anything
        return TFOptimSubWorkspace(θ, obj, opt, Optim.initial_state(opt, Optim.Options(), obj, p0), p0)
    end
    function TFOptimSubWorkspace(tfsm::TFSubmodel, loss::Function, only_s::Bool)
        if only_s
            θ = Flux.params(tfsm.lm.s)
        else
            θ = Flux.params(tfsm.lm.M, tfsm.lm.s, tfsm.lm.μ)
        end
        return TFOptimSubWorkspace(θ, loss)
    end
    function TFOptimSubWorkspace(tfsm1::TFSubmodel, tfsm2::TFSubmodel, loss::Function, only_s::Bool)
        if only_s
            θ = Flux.params(tfsm1.lm.s, tfsm2.lm.s)
        else
            θ = Flux.params(tfsm1.lm.M, tfsm1.lm.s, tfsm1.lm.μ, tfsm2.lm.M, tfsm2.lm.s, tfsm2.lm.μ)
        end
        return TFOptimSubWorkspace(θ, loss)
    end
    function TFOptimSubWorkspace(θ, obj, opt, optstate, p0)
        len = 0
        for i in 1:length(θ)
            len += length(θ[i])
        end
        @assert len == length(p0)
        new(θ, obj, opt, optstate, p0)
    end
end

abstract type TFOptimWorkspace end

struct TFWorkspace <: TFOptimWorkspace
    tel::TFOptimSubWorkspace
    star::TFOptimSubWorkspace
    rv::TFOptimSubWorkspace
    tfm::TFModel
    tfo::TFOutput
    tfd::TFData
    function TFWorkspace(tfm::TFModel, tfo::TFOutput, tfd::TFData; return_loss_f::Bool=false, only_s::Bool=false)
        loss, loss_tel, loss_star, _, loss_rv = loss_funcs(tfo, tfm, tfd)
        tel = TFOptimSubWorkspace(tfm.tel, loss_tel, only_s)
        star = TFOptimSubWorkspace(tfm.star, loss_star, only_s)
        rv = TFOptimSubWorkspace(tfm.rv, loss_rv, true)
        tfow = TFWorkspace(tel, star, rv, tfm, tfo, tfd)
        if return_loss_f
            return tfow, loss
        else
            return tfow
        end
    end
    TFWorkspace(tfm::TFModel, tfd::TFData, inds::AbstractVecOrMat; kwargs...) =
        TFWorkspace(tfm(inds), tfd(inds); kwargs...)
    TFWorkspace(tfm::TFModel, tfd::TFData; kwargs...) =
        TFWorkspace(tfm, TFOutput(tfm), tfd; kwargs...)
    function TFWorkspace(tel, star, rv, tfm, tfo, tfd)
        @assert length(tel.θ) == length(star.θ)
        @assert (length(tel.θ) == 1) || (length(tel.θ) == 3)
        @assert length(rv.θ) == 1
        new(tel, star, rv, tfm, tfo, tfd)
    end
end

struct TFWorkspaceTelStar <: TFOptimWorkspace
    telstar::TFOptimSubWorkspace
    rv::TFOptimSubWorkspace
    tfm::TFModel
    tfo::TFOutput
    tfd::TFData
    function TFWorkspaceTelStar(tfm::TFModel, tfo::TFOutput, tfd::TFData; return_loss_f::Bool=false, only_s::Bool=false)
        loss, _, _, loss_telstar, loss_rv = loss_funcs(tfo, tfm, tfd)
        telstar = TFOptimSubWorkspace(tfm.tel, tfm.star, loss_telstar, only_s)
        rv = TFOptimSubWorkspace(tfm.rv, loss_rv, true)
        tfow = TFWorkspaceTelStar(telstar, rv, tfm, tfo, tfd)
        if return_loss_f
            return tfow, loss
        else
            return tfow
        end
    end
    TFWorkspaceTelStar(tfm::TFModel, tfd::TFData, inds::AbstractVecOrMat; kwargs...) =
        TFWorkspaceTelStar(tfm(inds), tfd(inds); kwargs...)
    TFWorkspaceTelStar(tfm::TFModel, tfd::TFData; kwargs...) =
        TFWorkspaceTelStar(tfm, TFOutput(tfm), tfd; kwargs...)
    function TFWorkspaceTelStar(telstar, rv, tfm, tfo, tfd)
        @assert (length(telstar.θ) == 2) || (length(telstar.θ) == 6)
        @assert length(rv.θ) == 1
        new(telstar, rv, tfm, tfo, tfd)
    end
end

function _Flux_optimize!(θ::Flux.Params, obj::OnceDifferentiable, p0::Vector,
    opt::Optim.FirstOrderOptimizer, optstate::Optim.AbstractOptimizerState,
    options::Optim.Options)

    # Optim.optimize(obj, p0, LBFGS(); options)
    Optim.optimize(obj, p0, opt, options, optstate)
    copyto!(p0, θ)
end
_Flux_optimize!(tfosw::TFOptimSubWorkspace, options) =
    _Flux_optimize!(tfosw.θ, tfosw.obj, tfosw.p0, tfosw.opt, tfosw.optstate, options)

function train_TFModel!(tfow::TFWorkspace; options::Optim.Options=Optim.Options(iterations=10, f_tol=1e-3, g_tol=1e5))
    # optimize star
    _Flux_optimize!(tfow.star, options)
    tfow.tfo.star[:, :] = star_model(tfow.tfm)

    # optimize RVs
    tfow.tfm.rv.lm.M[:] = calc_doppler_component_RVSKL(tfow.tfm.star.λ, tfow.tfm.star.lm.μ)
    _Flux_optimize!(tfow.rv, options)
    tfow.tfo.rv[:, :] = rv_model(tfow.tfm)

    # optimize tellurics
    _Flux_optimize!(tfow.tel, options)
    tfow.tfo.tel[:, :] = tel_model(tfow.tfm)
end

function train_TFModel!(tfow::TFWorkspaceTelStar; options::Optim.Options=Optim.Options(iterations=10, f_tol=1e-3, g_tol=1e5))
    # optimize tellurics and star
    _Flux_optimize!(tfow.telstar, options)
    tfow.tfo.star[:, :] = star_model(tfow.tfm)
    tfow.tfo.tel[:, :] = tel_model(tfow.tfm)

    # optimize RVs
    tfow.tfm.rv.lm.M[:] = calc_doppler_component_RVSKL(tfow.tfm.star.λ, tfow.tfm.star.lm.μ)
    _Flux_optimize!(tfow.rv, options)
    tfow.tfo.rv[:, :] = rv_model(tfow.tfm)
end

function train_TFModel!(tfow::TFOptimWorkspace, n::Int; options::Optim.Options=Optim.Options(iterations=10, f_tol=1e-3, g_tol=1e5))
    for i in 1:n
        train_TFModel!(tfow; options=options)
    end
end
