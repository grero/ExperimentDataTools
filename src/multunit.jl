struct MultiUnitArgs <: DPHT.DPHDataArgs
    snr::Float64  # signal to noise ratio
end

struct MultiUnit <: DPHT.DPHData
    timestamps::Vector{Int64}
    sampling_rate::Float64
    args::MultiUnitArgs
end

DPHT.level(::Type{MultiUnit}) = "channel"
DPHT.filename(::Type{MultiUnit}) = "multiunit.mat"
DPHT.datatype(::Type{MultiUnitArgs}) = MultiUnit

function DPHT.load(::Type{MultiUnit}, fname=DPHT.filename(MultiUnit))
    Q = MAT.matread(fname)
    args = MultiUnitArgs(Q["snr"][1])
    MultiUnit(Q["timestamps"], Q["sampling_rate"][1], args)
end

function DPHT.save(X::MultiUnit)
    fname = DPHT.filename(X.args)
    Q = Dict()
    Q["snr"] = X.args.snr
    Q["timestamps"] = X.timestamps
    Q["sampling_rate"] = X.sampling_rate
    MAT.matwrite(fname, Q)
    nothing
end

function MultiUnit(args::MultiUnitArgs;do_save=true, force_redo=false, kvs...)
    fname = DPHT.filename(args)
    redo = !isfile(fname) || force_redo
    if !redo
        X = DPHT.load(args)
    else
        HH = HighpassData()
        data = HH.data
        μ,σ = SpikeExtraction.get_threshold(data, args.snr)
        #get 1.5 ms duration in data points
        nn = div(round(Int64,HH.sampling_rate),1000)
        pre_spike = div(nn,2)
        post_spike = nn
        timestamps,spikeshapes = SpikeExtraction.extract_spikes(data;μ=μ,σ=σ,nq=(pre_spike, post_spike),θ=args.snr)
        X = MultiUnit(timestamps, HH.sampling_rate, args)
        if do_save
            DPHT.save(X)
        end
    end
    return X
end

