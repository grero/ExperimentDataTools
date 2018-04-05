abstract type RawData <: DPHData end

import Base.zero, Base.hcat, Base.append!

struct Trials end
DataProcessingHierarchyTools.level(::Type{Trials}) = "session"
DataProcessingHierarchyTools.filename(::Type{Trials}) = "event_markers.csv"

function Trials()
    ndir = process_level(Trials)
    trials = cd(ndir) do
        fname = filename(Trials)
        if isfile(fname)
            _ddf = readtable(fname;eltypes=[String, Float64])
            words = Array(_ddf[:words])
            timestamps = Array(_ddf[:timestamps])
            trials = parse(Stimulus.NewTrial, words, timestamps)
        else
            #look for an nev file
            fname = convert(String, (first(split(readchomp(`find . -name "*nev"`),'\n'))))
            if !isempty(fname)
                trials = parse(Stimulus.NewTrial, File(format"NSHR", fname))
            else
                trials = NewTrial[]
            end
        end
        trials
    end
    return trials
end

struct OldTrials <: DPHData
    data::Vector{Stimulus.Trial}
    setid::AbstractVector{Int64}
end
DataProcessingHierarchyTools.level(::Type{OldTrials}) = "session"
DataProcessingHierarchyTools.filename(::Type{OldTrials}) = "event_data.mat"

function OldTrials()
    trials = Stimulus.loadTrialInfo("event_data.mat")
    OldTrials(trials, fill(1, length(trials)))
end

struct BroadbandData <: DPHData end
level(::Type{BroadbandData}) = "day"

function BroadbandData()
    # check for ripple files
    datafiles = split(readchomp(`find . -name "*.ns*" -d 1`), '\n')
    if isempty(datafiles)
        datafiles = split(readchomp(`find . -name "*.pl2*" -d 1`), '\n')
        datafile = convert(String, datafiles[1])

    else
        datafile = convert(String, datafiles[1])

    end
    if isempty(datafiles)
    end
end

mutable struct HighpassData{T1<:Real, T2<:Real} <: RawData
    data::Array{T1,1}
    channel::Int64
    sampling_rate::T2
    filter_coefs::FilterCoefficients
    filter_name::String
    filter_order::Int64
    low_freq::Float64
    high_freq::Float64
end

filename(X::HighpassData) = "highpass.mat"
DPHT.filename(::Type{HighpassData}) = "highpass.mat"
DPHT.level(::Type{HighpassData}) = "channel"
matname(X::HighpassData) = "highpassdata"
matname(::Type{HighpassData}) = "highpassdata"

HighpassData{T2<:Real, T3<:Real}(sampling_rate::T2, filter_coefs::FilterCoefficients, filter_name::String, cutoff::T3) = HighpassData(Float64[], 0, sampling_rate, filter_coefs, filter_name, cutoff)

function HighpassData{T2<:Real}(sampling_rate::T2, low_freq::Float64, high_freq::Float64, filter_method::Function, filter_order::Int64)
    filter_coefs = digitalfilter(Highpass(cutoff;fs=sampling_rate),filter_method(filter_order))
    filter_name = convert(String, split(string(filter_method), ".")[end])
    filter_name = "$(filter_name)($(filter_order))"
    HighpassData(sampling_rate, filter_coefs, filter_name, cutoff)
end

function Base.zero(::Type{HighpassData{T1, T2}}) where T1 <: Real where T2 <: Real
    HighpassData(T1[], 0, zero(T2), ZeroPoleGain([0.0], [0.0], 0.0),"", 0, 0.0, 0.0)
end

function HighpassData()
    fname = filename(HighpassData)
    if isfile(fname)
        hh = load_data(HighpassData, fname)
    else
        hh = zero(HighpassData)
    end
    hh
end

mutable struct LowpassData{T1<:Real, T2<:Real} <: RawData
    data::Array{T1,1}
    channel::Int64
    sampling_rate::T2
    filter_coefs::FilterCoefficients
    filter_name::String
    filter_order::Int64
    low_freq::Float64
    high_freq::Float64
end

function zero(::Type{LowpassData{T1, T2}}) where T1 <: Real where T2 <: Real
    LowpassData(T1[], 0, zero(T2), ZeroPoleGain([0.0], [0.0], 0.0),"", 0, 0.0, 0.0)
end

zero(::Type{LowpassData}) = zero(LowpassData{Float64, Float64})

filename(X::LowpassData) = "lowpass.mat"
filename(::Type{LowpassData}) = "lowpass.mat"
matname(X::LowpassData) = "lowpassdata"
matname(X::Type{LowpassData}) = "lowpassdata"
level(::Type{LowpassData}) = "channel"

struct AlignedLFP
    X::Matrix{Float64}
    setid::AbstractVector{Int64}
    t::AbstractVector{Float64}
end

function Base.hcat(al1::AlignedLFP, al2::AlignedLFP)
    if !(al1.t == al2.t)
        throw(ArgumentError("Both arguments should be aligned to the same time"))
    end
    X = hcat(al1.X, al2.X)
    setid = vcat(al1.setid, al2.setid + al1.setid[end] )
    AlignedLFP(X, setid, al1.t)
end

function Base.append!(al1::AlignedLFP, al2::AlignedLFP)
    if !(al1.t == al2.t)
        throw(ArgumentError("Both arguments should be aligned to the same time"))
    end
end

function splitsets(al::AlignedLFP)
    cc = countmap(al.setid)
    nsets = unique(values(cc))
    if length(nsets) != 1
        throw(ArgumentError("Unequal number of trials detected"))
    end
    X = reshape(al.X, size(al.X,1), first(nsets), length(cc))
    return X
end

function AlignedLFP(ldata::LowpassData, ttime::Vector{Int64}, window::Tuple{Int64,Int64})
    X,t = LFPTools.align_lfp(ldata.data, ttime, ldata.sampling_rate, window)
    AlignedLFP(X,fill(1,size(X,2)), t)
end

function AlignedLFP(start_time::Vector{Int64},window::Tuple{Int64,Int64})
    ldata = LowpassData()
    AlignedLFP(ldata, start_time, window)
end

function getpath(session::String, channel::Int)
    _array = div(channel-1,32) + 1
    chs = @sprintf "channel%03d" channel
    ahs = @sprintf "array%02d" _array
    _pth = joinpath(session,ahs, chs)
    _pth
end

function save_data(X::T, session::String) where T <: RawData
    _pth = getpath(session, X.channel)
    mkpath(_pth)
    fname = joinpath(_pth, filename(X))
    mat_dict = Dict()
    for k in fieldnames(X)
        if k == :filter_coefficients
            bb = getfield(X,bb)
            mat_dict[string(k)] = Dict("z" => bb.z, "p" => bb.p, "k" => bb.k)
        else
            mat_dict[string(k)] = getfield(X, k)
        end
    end
    MAT.matwrite(fname, Dict(matname(X) => Dict("data" => mat_dict)))
end

function load_data(::Type{T}, fname::String) where T <: RawData
    if ishdf5(fname)
        X = h5open(fname,"r") do ff
            _fn = read(ff, "highpassdata/data/filter_name")
            #convoluted way of reading a string from hdf5
            readbuf = IOBuffer(reinterpret(UInt8, _fn[:]))
            fn = String(read(readbuf))
            fn = replace(fn, "\0","")
            fo = read(ff, "highpassdata/data/filter_order")[1]
            bb = read(ff, "highpassdata/data/filter_coefs")
            low_freq = read(ff, "highpassdata/data/low_freq")[1]
            high_freq = read(ff, "highpassdata/data/high_freq")[1]
            channel = read(ff, "highpassdata/data/channel")[1]
            sampling_rate = read(ff, "highpassdata/data/sampling_rate")[1]
            bb["k"] = bb["k"][1]
            bb["z"] = [r.data[1] + r.data[2]*1im for r in bb["z"]]
            bb["p"] = [r.data[1] + r.data[2]*1im for r in bb["p"]]
            _filter = ZeroPoleGain(bb["z"], bb["p"], bb["k"])
            if ismmappable(ff["highpassdata/data/data"])
                _data = readmmap(ff["highpassdata/data/data"])
            else
                _data = read(ff,"highpassdata/data/data")
            end
            T(_data,channel, sampling_rate, _filter, fn, fo, low_freq, high_freq )
        end
    else
        mat_dict = MAT.matread(fname)
        _data = mat_dict[matname(T)]["data"]
        fn = _data["filter_name"]
        fo = _data["filter_order"]
        bb = _data["filter_coefs"]
        ff = ZeroPoleGain(bb["z"], bb["p"], bb["k"])
        X = T(_data["data"], _data["channel"], _data["sampling_rate"],
                     ff, fn,fo, _data["low_freq"], _data["high_freq"])
    end
    X
end

function LowpassData()
    fname = filename(LowpassData)
    if isfile(fname)
        return load_data(LowpassData, fname)
    end
    return zero(LowpassData)
end

struct ChannelConfig <: DPHT.DPHData
    config::Dict{String, UnitRange{Int64}}
end

DPHT.level(::Type{ChannelConfig}) = "subject"
DPHT.filename(::Type{ChannelConfig}) = "channel_config.csv"

function ChannelConfig()
    _config = cd(DPHT.process_level(level(ChannelConfig))) do
        ExperimentDataTools.channel_config(filename(ChannelConfig))
    end
    ChannelConfig(_config)
end
