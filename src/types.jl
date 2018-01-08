abstract type RawData end

import Base.zero

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
filename(::Type{HighpassData}) = "highpass.mat"
matname(X::HighpassData) = "highpassdata"
matname(::Type{HighpassData}) = "highpassdata"

HighpassData{T2<:Real, T3<:Real}(sampling_rate::T2, filter_coefs::FilterCoefficients, filter_name::String, cutoff::T3) = HighpassData(Float64[], 0, sampling_rate, filter_coefs, filter_name, cutoff)

function HighpassData{T2<:Real}(sampling_rate::T2, low_freq::Float64, high_freq::Float64, filter_method::Function, filter_order::Int64)
    filter_coefs = digitalfilter(Highpass(cutoff;fs=sampling_rate),filter_method(filter_order))
    filter_name = convert(String, split(string(filter_method), ".")[end])
    filter_name = "$(filter_name)($(filter_order))"
    HighpassData(sampling_rate, filter_coefs, filter_name, cutoff)
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

function getpath(session::String, channel::Int) 
    _array = div(channel,32) + 1
    chs = @sprintf "channel%03d" channel 
    ahs = @sprintf "array%02d" _array
    _pth = joinpath(session,ahs, chs)
    _pth
end

function save_data(X::RawData, session::String)
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
    mat_dict = MAT.matread(fname)
    _data = mat_dict[matname(T)]["data"]
    fn = _data["filter_name"]
    fo = _data["filter_order"]
    bb = _data["filter_coefs"]
    ff = ZeroPoleGain(bb["z"], bb["p"], bb["k"])
    T(_data["data"], _data["channel"], _data["sampling_rate"],
                 ff, fn,fo, _data["low_freq"], _data["high_freq"])
end

function LowpassData()
    fname = filename(LowpassData)
    if isfile(fname)
        return load_data(LowpassData, fname)
    end
    return zero(LowpassData)
end
