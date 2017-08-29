abstract type RawData end

mutable struct HighpassData{T1<:Real, T2<:Real,T3<:Real} <: RawData
    data::Array{T1,1}
    channel::Int64
    sampling_rate::T2
    filter_coefs::FilterCoefficients
    filter_name::String
    cutoff::T3
end

function save_data(X::HighpassData, session::String)
    _array = div(X.channel,32) + 1
    ch = X.channel - (_array-1)*32
    chs = @sprintf "channel%03d" ch
    ahs = @sprintf "array%02d" _array
    _pth = joinpath(session,ahs, chs)
    mkpath(_pth)
    fname = joinpath(_pth, "highpassdata.mat")
    mat_dict = Dict()
    for k in fieldnames(X)
        if k == :filter_coefficients
            bb = getfield(X,bb)
            mat_dict[string(k)] = Dict("z" => bb.z, "p" => bb.p, "k" => bb.k)
        else
            mat_dict[string(k)] = getfield(X, k)
        end
    end
    MAT.matwrite(fname, Dict("data" => mat_dict))
end

function load_data(::Type{HighpassData}, fname::String)
    mat_dict = MAT.matread(fname)
    _data = mat_dict["data"]
    fn = _data["filter_name"]
    bb = _data["filter_coefs"]
    ff = ZeroPoleGain(bb["z"], bb["p"], bb["k"])
    HighpassData(_data["data"], _data["channel"], _data["sampling_rate"],
                 ff, fn, _data["cutoff"])
end
