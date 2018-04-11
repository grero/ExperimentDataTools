module ExperimentDataTools
using ProgressMeter
using Eyelink
using Stimulus
using FileIO
using MAT
using DataFrames
using SpikeSorter
using Spiketrains
using DSP
using LFPTools
using RippleTools
using DataFrames
using Glob
using MAT
using LegacyFileReaders
using HDF5
import Base.parse
import LFPTools.align_lfp
using DataProcessingHierarchyTools
const DPHT = DataProcessingHierarchyTools
import DataProcessingHierarchyTools: filename, level

include("$(Pkg.dir("LFPTools"))/src/plots.jl")
include("types.jl")
include("sessions.jl")
include("behaviour.jl")
#include("spiketrains.jl")
include("spikesorting.jl")

export NPTData, HighpassData, LowpassData, OldTrials

const levels = ["days", "day", "session", "array", "channel", "cell"]

"""
Get the spike times from `wf` that falls witin the session boundaries of `eyelinkdata`.
"""
function get_session_spiketimes(wf::SpikeWaveforms, session_markers::Dict, session_marker_timestamps::Dict)
    get_session_spiketimes(wf.timestamps, session_markers,session_marker_timestamps)
end

function get_session_spiketimes(spiketimes::Array{Float64,1}, session_markers::Dict, session_marker_timestamps::Dict)
    #identity strobe markers for sessions
    session_timestamps = Dict{Int64, Array{Float64,1}}()
    for k in keys(session_markers)
        session_timestamps[k] = filter(t->session_marker_timestamps[k][1] < t < session_marker_timestamps[k][end], spiketimes) - session_marker_timestamps[k][1]
    end
    session_timestamps
end

function get_session_markers()
    return get_session_markers(get_markers()...)
end

function get_markers()
    dframe = readtable("event_markers.txt";eltypes=[String, Float64])
    markers = [string(m) for m in dframe[:markers]]
    return markers, Array(dframe[:timestamps])
end

function get_session_markers(pl2_file::String)
    markers, marker_timestamps = PlexonTools.extract_markers(pl2_file)
    get_session_markers(markers, marker_timestamps)
end

function get_session_markers(markers, marker_timestamps)
    session_markers = Dict{Int64, Array{String,1}}()
    session_marker_timestamps = Dict{Int64, Array{Float64,1}}()
    #FIXME: This is a complete hack. Figure out why we sometimes get 11111111 strobes
    sidx = find(m->m[1:3] == "110", markers)
    #END FIXME
    push!(sidx, length(markers)+1)
    for i in 1:length(sidx)-1
        k = parse(Int64,markers[sidx[i]][3:end],2)
        session_markers[k] = markers[sidx[i]:sidx[i+1]-1]
        session_marker_timestamps[k] = marker_timestamps[sidx[i]:sidx[i+1]-1]
    end
    session_markers, session_marker_timestamps
end

function get_session_spiketimes(wf::SpikeWaveforms, cwd=pwd())
    markers, marker_timestamps = PlexonTools.extract_markers()
    get_session_spiketimes(wf, markers, marker_timestamps)
end

function get_session_spiketimes(ff::File{format"MAT"}, cwd=pwd())
    data = matread(ff.filename)
    get_session_spiketimes(ff.filename, data)
end

function get_session_spiketimes(unit_name::String, unit_data::Dict, cwd=pwd())
    waveform = unit_data["waveform"]
    spikeidx = unit_data["timestamps"]
    sampling_rate = unit_data["sampling_rate"]
    markers, marker_timestamps = get_session_markers()
    session_timestamps = get_session_spiketimes(spikeidx./sampling_rate, markers, marker_timestamps)
    for (session,timestamps) in session_timestamps
        sn = @sprintf "session%02d" session
        mm = match(r"g([0-9]*)c([0-9]*)", unit_name)
        channel = parse(Int64,mm[1])
        _array = div(channel,32) + 1
        channel -= (_array-1)*32
        an = @sprintf "array%02d" _array
        cell = parse(Int64, mm[2])
        chn = @sprintf "channel%03d" channel
        cn = @sprintf("cell%02d", cell)
        dirn = joinpath(sn,an,chn,cn)
        if !isdir(dirn)
            mkpath(dirn)
        end
        fn = joinpath(sn,an,chn,cn,"unit.mat")
        MAT.matwrite(fn, Dict("timestamps" => timestamps*sampling_rate,
                              "waveform" => waveform,
                              "sampling_rate" => sampling_rate))
    end
end

"""
Transfer data from `src` to `dest`
"""
function get_data(src::String,dest=pwd())
    #transfer pl2 file
    pl2_files = split(readchomp(`find $src -name "*.pl2"`),"\n")
    for f in pl2_files
        cp($f, "$(dest)/$f")
        #extract markers
    end
    edf_files = split(readchomp(`find $src -name "*.edf"`), "\n")
    for f in edf_files
        #transfer
    end
end

function HighpassData(X::Array{Float64,1}, channel::Int64, sampling_rate::Float64, cutoff::Float64=300.0, filter_method=Butterworth, filter_order=4)
    ff = digitalfilter(Highpass(cutoff;fs=sampling_rate),filter_method(filter_order))
    Y = filtfilt(ff, X)
    #clunky way of getting the filter name
    filter_name = convert(String, split(string(filter_method), ".")[end])
    filter_name = "$(filter_name)($(filter_order))"
    HighpassData(Y, channel, sampling_rate, ff, filter_name, filter_order, cutoff, 10000.0)
end

function recompute!(H::HighpassData, data::Array{Float64,1}, channel::Int64)
    H.data = filtfilt(H.filter_coefs, data)
    H.channel = channel
    nothing
end

function get_triggers(rfile::File{format"NSHR"})
    #extract triggers
    dd, bn = splitdir(rfile.filename)
    if isempty(dd)
        dd = "."
    end
    marker_file = "$(dd)/event_markers.csv"
    if isfile(marker_file)
        _ddf = readtable(marker_file;eltypes=[String, Float64])
        words = Array(_ddf[:words])
        timestamps = Array(_ddf[:timestamps])
    else
        strobes, timestamps = RippleTools.extract_markers(rfile.filename)
        words = RippleTools.parse_strobe.(strobes)
        writetable(marker_file, DataFrame(words=words, timestamps=timestamps))
    end
    words, timestamps
end

function Base.parse(::Type{Stimulus.NewTrial}, rfile::File{format"NSHR"})
    words, timestamps = get_triggers(rfile)
    trials = parse(Stimulus.NewTrial, words, timestamps)
    trials
end

function get_session_starts()
    fname = filename(Trials)
    session_start = Dict()
    if isfile(fname)
        _ddf = readtable(fname;eltypes=[String, Float64])
        for r in eachrow(_ddf)
            w = r[:words]
            if w[1:2] == "11"
                sid = parse(Int64,w[3:end], 2)
                session_start[sid] = r[:timestamps]
            end
        end
    end
    # if we did not record the start of the first session, simply set it to 0
    if !(1 in keys(session_start)) && (2 in keys(session_start))
        session_start[1] = 0.0
    end
    session_start
end

function process_rawdata(rfile::File{format"NSHR"}, channels=1:128, fs=30_000)
    #get the trial structure
    session_start = get_session_starts()
    _dd, bn = splitdir(rfile.filename)
    if isempty(_dd)
        _dd = "."
    end
    session_name, ~ = splitext(bn)
    println("Processing data...")
    open(rfile.filename, "r") do ff
        dd = RippleTools.DataPacket(ff)
        @showprogress 1 "Processing channels... " for ch in channels
            #separate into sessions
            idx1 = 1
            sessions = collect(keys(session_start))
            sort!(sessions)
            for session in sessions
                session_name = @sprintf "session%02d" session
                _start = max(session_start[session],1.0)
                _end = get(session_start, session+1, size(dd.data,2)/fs)
                idx1 = round(Int64, _start*fs)
                idx2 = round(Int64, _end*fs)
                # lowpass data
                filepath = "$(getpath(session_name, ch))/$(filename(LowpassData))"
                if !isfile(filepath)
                    ldata,ff = LFPTools.lowpass_filter(float(dd.data[ch,idx1:idx2]),0.1, 250.0,fs)
                    lldata = LowpassData(ldata, ch, 1000.0, ff, "Butterworth", 4, 0.1, 250.0)
                    save_data(lldata, session_name)
                end
                filepath = "$(getpath(session_name, ch))/$(filename(HighpassData))"
                if !isfile(filepath)
                # highpass data
                    hdata, ffh = LFPTools.bandpass_filter(float(dd.data[ch, idx1:idx2]), 300.0, 10_000.0)
                    hhdata = HighpassData(hdata, ch, fs, ffh, "Butterworth", 4, 300.0, 10_000.0)
                    save_data(hhdata, session_name)
                end
            end
        end
    end
end

"""
Process the current directory, looking for data of type `T`.
"""
function process_dir(::Type{T}, dir=pwd();redo=false, save=true) where T <: Any
    cd(dir) do
    end
end

function process_level(::Type{T}, dir=pwd();kvs...) where T <: Any
    target_level = level(T)
    process_level(target_level, dir;kvs...)
end

function process_level(target_level::String, dir=pwd();kvs...)
    # get the current level
    this_level = level(dir)
    this_idx = findfirst(l->this_level==l, levels)
    target_idx = findfirst(l->target_level==l, levels)
    for lidx in [this_idx, target_idx]
        if !(0 < lidx <= length(levels))
            throw(ArgumentError("Unknown level"))
        end
    end
    pl = ["."]
    append!(pl, [".." for i in 1:(this_idx - target_idx)])
    dirstring = joinpath(pl...)
end

"""
Get the name of the requested level
"""
function get_level_name(target_level::String, dir=pwd())
    this_level = level(dir)
    this_idx = findfirst(l->this_level==l, levels)
    target_idx = findfirst(l->target_level==l, levels)
    i = this_idx
    cwd = dir
    pp = ""
    while i >= target_idx
        cwd, pp = splitdir(cwd)
        i -= 1
    end
    pp
end

function process_dirs(::Type{T}, dirs::Vector{String}, args...;kvs...) where T <: NPTData
    pp = cd(dirs[1]) do
        T(args...;kvs...)
    end
    @showprogress 1 "Processing dirs..." for d in dirs[2:end]
        _pp = cd(d) do
            T(args...;kvs...)
        end
        pp = hcat(pp, _pp)
    end
    return pp
end

"""
Convert old data to the new format. Basically, old data were split into chunks, and all channels for a particular chunk was stored int eh same file. 
"""
function process_old_data(channels::AbstractVector{Int64})
    files = glob("highpass/*highpass.*")
    sort!(files)
    data = LegacyFileReaders.load(File(format"NPTD", files[1]))
    sampling_rate = data.header.samplingrate
    if data.header.transpose
        nchannels = size(data.data,2)
    else
        nchannels = size(data.data,1)
    end
    filter_coefs = digitalfilter(Bandpass(250.0, 10000.0;fs=sampling_rate),Butterworth(4))
    @showprogress 1 "Processing channels..." for ch in intersect(channels,1:nchannels)
        hdata = Array{eltype(data.data),1}(0)
        for f in files
            data = LegacyFileReaders.load(File(format"NPTD", f))
            if data.header.transpose
                append!(hdata, data.data[:,ch])
            else
                append!(hdata, data.data[ch,:])
            end
        end
        H = HighpassData(hdata, ch,sampling_rate, filter_coefs, "Butterworth", 4, 250.0, 10000.0) 
        save_data(H, ".")
    end
end

function process_old_data()
    files = glob("*highpass.*")
end

end#module
