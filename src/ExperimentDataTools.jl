module ExperimentDataTools
using SpikeSorter
using FileIO
using MAT
using DataFrames
using DSP

include("types.jl")

export HighpassData

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
    dframe = readtable("event_markers.txt";eltypes=[String, Float64])
    markers = [string(m) for m in dframe[:markers]]
    return get_session_markers(markers, dframe[:timestamps])
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
    HighpassData(Y, channel, sampling_rate, ff, filter_name, cutoff)
end

function recompute!(H::HighpassData, data::Array{Float64,1}, channel::Int64)
    H.data = filtfilt(H.filter_coefs, data)
    H.channel = channel
    nothing
end

end#module
