module ExperimentDataTools
using SpikeSorter
using PlexonTools
using FileIO
using MAT

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
    if isfile("event_markers.txt")
        dframe = readtable("event_markers.txt";eltypes=[String, Float64])
        markers = [string(m) for m in dframe[:markers]]
        return get_session_markers(markers, dframe[:timestamps])
    else
        pl2_file = split(readchomp(`find . -name "*.pl2"`))
        if isempty(pl2_file)
            throw(ArgumentError("No pl2 file found"))
        end
        return get_session_markers(convert(String,first(pl2_file)))
    end
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
    get_session_spiketimes(data)
end

function get_session_spiketimes(unit_name::String, unit_data::Dict, cwd=pwd())
    waveform = unit_data["waveform"]
    spikeidx = unit_data["timestamps"]
    sampling_rate = unit_data["sampling_rate"]
    markers, marker_timestamps = get_session_markers()
    session_timestamps = get_session_spiketimes(spikeidx./sampling_rate, markers, marker_timestamps)
    for (session,timestamps) in session_timestamps
        sn = @sprintf "session%02d" session
        if !isdir(sn)
            mkdir(sn)
        end
        mm = match(r"g([0-9]*)c([0-9]*)", unit_name)
        channel = parse(Int64,mm[1])
        cell = parse(Int64, mm[2])
        chn = @sprintf "channel%03d" channel
        dirn = joinpath(sn,chn)
        if !isdir(dirn)
            mkdir(dirn)
        end
        cn = @sprintf("cell%02d", cell)
        dirn = joinpath(sn,chn, cn)
        if !isdir(dirn)
            mkdir(dirn)
        end
        fn = joinpath(sn,chn,cn,"unit.mat")
        MAT.matwrite(fn, Dict("timestamps" => timestamps,
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

end#module
