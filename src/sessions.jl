

"""
Break up the continuous data according to the session markers
"""
function get_session_idx{T<:AbstractArray{Int64,1}}(full_idx::T, session_markers::Array{String,1}, session_timestamps::Array{Float64,1}, sampling_rate::Real)
    session_data = Dict{Int64,T}()
    sidx = find(m->m[1:3] == "110", session_markers)
    for i in 1:length(sidx)
        k = parse(Int64,session_markers[sidx[i]][3:end],2)
        t0 = round(Int64,session_timestamps[sidx[i]]*sampling_rate)
        if i < length(sidx)
            t1 = round(Int64,session_timestamps[sidx[i+1]]*sampling_rate)
        else
            t1 = length(full_idx)
        end
        session_data[k] = full_idx[t0:t1]
    end
    return session_data
end

"""
Split `data` into sessions and save as typeof `T`.
"""
function save_session_data{T<:RawData}(H::T, data::Array{Float64,1}, session_markers::Array{String,1}, session_timestamps::Array{Float64,1}, channel, sampling_rate,cutoff, filter_method, filter_order)
    idx = get_session_idx(1:length(data), session_markers, session_timestamps,sampling_rate)
    for (k,v) in idx
        recompute!(H, data[v], channel)
        session = @sprintf "session%02d" k
        save_data(H, session)
    end
end
