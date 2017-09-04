

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

"""
Discovery git repository for the current path
"""
function discover_repo(path)
    if isdir("$(path)/.git")
        return LibGit2.GitRepo(path)
    else
        parts = split(path, "/")
        if length(parts) > 1
            return discover_repo(join(parts[1:end-1], "/"))
        else
            return nothing
        end
    end
end

"""
Re-organize the current directory by moving session related files to their dedicated directories, e.g.
w7_11_1.edf -> session01/w7_11_1.edf
If files are under revision control, use git mv instead and opens the default editor asking for a commit message.
"""
function process()
    repo = discover_repo(pwd())
    repo_path = LibGit2.path(repo)*"/"
    files = split(chomp(readstring(`find . -name "*_settings.txt"`)))
    index_updated = false
    for f in files
        mm = match(r"([a-zA-Z]*)([0-9])_([0-9])_([0-9])",f)
        animal,month,day,session = mm.captures
        session = parse(Int64,session)
        # find all files belonging to this session
        sfiles = split(chomp(readstring(`find $(pwd()) -name "$(animal)$(month)_$(day)_$(session)*"`)))
        session_name = @sprintf "session%02d" session
        mkpath(session_name)
        for sf in sfiles
            if repo != nothing
                qq = replace(sf, repo_path, "") # get the path relative to the repository
                ss = LibGit2.status(repo, qq)
            else
                ss = Nullable{UInt32}()
            end
            if isnull(ss) # not under revision control
                mv(sf, "$(session_name)/$sf")
            else # under revision control
                vv = replace(sf, pwd()*"/", "")
                run(`git mv $(vv) $(session_name)/$vv`)
                index_updated = true
            end
        end
    end
    if index_updated
        run(`git commit`)
    end
end
