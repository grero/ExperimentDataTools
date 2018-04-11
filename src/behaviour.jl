struct ExpResults <: DPHT.DPHData
    data::DataFrame
    setid::AbstractVector{Int64}
end

Base.zero(::Type{ExpResults}) = ExpResults(DataFrame(), Int64[])

function Base.hcat(r1::ExpResults, r2::ExpResults)
    if names(r1.data) != names(r2.data)
        throw(ArgumentError("Column names do not match"))
    end
    data = [r1.data; r2.data]
    setid = [r1.setid; r1.setid[end] + r2.setid]
    return ExpResults(data, setid)
end

level(::Type{ExpResults}) = "session"


function ExpResults()
    fname = split(readchomp(`find . -d 1 -name "*results.txt"`), '\n')
    if !isempty(fname)
        data = readtable(convert(String, fname[1]))
        return ExpResults(data, fill(1, size(data,1)))
    end
    return ExpResults()
end

struct EyeTrials <: DPHT.DPHData
    data::Vector{Stimulus.EyeTraceTrial}
    setid::Vector{Int64}
end
filename(::Type{EyeTrials}) = "eyetrials.mat"
matname(::Type{EyeTrials}) = "et"

function get_value(data, ff, i)
    v = 0.0
    if !isempty(data[ff][i])
        v = data[ff][i]
    end
    v
end

function EyeTrials()
    fname = filename(EyeTrials)
    ff = MAT.matread(fname)
    data = ff["et"]["data"]["trials"]
    ntrials = length(data["start"])
    trials = Vector{Stimulus.EyeTraceTrial}(ntrials)
    for i in 1:ntrials
        trial_start = float(data["start"][i])
        prestim = get_value(data, "fixation_start", i)
        stimstart = 0.0
        trial_response = get_value(data,"response_cue",i)
        reward = get_value(data,"reward",i)
        reward_duration = 0.0
        failure = get_value(data,"failure",i)
        delay = get_value(data,"delay",i)
        gazex = data["gazex"][i][:]
        gazey = data["gazey"][i][:]
        pupil = data["pupil"][i][:]
        if !("gtime" in keys(data))
            gtime = 0:length(gazex)-1
        else
            gtime = data["gtime"][i][:]
        end
        _target = data["target"][i]
        if isa(_target, Dict)
            target = Stimulus.Target(_target["row"], _target["column"], _target["onset"])
        else
            target = Stimulus.Target(0.0, 0.0, NaN)
        end
        _distr = data["distractor"][i]
        if isa(_distr, Dict)
            distractors = [Stimulus.Distractor(_distr["row"], _distr["column"], _distr["onset"])]
        else
            distractors = [Stimulus.Distractor(0.0, 0.0, NaN)]
        end
        _saccade = data["saccade"][i]
        if isempty(_saccade)
            nsaccade = 0
        else
            nsaccade = length(_saccade["endx"])
        end
        saccade = Vector{Eyelink.Saccade}(nsaccade)
        for j in 1:length(saccade)
            endx = _saccade["endx"][j]
            startx = _saccade["startx"][j]
            endy = _saccade["endy"][j]
            starty = _saccade["starty"][j]
            onset = _saccade["onset"][j]
            offset = _saccade["offset"][j]
            saccade[j] = Eyelink.Saccade(onset, offset, startx, starty, endx, endy, i)
        end
        trial_end = data["end"][i]
        trials[i] = Stimulus.EyeTraceTrial(target, distractors, trial_start, trial_end, prestim,
                                           delay, stimstart, trial_response, reward, reward_duration,
                                           failure, gazex, gazey, pupil, gtime, saccade)
    end
    EyeTrials(trials, fill(1, length(trials)))
end
