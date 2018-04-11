struct NewSpiketrains end

level(::Type{NewSpiketrains}) = "cell"

abstract type AbstractSortedSpikes end

struct HMMSortedSpikes <: AbstractSortedSpikes
    channel::Int64
    timestamps::Vector{Float64}
    waveforms::Vector{Float64}
end

function level(::Type{T}) where T<: AbstractSortedSpikes
   "channel"
end

filename(::Type{HMMSortedSpikes}) = "hmmsortedspikes.jld"

function load(::Type{T}, args...) where T <: AbstractSortedSpikes
    dir = process_level(T)
    qq = cd(dir) do
        fname = filename(T)
        if isfile(fname)
            spike_model = T(fname)
        else
           # code to initiate spike sorting for this channel 
           X = load(HighpassData)
           spike_model = fit(HMMSpikeSorter.HMMSpikingModel, X.data, args...)
        end
        spike_model
    end
    qq
end

function HMMSortedSpikes(fname::String)
    JLD.load(fname)
end

HMMSortedSpikes() = load(HMMSortedSpikes)



