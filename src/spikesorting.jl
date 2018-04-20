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

struct SpikeTemplateArgs <: DPHT.DPHDataArgs
end

struct SpikeTemplates <: DPHT.DPHData
    spikeshapes::Array{Float64,3}
    cinv::Matrix{Float64}
    pp::Vector{Float64}
    channel::Int64
    args::SpikeTemplateArgs
end

SpikeTemplates() = SpikeTemplates(Array{Float64}(0,0,0), Matrix{Float64}(0,0), Vector{Float64}(0), -1, SpikeTemplateArgs())

DPHT.level(::Type{SpikeTemplates}) = "channel"
DPHT.filename(::Type{SpikeTemplates}) = "spike_templates.mat"
DPHT.datatype(::Type{SpikeTemplateArgs}) = SpikeTemplates


function SpikeTemplates(args::SpikeTemplateArgs;kvs...)
    fnames = glob("spike_templates.hdf5")
    if isempty(fnames)
        return SpikeTemplates()
    end
    fname = first(fnames)
    channel = parse(Int64, DPHT.get_numbers(DPHT.get_level_name(DPHT.level(SpikeTemplates))))
    HDF5.h5open(fname, "r") do ff
        spikeshapes = read(ff, "spikeForms") 
        cinv = read(ff, "cinv")
        if ndims(cinv) == 1
            cinv = cinv[:,1:1]
        end
        pp = read(ff,"p") 
        SpikeTemplates(spikeshapes, cinv, pp, channel,args)
    end
end
