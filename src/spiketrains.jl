#struct SaccadeAlignedSpiketrain <: NPTData
#    aligned_spikes::Spiketrains.AlignedSpiketrains
#end

struct SpikeTrain
    timestamps::Vector{Float64}
    spikeshape::Vector{Float64}
end

struct SpikeTrainData <: DPHData
    data::SpikeTrain
    setid::Vector{Int64}
end

DPHT.level(::Type{SpikeTrainData}) = "cell"
DPHT.filename(::Type{SpikeTrainData}) = "spiketrain.mat"

function SpikeTrainData()
    ss = MAT.matread(filename(SpikeTrainData))
    _data = SpikeTrain(ss["timestamps"][:], ss["spikeForm"][:])
    SpikeTrainData(_data, [1])
end
