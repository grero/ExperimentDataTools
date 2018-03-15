module ExperimentDataPlots
using ExperimentDataTools
const EDT = ExperimentDataTools
using PyPlot

function plot_data(::Type{EDT.HighpassData}, window::AbstractVector{Int64})
    hdata = HighpassData()
    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:plot](hdata.data[window])
    fig
end

function plot_data(::Type{EDT.SpikeTrainData},window::AbstractVector{Float64})
    sptrain = EDT.SpikeTrainData()
    tt = sptrain.data.timestamps
    idx1 = findfirst(t->t>window[1], tt)
    idx2 = findlast(t->t<window[end], tt)
    _tt = tt[idx1:idx2]
    y = [-1.0 1.0].*ones(idx2-idx1+1,2)
    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ax[:plot]([_tt _tt]',y';color="k") 
    fig
end

function plot_data(::Type{EDT.HighpassData}, window::AbstractVector{Int64}, dirs::Vector{String})
    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    offset = 0.0
    for d in dirs
        cd(d) do
            hdata = EDT.HighpassData()
            _data = hdata.data[window]
            l,h = extrema(_data)
            offset += -l
            ax[:plot](_data + offset)
            offset += h 
        end
    end
    ax[:spines]["top"][:set_visible](false)
    ax[:spines]["right"][:set_visible](false)
    ax[:spines]["left"][:set_visible](false)
    ax[:set_yticks]([])
    ax[:set_yticklabels]([])
    fig
end

end #module
