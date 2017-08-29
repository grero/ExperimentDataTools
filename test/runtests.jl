using ExperimentDataTools
using Distributions
using DataFrames
using MAT
using DSP
using Base.Test

function create_data()
    sampling_rate = 40_000.0
    timestamps = cumsum(rand(Exponential(0.1),10000))
    nsessions = 1
    session_timestamps = rand(0.0:maximum(timestamps),nsessions)
    bb = 192
    session_markers = [bin(bb+i) for i in 1:nsessions]
    push!(session_markers, "00110000")
    push!(session_timestamps, timestamps[end] + 1.0)
    t = linspace(0,59,60)
    waveform = sin.(2pi*t/60)
    data = Dict("timestamps" => timestamps.*sampling_rate,
         "waveform" => waveform,
         "sampling_rate" => sampling_rate)
    return Dict("g003c01_spiketrains.mat"=> data),session_markers, session_timestamps
 end

srand(1234)
data, session_markers, session_timestamps = create_data()
dframe = DataFrame(markers=session_markers, timestamps=session_timestamps)
cwd = pwd()
dd = tempdir()
cd(dd)
writetable("event_markers.txt", dframe)
for (k,v) in data
    ExperimentDataTools.get_session_spiketimes(k,v)
end
new_data = MAT.matread("session01/array01/channel003/cell01/unit.mat")
sm,st = ExperimentDataTools.get_session_markers()

cd(cwd)

@test data["g003c01_spiketrains.mat"]["waveform"] == new_data["waveform"]
t1 = filter(t->t/40000.0>=session_timestamps[1], data["g003c01_spiketrains.mat"]["timestamps"])
t2 = new_data["timestamps"] + session_timestamps[1]*40_000
@test t1 ≈ t2

function test_highpass()
    cwd = pwd()
    dd = tempdir()
    cd(dd)
    srand(1234)
    X = zeros(100_000);
    for i in 2:length(X)
        X[i] = X[i-1] + randn()
    end
    H = HighpassData(X, 1, 40_000.0, 300.0, Butterworth, 4)
    ExperimentDataTools.save_data(H, "session01")
    H2 = ExperimentDataTools.load_data(ExperimentDataTools.HighpassData, "session01/array01/channel001/highpassdata.mat")
    @test H.channel == H2.channel
    @test H.sampling_rate ≈ H2.sampling_rate
    @test H.cutoff ≈ H2.cutoff
    @test H.data ≈ H2.data
    @test H.filter_name == H2.filter_name
    @test H.filter_coefs.p ≈ H2.filter_coefs.p
    @test H.filter_coefs.k ≈ H2.filter_coefs.k
    @test H.filter_coefs.z ≈ H2.filter_coefs.z
    cd(cwd)
end

test_highpass()
