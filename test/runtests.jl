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
    H2 = cd("session01/array01/channel001") do
        ExperimentDataTools.HighpassData()
    end
    #cleanup
    rm("session01/array01/channel001";recursive=true)
    @test H.channel == H2.channel
    @test H.sampling_rate ≈ H2.sampling_rate
    @test H.low_freq ≈ H2.low_freq
    @test H.high_freq ≈ H2.high_freq
    @test H.data ≈ H2.data
    @test H.filter_name == H2.filter_name
    @test H.filter_coefs.p ≈ H2.filter_coefs.p
    @test H.filter_coefs.k ≈ H2.filter_coefs.k
    @test H.filter_coefs.z ≈ H2.filter_coefs.z
    cd(cwd)
end

function test_reorganising_sessions()
    cwd = pwd()
    dd = tempdir()
    cd(dd) do
        repo_path = "test"
        repo = LibGit2.init(repo_path)
        files = ["w7_11_1.edf", "w7_11_1_settings.txt", "w7_11_1_results.txt"]
        for f in files
            touch("$(repo_path)/$f")
            LibGit2.add!(repo, f)
        end
        LibGit2.commit(repo, "Adds files")
        cd(repo_path) do
            ExperimentDataTools.process(;commit_message="Cleanup")
            for f in files
                nf = "session01/$f"
                @test isfile(nf)
            end
        end
        rm(repo_path;recursive=true)
    end
end

@testset "Utils" begin
    @testset "Cellnames" begin
        _config = Dict("vFEF" => 1:32, "dFEF" => 33:64)
        cellname = "dFEF_g30c1"
        array, channel, cc = ExperimentDataTools.get_flat_channel(cellname, _config)
        @test cc == 1
        @test array == 2
        @test channel == 62
        pth = ExperimentDataTools.get_cell_path(cellname, _config)
        @test pth == "array02/channel062/cell01"
    end
end


test_highpass()
test_reorganising_sessions()
