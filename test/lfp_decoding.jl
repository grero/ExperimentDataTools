using PyPlot
using Colors

include("$(homedir())/Documents/programming/julia/Decoding/src/Decoding.jl")
function test()
    cd("/Users/roger/Documents/research/monkey/newWorkingMemory/Wiesel/20171128/session01") do
        trials = ExperimentDataTools.parse(Stimulus.NewTrial, File(format"NSHR", "w11_28_1.nev"));
        ctrials = Stimulus.getTrialType(trials, :reward_on)
        rtime = [round(Int64, (trial.trial_start + trial.response_on)*30_000) for trial in ctrials]
        target_label = [trial.stimulus[1].locidx for trial in ctrials]

        cd("array01")

        Z = zeros(701,595,31)
        t = -0.1:0.001:0.6
        for ch in 1:31
            ss = @sprintf "channel%03d" ch
            cd(ss) do
                 ldata = ExperimentDataTools.load_data(ExperimentDataTools.LowpassData, "lowpass.mat")
                 Z[:,:,ch],t = LFPTools.align_lfp(ldata.data, div.(rtime,30),1000, (-100,600));
            end
        end

        Q = zeros(size(Z,1)-50, size(Z,2),size(Z,3))

        for k in 1:size(Z,3)
            for j in 1:size(Z,2)
                for i in 1:size(Z,1)-50
                    Q[i,j,k] = sum(Z[i:i+49, j, k])
                end
            end
        end

        q_rr,test_idx_rr = Decoding.decode(Decoding.MV.MulticlassLDA{Float64}, permutedims(Q,[3,2,1]),trial_label,100)

        AA = squeeze(sum(q_rr .== trial_label[test_idx_rr],1),1)./size(q_rr,1);
        μ = mean(AA,2)[:]
        σ = std(AA,2)[:]

        plt[:figure]()
        plot(t[1:end-50],μ) 
        fill_between(t[1:end-50], μ-σ, μ+σ;color=(0.8, 0.8, 1.0),zorder=-1)
        plt[:xlabel]("Time from response_cue [s]")
        plt[:ylabel]("Decoding performance")
        cc = parse(Colorant, "tomato")
        plt[:axvline](0.2;label="Minimum reward time",color=(cc.r, cc.g, cc.b))
        plt[:axvline](0.0;label="Response on")
        plt[:legend]()

        plt[:savefig]("/tmp/w11_28_lfp_target_deccoding_response_ch_1-32.png";dpi=300)
    end
end

function test2(alignment::Symbol=:target_on)
    dirs = convert(Vector{String}, split(readchomp(`find . -name "channel*" -type d`), '\n'))
    sort!(dirs)
    trials = ExperimentDataTools.Trials()
    ctrials = Stimulus.getTrialType(trials, :reward_on)
    trial_label = Stimulus.get_label.(ctrials)
    if alignment == :target_on
        rtime = [round(Int64, (trial.trial_start + trial.stimulus[1].onset)*30_000) for trial in ctrials]
        window = (-300, 1000)
    else
        rtime = [round(Int64, (trial.trial_start + trial.response_on)*30_000) for trial in ctrials]
        window = (-100, 600)
    end
    al = ExperimentDataTools.process_dirs(ExperimentDataTools.AlignedLFP, dirs, div.(rtime,30), window)
    Z = ExperimentDataTools.splitsets(al)
    Q = zeros(size(Z,1)-50, size(Z,2),size(Z,3))

    for k in 1:size(Z,3)
        for j in 1:size(Z,2)
            for i in 1:size(Z,1)-50
                Q[i,j,k] = sum(Z[i:i+49, j, k])
            end
        end
    end

    q_rr,test_idx_rr = Decoding.decode(Decoding.MV.MulticlassLDA{Float64}, permutedims(Q,[3,2,1]),trial_label,100)

    AA = squeeze(sum(q_rr .== trial_label[test_idx_rr],1),1)./size(q_rr,1);
    μ = mean(AA,2)[:]
    σ = std(AA,2)[:]

    fig = plt[:figure]()
    plot(al.t[1:end-50],μ) 
    fill_between(al.t[1:end-50], μ-σ, μ+σ;color=(0.8, 0.8, 1.0),zorder=-1)
    plt[:xlabel]("Time from response_cue [s]")
    plt[:ylabel]("Decoding performance")
    cc = parse(Colorant, "tomato")
    if alignment == :target_on
        plt[:axvline](0.0;label="Target on")
        plt[:axvline](0.35;label="Target off", color=(cc.r, cc.g, cc.b))
    else
        plt[:axvline](0.2;label="Minimum reward time",color=(cc.r, cc.g, cc.b))
        plt[:axvline](0.0;label="Response on")
    end
    plt[:legend]()
    fig
end
