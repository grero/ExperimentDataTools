using JSON
using PyPlot
fig_save_path = "/Users/roger/Documents/research/monkey/training"

function plot_saccade_endpoints(eyetrials::Vector{Stimulus.NewEyelinkTrial}, screen_width=1920.0, screen_height=1200.0)
    rtrials = Stimulus.getTrialType(eyetrials, :response_on)
    rsaccades = Stimulus.get_response_saccade(rtrials)
    delay = [t.response_on for t in rtrials] - [t.stimulus[1].offset for t in rtrials]
    sacpos = [(s.end_x, s.end_y) for s in rsaccades]

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ss = ax[:scatter]([sc[1] for sc in sacpos], [screen_height - sc[2] for sc in sacpos], 5.0, delay)
    ax[:set_xlim](0, screen_width-1)
    ax[:set_ylim](0, screen_height-1)
    plt[:colorbar](ss;ax=ax,label="Delay duration [s]")
    fig[:set_size_inches](7.8, 4.0)
    fig[:tight_layout]()
    fig
end

function plot_saccade_endpoints(fname::String)
    subject_prefix = fname[1:1]
    if subject_prefix == "w"
        subject = "Wiesel"
    end
    figbase = replace(fname, ".edf", "_response_saccades_delay.pdf")
    figfname = "$(fig_save_path)/$(subject)/$(figbase)"
    eyedata = Eyelink.load(fname)
    eyetrials = Stimulus.parse(Stimulus.NewEyelinkTrial, eyedata.events)
    settings = JSON.parsefile(replace(fname, ".edf", "_settings.txt"))
    fig = plot_saccade_endpoints(eyetrials,settings["screen_width"], settings["screen_height"])
    fig[:savefig](figfname)
    fig
end

function get_reward_duration()
    files = split(readchomp(`find . -name "*results.txt"`), "\n")
    reward_duration = Float64[]
    delay = Float64[]
    for f in files
        results = readtable(f)
        ridx  = find(results[:reward_on] .> 0)
        append!(reward_duration, (results[:reward_off] - results[:reward_on])[ridx])
        append!(delay, (results[:response_on] - results[:target_off])[ridx])
    end
    reward_duration, delay
end

function plot_reward_duration_distr()
    reward_duration, delay = get_reward_duration()
    plot_reward_duration(reward_duration, delay)
end

function plot_reward_duration_distr(reward_duration::Vector{Float64}, delay::Vector{Float64})
    fig = plot_marginal_hist(reward_duration, delay)
    ax = fig[:axes][1]
    ax[:set_xlabel]("Reward duration [s]")
    ax[:set_ylabel]("Delay duration [s]")
    fig
end

function plot_marginal_hist(y1::Vector{Float64}, y2::Vector{Float64}) 
    fig = plt[:figure]()
    plot_marginal_hist(fig, y1, y2)
    fig
end

function plot_marginal_hist(fig, y1::Vector{Float64}, y2::Vector{Float64}) 
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histy1 = [left, bottom_h, width, 0.2]
    rect_histy2 = [left_h, bottom, 0.2, height]
    
    ax_scatter = fig[:add_axes](rect_scatter)
    ax_histy1 = fig[:add_axes](rect_histy1)
    ax_histy2 = fig[:add_axes](rect_histy2)
    
    h1 = KernelDensity.kde_lscv(y1)
    h2 = KernelDensity.kde_lscv(y2)
    ax_histy1[:plot](h1.x, h1.density)
    ax_histy2[:plot](h2.density, h2.x)
    ax_scatter[:scatter](y1,y2)
    ax_histy1[:set_xlim](extrema(h1.x))
    ax_histy1[:set_xticklabels]([])
    ax_scatter[:set_xlim](extrema(h1.x))
    ax_histy2[:set_ylim](extrema(h2.x))
    ax_histy2[:set_yticklabels]([])
    ax_scatter[:set_ylim](extrema(h2.x))
end
