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
