function plot_sacacde_endpoints(eyetrials::Vector{Stimulus.NewEyelinkTrial}, screen_width=1920.0, screen_height=1200.0)
    rtrials = Stimulus.getTrialType(eyetrials, :response_on)
    rsaccades = Stimulus.get_response_saccade(rtrials)
    delay = [t.response_on for t in rtrials] - [t.stimulus[1].offset for t in rtrials]
    sacpos = [(s.end_x, s.end_y) for s in rsaccades]

    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    ss = ax[:scatter]([sc[1] for sc in sacpos], [settings["screen_height"] - sc[2] for sc in sacpos], 5.0, delay)
    ax[:set_xlim](0, screen_width-1)
    ax[:set_ylim](0, screen_height-1)
    plt[:colorbar](ss;ax=ax,label="Delay duration [s]")
    fig[:set_size_inches](7.8, 4.0)
    fig[:tight_layout]()
    fig 
end
