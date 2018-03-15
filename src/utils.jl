import Base.parse

function dissect_cell(cell::String)
	m = match(r"([A-Za-z0-9]*)[_]*g([0-9]*)c([0-9]*)",cell) 
	a,g,c = m.captures
	return a,parse(Int64,g),parse(Int64,c)
end

function parse(::Type{T},a::String) where  T <: Range{T2} where T2 <: Integer 
	m = match(r"([0-9]*)[-:]([0-9]*)",a)
    mc = m.captures
    isempty(mc) && return nothing
    return T(map(s->parse(T2,String(s)),mc)...)
end

function channel_config(fname::String)
	dd = CSV.read(fname)
    Q = Dict{String, UnitRange{Int64}}()
    for rr in CSV.eachrow(dd)
        v = parse(UnitRange{Int64}, convert(String,rr[:channels].value))
        k = convert(String, rr[:area].value)
        Q[k] = v
    end
    return Q
end

function get_flat_channel(cellname::String,config::Dict{String, T}) where T <: Range{T2} where T2 <: Integer
    aa,gg,cc = dissect_cell(cellname)
    channel = 0
    array = 0
    i = 0
    #sort by first channel 
    channel_ranges = collect(values(config))
    areas = collect(keys(config))
    sidx = sortperm(channel_ranges;by=v->v[1])
    for k in areas[sidx]
        v = config[k]
        i += 1
        if aa == k
            channel = v[gg]
            array = i
            break
        end
    end
    array, channel 
end

function get_cell_path(cellname, config)
    array,channel = get_flat_channel(cellname,config)
    arrayname = @sprintf "array%02d" array
    channelname = @sprintf "channel%03d" channel
    joinpath(arrayname, channelname)
end
