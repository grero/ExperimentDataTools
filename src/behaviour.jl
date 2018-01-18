struct ExpResults <: NPTData
    data::DataFrame
    setid::AbstractVector{Int64}
end

Base.zero(::Type{ExpResults}) = ExpResults(DataFrame(), Int64[])

function Base.hcat(r1::ExpResults, r2::ExpResults)
    if names(r1.data) != names(r2.data)
        throw(ArgumentError("Column names do not match"))
    end
    data = [r1.data; r2.data]
    setid = [r1.setid; r1.setid[end] + r2.setid]
    return ExpResults(data, setid)
end

level(::Type{ExpResults}) = "session"


function ExpResults()
    fname = split(readchomp(`find . -d 1 -name "*results.txt"`), '\n')
    if !isempty(fname)
        data = readtable(convert(String, fname[1]))
        return ExpResults(data, fill(1, size(data,1)))
    end
    return ExpResults()
end
