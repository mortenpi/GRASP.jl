function read_isodata(filename)
    lines = readlines(filename)
    values = map(1:div(length(lines), 2)) do i
        # First line: drop the : and return the rest as the "key"
        m = match(r"^(.*):", lines[2i-1])
        isnothing(m) && return nothing
        key = strip(m[1])
        value = if key == "NNNP"
            parse(Int, lines[2i])
        else
            parse(Float64, lines[2i])
        end
        string(key) => value
    end |> Dict{Union{String,Symbol},Real}
    # We will now generate a few normalized entries in values, based on the existing values
    isodata_rename_copy!(values, "Atomic number", :Z, filename=filename)
    isodata_rename_copy!(values, "Mass number (integer)", :A, filename=filename)
    isodata_rename_copy!(values, "Mass of nucleus (in amu)", :m, filename=filename)
    isodata_rename_copy!(values, "Fermi distribution parameter a", :fermi_a, filename=filename)
    isodata_rename_copy!(values, "Fermi distribution parameter c", :fermi_c, filename=filename)
    return values
end

function isodata_rename_copy!(d::Dict, key::AbstractString, newkey::Symbol; filename)
    if !haskey(d, key)
        @warn "isodata is missing key '$key'" filename
        return
    end
    if haskey(d, newkey)
        @warn "isodata file already contains key '$newkey'" filename
        return
    end
    d[newkey] = d[key]
    return
end
