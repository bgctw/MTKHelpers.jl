"""
Construct a `Dict{Symbol => Num}` for all properties in `sys`.
All Symbols are prefixed with `<string_sys>₊`
"""
function get_system_symbol_dict(sys::AbstractSystem, string_sys=string(sys.name))
    Dict(Symbol(string_sys*"₊"*string(p)) => getproperty(sys,p) for p in propertynames(sys))
end


"""
    system_num_dict(d, sys::AbstractSystem, string_sys=sys.name)
    system_num_dict(d, symbol_dict::AbstractDict)

Create a Dictionary Num=>value from symbolic Dictionary or ComponentVector.

Omit pairs where no Num was found.
"""
function system_num_dict(d, sys::AbstractSystem, string_sys=string(sys.name))
    system_num_dict(d, get_system_symbol_dict(sys, string_sys))
end
function system_num_dict(d, symbol_dict::AbstractDict)
    Dict([get(symbol_dict, k, missing) => v for (k,v) in d if !ismissing(get(symbol_dict, k, missing))])
end
function system_num_dict(ca::CA.ComponentVector, symbol_dict::AbstractDict)
    system_num_dict(Dict(propertynames(ca) .=> values(ca)), symbol_dict)
end

