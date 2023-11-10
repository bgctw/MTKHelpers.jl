"""
    get_system_symbol_dict(sys::AbstractSystem, string_sys::String=string(nameof(sys)))
    get_system_symbol_dict(systems...)

Construct a `Dict{Symbol => Num}` for all properties in `sys`.
All Symbols are prefixed with `<string_sys>₊`

The second variant merges the dictionaries obtained from several systems.
"""
function get_system_symbol_dict(sys::AbstractSystem,
    string_sys::String = string(nameof(sys)))
    prefix = isempty(string_sys) ? "" : string_sys * "₊"
    Dict(Symbol(prefix * string(p)) => getproperty(sys, p) for
         p in propertynames(sys))
end
function get_system_symbol_dict(systems...)
    dicts = map(systems) do sys
        get_system_symbol_dict(sys)
    end
    merge(dicts...)
end

"""
    system_num_dict(d, sys::AbstractSystem, string_sys=nameof(sys))
    system_num_dict(d, symbol_dict::AbstractDict)
    system_num_dict(d, systems::NTuple) 

Create a Dictionary Num=>value from symbolic Dictionary or ComponentVector.

Omit pairs where no Num was found.

In the third variant, a tuple of AbstractSystems can be specified to 
replace Nums of several (sub-)systems in the dictionary.
"""
function system_num_dict(d, sys::AbstractSystem, string_sys = string(nameof(sys)))
    system_num_dict(d, get_system_symbol_dict(sys, string_sys))
end
function system_num_dict(d, symbol_dict::AbstractDict)
    Dict([get(symbol_dict, k, missing) => v for
          (k, v) in d if !ismissing(get(symbol_dict, k, missing))])
end
function system_num_dict(ca::CA.ComponentVector, symbol_dict::AbstractDict)
    system_num_dict(Dict(propertynames(ca) .=> values(ca)), symbol_dict)
end
function system_num_dict(d, systems::NTuple)
    system_num_dict(d, get_system_symbol_dict(systems...))
end
