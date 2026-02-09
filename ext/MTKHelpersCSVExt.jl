module MTKHelpersCSVExt

# function __init__()
#     @info "MTKHelpers: loading MTKHelpersCSVExt"
# end

isdefined(Base, :get_extension) ? (using CSV) : (using ..CSV)
using MTKHelpers

using ComponentArrays: ComponentArrays as CA
using ComponentArrays # when parsing Axis needs to be defined
using DataFrames
using CSV
using Chain

"""
    write_csv_cv(con::IO, df0::DataFrame)
    read_csv_cv(con::IO)::DataFrame

Write/Read a DataFrame that contains columns of type ComponentArrays   
to a CSV file.
Each ComponentArray column is replaced by length(cv) columns,
and information on the axis is generated as a comment in front of the CSV.

When re-reading the CSV, the respective ComponentArrays are recreated 
usingn this information and are replaced for the read plain columns.

Make sure to not duplicate names. If there is a ComponentArray column
`u0` containing three Floats, the following three columns are added 
(and should not exist before): `u0_1, u0_2, u0_3`.
"""
function MTKHelpers.write_csv_cv(con::IO, df0::DataFrame)
    cv_info = get_cv_info(df0)
    df = spread_cv!(copy(df0); cv_info)
    print_cv_info(con, cv_info)
    CSV.write(con, df; append = true, header = true)
    cv_info
end

"""
    get_cv_info(df::DataFrame)

Get the information on axes of all ComponentArray columns in DataFrame.
Returns a DataFrame with columns
- col: Symbol of the column name
- axis: String representation of the axis    
"""
function get_cv_info(df::DataFrame)
    cv_info = DataFrame(col = Symbol[], axis = String[])
    io = IOBuffer()
    map(propertynames(df)) do col
        val1 = df[1, col]
        if val1 isa CA.ComponentArray
            ax = CA.getaxes(val1)
            show(io, ax)
            s1 = String(take!(io))
            push!(cv_info, (col, s1))
        end
    end
    return (cv_info)
end

function spread_cv!(df::DataFrame; cv_info = get_cv_info(df))
    for col in cv_info.col
        spread_cv!(df, col)
    end
    df
end

function spread_cv!(df::DataFrame, col::Symbol)
    l = length(df[1, col])
    newcols = ["$(col)_$i" for i in 1:l]
    transform!(df, col => ByRow(collect) => newcols)
    select!(df, Not(col))
end

"""
print the information of ComponentVector columns to IO connection con
"""
function print_cv_info(con::IO, cv_info::DataFrame)
    for r in eachrow(cv_info)
        print(con, "# $(r.col)=$(r.axis)\n")
    end
    return (nothing)
end

function MTKHelpers.read_csv_cv(con)
    cv_info = read_cv_info(con)
    dfs = CSV.read(con, DataFrame; comment = "#", ignoreemptyrows = true)
    df = gather_cv!(copy(dfs); cv_info)
end

function read_cv_info(con::IO)
    cv_info = DataFrame(col = Symbol[], axis = String[])
    while (true)
        mark(con)
        s = readline(con)
        _rm = match(r"#\s*([[:alnum:]]+)=(\(Axis\(.*)$", s)
        isnothing(_rm) && break
        push!(cv_info, (Symbol(_rm[1]), _rm[2]))
    end
    reset(con) # reset to before reading the non-matching line
    return (cv_info)
end

function eval_axis_dict(cv_info)
    _parse_axis = (axis) -> begin
        ex = Meta.parse(axis)
        valex = validate_axes_expression(ex)
        !valex.valid && error(valex.msg)
        eval(ex)
    end
    _df = transform(cv_info, :axis => ByRow(_parse_axis) => :ax)
    Dict(_df.col .=> _df.ax)
end

function validate_axes_expression(ex)
    ex.head == :tuple ||
        return (; valid = false,
            msg = "Expected expression of tuple, but was $(ex.head)")
    for i in 1:length(ex.args)
        (ex.args[i].head == :call &&
         ex.args[i].args[1] == :Axis) ||
            return (;
                valid = false,
                msg = "Expected tuple entry $i to be an Axis, but was not.")
    end
    return (valid = true, msg = "")
end

function gather_cv!(df::DataFrame; cv_info::DataFrame)
    axes_dict = eval_axis_dict(cv_info)
    col = cv_info.col[1]
    for col in cv_info.col
        ax = get(axes_dict, col, nothing)
        if !isnothing(ax)
            gather_cv!(df, col, ax)
        end
    end
    return (df)
end

function gather_cv!(df::DataFrame, col::Symbol, ax)
    oldcols = filter(s -> occursin(Regex("^$(col)_\\d+\$"), s), names(df))
    # make sure to supply columns in ascending order
    ind = map(cn -> match(r"(\d+)$", cn)[1], oldcols)
    oldcols_ord = oldcols[sortperm(ind)]
    _make_cv = (comps...) -> CA.ComponentVector(collect(comps), ax)
    df2 = @chain df begin
        transform!(oldcols_ord => ByRow(_make_cv) => col)
        select!(Not(oldcols_ord))
    end
end

end # module
