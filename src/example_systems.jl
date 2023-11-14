function samplesystem_vec(; name, τ = 3.0, i=0.1, p = [1.1, 1.2, 1.3])
    n_comp = 2
    @parameters t
    D = Differential(t)
    @variables x(..)[1:n_comp] #dx(t)[1:2]  # observed dx now can be accessed
    #sts = @variables x[1:n_comp](t) 
    #ps = @parameters τ=τ p[1:n_comp]=p i=i       # parameters
    ps = @parameters τ=τ i=i p[1:3]=p 
    sts = [x(t)[i] for i in 1:n_comp]
    eq = [
        D(x(t)[1]) ~ i - p[1] * x(t)[1] + (p[2] - x(t)[1]^2) / τ, 
        D(x(t)[2]) ~ i - p[3] * x(t)[2], 
     ]
     #ODESystem(eq, t; name)
     #ODESystem(eq, t, sts, [τ, p[1], p[2], i]; name)
     sys = ODESystem(eq, t, sts, vcat(ps...); name)
     #sys = ODESystem(eq, t, sts, ps; name)
     return sys
end

function indices_of_nums(nums)
    op_syms = symbol_op.(nums)
    #vec_of_pairs = Vector{Pair{Symbol, Union{Int,UnitRange{Int}}}}()
    vec_of_pairs = Vector{Pair{Symbol, UnitRange{Int}}}()
    for (pos, sym) in enumerate(op_syms) 
        #Main.@infiltrate_main
        if length(vec_of_pairs) == 0 
            push!(vec_of_pairs, sym => pos:pos)
        else 
            current_pair = vec_of_pairs[end]
            if sym != first(current_pair)
                push!(vec_of_pairs, sym => pos:pos)
            else
                current_range = Base.last(current_pair)
                vec_of_pairs[end] = sym => first(current_range):pos
            end
        end
    end
    return vec_of_pairs
end

function axis_of_nums(nums)
    pos_nums = indices_of_nums(nums)
    pos_nums_scalar = [first(pr) => length(last(pr)) == 1 ? first(last(pr)) : last(pr) for pr in pos_nums]
    #CA._component_axis(first(axes(CA.ComponentArray(;pos_nums_scalar...))))
    first(getaxes(CA.ComponentArray(;pos_nums_scalar...)))
end


