@testset "_update_cv_top" begin
    p1 = ComponentVector(a=1, b=[2,3,4], c=5)
    ptmp = p1[(:c,:b)] .* 10
    pup = MTKHelpers._update_cv_top(p1, ptmp)
    MTKHelpers._get_axis(pup) == MTKHelpers._get_axis(p1)
    @test typeof(getdata(pup)) == typeof(getdata(p1))
    @test pup.a == p1.a
    @test pup.b == ptmp.b
    @test pup.c == ptmp.c
    #
    # empty updater - nothing to update
    ptmp = ComponentVector()
    pup = MTKHelpers._update_cv_top(p1, ptmp)
    @test pup == p1
    @test !(pup === p1) # is not a reference to p1
end;

benchmark_update_cv_top = () -> begin
    ptmp = p1.*2 
    is_updated = [true,true,false]
    tmp = MTKHelpers._update_cv_top(p1, ptmp, is_updated)
    #using BenchmarkTools
    #@btime MTKHelpers._update_cv_top($p1, $ptmp, $is_updated) 
    #
    tmp = MTKHelpers._update_cv_top(u0, u0, is_updated)

end

# @testset "_update_cv" begin
#     cv = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     s = ComponentVector(a=(a1=11,a2=2,a3=3),b=20)
#     cv1 = _update_cv(cv, s)
#     @test cv1 == s
#     s = ComponentVector(b=21)
#     cv1 = _update_cv(cv, s)
#     @test (cv1.a == cv.a) & (cv1.b == s.b)
#     s = ComponentVector(a=(a1=11,))
#     cv1 = _update_cv(cv, s)
#     @test (cv1.a.a1 == 11) & (cv1.a.a2 == cv.a.a2) & (cv1.b == cv.b)
#     s = ComponentVector(a=(aTypo=11,))
#     @test_throws ErrorException _update_cv(cv, s)
#     s = ComponentVector(a=11:13)
#     cv1 = _update_cv(cv, s)
#     @test cv1.a.a3 == 13
#     s = ComponentVector(a=1)
#     @test_throws ErrorException _update_cv(cv, s) # wrong length ComponentVector
#     s = ComponentVector(b=20:21) # wrong length flat axis
#     @test_throws ErrorException _update_cv(cv, s) 
#     ax = first(getaxes(cv))
#     CA.last_index(ax)
#     lastindex(ax)
# end;

# @testset "_get_index_axis" begin
#     # TODO update src and test in ComponentArrays.jl bgctw branch
#     cv = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     #s = ComponentVector(b=1)
#     s = ComponentVector(a=1)
#     ax = first(getaxes(s))
#     @test_throws ErrorException _get_index_axis(cv, ax)
# end;

# @testset "getindex_axis ComponentVector" begin
#     cv = ComponentVector(a=(a1=100,a2=(a21=210, a22=220)), b=2, c = (c1=reshape(1:4,(2,2)),))
#     cr = cv[first(axes(cv))]
#     @test  cr == cv
#     cr.a.a1 = 3100;  @test cv.a.a1 == 100 # does not modify original cv
#     # extract a single nested component
#     # construct Axis by template ComponentVector
#     cs = ComponentVector(a=(a2=(a22=1,),)) # the "," is important to make it a tuple
#     cr = cv[first(axes(cs))]
#     @test axes(cr) == axes(cs)
#     @test cr.a.a2.a22 == cv.a.a2.a22
#     # extract shaped component
#     @test cv[Axis(:c)].c == cv.c
#     # extract two upper-level components with structure
#     cr = cv[Axis(:a,:c)] # :b not include
#     @test keys(cr) == (:a,:c)
#     @test first(axes(cr.a)) == first(axes(cv.a))
#     @test cr.a == cv.a
#     @test cr.a.a1 == cv.a.a1
#     # TODO: documentation do not forget comma for tuple
#     cv[first(axes(ComponentVector(a=(a1=1,))))] # the "," is important to make it a tuple
#     # wrongly extracts full a including a1 including a2, because there is no a1 in axis
#     cv[first(axes(ComponentVector(a=(a1=1))))] 
#     #
#     # ugly error on not matching template
#     cs = ComponentVector(a=(aTypo=1,)) # the "," is important to make it a tuple
#     @test_throws ErrorException cv[first(axes(cs))] # has aTypo in error message

#     cv = ComponentVector(a=(a1=100,a2=(a21=210, a22=reshape(1:4,(2,2))),a3=300), b=3, c=4)
#     # extracting some of the nested a-components and the b component
#     # note that I do not specify the structure of a2
#     cv_indextemplate = ComponentVector(a=(a2=1, a3=1), b=1)
#     ax = first(getaxes(cv_indextemplate)) # todo specify without template
#     cv_sub = cv[ax]

# end

# @testset "_get_index_axis AxisArray" begin
#     u1 = ComponentVector(L = 10.0,)
#     p1 = ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
#     popt = ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
#     u1a = AxisArray(collect(u1); row = keys(u1))
#     p1a = AxisArray(collect(p1); row = keys(p1))
#     vcat(u1a, p1a) # does not preserve type - returns plain vector
#     # hence also ComponentArray cannot do anything
#     u1av = ComponentVector(u1a, (Axis(keys(u1)),))
#     p1av = ComponentVector(p1a, (Axis(keys(p1)),))
#     tmp = vcat(u1av, p1av) # skips underlying type
#     typeof(tmp)
#     # try with different AbstractArray than AxisArray
# end;

# @testset "_get_index_axis NamedArray" begin
#     u1 = ComponentVector(L = 10.0,)
#     p1 = ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
#     popt = ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
#     u1a = NamedArray(collect(u1), (collect(keys(u1)),))
#     p1a = NamedArray(collect(p1), (collect(keys(p1)),))
#     vcat(u1a, p1a) # does preserve type 
#     # hence also ComponentArray cannot do anything
#     u1av = ComponentVector(u1a, (Axis(keys(u1)),))
#     p1av = ComponentVector(p1a, (Axis(keys(p1)),))
#     tmp = vcat(u1av, p1av) # skips underlying type
#     tmp2 = p1av.k_L
#     _get_index_axis(p1av, Axis(:k_L, :k_R))
#     # try with different AbstractArray than AxisArray
# end;

@testset "_labels" begin
    x = first(getaxes(ComponentArray(a = 1, b = [1, 2])))
    @test MTKHelpers._labels(x) == [".a", ".b[1]", ".b[2]"]
    x = first(getaxes(ComponentArray(c = (a = 1, b = [1, 2]))))
    @test MTKHelpers._labels(x) == [".c.a", ".c.b[1]", ".c.b[2]"]
    x = (a = 1, b = [1, 2])
    @test MTKHelpers._labels(x) == [".a", ".b[1]", ".b[2]"]
    x = first(getaxes(ComponentArray(c = [(a = [1, 2],), (a = [2, 3],), (a = [3, 4],)])))
    @test MTKHelpers._labels(x) == [
        ".c[1].a[1]",
        ".c[1].a[2]",
        ".c[2].a[1]",
        ".c[2].a[2]",
        ".c[3].a[1]",
        ".c[3].a[2]",
    ]
    x = first(getaxes(ComponentArray(c = (b = [1 2; 5 6]))))
    @test MTKHelpers._labels(x) == [".c[1,1]", ".c[2,1]", ".c[1,2]", ".c[2,2]"]
    nt2 = (a = 5,
        b = [(a = (a = 20, b = 1), b = 0), (a = (a = 33, b = 1), b = 0)],
        c = (a = (a = 2, b = [1, 2]), b = [1.0 2.0; 5 6]))
    ca2 = ComponentArray(nt2)
    x = first(getaxes(ca2))
    lab = MTKHelpers._labels(x)
    @test lab == [
        ".a",
        ".b[1].a.a",
        ".b[1].a.b",
        ".b[1].b",
        ".b[2].a.a",
        ".b[2].a.b",
        ".b[2].b",
        ".c.a.a",
        ".c.a.b[1]",
        ".c.a.b[2]",
        ".c.b[1,1]",
        ".c.b[2,1]",
        ".c.b[1,2]",
        ".c.b[2,2]",
    ]
    #MTKHelpers.labels_noprefix(x)
    syms = MTKHelpers._ax_symbols_tuple(x)
    @test syms isa NTuple
    @test syms[2] == Symbol("b[1]₊a₊a")
end;

# @testset "subaxis" begin
#     cv = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     CP.subaxis(cv, :a)    
# end;
