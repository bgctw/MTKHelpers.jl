using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
#using StaticArrays: StaticArrays as SA

#include("test/testset_utils.jl") 
include("testset_utils.jl") # @testset_skip, get_pkg_version

@testset "_get_axis" begin
    p1 = CA.ComponentVector(CA.ComponentVector(x = [1, 2, 3]),
        par = CA.ComponentVector(a = 1, b = [2, 3, 4], c = 5))
    # axes returns a CombinedAxis, whld CA.getaxes returns ComponentAxis directly
    cax = first(axes(p1))
    res = MTKHelpers._get_axis(cax)
    @test res == MTKHelpers._get_axis(p1)
    # 
    u0 = (x = 2,)
    res = MTKHelpers._get_axis(keys(u0))
    @test res == MTKHelpers._get_axis(CA.ComponentVector(u0))
end;

@testset "attach_axis" begin
    c = (a = 2, b = [1, 2])
    x = CA.ComponentArray(a = 1, b = [2, 1, 4.0], c = c)
    x2 = x .* x'
    x2bc = x2[:b, :c]
    tmp = CA.getaxes(x2bc)
    # b has no axis, create one
    ax = CA.Axis(b1 = 1, b2 = 2, b3 = 3)
    res = MTKHelpers.attach_x_axis(x2bc, ax)
    @test res[:b1, :] == x2bc[1, :]
end;

@testset "_update_cv_top" begin
    p1 = CA.ComponentVector(a = 1, b = [2, 3, 4], c = 5)
    ptmp = p1[(:c, :b)] .* 10
    pup = MTKHelpers._update_cv_top(p1, ptmp)
    MTKHelpers._get_axis(pup) == MTKHelpers._get_axis(p1)
    @test typeof(CA.getdata(pup)) == typeof(CA.getdata(p1))
    @test pup.a == p1.a
    @test pup.b == ptmp.b
    @test pup.c == ptmp.c
    #
    # empty updater - nothing to update
    ptmp = CA.ComponentVector()
    pup = MTKHelpers._update_cv_top(p1, ptmp)
    @test pup == p1
    @test !(pup === p1) # is not a reference to p1
end;

benchmark_update_cv_top = () -> begin
    ptmp = p1 .* 2
    is_updated = [true, true, false]
    tmp = MTKHelpers._update_cv_top(p1, ptmp, is_updated)
    #using BenchmarkTools
    #@btime MTKHelpers._update_cv_top($p1, $ptmp, $is_updated) 
    #
    tmp = MTKHelpers._update_cv_top(u0, u0, is_updated)
end

# @testset "_update_cv" begin
#     cv = CA.ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     s = CA.ComponentVector(a=(a1=11,a2=2,a3=3),b=20)
#     cv1 = _update_cv(cv, s)
#     @test cv1 == s
#     s = CA.ComponentVector(b=21)
#     cv1 = _update_cv(cv, s)
#     @test (cv1.a == cv.a) & (cv1.b == s.b)
#     s = CA.ComponentVector(a=(a1=11,))
#     cv1 = _update_cv(cv, s)
#     @test (cv1.a.a1 == 11) & (cv1.a.a2 == cv.a.a2) & (cv1.b == cv.b)
#     s = CA.ComponentVector(a=(aTypo=11,))
#     @test_throws ErrorException _update_cv(cv, s)
#     s = CA.ComponentVector(a=11:13)
#     cv1 = _update_cv(cv, s)
#     @test cv1.a.a3 == 13
#     s = CA.ComponentVector(a=1)
#     @test_throws ErrorException _update_cv(cv, s) # wrong length CA.ComponentVector
#     s = CA.ComponentVector(b=20:21) # wrong length flat axis
#     @test_throws ErrorException _update_cv(cv, s) 
#     ax = first(CA.getaxes(cv))
#     CA.last_index(ax)
#     lastindex(ax)
# end;

# @testset "_get_index_axis" begin
#     # TODO update src and test in ComponentArrays.jl bgctw branch
#     cv = CA.ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     #s = CA.ComponentVector(b=1)
#     s = CA.ComponentVector(a=1)
#     ax = first(CA.getaxes(s))
#     @test_throws ErrorException _get_index_axis(cv, ax)
# end;

# @testset "getindex_axis CA.ComponentVector" begin
#     cv = CA.ComponentVector(a=(a1=100,a2=(a21=210, a22=220)), b=2, c = (c1=reshape(1:4,(2,2)),))
#     cr = cv[first(axes(cv))]
#     @test  cr == cv
#     cr.a.a1 = 3100;  @test cv.a.a1 == 100 # does not modify original cv
#     # extract a single nested component
#     # construct Axis by template CA.ComponentVector
#     cs = CA.ComponentVector(a=(a2=(a22=1,),)) # the "," is important to make it a tuple
#     cr = cv[first(axes(cs))]
#     @test axes(cr) == axes(cs)
#     @test cr.a.a2.a22 == cv.a.a2.a22
#     # extract shaped component
#     @test cv[CA.Axis(:c)].c == cv.c
#     # extract two upper-level components with structure
#     cr = cv[CA.Axis(:a,:c)] # :b not include
#     @test keys(cr) == (:a,:c)
#     @test first(axes(cr.a)) == first(axes(cv.a))
#     @test cr.a == cv.a
#     @test cr.a.a1 == cv.a.a1
#     # TODO: documentation do not forget comma for tuple
#     cv[first(axes(CA.ComponentVector(a=(a1=1,))))] # the "," is important to make it a tuple
#     # wrongly extracts full a including a1 including a2, because there is no a1 in axis
#     cv[first(axes(CA.ComponentVector(a=(a1=1))))] 
#     #
#     # ugly error on not matching template
#     cs = CA.ComponentVector(a=(aTypo=1,)) # the "," is important to make it a tuple
#     @test_throws ErrorException cv[first(axes(cs))] # has aTypo in error message

#     cv = CA.ComponentVector(a=(a1=100,a2=(a21=210, a22=reshape(1:4,(2,2))),a3=300), b=3, c=4)
#     # extracting some of the nested a-components and the b component
#     # note that I do not specify the structure of a2
#     cv_indextemplate = CA.ComponentVector(a=(a2=1, a3=1), b=1)
#     ax = first(CA.getaxes(cv_indextemplate)) # todo specify without template
#     cv_sub = cv[ax]

# end

# @testset "_get_index_axis AxisArray" begin
#     u1 = CA.ComponentVector(L = 10.0,)
#     p1 = CA.ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
#     popt = CA.ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
#     u1a = AxisArray(collect(u1); row = keys(u1))
#     p1a = AxisArray(collect(p1); row = keys(p1))
#     vcat(u1a, p1a) # does not preserve type - returns plain vector
#     # hence also CA.ComponentArray cannot do anything
#     u1av = CA.ComponentVector(u1a, (CA.Axis(keys(u1)),))
#     p1av = CA.ComponentVector(p1a, (CA.Axis(keys(p1)),))
#     tmp = vcat(u1av, p1av) # skips underlying type
#     typeof(tmp)
#     # try with different AbstractArray than AxisArray
# end;

# @testset "_get_index_axis NamedArray" begin
#     u1 = CA.ComponentVector(L = 10.0,)
#     p1 = CA.ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
#     popt = CA.ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
#     u1a = NamedArray(collect(u1), (collect(keys(u1)),))
#     p1a = NamedArray(collect(p1), (collect(keys(p1)),))
#     vcat(u1a, p1a) # does preserve type 
#     # hence also ComponentArray cannot do anything
#     u1av = CA.ComponentVector(u1a, (CA.Axis(keys(u1)),))
#     p1av = CA.ComponentVector(p1a, (CA.Axis(keys(p1)),))
#     tmp = vcat(u1av, p1av) # skips underlying type
#     tmp2 = p1av.k_L
#     _get_index_axis(p1av, CA.Axis(:k_L, :k_R))
#     # try with different AbstractArray than AxisArray
# end;

@testset "_labels" begin
    x = first(CA.getaxes(CA.ComponentArray(a = 1, b = [1, 2])))
    @test MTKHelpers._labels(x) == [".a", ".b[1]", ".b[2]"]
    x = first(CA.getaxes(CA.ComponentArray(c = (a = 1, b = [1, 2]))))
    @test MTKHelpers._labels(x) == [".c.a", ".c.b[1]", ".c.b[2]"]
    x = (a = 1, b = [1, 2])
    @test MTKHelpers._labels(x) == [".a", ".b[1]", ".b[2]"]
    x = first(CA.getaxes(CA.ComponentArray(c = [
        (a = [1, 2],),
        (a = [2, 3],),
        (a = [3, 4],),
    ])))
    @test MTKHelpers._labels(x) == [
        ".c[1].a[1]",
        ".c[1].a[2]",
        ".c[2].a[1]",
        ".c[2].a[2]",
        ".c[3].a[1]",
        ".c[3].a[2]",
    ]
    x = first(CA.getaxes(CA.ComponentArray(c = (b = [1 2; 5 6]))))
    @test MTKHelpers._labels(x) == [".c[1,1]", ".c[2,1]", ".c[1,2]", ".c[2,2]"]
    nt2 = (a = 5,
        b = [(a = (a = 20, b = 1), b = 0), (a = (a = 33, b = 1), b = 0)],
        c = (a = (a = 2, b = [1, 2]), b = [1.0 2.0; 5 6]))
    ca2 = CA.ComponentArray(nt2)
    x = first(CA.getaxes(ca2))
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
#     cv = CA.ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     CP.subaxis(cv, :a)    
# end;

@testset_skip "_ax_symbols_tuple" begin
    #cv = CA.ComponentVector(a = CA.ComponentVector()) # fail in CA 0.13.8
    ax = first(CA.getaxes(cv))
    @test MTKHelpers._ax_symbols_tuple(ax) == ()
end;

@testset "flatten1" begin
    cv = CA.ComponentVector(state=(x=1,y=2), par=(k=[3,4],))
    # https://discourse.julialang.org/t/how-to-execute-code-depending-on-julia-version/75029/3?u=progtw1
    @static if get_pkg_version("ComponentArrays") >= VersionNumber("0.15") 
        # test for empty subarray given sufficient version of StaticArrays 
        cv = vcat(cv, CA.ComponentVector(empty = []))
    end
    cvf = flatten1(cv)
    @test keys(cvf) == (:x, :y, :k)
    @test cvf.x == cv.state.x
    @test cvf.y == cv.state.y
    @test cvf.k == cv.par.k
    #
    cv = CA.ComponentVector(par=(k=3,))
    cvf = flatten1(cv)
    @test keys(cvf) == (:k,)
    #
    cv = CA.ComponentVector(k=3)
    cvf = flatten1(cv) 
    cvf == 3
    cvf isa Int
    #CA.getaxes(cvf)
    #@test keys(cvf) == ()
    #
    #no method @test flatten1(cvf) == cvf

end;
