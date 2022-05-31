@testset "_update_cv" begin
    cv = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
    s = ComponentVector(a=(a1=11,a2=2,a3=3),b=20)
    cv1 = _update_cv(cv, s)
    @test cv1 == s
    s = ComponentVector(b=21)
    cv1 = _update_cv(cv, s)
    @test (cv1.a == cv.a) & (cv1.b == s.b)
    s = ComponentVector(a=(a1=11,))
    cv1 = _update_cv(cv, s)
    @test (cv1.a.a1 == 11) & (cv1.a.a2 == cv.a.a2) & (cv1.b == cv.b)
    s = ComponentVector(a=(aTypo=11,))
    @test_throws ErrorException _update_cv(cv, s)
    s = ComponentVector(a=11:13)
    cv1 = _update_cv(cv, s)
    @test cv1.a.a3 == 13
    s = ComponentVector(a=1)
    @test_throws ErrorException _update_cv(cv, s) # wrong length ComponentVector
    s = ComponentVector(b=20:21) # wrong length flat axis
    @test_throws ErrorException _update_cv(cv, s) 
    ax = first(getaxes(cv))
    CA.last_index(ax)
    lastindex(ax)
end;

@testset "_get_index_axis" begin
    # TODO update src and test in ComponentArrays.jl bgctw branch
    cv = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
    #s = ComponentVector(b=1)
    s = ComponentVector(a=1)
    ax = first(getaxes(s))
    @test_throws ErrorException _get_index_axis(cv, ax)
end;




