using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using DataFrames, CSV
using ComponentArrays: ComponentArrays as CA

df0 = DataFrame(
    indiv_id = ["A", "B", "C"],
    u0 = [
        CA.ComponentVector(sv₊x = [2.0383292042153554, 2.0383292042153554]),
        CA.ComponentVector(sv₊x = [2.106525817516089, 2.038672471649886]),
        CA.ComponentVector(sv₊x = [2.010654503237803, 2.0510192980037196])
    ],
    p = [
        CA.ComponentVector(
            sv₊τ = 1.4009482259635606, sv₊i = 1.518711604434893, sv₊i2 = 0.1,
            sv₊p = [2.400789101642099, 2.400789101642099, 2.400789101642099]),
        CA.ComponentVector(
            sv₊τ = 1.4752043120005407, sv₊i = 1.518711604434893, sv₊i2 = 0.1,
            sv₊p = [2.400789101642099, 2.400789101642099, 2.400789101642099]),
        CA.ComponentVector(
            sv₊τ = 1.4034321912259409, sv₊i = 1.518711604434893, sv₊i2 = 0.1,
            sv₊p = [2.400789101642099, 2.400789101642099, 2.400789101642099])]
)

@testset "write_csv_cv" begin
    con = IOBuffer()
    write_csv_cv(con, df0)
    s = String(take!(con)); 
    # information on axes in comments and cv spread across several columns
    # print(s)
    # can read the the csv-data by usual csv-tools albeit in scalar columns
    df_plain = CSV.read(IOBuffer(s), DataFrame; comment="#", ignoreemptyrows=true)
    # will reconstruct ComponentVectors
    df2 = read_csv_cv(IOBuffer(s))
    @test df2 == df0
end;