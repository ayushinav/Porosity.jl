@testitem "types" tags=[:types] begin
    m = MTModel([100.0, 100.0], [100.0])
    @test typeof(m) <: AbstractGeophyModel
    @test typeof(zero(m)) <: AbstractGeophyModel
    @test typeof(copy(m)) <: AbstractGeophyModel

    ω = [1.0, 10.0]
    resp = forward(m, ω)

    @test typeof(resp) <: AbstractGeophyResponse
    @test typeof(zero(resp)) <: AbstractGeophyResponse
    @test typeof(copy(resp)) <: AbstractGeophyResponse
end
