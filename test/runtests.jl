using Test

@info "Running tests with $(Base.Threads.nthreads()) Julia threads active."

Test.@testset "Package MCBench" begin
    a = 42
    @test a == 42
end