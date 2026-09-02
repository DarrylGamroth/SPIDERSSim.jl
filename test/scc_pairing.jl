@testset "C-RED 2 SCC phase-tagged pairing" begin
    plan, state, products = prepare_scc_pairing(Float32; frame_shape=(2, 3))
    fringed = fill(4f0, 2, 3)
    unfringed = fill(1f0, 2, 3)

    @test accept_scc_frame!(
        products,
        state,
        plan,
        fringed,
        UInt64(1),
        SpidersFringed,
    ) === products
    @test !products.valid
    @test all(iszero, products.difference)

    accept_scc_frame!(
        products,
        state,
        plan,
        unfringed,
        UInt64(2),
        SpidersUnfringed,
    )
    @test products.valid
    @test products.difference == fill(3f0, 2, 3)
    @test products.fringed_sequence == UInt64(1)
    @test products.unfringed_sequence == UInt64(2)

    fill!(fringed, 5f0)
    accept_scc_frame!(
        products,
        state,
        plan,
        fringed,
        UInt64(3),
        SpidersFringed,
    )
    @test products.valid
    @test products.difference == fill(4f0, 2, 3)

    accept_scc_frame!(
        products,
        state,
        plan,
        unfringed,
        UInt64(5),
        SpidersUnfringed,
    )
    @test !products.valid
    @test all(iszero, products.difference)

    accept_scc_frame!(
        products,
        state,
        plan,
        fringed,
        UInt64(6),
        SpidersFringed,
    )
    @test products.valid
    @test products.difference == fill(4f0, 2, 3)
    @test @allocated(accept_scc_frame!(
        products,
        state,
        plan,
        unfringed,
        UInt64(7),
        SpidersUnfringed,
    )) == 0
    @test_throws ArgumentError accept_scc_frame!(
        products,
        state,
        plan,
        fringed,
        UInt64(7),
        SpidersFringed,
    )

    @test reset_scc_pairing!(state, products) === products
    @test !products.valid
    @test state.last_sequence == UInt64(0)
    @test all(iszero, state.fringed)
    @test_throws ArgumentError prepare_scc_pairing(
        Float32;
        frame_shape=(0, 3),
    )
end
