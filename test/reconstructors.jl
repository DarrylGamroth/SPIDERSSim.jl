@testset "GoldEye LLOWFS processor" begin
    reference = fill(0.25f0, 2, 2)
    image_mask = Bool[1 0; 1 1]
    centroid_mask = trues(2, 2)
    image_to_modes = Float32[1 0; 0 1; 1 1]
    modes_to_actuators = Float32[1 2; 3 4]
    slopes_to_tiptilt = Matrix{Float32}(I, 2, 2)
    plan, workspace, products = prepare_llowfs(
        reference,
        image_mask,
        centroid_mask,
        image_to_modes,
        modes_to_actuators,
        slopes_to_tiptilt,
        (1.5f0, 1.5f0),
    )

    image = Float32[2 2; 2 2]
    @test process_llowfs!(products, workspace, plan, image) === products
    @test products.valid
    @test products.total_flux == 8f0
    @test products.centroid_flux == 8f0
    @test products.modes == Float32[0, 0]
    @test products.tiptilt == Float32[0, 0]
    @test products.actuator_coefficients == Float32[0, 0]

    image[1, 1] = 4f0
    process_llowfs!(products, workspace, plan, image)
    @test products.valid
    @test products.modes ≈ Float32[0.1, -0.1]
    @test products.actuator_coefficients ≈ Float32[-0.2, -0.2]
    @test products.tiptilt ≈ Float32[-0.1, -0.1]
    @test @allocated(process_llowfs!(products, workspace, plan, image)) == 0
    @test @inferred(process_llowfs!(products, workspace, plan, image)) === products

    strided_parent = repeat(image, inner=(2, 2))
    strided_image = @view strided_parent[1:2:end, 1:2:end]
    process_llowfs!(products, workspace, plan, strided_image)
    @test products.valid
    @test products.actuator_coefficients ≈ Float32[-0.2, -0.2]

    fill!(image, 0f0)
    process_llowfs!(products, workspace, plan, image)
    @test !products.valid
    @test all(iszero, products.actuator_coefficients)

    @test_throws DimensionMismatch prepare_llowfs(
        reference,
        trues(3, 3),
        centroid_mask,
        image_to_modes,
        modes_to_actuators,
        slopes_to_tiptilt,
        (1.5f0, 1.5f0),
    )
end

@testset "C-RED 2 SCC processor" begin
    reference = zeros(Float32, 2, 3)
    image_mask = Bool[1 1 0; 1 0 1]
    correction_mask = Bool[1 0 0; 1 0 1]
    image_to_actuators = Float32[1 2; 3 4; 5 6; 7 8]
    plan, workspace, products = prepare_scc(
        reference,
        image_mask,
        correction_mask,
        image_to_actuators,
    )
    difference = Float32[1 10 0; 2 0 3]
    @test process_scc!(products, workspace, plan, difference) === products
    @test products.valid
    @test workspace.selected_residual == Float32[1, 2, 0, 3]
    @test products.actuator_coefficients == Float32[28, 34]
    @test @allocated(process_scc!(products, workspace, plan, difference)) == 0
    @test @inferred(process_scc!(products, workspace, plan, difference)) === products

    strided_parent = repeat(difference, inner=(2, 2))
    strided_difference = @view strided_parent[1:2:end, 1:2:end]
    process_scc!(products, workspace, plan, strided_difference)
    @test products.valid
    @test products.actuator_coefficients == Float32[28, 34]

    difference[1, 1] = Float32(Inf)
    process_scc!(products, workspace, plan, difference)
    @test !products.valid
    @test all(iszero, products.actuator_coefficients)

    @test_throws DimensionMismatch prepare_scc(
        reference,
        trues(2, 2),
        correction_mask,
        image_to_actuators,
    )
end
