function _spiders_tilt!(opd::AbstractMatrix{T}, amplitude::T) where {T}
    n = size(opd, 2)
    center = T(n + 1) / T(2)
    @inbounds for column in axes(opd, 2)
        value = amplitude * (T(column) - center) / T(n)
        for row in axes(opd, 1)
            opd[row, column] = value
        end
    end
    return opd
end

function _spiders_graph(configuration, pupil_opd, pupil_amplitude)
    return prepare_algorithm_graph(
        _spiders_test_graph(
            configuration,
            pupil_opd,
            pupil_amplitude,
        );
        target=AdaptiveOpticsSim.Backends.HostComputeDevice(),
    )
end

function _spiders_test_graph(configuration, pupil_opd, pupil_amplitude)
    return algorithm_graph(
        (proper_propagation_node(:proper, configuration),);
        name=:spiders_test,
        inputs=(
            graph_input(:pupil_opd, :proper => :pupil_opd, pupil_opd),
            graph_input(
                :pupil_amplitude,
                :proper => :pupil_amplitude,
                pupil_amplitude,
            ),
        ),
        outputs=(graph_output(:intensity, :proper => :output),),
    )
end

@testset "provisional SPIDERS Proper prescriptions" begin
    profile = provisional_spiders_h_regular_profile(Float32)
    llowfs_configuration = provisional_spiders_llowfs_configuration(profile)
    scc_configuration = provisional_spiders_scc_configuration(profile)

    @test llowfs_configuration.resolution == 512
    @test llowfs_configuration.pupil_resolution == 128
    @test llowfs_configuration.output_rows == 34
    @test llowfs_configuration.output_columns == 47
    @test llowfs_configuration.output_schema ==
        PROVISIONAL_SPIDERS_LLOWFS_SCHEMA
    @test scc_configuration.resolution == 512
    @test scc_configuration.pupil_resolution == 128
    @test scc_configuration.output_rows == 416
    @test scc_configuration.output_columns == 380
    @test scc_configuration.output_schema == PROVISIONAL_SPIDERS_SCC_SCHEMA
    @test length(spiders_prescription_claims(llowfs_configuration)) == 4
    @test length(spiders_prescription_claims(scc_configuration)) == 4
    @test length(spiders_configuration_claims(llowfs_configuration)) == 12
    @test all(
        claim -> claim.evidence === SpidersPlaceholder,
        spiders_prescription_claims(llowfs_configuration),
    )
    @test Set(
        claim.name for claim in spiders_prescription_claims(scc_configuration)
    ) == Set((
        :pupil_embedding,
        :tilt_gaussian_vortex_surrogate,
        :proper_relay,
        :provisional_detector_model,
    ))

    pupil_opd = zeros(Float32, 128, 128)
    pupil_amplitude = ones(Float32, 128, 128)
    llowfs_graph = _spiders_graph(
        llowfs_configuration,
        pupil_opd,
        pupil_amplitude,
    )
    scc_graph = _spiders_graph(
        scc_configuration,
        pupil_opd,
        pupil_amplitude,
    )

    @test step_graph!(llowfs_graph) === llowfs_graph
    @test step_graph!(scc_graph) === scc_graph
    llowfs_reference = copy(graph_output(llowfs_graph, :intensity))
    scc_reference = copy(graph_output(scc_graph, :intensity))
    @test size(llowfs_reference) == profile.llowfs_output_shape
    @test size(scc_reference) == profile.scc_output_shape
    @test all(isfinite, llowfs_reference)
    @test all(isfinite, scc_reference)
    @test minimum(llowfs_reference) >= -10eps(Float32)
    @test minimum(scc_reference) >= -10eps(Float32)
    @test sum(llowfs_reference) > 0
    @test sum(scc_reference) > 0

    _spiders_tilt!(pupil_opd, 2f-6)
    @test step_graph!(llowfs_graph) === llowfs_graph
    @test step_graph!(scc_graph) === scc_graph
    @test maximum(abs, graph_output(llowfs_graph, :intensity) .-
        llowfs_reference) > 1f-4
    @test maximum(abs, graph_output(scc_graph, :intensity) .-
        scc_reference) > 1f-4
    @test @allocated(step_graph!(llowfs_graph)) == 0
    @test @allocated(step_graph!(scc_graph)) == 0
    @test @inferred(step_graph!(llowfs_graph)) === llowfs_graph
    @test @inferred(step_graph!(scc_graph)) === scc_graph

    llowfs_owner = AlgorithmGraphs.prepared_graph_node(
        llowfs_graph,
        :proper,
    )
    scc_owner = AlgorithmGraphs.prepared_graph_node(scc_graph, :proper)
    @test AlgorithmGraphs.graph_node_capture_capability(llowfs_owner) isa
        AlgorithmGraphs.GraphNodeCaptureSafe
    @test AlgorithmGraphs.graph_node_capture_capability(scc_owner) isa
        AlgorithmGraphs.GraphNodeCaptureSafe
    @test llowfs_owner.pupil_opd === pupil_opd
    @test scc_owner.pupil_opd === pupil_opd
    @test size(llowfs_owner.workspace.assets.padded_pupil_opd) == (512, 512)
    @test size(scc_owner.workspace.assets.focal_plane_mask_internal) ==
        (512, 512)
    @test all(isfinite, scc_owner.workspace.assets.focal_plane_mask_internal)
    @test sum(scc_owner.workspace.assets.lyot_mask) > 0

    closed_profile = provisional_spiders_h_regular_profile(
        Float32;
        scc_reference_beam=:closed,
    )
    closed_configuration = provisional_spiders_scc_configuration(
        closed_profile,
    )
    closed_graph = _spiders_graph(
        closed_configuration,
        pupil_opd,
        pupil_amplitude,
    )
    @test step_graph!(closed_graph) === closed_graph
    @test maximum(abs, graph_output(closed_graph, :intensity) .-
        graph_output(scc_graph, :intensity)) > 1f-6

    @test_throws ArgumentError provisional_spiders_scc_configuration(
        profile;
        resolution=256,
    )
    @test_throws ArgumentError provisional_spiders_llowfs_configuration(
        profile;
        resolution=512,
        pupil_resolution=0,
    )
    @test_throws ArgumentError provisional_spiders_scc_configuration(
        profile;
        beam_diameter_fraction=0.5,
    )
end
