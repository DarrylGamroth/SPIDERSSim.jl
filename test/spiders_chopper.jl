function _step_chopped_scc!(graph)
    step_graph!(graph)
    return spiders_chopper_frame_status(graph, Val(:proper))
end

@testset "frame-synchronous SPIDERS chopper" begin
    profile = provisional_spiders_h_regular_profile(Float32)
    configuration = provisional_spiders_chopped_scc_configuration(profile)

    @test configuration.output_schema ==
        PROVISIONAL_SPIDERS_CHOPPED_SCC_SCHEMA
    @test configuration.prescription.chopper.initial_phase === SpidersFringed
    @test length(spiders_prescription_claims(configuration)) == 5
    @test length(spiders_configuration_claims(configuration)) == 14
    @test last(spiders_prescription_claims(configuration)).name ===
        :frame_synchronous_chopper

    pupil_opd = zeros(Float32, 128, 128)
    _spiders_tilt!(pupil_opd, 2f-7)
    pupil_amplitude = ones(Float32, 128, 128)
    chopped_graph = _spiders_graph(
        configuration,
        pupil_opd,
        pupil_amplitude,
    )
    @test_throws ArgumentError spiders_chopper_frame_status(chopped_graph)

    first_status = _step_chopped_scc!(chopped_graph)
    fringed = copy(graph_output(chopped_graph, :intensity))
    @test first_status == SpidersChopperFrameStatus(
        UInt64(1),
        SpidersFringed,
    )
    @test spiders_chopper_sequence(first_status) == UInt64(1)
    @test spiders_chopper_phase(first_status) === SpidersFringed
    @test spiders_chopper_frame_status(chopped_graph, :proper) == first_status

    second_status = _step_chopped_scc!(chopped_graph)
    unfringed = copy(graph_output(chopped_graph, :intensity))
    @test second_status == SpidersChopperFrameStatus(
        UInt64(2),
        SpidersUnfringed,
    )
    @test maximum(abs, fringed .- unfringed) > 1f-6

    open_graph = _spiders_graph(
        provisional_spiders_scc_configuration(profile),
        pupil_opd,
        pupil_amplitude,
    )
    closed_profile = provisional_spiders_h_regular_profile(
        Float32;
        scc_reference_beam=:closed,
    )
    closed_graph = _spiders_graph(
        provisional_spiders_scc_configuration(closed_profile),
        pupil_opd,
        pupil_amplitude,
    )
    step_graph!(open_graph)
    step_graph!(closed_graph)
    @test fringed == graph_output(open_graph, :intensity)
    @test unfringed == graph_output(closed_graph, :intensity)

    third_status = _step_chopped_scc!(chopped_graph)
    @test third_status == SpidersChopperFrameStatus(
        UInt64(3),
        SpidersFringed,
    )
    @test graph_output(chopped_graph, :intensity) == fringed
    @test @allocated(_step_chopped_scc!(chopped_graph)) == 0
    @test @inferred(_step_chopped_scc!(chopped_graph)) ==
        SpidersChopperFrameStatus(UInt64(5), SpidersFringed)

    @test reset_graph!(chopped_graph) === chopped_graph
    @test_throws ArgumentError spiders_chopper_frame_status(chopped_graph)
    reset_status = _step_chopped_scc!(chopped_graph)
    @test reset_status == SpidersChopperFrameStatus(
        UInt64(1),
        SpidersFringed,
    )
    @test graph_output(chopped_graph, :intensity) == fringed

    shared_chopped_graph = _shared_spiders_graph(
        provisional_spiders_llowfs_configuration(profile),
        configuration,
        pupil_opd,
        pupil_amplitude,
    )
    @test step_graph!(shared_chopped_graph) === shared_chopped_graph
    @test spiders_chopper_frame_status(
        shared_chopped_graph,
        Val(:spiders_optics),
    ) == SpidersChopperFrameStatus(UInt64(1), SpidersFringed)
    @test graph_output(shared_chopped_graph, :scc_intensity) == fringed
    @test step_graph!(shared_chopped_graph) === shared_chopped_graph
    @test spiders_chopper_frame_status(
        shared_chopped_graph,
        Val(:spiders_optics),
    ) == SpidersChopperFrameStatus(UInt64(2), SpidersUnfringed)
    @test graph_output(shared_chopped_graph, :scc_intensity) == unfringed
    shared_chopped_owner = AlgorithmGraphs.prepared_graph_node(
        shared_chopped_graph,
        :spiders_optics,
    )
    @test AlgorithmGraphs.graph_node_capture_capability(
        shared_chopped_owner,
    ) isa AlgorithmGraphs.GraphNodeCaptureUnsupported
    @test @allocated(step_graph!(shared_chopped_graph)) == 0
    @test reset_graph!(shared_chopped_graph) === shared_chopped_graph
    @test_throws ArgumentError spiders_chopper_frame_status(
        shared_chopped_graph,
        Val(:spiders_optics),
    )
    @test step_graph!(shared_chopped_graph) === shared_chopped_graph
    @test spiders_chopper_frame_status(
        shared_chopped_graph,
        Val(:spiders_optics),
    ) == SpidersChopperFrameStatus(UInt64(1), SpidersFringed)

    reset_graph!(chopped_graph)
    driver = FixedStepModelTimeDriver(
        AlgorithmGraphs.PeriodicSchedule(period_ns=1_000_000),
    )
    @test step_graph_at!(chopped_graph, driver) ==
        ModelTimestamp(0)
    @test spiders_chopper_frame_status(chopped_graph) ==
        SpidersChopperFrameStatus(UInt64(1), SpidersFringed)
    @test step_graph_at!(chopped_graph, driver) ==
        ModelTimestamp(1_000_000)
    @test spiders_chopper_frame_status(chopped_graph) ==
        SpidersChopperFrameStatus(UInt64(2), SpidersUnfringed)

    unfringed_first = provisional_spiders_chopped_scc_configuration(
        profile;
        initial_phase=SpidersUnfringed,
    )
    unfringed_first_graph = _spiders_graph(
        unfringed_first,
        pupil_opd,
        pupil_amplitude,
    )
    @test _step_chopped_scc!(unfringed_first_graph) ==
        SpidersChopperFrameStatus(UInt64(1), SpidersUnfringed)
    @test graph_output(unfringed_first_graph, :intensity) == unfringed
    @test _step_chopped_scc!(unfringed_first_graph) ==
        SpidersChopperFrameStatus(UInt64(2), SpidersFringed)
    @test graph_output(unfringed_first_graph, :intensity) == fringed

    reset_graph!(chopped_graph)
    owner = AlgorithmGraphs.prepared_graph_node(chopped_graph, :proper)
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa
        AlgorithmGraphs.GraphNodeCaptureUnsupported
    owner.state.prescription.completed_sequence = typemax(UInt64)
    prior_output = copy(graph_output(chopped_graph, :intensity))
    @test_throws ArgumentError step_graph!(chopped_graph)
    @test graph_failed(chopped_graph)
    @test graph_output(chopped_graph, :intensity) == prior_output
    @test_throws ArgumentError spiders_chopper_frame_status(chopped_graph)
    @test reset_graph!(chopped_graph) === chopped_graph
    @test _step_chopped_scc!(chopped_graph) ==
        SpidersChopperFrameStatus(UInt64(1), SpidersFringed)

    node_types = merge(
        builtin_graph_node_types(),
        (
            spiders_scc_proper=
                proper_propagation_node_factory(configuration),
        ),
    )
    file_graph = prepare_algorithm_graph(
        load_algorithm_graph(
            joinpath(
                dirname(@__DIR__),
                "graphs",
                "spiders_scc_hil.toml",
            );
            node_types,
            bindings=(; pupil_opd, pupil_amplitude),
        );
        target=AdaptiveOpticsSim.Backends.HostComputeDevice(),
    )
    @test graph_name(file_graph) === :spiders_scc_hil
    @test_throws ArgumentError spiders_chopper_frame_status(
        file_graph,
        Val(:scc_optics),
    )
    file_driver = FixedStepModelTimeDriver(
        AlgorithmGraphs.PeriodicSchedule(period_ns=1_000_000),
    )
    @test step_graph_at!(file_graph, file_driver) ==
        ModelTimestamp(0)
    @test spiders_chopper_frame_status(file_graph, Val(:scc_optics)) ==
        SpidersChopperFrameStatus(UInt64(1), SpidersFringed)
    @test graph_output(file_graph, :scc_observation) == fringed

    @test_throws ArgumentError provisional_spiders_chopped_scc_configuration(
        closed_profile,
    )
    @test_throws ArgumentError spiders_chopper_frame_status(open_graph)
end
