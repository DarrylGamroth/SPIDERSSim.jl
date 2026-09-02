function synthetic_bax307_fixture(::Type{T}=Float32) where {T}
    actuator_map = trues(BAX307_GRID_SIZE, BAX307_GRID_SIZE)
    actuator_map[(BAX307_ACTUATOR_COUNT + 1):end] .= false
    valid_map = copy(actuator_map)
    actuator_indices = findall(actuator_map)
    for command in (38, 57, 240)
        valid_map[actuator_indices[command]] = false
    end
    support = falses(4, 4)
    support[[1, 3, 8, 12, 16]] .= true
    samples = zeros(T, count(support), BAX307_ACTUATOR_COUNT)
    samples[:, 1] .= T.(1:count(support))
    samples[:, 38] .= T(10)
    return (; actuator_map, valid_map, support, samples)
end

function warmed_bax307_graph_step_allocation_bytes(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

@testset "BAX307 calibration plan" begin
    fixture = synthetic_bax307_fixture()
    plan = prepare_bax307_calibration(
        fixture.actuator_map,
        fixture.valid_map,
        fixture.support,
        transpose(fixture.samples);
        influence_sample_unit_m=1.0f-6,
    )

    @test count(plan.actuator_map) == BAX307_ACTUATOR_COUNT
    @test count(plan.valid_actuator_map) == 465
    @test size(plan.influence_samples) == (5, BAX307_ACTUATOR_COUNT)
    @test findall(!, bax307_valid_command_mask(plan)) == [38, 57, 240]
    @test plan.actuator_map !== fixture.actuator_map
    @test plan.valid_actuator_map !== fixture.valid_map

    coordinates = bax307_actuator_coordinates(plan)
    @test size(coordinates) == (2, BAX307_ACTUATOR_COUNT)
    @test coordinates[:, 1] == Float32[-1, -1]
    @test coordinates[:, 24] == Float32[1, -1]

    modes = bax307_influence_functions_opd(plan)
    @test size(modes) == (16, BAX307_ACTUATOR_COUNT)
    @test modes[vec(plan.influence_support), 1] ==
        Float32[2, 4, 6, 8, 10] .* 1.0f-6
    @test all(iszero, modes[.!vec(plan.influence_support), :])

    command = zeros(Float32, BAX307_ACTUATOR_COUNT)
    command[1] = 1
    command[38] = 0.5
    surface = bax307_surface_opd(plan, command)
    @test surface[plan.influence_support] ==
        Float32[1, 2, 3, 4, 5] .* 1.0f-6

    fixture.samples[1, 1] = 99
    @test plan.influence_samples[1, 1] == 1

    @test_throws ArgumentError prepare_bax307_calibration(
        falses(BAX307_GRID_SIZE, BAX307_GRID_SIZE),
        fixture.valid_map,
        fixture.support,
        fixture.samples;
        influence_sample_unit_m=1.0f-6,
    )
    @test_throws DimensionMismatch prepare_bax307_calibration(
        fixture.actuator_map,
        fixture.valid_map,
        fixture.support,
        zeros(Float32, 4, BAX307_ACTUATOR_COUNT);
        influence_sample_unit_m=1.0f-6,
    )
    @test_throws ArgumentError prepare_bax307_calibration(
        fixture.actuator_map,
        fixture.valid_map,
        fixture.support,
        fixture.samples;
        influence_sample_unit_m=0,
    )
end

@testset "BAX307 measured graph node" begin
    fixture = synthetic_bax307_fixture()
    plan = prepare_bax307_calibration(
        fixture.actuator_map,
        fixture.valid_map,
        fixture.support,
        fixture.samples;
        influence_sample_unit_m=1.0f-6,
    )
    command = zeros(Float32, BAX307_ACTUATOR_COUNT)
    command[1] = 1
    command[38] = 0.5
    node = bax307_deformable_mirror_node(
        :bax307,
        plan;
        pupil_diameter_m=0.033,
    )
    definition = algorithm_graph(
        (node,);
        name=:bax307_measured_surface,
        inputs=(
            graph_input(:command, :bax307 => :pdm_command, command),
        ),
        outputs=(
            graph_output(:surface_opd, :bax307 => :surface_opd),
        ),
        parameters=bax307_graph_parameters(
            :bax307,
            plan,
            AdaptiveOpticsSim.Backends.HostComputeDevice(),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    step_graph!(graph)
    surface = graph_output(graph, Val(:surface_opd))
    expected = zeros(Float32, 4, 4)
    expected[plan.influence_support] .=
        Float32[1, 2, 3, 4, 5] .* 1.0f-6
    @test surface == expected
    @test warmed_bax307_graph_step_allocation_bytes(graph) == 0
end

@testset "BAX307 FITS calibration extension" begin
    fixture = synthetic_bax307_fixture(Float64)
    mktempdir() do root
        dm_root = joinpath(root, "dm")
        mkpath(dm_root)
        FITS(joinpath(dm_root, "BAX307-actu-map.fits"), "w") do file
            write(file, UInt8.(fixture.actuator_map))
        end
        FITS(joinpath(dm_root, "BAX307-valid-actu-map.fits"), "w") do file
            write(file, UInt8.(fixture.valid_map))
        end
        FITS(joinpath(dm_root, "AX307_Influences.fits"), "w") do file
            write(file, Matrix(transpose(fixture.samples)))
            write(file, UInt8.(fixture.support))
        end
        best_flat = collect(range(-0.25f0, 0.25f0;
            length=BAX307_ACTUATOR_COUNT))
        FITS(joinpath(dm_root, "bestflat.fits"), "w") do file
            write(file, best_flat)
        end

        plan = load_bax307_calibration(
            root;
            influence_sample_unit_m=1.0e-6,
            T=Float64,
        )
        @test plan.influence_samples == fixture.samples
        @test load_bax307_best_flat(root) == best_flat
    end
end
