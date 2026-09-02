function _scc_propagation_allocation_bytes(prepared)
    propagate_scc!(prepared)
    return @allocated propagate_scc!(prepared)
end

function _llowfs_propagation_allocation_bytes(prepared)
    propagate_llowfs!(prepared)
    return @allocated propagate_llowfs!(prepared)
end

@testset "SPIDERS Proper optical configuration" begin
    configuration = SCCPropagationConfiguration()
    @test configuration.pupil_mode === :analytic
    @test configuration.apodizer_mode === :polynomial
    @test SPIDERSSim.ProperOptics.validate(configuration) === configuration

    @test_throws ArgumentError SPIDERSSim.ProperOptics.validate(
        SCCPropagationConfiguration(pupil_mode=:fits),
    )
    @test_throws ArgumentError SPIDERSSim.ProperOptics.validate(
        SCCPropagationConfiguration(fpm_mode=:map),
    )

    relay = LLOWFSRelayConfiguration()
    @test SPIDERSSim.ProperOptics.validate(relay) === relay
    @test_throws ArgumentError SPIDERSSim.ProperOptics.validate(
        LLOWFSRelayConfiguration(detector_pixel_pitch_m=14e-6),
    )
end

@testset "SPIDERS reflective Lyot stop" begin
    wavefront = prop_begin(4.068e-3, 1.25e-6, 128;
        beam_diam_fraction=0.4)
    prop_circular_aperture(wavefront, 4.068e-3 / 2)
    prop_define_entrance(wavefront)
    stop = scc_lyot_stop(
        wavefront,
        4.068e-3;
        pinhole_diameter_m=4.068e-3 / 20,
        obscuration_diameter_m=0.29 * 4.068e-3,
    )

    @test size(stop.amplitude_transmission) == (128, 128)
    @test size(stop.amplitude_reflection) == (128, 128)
    @test all(isapprox.(
        abs2.(stop.amplitude_transmission) .+
        abs2.(stop.amplitude_reflection),
        1;
        atol=eps(Float64),
    ))
    transmitted = transmitted_lyot_field(wavefront, stop)
    reflected = reflected_lyot_field(wavefront, stop)
    incident = sum(abs2, SPIDERSSim.ProperOptics.centered_field(wavefront))
    @test sum(abs2, transmitted) + sum(abs2, reflected) ≈ incident
end

@testset "SPIDERS Zemax LLOWFS relay" begin
    indices = llowfs_refractive_indices(1.55e-6)
    @test indices.s_lah71 ≈ 1.8134701603020393 atol=1e-14
    @test indices.s_nph5 ≈ 1.8090114983950631 atol=1e-14
    @test indices.fused_silica ≈ 1.4440236216697244 atol=1e-14

    relay = llowfs_relay_summary(1.55e-6)
    @test relay.effective_focal_length_m ≈ 157.691076323e-3 atol=1e-12
    @test relay.back_focal_length_m ≈ 153.859769413e-3 atol=1e-12
    @test relay.conjugate_distance_m ≈ 223.444722793e-3 atol=1e-12
    @test relay.detector_distance_m ≈ 223.260958398e-3 atol=1e-12
    @test relay.detector_defocus_m ≈ -183.764395e-6 atol=1e-12
    @test relay.lateral_magnification ≈ -0.441273881834 atol=1e-12
end

@testset "Prepared SCC and LLOWFS Proper propagation" begin
    pupil_amplitude = ones(Float32, 64, 64)
    pupil_opd = zeros(Float32, 64, 64)
    configuration = SCCPropagationConfiguration(reference_pinhole=true)
    pupil_sampling_m = configuration.telescope_diameter_m /
        size(pupil_amplitude, 1)

    scc = prepare_scc_propagation(
        1.55e-6,
        128;
        configuration,
        pupil_amplitude,
        pupil_sampling_m,
        opd_m=pupil_opd,
        T=Float32,
    )
    llowfs = prepare_llowfs_propagation(
        1.55e-6,
        128;
        configuration,
        pupil_amplitude,
        pupil_sampling_m,
        opd_m=pupil_opd,
        T=Float32,
    )

    scc_output = propagate_scc!(scc)
    llowfs_output = propagate_llowfs!(llowfs)
    @test scc_output === scc_intensity(scc)
    @test llowfs_output === llowfs_intensity(llowfs)
    @test size(scc_output) == (128, 128)
    @test size(llowfs_output) == (128, 128)
    @test all(isfinite, scc_output)
    @test all(isfinite, llowfs_output)
    @test minimum(scc_output) >= 0
    @test minimum(llowfs_output) >= 0
    @test maximum(scc_output) > 0
    @test maximum(llowfs_output) > 0
    @test _scc_propagation_allocation_bytes(scc) == 0
    @test _llowfs_propagation_allocation_bytes(llowfs) == 0

    scc_baseline = copy(scc_output)
    llowfs_baseline = copy(llowfs_output)
    for column in axes(pupil_opd, 2)
        @views pupil_opd[:, column] .= column * 1f-9
    end
    @test propagate_scc!(scc) === scc_output
    @test propagate_llowfs!(llowfs) === llowfs_output
    @test maximum(abs, scc_output .- scc_baseline) > 0
    @test maximum(abs, llowfs_output .- llowfs_baseline) > 0
    @test _scc_propagation_allocation_bytes(scc) == 0
    @test _llowfs_propagation_allocation_bytes(llowfs) == 0
end
