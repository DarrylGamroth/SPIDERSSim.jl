using SpidersProper
using Proper
using Random
using Test
using Aqua

@testset "Package quality" begin
    Aqua.test_all(SpidersProper)
end

function prepared_allocation_bytes(prepared)
    spiders_propagate!(prepared)
    return @allocated spiders_propagate!(prepared)
end

@testset "SPIDERS masks" begin
    wf = prop_begin(7.92, 1.25e-6, 128; beam_diam_fraction=0.1)
    prop_circular_aperture(wf, 7.92 / 2)
    prop_define_entrance(wf)

    apodizer = radial_apodizer(
        SpidersProper.DEFAULT_APODIZER_POLYNOMIAL,
        128,
        prop_get_beamradius(wf),
        0.29,
        prop_get_sampling(wf),
    )
    @test size(apodizer) == (128, 128)
    @test all(isfinite, apodizer)
    @test minimum(apodizer) >= 0
    @test maximum(apodizer) ≈ 1

    prop_lens(wf, 64 * 7.92, "test focus")
    prop_propagate(wf, 64 * 7.92, "test focus")
    fpm = tilt_gaussian_fpm(
        wf;
        radius_m=459.62e-6 / 2,
        tilt_rad=1.6 * atan(1 / 64),
    )
    @test size(fpm) == (128, 128)
    @test all(isfinite, fpm)
    @test length(unique(fpm)) > 16 # centering tilt is applied after quantization
end

@testset "Lyot stop" begin
    wf = prop_begin(4.068e-3, 1.25e-6, 128; beam_diam_fraction=0.4)
    prop_circular_aperture(wf, 4.068e-3 / 2)
    prop_define_entrance(wf)
    lyot = scc_lyot_stop(
        wf,
        4.068e-3;
        pinhole_diameter_m=4.068e-3 / 20,
        obscuration_diameter_m=0.29 * 4.068e-3,
    )
    @test size(lyot.amplitude_transmission) == (128, 128)
    @test size(lyot.amplitude_reflection) == (128, 128)
    @test all(isapprox.(
        abs2.(lyot.amplitude_transmission) .+ abs2.(lyot.amplitude_reflection),
        1;
        atol=eps(Float64),
    ))
    transmitted = transmitted_lyot_field(wf, lyot)
    reflected = reflected_lyot_field(wf, lyot)
    incident_power = sum(abs2, SpidersProper.centered_field(wf))
    @test sum(abs2, transmitted) + sum(abs2, reflected) ≈ incident_power
    @test 0 < lyot.main_pupil_transmission < 1
    @test lyot.pinhole_transmission >= 0
end

@testset "Configuration validation" begin
    @test_throws ArgumentError SpidersProper.validate(SpidersConfig(pupil_mode=:bad))
    @test_throws ArgumentError SpidersProper.validate(SpidersConfig(fpm_mode=:map))
    @test SpidersProper.validate(SpidersConfig()) isa SpidersConfig
end

@testset "End-to-end prescription" begin
    config = SpidersConfig(reference_pinhole=true)
    result = spiders_proper(1.25e-6, 128; config=config, rng=MersenneTwister(42))
    intensity = spiders_intensity(result)
    @test size(intensity) == (128, 128)
    @test all(isfinite, intensity)
    @test all(>=(0), intensity)
    @test maximum(intensity) > 0
    @test prop_get_sampling(result.wavefront) > 0
    @test 0 <= result.lyot_stop.main_pupil_transmission <= 1
    @test result.lyot_stop.pinhole_transmission >= 0
    @test all(x -> 0 <= x <= 1, values(result.aperture_throughput))

    contrast = psf_contrast(result)
    @test size(contrast.intensity) == (3, 64)
    @test all(isfinite, contrast.intensity)

    detector = simulate_cred2(
        result;
        pixels=16,
        exposure_s=1e-3,
        rng=MersenneTwister(7),
    )
    @test size(detector) == (16, 16)
    @test eltype(detector) === Int

    cdi = cdi_otf(result)
    @test size(cdi.mtf) == (128, 128)
    @test all(isfinite, cdi.mtf)
    @test isfinite(cdi.visibility)
end

@testset "Alternative analytic configuration" begin
    config = SpidersConfig(
        pupil_mode=:analytic,
        apodizer_mode=:polynomial,
        fpm_type=:TGV,
        reference_pinhole=false,
    )
    result = spiders_proper(1.55e-6, 128; config=config)
    @test all(isfinite, spiders_intensity(result))
    @test iszero(result.lyot_stop.pinhole_transmission)
end

@testset "Bundled focal-plane map" begin
    config = SpidersConfig(
        fpm_mode=:map,
        fpm_map_path=joinpath(pkgdir(SpidersProper), "FPM-J-1250nm-Hlevels-500nmPixel.fits"),
        fpm_map_sampling_m=0.5e-6,
        fpm_map_type=:phase,
    )
    result = spiders_proper(1.25e-6, 128; config=config)
    @test all(isfinite, spiders_intensity(result))
end

@testset "Bundled AO residual and seeded scintillation" begin
    config = SpidersConfig(
        turbulence=true,
        scintillation=true,
        coronagraph=false,
        apodizer_mode=:none,
    )
    result = spiders_proper(1.25e-6, 128; config=config, rng=MersenneTwister(11))
    @test all(isfinite, spiders_intensity(result))
end

@testset "AO movie accepts in-memory frames" begin
    frames = zeros(Float64, 16, 16, 1)
    config = SpidersConfig(reference_pinhole=true)
    cube = spiders_movie(
        1.25e-6,
        frames;
        gridsize=128,
        output_size=128,
        final_sampling_m=nothing,
        config=config,
        rng=MersenneTwister(1),
    )
    @test size(cube) == (128, 128, 1)
    @test all(isfinite, cube)
end

@testset "Prepared repeated propagation" begin
    residual_nm = zeros(Float32, 128, 128)
    config = SpidersConfig(reference_pinhole=true)
    reference = spiders_proper(
        1.25e-6,
        128;
        config,
        ao_residual_nm=residual_nm,
    )
    prepared = prepare_spiders(
        1.25e-6,
        128;
        config,
        ao_residual_nm=residual_nm,
    )
    output = spiders_propagate!(prepared)
    @test output === spiders_intensity(prepared)
    @test output == spiders_intensity(reference)
    @test prepared_allocation_bytes(prepared) == 0

    baseline = copy(output)
    for column in axes(residual_nm, 2)
        @views residual_nm[:, column] .= column
    end
    @test spiders_propagate!(prepared) === output
    @test maximum(abs, output .- baseline) > 0
    @test prepared_allocation_bytes(prepared) == 0

    pupil = ones(Float32, 64, 64)
    opd_m = zeros(Float32, 64, 64)
    external = prepare_spiders(
        1.25e-6,
        128;
        config,
        pupil_amplitude=pupil,
        pupil_sampling_m=config.telescope_diameter_m / size(pupil, 1),
        opd_m,
        T=Float32,
    )
    external_output = spiders_propagate!(external)
    @test all(isfinite, external_output)
    @test all(>=(0), external_output)
    @test prepared_allocation_bytes(external) == 0
    external_baseline = copy(external_output)
    for column in axes(opd_m, 2)
        @views opd_m[:, column] .= column * 1e-9
    end
    spiders_propagate!(external)
    @test maximum(abs, external_output .- external_baseline) > 0
end
