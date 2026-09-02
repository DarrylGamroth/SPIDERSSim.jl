@testset "provisional SPIDERS profile provenance" begin
    profile = provisional_spiders_h_regular_profile(Float32)

    @test spiders_profile_is_provisional(profile)
    @test profile.mode === :h_regular
    @test profile.wavelength_um === 1.55f0
    @test profile.focal_plane_mask.stage_serial == "27260961"
    @test profile.focal_plane_mask.stage_counts == 424_960
    @test profile.focal_plane_mask.diameter_lambda_over_d === 5.7f0
    @test profile.focal_plane_mask.vortex_charge == 4
    @test profile.galil_alignment.fpm_x_mm === 6.91f0
    @test profile.galil_alignment.dm_y_mm === 6.827f0

    @test profile.llowfs_filter.wheel == 2
    @test profile.llowfs_filter.serial == "TP02692146-23357"
    @test profile.scc_filter.wheel == 3
    @test profile.scc_filter.serial == "TP02659529-22864"
    @test profile.llowfs_filter.slot == profile.scc_filter.slot == 5
    @test profile.llowfs_filter.center_wavelength_um === 1.55f0
    @test profile.llowfs_filter.full_width_um === 0.025f0
    @test profile.llowfs_output_shape == (34, 47)
    @test profile.scc_output_shape == (416, 380)

    optics = profile.optics
    @test optics.apodizer_outer_diameter_ratio ≈ 7840.8f0 / 7920f0
    @test optics.lyot_outer_diameter_ratio ≈ 7527.2f0 / 7920f0
    @test optics.reference_pinhole_diameter_m === 203f-6
    @test optics.reference_pinhole_separation_pupils === 1.6f0
    @test optics.camera_rotation_deg === 31f0

    claims = spiders_profile_claims(profile)
    @test length(claims) == 9
    @test all(claim -> !isempty(claim.current_choice), claims)
    @test all(claim -> !isempty(claim.source), claims)
    @test all(claim -> !isempty(claim.qualification_question), claims)
    @test Set(claim.name for claim in claims) == Set((
        :science_wheel_ownership,
        :science_filter_slot,
        :llowfs_selector_state,
        :scc_selector_state,
        :reference_pinhole,
        :focal_plane_mask_map,
        :entrance_pupil_map,
        :apodizer_map,
        :detector_mapping,
    ))
    @test only(filter(
        claim -> claim.name === :science_wheel_ownership,
        claims,
    )).evidence === SpidersPlaceholder
    @test only(filter(
        claim -> claim.name === :science_filter_slot,
        claims,
    )).evidence === SpidersInferred
    detector_claim = only(filter(
        claim -> claim.name === :detector_mapping,
        claims,
    ))
    @test detector_claim.evidence === SpidersDeploymentConfigured
    @test occursin("GoldEye", detector_claim.current_choice)
    @test occursin("C-RED 2", detector_claim.current_choice)

    entrance_pupil = provisional_spiders_entrance_pupil_amplitude(
        profile,
        128,
    )
    @test size(entrance_pupil) == (128, 128)
    @test eltype(entrance_pupil) === Float32
    @test Set(entrance_pupil) == Set((0.0f0, 1.0f0))
    @test entrance_pupil == reverse(entrance_pupil; dims=(1, 2))
    @test iszero(entrance_pupil[64, 64])
    @test iszero(entrance_pupil[1, 1])
    @test sum(entrance_pupil) > 0
    annular_count = count(
        CartesianIndices(entrance_pupil),
    ) do index
        row, column = Tuple(index)
        x = (Float32(column) - 64.5f0) / 64.0f0
        y = (Float32(row) - 64.5f0) / 64.0f0
        radius = hypot(x, y)
        0.29f0 <= radius <= 1.0f0
    end
    @test count(!iszero, entrance_pupil) < annular_count
    @test_throws ArgumentError provisional_spiders_entrance_pupil_amplitude(
        profile,
        0,
    )

    reversed = provisional_spiders_h_regular_profile(
        Float64;
        llowfs_wheel=3,
        scc_wheel=2,
        reference_pinhole_diameter_m=291e-6,
        reference_pinhole_position_angle_deg=45,
    )
    @test reversed.llowfs_filter.wheel == 3
    @test reversed.scc_filter.wheel == 2
    @test reversed.optics.reference_pinhole_diameter_m === 291e-6
    @test reversed.optics.reference_pinhole_position_angle_deg === 45.0

    @test_throws ArgumentError provisional_spiders_h_regular_profile(
        Float32;
        llowfs_wheel=2,
        scc_wheel=2,
    )
    @test_throws ArgumentError provisional_spiders_h_regular_profile(
        Float32;
        llowfs_plane=:unknown,
    )
    @test_throws ArgumentError provisional_spiders_h_regular_profile(
        Float32;
        reference_pinhole_diameter_m=0,
    )
end
