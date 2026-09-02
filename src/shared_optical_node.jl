"""
    provisional_spiders_optical_node(name, llowfs_configuration, scc_configuration)

Declare one provisional SPIDERS graph node that shares coherent propagation
through the Lyot plane before producing the LLOWFS and SCC detector products.

Both configurations must use the same propagation and external-pupil grids,
entrance schemas, wavelength, optical geometry, focal-plane mask, and common
surrogate-mask parameters. Detector defocus, output shape, output schema, and
camera magnification remain branch-specific. The SCC configuration may use the
fixed or camera-synchronous chopped prescription.
"""
function provisional_spiders_optical_node(
    name::Symbol,
    llowfs_configuration::ProperPropagationConfiguration{
        T,
        LLOWFS,
    },
    scc_configuration::ProperPropagationConfiguration{T,SCC},
) where {
    T,
    LLOWFS<:ProvisionalSpidersLLOWFSPrescription{T},
    SCC<:Union{
        ProvisionalSpidersSCCPrescription{T},
        ProvisionalSpidersChoppedSCCPrescription{T},
    },
}
    configuration = _validate_shared_optical_configuration(
        llowfs_configuration,
        scc_configuration,
    )
    return _AOG.algorithm_node(
        name,
        _ProvisionalSpidersOpticalNode{T,LLOWFS,SCC},
        configuration;
        props=NamedTuple(),
    )
end

struct _ProvisionalSpidersOpticalConfiguration{LLOWFS,SCC}
    llowfs::LLOWFS
    scc::SCC
end

struct _ProvisionalSpidersOpticalNode{T,LLOWFS,SCC} end

struct _ProvisionalSpidersOpticalPlan{Configuration,Target,SCC}
    configuration::Configuration
    target::Target
    prescription::SCC
end

mutable struct _ProvisionalSpidersOpticalState{RNG,SCCState}
    rng::RNG
    initial_rng::RNG
    prescription::SCCState
end

struct _ProvisionalSpidersOpticalWorkspace{
    CommonContext,
    SCCContext,
    CommonField,
    SCCField,
    CommonWavefront,
    SCCWavefront,
    LLOWFSAssets,
    SCCAssets,
}
    common_context::CommonContext
    scc_context::SCCContext
    common_field::CommonField
    scc_field::SCCField
    common_wavefront::CommonWavefront
    scc_wavefront::SCCWavefront
    llowfs_assets::LLOWFSAssets
    scc_assets::SCCAssets
end

struct _ProvisionalSpidersOpticalOwner{
    Plan,
    State,
    Workspace,
    OPD,
    Amplitude,
    LLOWFSOutput,
    SCCOutput,
}
    plan::Plan
    state::State
    workspace::Workspace
    pupil_opd::OPD
    pupil_amplitude::Amplitude
    llowfs_output::LLOWFSOutput
    scc_output::SCCOutput
end

function _require_shared_optical_value(
    name::AbstractString,
    llowfs_value,
    scc_value,
)
    isequal(llowfs_value, scc_value) || throw(ArgumentError(
        "shared SPIDERS $name differs between the LLOWFS and SCC configurations",
    ))
    return llowfs_value
end

function _validate_shared_optical_configuration(llowfs, scc)
    _require_shared_optical_value(
        "propagation resolution",
        llowfs.resolution,
        scc.resolution,
    )
    _require_shared_optical_value(
        "pupil resolution",
        llowfs.pupil_resolution,
        scc.pupil_resolution,
    )
    _require_shared_optical_value(
        "pupil diameter",
        llowfs.diameter_m,
        scc.diameter_m,
    )
    _require_shared_optical_value(
        "wavelength",
        llowfs.wavelength_um,
        scc.wavelength_um,
    )
    _require_shared_optical_value(
        "RNG seed",
        llowfs.rng_seed,
        scc.rng_seed,
    )
    _require_shared_optical_value(
        "pupil OPD schema",
        llowfs.pupil_opd_schema,
        scc.pupil_opd_schema,
    )
    _require_shared_optical_value(
        "pupil amplitude schema",
        llowfs.pupil_amplitude_schema,
        scc.pupil_amplitude_schema,
    )

    llowfs_prescription = llowfs.prescription
    scc_prescription = scc.prescription
    _require_shared_optical_value(
        "optical geometry",
        llowfs_prescription.profile.optics,
        scc_prescription.profile.optics,
    )
    _require_shared_optical_value(
        "focal-plane mask",
        llowfs_prescription.profile.focal_plane_mask,
        scc_prescription.profile.focal_plane_mask,
    )
    llowfs_parameters = llowfs_prescription.parameters
    scc_parameters = scc_prescription.parameters
    for name in (
        :beam_diameter_fraction,
        :central_tilt_cycles,
        :central_gaussian_phase_rad,
        :central_gaussian_sigma_ratio,
    )
        _require_shared_optical_value(
            String(name),
            getproperty(llowfs_parameters, name),
            getproperty(scc_parameters, name),
        )
    end
    return _ProvisionalSpidersOpticalConfiguration(llowfs, scc)
end

function _AOG.graph_node_ports(
    ::Type{_ProvisionalSpidersOpticalNode{T,LLOWFS,SCC}},
    configuration::_ProvisionalSpidersOpticalConfiguration,
) where {T,LLOWFS,SCC}
    llowfs = configuration.llowfs
    scc = configuration.scc
    pupil_shape = (llowfs.pupil_resolution, llowfs.pupil_resolution)
    return (
        _AOG.graph_port_contract(
            :pupil_opd,
            :input,
            :data,
            T,
            pupil_shape,
            llowfs.pupil_opd_schema,
            :column_major,
        ),
        _AOG.graph_port_contract(
            :pupil_amplitude,
            :input,
            :data,
            T,
            pupil_shape,
            llowfs.pupil_amplitude_schema,
            :column_major,
        ),
        _AOG.graph_port_contract(
            :llowfs_output,
            :output,
            :data,
            T,
            (llowfs.output_rows, llowfs.output_columns),
            llowfs.output_schema,
            :column_major,
        ),
        _AOG.graph_port_contract(
            :scc_output,
            :output,
            :data,
            T,
            (scc.output_rows, scc.output_columns),
            scc.output_schema,
            :column_major,
        ),
    )
end

function _validate_shared_optical_buffer(
    target,
    values,
    ::Type{T},
    shape::Tuple{Int,Int},
    role::AbstractString,
) where {T}
    eltype(values) === T || throw(ArgumentError(
        "$role element type $(eltype(values)) does not match $T",
    ))
    size(values) == shape || throw(DimensionMismatch(
        "$role shape $(size(values)) does not match $shape",
    ))
    _Backends.compute_device(values) == target || throw(ArgumentError(
        "$role target $(_Backends.compute_device(values)) does not match $target",
    ))
    return values
end

function _allocate_shared_optical_field(target, ::Type{T}, resolution) where {T}
    field = _Backends.allocate_device_array(
        target,
        Complex{T},
        resolution,
        resolution,
    )
    fill!(field, zero(Complex{T}))
    return field
end

function _AOG.prepare_graph_node(
    ::Type{_ProvisionalSpidersOpticalNode{T,LLOWFS,SCC}},
    configuration::_ProvisionalSpidersOpticalConfiguration,
    ::NamedTuple{()},
    inputs::NamedTuple{(:pupil_opd,:pupil_amplitude)},
    outputs::NamedTuple{(:llowfs_output,:scc_output)},
    ::NamedTuple{()},
    target::_Backends.AbstractComputeDevice,
) where {T,LLOWFS,SCC}
    llowfs = configuration.llowfs
    scc = configuration.scc
    pupil_shape = (llowfs.pupil_resolution, llowfs.pupil_resolution)
    _validate_shared_optical_buffer(
        target,
        inputs.pupil_opd,
        T,
        pupil_shape,
        "pupil OPD",
    )
    _validate_shared_optical_buffer(
        target,
        inputs.pupil_amplitude,
        T,
        pupil_shape,
        "pupil amplitude",
    )
    _validate_shared_optical_buffer(
        target,
        outputs.llowfs_output,
        T,
        (llowfs.output_rows, llowfs.output_columns),
        "LLOWFS output",
    )
    _validate_shared_optical_buffer(
        target,
        outputs.scc_output,
        T,
        (scc.output_rows, scc.output_columns),
        "SCC output",
    )

    common_field = _allocate_shared_optical_field(
        target,
        T,
        llowfs.resolution,
    )
    scc_field = _allocate_shared_optical_field(target, T, llowfs.resolution)
    rng = AdaptiveOpticsSim.runtime_rng(llowfs.rng_seed)
    common_context = Proper.RunContext(typeof(common_field); rng)
    scc_context = Proper.RunContext(typeof(scc_field); rng)
    common_wavefront = Proper.prop_begin!(
        common_field,
        llowfs.diameter_m,
        llowfs.wavelength_um * T(1e-6);
        beam_diam_fraction=one(T),
        context=common_context,
    )
    scc_wavefront = Proper.prop_begin!(
        scc_field,
        scc.diameter_m,
        scc.wavelength_um * T(1e-6);
        beam_diam_fraction=one(T),
        context=scc_context,
    )
    llowfs_assets = prepare_proper_assets(
        llowfs.prescription,
        target,
        T,
        llowfs.resolution,
        (llowfs.output_rows, llowfs.output_columns),
    )
    scc_assets = prepare_proper_assets(
        scc.prescription,
        target,
        T,
        scc.resolution,
        (scc.output_rows, scc.output_columns),
    )
    plan = _ProvisionalSpidersOpticalPlan(
        configuration,
        target,
        scc.prescription,
    )
    state = _ProvisionalSpidersOpticalState(
        rng,
        copy(rng),
        prepare_proper_state(scc.prescription, T),
    )
    workspace = _ProvisionalSpidersOpticalWorkspace(
        common_context,
        scc_context,
        common_field,
        scc_field,
        common_wavefront,
        scc_wavefront,
        llowfs_assets,
        scc_assets,
    )
    return _ProvisionalSpidersOpticalOwner(
        plan,
        state,
        workspace,
        inputs.pupil_opd,
        inputs.pupil_amplitude,
        outputs.llowfs_output,
        outputs.scc_output,
    )
end

@inline function _copy_spiders_wavefront!(destination, source)
    copyto!(destination.field, source.field)
    destination.wavelength_m = source.wavelength_m
    destination.sampling_m = source.sampling_m
    destination.z_m = source.z_m
    destination.beam_diameter_m = source.beam_diameter_m
    destination.z_w0_m = source.z_w0_m
    destination.w0_m = source.w0_m
    destination.z_rayleigh_m = source.z_rayleigh_m
    destination.current_fratio = source.current_fratio
    destination.reference_surface = source.reference_surface
    destination.beam_type_old = source.beam_type_old
    destination.propagator_type = source.propagator_type
    destination.rayleigh_factor = source.rayleigh_factor
    return destination
end

@inline _shared_scc_lyot_mask(
    ::ProvisionalSpidersSCCPrescription,
    state,
    assets,
) = assets.lyot_mask

@inline function _shared_scc_lyot_mask(
    prescription::ProvisionalSpidersChoppedSCCPrescription,
    state::_SpidersChopperState,
    assets,
)
    phase = _spiders_next_chopper_phase(prescription, state)
    return _spiders_chopper_lyot_mask(
        phase,
        assets.fringed_lyot_mask,
        assets.unfringed_lyot_mask,
    )
end

@inline _commit_shared_scc_frame!(
    ::ProvisionalSpidersSCCPrescription,
    state,
) = nothing

@inline _commit_shared_scc_frame!(
    ::ProvisionalSpidersChoppedSCCPrescription,
    state::_SpidersChopperState,
) = _spiders_commit_chopper_frame!(state)

@inline function _AOG.step_graph_node!(
    owner::_ProvisionalSpidersOpticalOwner,
)
    plan = owner.plan
    state = owner.state
    workspace = owner.workspace
    llowfs = plan.configuration.llowfs
    scc = plan.configuration.scc
    llowfs_assets = workspace.llowfs_assets
    scc_assets = workspace.scc_assets
    pupil_magnification = _spiders_propagate_to_lyot!(
        llowfs.prescription,
        llowfs.wavelength_um * typeof(llowfs.wavelength_um)(1e-6),
        llowfs.resolution;
        pupil_opd=owner.pupil_opd,
        pupil_amplitude=owner.pupil_amplitude,
        diameter_m=llowfs.diameter_m,
        field=workspace.common_field,
        wavefront=workspace.common_wavefront,
        run_context=workspace.common_context,
        apodizer_mask=llowfs_assets.apodizer_mask,
        focal_plane_mask_internal=llowfs_assets.focal_plane_mask_internal,
        padded_pupil_opd=llowfs_assets.padded_pupil_opd,
        padded_pupil_amplitude=llowfs_assets.padded_pupil_amplitude,
        pupil_resample_context=llowfs_assets.pupil_resample_context,
    )
    _copy_spiders_wavefront!(
        workspace.scc_wavefront,
        workspace.common_wavefront,
    )
    _spiders_propagate_from_lyot!(
        llowfs.prescription,
        pupil_magnification,
        workspace.common_wavefront;
        output=owner.llowfs_output,
        run_context=workspace.common_context,
        lyot_mask=llowfs_assets.lyot_mask,
        square_intensity=llowfs_assets.square_intensity,
        rotated_intensity=llowfs_assets.rotated_intensity,
        centered_field=llowfs_assets.centered_field,
        detector_resample_context=llowfs_assets.detector_resample_context,
    )
    scc_lyot_mask = _shared_scc_lyot_mask(
        scc.prescription,
        state.prescription,
        scc_assets,
    )
    _spiders_propagate_from_lyot!(
        scc.prescription,
        pupil_magnification,
        workspace.scc_wavefront;
        output=owner.scc_output,
        run_context=workspace.scc_context,
        lyot_mask=scc_lyot_mask,
        square_intensity=scc_assets.square_intensity,
        rotated_intensity=scc_assets.rotated_intensity,
        centered_field=scc_assets.centered_field,
        detector_resample_context=scc_assets.detector_resample_context,
    )
    _commit_shared_scc_frame!(scc.prescription, state.prescription)
    return nothing
end

function _AOG.reset_graph_node!(owner::_ProvisionalSpidersOpticalOwner)
    plan = owner.plan
    workspace = owner.workspace
    copy!(owner.state.rng, owner.state.initial_rng)
    reset_proper_state!(plan.prescription, owner.state.prescription)
    reset_proper_assets!(
        plan.configuration.llowfs.prescription,
        workspace.llowfs_assets,
    )
    reset_proper_assets!(plan.prescription, workspace.scc_assets)
    return nothing
end

@inline _AOG.graph_node_capture_capability(
    ::_ProvisionalSpidersOpticalOwner,
) = _AOG.GraphNodeCaptureUnsupported()
