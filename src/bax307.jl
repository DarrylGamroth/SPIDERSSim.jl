const BAX307_ACTUATOR_COUNT = 468
const BAX307_GRID_SIZE = 24

const BAX307_PDM_COMMAND_SCHEMA =
    "org.subaru.spiders.bax307-normalized-command/1"
const BAX307_SURFACE_OPD_SCHEMA =
    "org.subaru.spiders.bax307-surface-opd-m/1"
const BAX307_ACTUATOR_COORDINATES_SCHEMA =
    "org.subaru.spiders.bax307-nominal-grid-coordinates/1"
const BAX307_INFLUENCE_FUNCTIONS_SCHEMA =
    "org.subaru.spiders.bax307-measured-influence-opd-m/1"

"""
    BAX307CalibrationPlan

Run-immutable numerical calibration for the SPIDERS ALPAO BAX307 deformable
mirror. `influence_samples` stores one column per command coordinate and one
row per true element of `influence_support`, in the support mask's native
column-major order. Samples remain in their source calibration unit; the
explicit `influence_sample_unit_m` converts them to mirror-surface metres.

The plan snapshots caller arrays during preparation. It does not own a best
flat, live command, command contributors, device storage, or a FITS file.
"""
struct BAX307CalibrationPlan{T<:Union{Float32,Float64}}
    actuator_map::BitMatrix
    valid_actuator_map::BitMatrix
    influence_support::BitMatrix
    influence_samples::Matrix{T}
    influence_sample_unit_m::T
    surface_to_opd::T
    command_limit::T
end

function _bax307_map(
    values::AbstractMatrix{<:Real},
    role::AbstractString,
    expected_shape::Tuple{Int,Int},
)
    Base.require_one_based_indexing(values)
    size(values) == expected_shape || throw(DimensionMismatch(
        "$role shape $(size(values)) does not match $expected_shape",
    ))
    all(isfinite, values) || throw(ArgumentError("$role must be finite"))
    return BitMatrix(values .!= zero(eltype(values)))
end

function _bax307_influence_samples(
    values::AbstractMatrix{<:Real},
    support_count::Int,
    ::Type{T},
) where {T<:Union{Float32,Float64}}
    Base.require_one_based_indexing(values)
    expected = (support_count, BAX307_ACTUATOR_COUNT)
    samples = if size(values) == expected
        Matrix{T}(values)
    elseif size(values) == reverse(expected)
        Matrix{T}(transpose(values))
    else
        throw(DimensionMismatch(
            "BAX307 influence samples shape $(size(values)) must be " *
            "$expected or $(reverse(expected))",
        ))
    end
    all(isfinite, samples) || throw(ArgumentError(
        "BAX307 influence samples must be finite",
    ))
    return samples
end

"""
    prepare_bax307_calibration(actuator_map, valid_actuator_map,
                               influence_support, influence_samples;
                               influence_sample_unit_m, ...)

Validate and snapshot one BAX307 calibration bundle. The source influence unit
must be supplied explicitly because the available manufacturer FITS product
does not carry a trustworthy physical-unit declaration. For a reflective DM,
the default `surface_to_opd=2` converts mirror displacement to optical path
difference.
"""
function prepare_bax307_calibration(
    actuator_map::AbstractMatrix{<:Real},
    valid_actuator_map::AbstractMatrix{<:Real},
    influence_support::AbstractMatrix{<:Real},
    influence_samples::AbstractMatrix{<:Real};
    influence_sample_unit_m::Real,
    surface_to_opd::Real=2,
    command_limit::Real=0.5,
    T::Type{<:Union{Float32,Float64}}=Float32,
)
    Base.require_one_based_indexing(influence_support)
    actuators = _bax307_map(
        actuator_map,
        "BAX307 actuator map",
        (BAX307_GRID_SIZE, BAX307_GRID_SIZE),
    )
    count(actuators) == BAX307_ACTUATOR_COUNT || throw(ArgumentError(
        "BAX307 actuator map contains $(count(actuators)) commands; " *
        "expected $BAX307_ACTUATOR_COUNT",
    ))
    valid = _bax307_map(
        valid_actuator_map,
        "BAX307 valid-actuator map",
        size(actuators),
    )
    all(valid .<= actuators) || throw(ArgumentError(
        "BAX307 valid-actuator map includes cells outside the actuator map",
    ))
    any(valid) || throw(ArgumentError(
        "BAX307 valid-actuator map must contain at least one actuator",
    ))

    size(influence_support, 1) == size(influence_support, 2) ||
        throw(DimensionMismatch(
            "BAX307 influence support must be square, got " *
            "$(size(influence_support))",
        ))
    all(isfinite, influence_support) || throw(ArgumentError(
        "BAX307 influence support must be finite",
    ))
    support = BitMatrix(
        influence_support .!= zero(eltype(influence_support)),
    )
    support_count = count(support)
    support_count > 0 || throw(ArgumentError(
        "BAX307 influence support must contain at least one sample",
    ))
    samples = _bax307_influence_samples(
        influence_samples,
        support_count,
        T,
    )

    unit_m = T(influence_sample_unit_m)
    isfinite(unit_m) && unit_m > zero(T) || throw(ArgumentError(
        "influence_sample_unit_m must be finite and positive",
    ))
    opd_factor = T(surface_to_opd)
    isfinite(opd_factor) && opd_factor > zero(T) || throw(ArgumentError(
        "surface_to_opd must be finite and positive",
    ))
    limit = T(command_limit)
    isfinite(limit) && limit > zero(T) || throw(ArgumentError(
        "command_limit must be finite and positive",
    ))
    return BAX307CalibrationPlan(
        actuators,
        valid,
        support,
        samples,
        unit_m,
        opd_factor,
        limit,
    )
end

"""Return the current valid/dead policy in exact 468-command order."""
function bax307_valid_command_mask(plan::BAX307CalibrationPlan)
    return Vector{Bool}(plan.valid_actuator_map[plan.actuator_map])
end

"""
Return nominal coordinates in the BAX307 command order. Coordinates use the
first/fast map axis as x, the second axis as y, a centered 24-by-24 grid, and
23 actuator pitches across the normalized pupil diameter. They preserve
command topology but do not claim an as-built rotation, parity, or decenter.
"""
function bax307_actuator_coordinates(
    plan::BAX307CalibrationPlan{T},
) where {T}
    pitch = T(2) / T(BAX307_GRID_SIZE - 1)
    center = T(BAX307_GRID_SIZE + 1) / T(2)
    coordinates = Matrix{T}(undef, 2, BAX307_ACTUATOR_COUNT)
    command = 1
    for index in CartesianIndices(plan.actuator_map)
        plan.actuator_map[index] || continue
        coordinates[1, command] = (T(index[1]) - center) * pitch
        coordinates[2, command] = (T(index[2]) - center) * pitch
        command += 1
    end
    command == BAX307_ACTUATOR_COUNT + 1 || error(
        "validated BAX307 actuator ordering changed unexpectedly",
    )
    return coordinates
end

"""
Materialize the measured influence functions as full-grid reflected OPD in
metres per normalized BAX307 command unit. The returned matrix has shape
`(resolution^2, 468)` using the support grid's native column-major order.
"""
function bax307_influence_functions_opd(
    plan::BAX307CalibrationPlan{T},
) where {T}
    resolution = size(plan.influence_support, 1)
    modes = zeros(T, resolution * resolution, BAX307_ACTUATOR_COUNT)
    support = vec(plan.influence_support)
    scale = plan.influence_sample_unit_m * plan.surface_to_opd
    @views modes[support, :] .= plan.influence_samples .* scale
    return modes
end

"""Evaluate one complete BAX307 command into a host surface-OPD matrix."""
function bax307_surface_opd(
    plan::BAX307CalibrationPlan{T},
    command::AbstractVector{<:Real},
) where {T}
    Base.require_one_based_indexing(command)
    length(command) == BAX307_ACTUATOR_COUNT || throw(DimensionMismatch(
        "BAX307 command has $(length(command)) elements; expected " *
        "$BAX307_ACTUATOR_COUNT",
    ))
    all(isfinite, command) || throw(ArgumentError(
        "BAX307 command must be finite",
    ))
    applied = clamp.(T.(command), -plan.command_limit, plan.command_limit)
    applied .*= bax307_valid_command_mask(plan)
    resolution = size(plan.influence_support, 1)
    surface = bax307_influence_functions_opd(plan) * applied
    return reshape(surface, resolution, resolution)
end

"""Construct the BAX307 clipping and dead-actuator response model."""
function bax307_actuator_model(plan::BAX307CalibrationPlan{T}) where {T}
    health = T.(bax307_valid_command_mask(plan))
    return AdaptiveOpticsSim.Optics.CompositeDMActuatorModel(
        AdaptiveOpticsSim.Optics.ClippedActuators(
            -plan.command_limit,
            plan.command_limit,
        ),
        AdaptiveOpticsSim.Optics.ActuatorHealthMap(health),
    )
end

"""Declare the AOS measured-influence node for one BAX307 plan."""
function bax307_deformable_mirror_node(
    name::Symbol,
    plan::BAX307CalibrationPlan{T};
    pupil_diameter_m::Real,
    pdm_command_schema::AbstractString=BAX307_PDM_COMMAND_SCHEMA,
    surface_opd_schema::AbstractString=BAX307_SURFACE_OPD_SCHEMA,
) where {T}
    resolution = size(plan.influence_support, 1)
    return _AOG.deformable_mirror_surface_node(
        name;
        resolution,
        telescope_diameter_m=pupil_diameter_m,
        actuator_count=BAX307_ACTUATOR_COUNT,
        actuator_model=bax307_actuator_model(plan),
        pdm_command_schema,
        surface_opd_schema,
        actuator_coordinates_schema=BAX307_ACTUATOR_COORDINATES_SCHEMA,
        influence_functions_schema=BAX307_INFLUENCE_FUNCTIONS_SCHEMA,
        T,
    )
end

"""Return the startup sparse parameters for a named BAX307 graph node."""
function bax307_graph_parameters(
    name::Symbol,
    plan::BAX307CalibrationPlan,
)
    return (
        _AOG.sparse_parameter(
            name => :actuator_coordinates,
            bax307_actuator_coordinates(plan),
        ),
        _AOG.sparse_parameter(
            name => :influence_functions,
            bax307_influence_functions_opd(plan),
        ),
    )
end

function _bax307_target_array(target, source::Array{T,N}) where {T,N}
    destination = _Backends.allocate_device_array(target, T, size(source)...)
    copyto!(destination, source)
    return destination
end

"""Return BAX307 startup parameters copied to the selected graph target."""
function bax307_graph_parameters(
    name::Symbol,
    plan::BAX307CalibrationPlan,
    target,
)
    actuator_coordinates = _bax307_target_array(
        target,
        bax307_actuator_coordinates(plan),
    )
    influence_functions = _bax307_target_array(
        target,
        bax307_influence_functions_opd(plan),
    )
    return (
        _AOG.sparse_parameter(
            name => :actuator_coordinates,
            actuator_coordinates,
        ),
        _AOG.sparse_parameter(
            name => :influence_functions,
            influence_functions,
        ),
    )
end

"""
Load a BAX307 calibration from instrument files. Load `FITSIO` first to enable
the optional SPIDERSSim FITS extension.
"""
function load_bax307_calibration end

"""
Load one 468-element BAX307 best-flat command. Load `FITSIO` first to enable
the optional SPIDERSSim FITS extension.
"""
function load_bax307_best_flat end
