module SPIDERSSimFITSIOExt

using FITSIO
using SPIDERSSim

function _read_fits(path::AbstractString, hdu::Integer=1)
    isfile(path) || throw(ArgumentError("FITS file not found: $path"))
    file = FITS(path, "r")
    try
        return read(file[Int(hdu)])
    finally
        close(file)
    end
end

function SPIDERSSim.load_bax307_calibration(
    data_root::AbstractString;
    influence_sample_unit_m::Real,
    surface_to_opd::Real=2,
    command_limit::Real=0.5,
    T::Type{<:Union{Float32,Float64}}=Float32,
    actuator_map_path::AbstractString=joinpath(
        data_root,
        "dm",
        "BAX307-actu-map.fits",
    ),
    valid_actuator_map_path::AbstractString=joinpath(
        data_root,
        "dm",
        "BAX307-valid-actu-map.fits",
    ),
    influence_path::AbstractString=joinpath(
        data_root,
        "dm",
        "AX307_Influences.fits",
    ),
)
    actuator_map = _read_fits(actuator_map_path)
    valid_actuator_map = _read_fits(valid_actuator_map_path)
    influence_samples = _read_fits(influence_path, 1)
    influence_support = _read_fits(influence_path, 2)
    return SPIDERSSim.prepare_bax307_calibration(
        actuator_map,
        valid_actuator_map,
        influence_support,
        influence_samples;
        influence_sample_unit_m,
        surface_to_opd,
        command_limit,
        T,
    )
end

function SPIDERSSim.load_bax307_best_flat(
    data_root::AbstractString;
    path::AbstractString=joinpath(data_root, "dm", "bestflat.fits"),
    T::Type{<:Union{Float32,Float64}}=Float32,
)
    values = vec(_read_fits(path))
    length(values) == SPIDERSSim.BAX307_ACTUATOR_COUNT ||
        throw(DimensionMismatch(
            "BAX307 best flat has $(length(values)) commands; expected " *
            "$(SPIDERSSim.BAX307_ACTUATOR_COUNT)",
        ))
    all(isfinite, values) || throw(ArgumentError(
        "BAX307 best flat must be finite",
    ))
    converted = Vector{T}(values)
    all(isfinite, converted) || throw(ArgumentError(
        "BAX307 best flat is not finite after conversion to $T",
    ))
    return converted
end

end
