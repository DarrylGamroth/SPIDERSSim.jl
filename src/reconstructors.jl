"""Run-immutable calibration contract for the GoldEye LLOWFS processor."""
struct SpidersLLOWFSPlan{T,R<:Matrix{T},I<:Matrix{T},M<:Matrix{T},S<:Matrix{T}}
    reference_image::R
    image_mask::BitMatrix
    centroid_mask::BitMatrix
    image_to_modes::I
    modes_to_actuators::M
    slopes_to_tiptilt::S
    centroid_reference::NTuple{2,T}
    minimum_total_flux::T
    minimum_centroid_flux::T
end

"""Replaceable scratch owned by one GoldEye LLOWFS execution context."""
struct SpidersLLOWFSWorkspace{T,V<:Vector{T}}
    normalized_residual::V
end

"""Caller-visible products from one GoldEye LLOWFS frame."""
mutable struct SpidersLLOWFSProducts{T,V<:Vector{T}}
    modes::V
    tiptilt::Vector{T}
    actuator_coefficients::Vector{T}
    total_flux::T
    centroid_flux::T
    valid::Bool
end

function prepare_llowfs(
    reference_image::AbstractMatrix{T},
    image_mask::AbstractMatrix{Bool},
    centroid_mask::AbstractMatrix{Bool},
    image_to_modes::AbstractMatrix{T},
    modes_to_actuators::AbstractMatrix{T},
    slopes_to_tiptilt::AbstractMatrix{T},
    centroid_reference::NTuple{2,T};
    minimum_total_flux::T=eps(T),
    minimum_centroid_flux::T=eps(T),
) where {T<:AbstractFloat}
    image_shape = size(reference_image)
    size(image_mask) == image_shape || throw(DimensionMismatch(
        "GoldEye image mask must match the reference image",
    ))
    size(centroid_mask) == image_shape || throw(DimensionMismatch(
        "GoldEye centroid mask must match the reference image",
    ))
    valid_pixels = count(image_mask)
    valid_pixels > 0 || throw(ArgumentError(
        "GoldEye image mask must select at least one pixel",
    ))
    count(centroid_mask) > 0 || throw(ArgumentError(
        "GoldEye centroid mask must select at least one pixel",
    ))
    size(image_to_modes, 1) == valid_pixels || throw(DimensionMismatch(
        "image_to_modes must have one row per selected GoldEye pixel",
    ))
    mode_count = size(image_to_modes, 2)
    size(modes_to_actuators, 1) == mode_count || throw(DimensionMismatch(
        "modes_to_actuators must have one row per reconstructed mode",
    ))
    size(slopes_to_tiptilt) == (2, 2) || throw(DimensionMismatch(
        "slopes_to_tiptilt must be a 2-by-2 matrix",
    ))
    isfinite(minimum_total_flux) && minimum_total_flux >= zero(T) ||
        throw(ArgumentError("minimum_total_flux must be finite and nonnegative"))
    isfinite(minimum_centroid_flux) && minimum_centroid_flux >= zero(T) ||
        throw(ArgumentError(
            "minimum_centroid_flux must be finite and nonnegative",
        ))
    all(isfinite, reference_image) || throw(ArgumentError(
        "GoldEye reference image must contain only finite values",
    ))
    all(isfinite, image_to_modes) || throw(ArgumentError(
        "GoldEye image-to-modes calibration must contain only finite values",
    ))
    all(isfinite, modes_to_actuators) || throw(ArgumentError(
        "GoldEye modes-to-actuators calibration must contain only finite values",
    ))
    all(isfinite, slopes_to_tiptilt) || throw(ArgumentError(
        "GoldEye slopes-to-tiptilt calibration must contain only finite values",
    ))
    all(isfinite, centroid_reference) || throw(ArgumentError(
        "GoldEye centroid reference must contain only finite values",
    ))

    plan = SpidersLLOWFSPlan(
        Matrix{T}(reference_image),
        BitMatrix(image_mask),
        BitMatrix(centroid_mask),
        Matrix{T}(image_to_modes),
        Matrix{T}(modes_to_actuators),
        Matrix{T}(slopes_to_tiptilt),
        centroid_reference,
        minimum_total_flux,
        minimum_centroid_flux,
    )
    workspace = SpidersLLOWFSWorkspace(zeros(T, valid_pixels))
    products = SpidersLLOWFSProducts(
        zeros(T, mode_count),
        zeros(T, 2),
        zeros(T, size(modes_to_actuators, 2)),
        zero(T),
        zero(T),
        false,
    )
    return plan, workspace, products
end

@inline function _invalidate_llowfs!(products::SpidersLLOWFSProducts{T}) where {T}
    fill!(products.modes, zero(T))
    fill!(products.tiptilt, zero(T))
    fill!(products.actuator_coefficients, zero(T))
    products.valid = false
    return products
end

function process_llowfs!(
    products::SpidersLLOWFSProducts{T},
    workspace::SpidersLLOWFSWorkspace{T},
    plan::SpidersLLOWFSPlan{T},
    image::AbstractMatrix{T},
) where {T}
    Base.require_one_based_indexing(image)
    size(image) == size(plan.reference_image) || throw(DimensionMismatch(
        "GoldEye LLOWFS frame shape changed after preparation",
    ))

    total_flux = zero(T)
    centroid_flux = zero(T)
    centroid_x = zero(T)
    centroid_y = zero(T)
    @inbounds for column in axes(image, 2), row in axes(image, 1)
        value = image[row, column]
        total_flux += value
        if plan.centroid_mask[row, column]
            centroid_flux += value
            centroid_x += T(row) * value
            centroid_y += T(column) * value
        end
    end
    products.total_flux = total_flux
    products.centroid_flux = centroid_flux
    if !(isfinite(total_flux) && total_flux > plan.minimum_total_flux &&
         isfinite(centroid_flux) &&
         centroid_flux > plan.minimum_centroid_flux)
        return _invalidate_llowfs!(products)
    end

    selected = 0
    @inbounds for index in eachindex(plan.image_mask, image, plan.reference_image)
        plan.image_mask[index] || continue
        selected += 1
        workspace.normalized_residual[selected] =
            image[index] / total_flux - plan.reference_image[index]
    end
    mul!(
        products.modes,
        transpose(plan.image_to_modes),
        workspace.normalized_residual,
    )

    slope_x = centroid_x / centroid_flux - plan.centroid_reference[1]
    slope_y = centroid_y / centroid_flux - plan.centroid_reference[2]
    products.tiptilt[1] =
        plan.slopes_to_tiptilt[1, 1] * slope_x +
        plan.slopes_to_tiptilt[1, 2] * slope_y
    products.tiptilt[2] =
        plan.slopes_to_tiptilt[2, 1] * slope_x +
        plan.slopes_to_tiptilt[2, 2] * slope_y
    mul!(
        products.actuator_coefficients,
        transpose(plan.modes_to_actuators),
        products.modes,
    )
    products.valid = all(isfinite, products.actuator_coefficients)
    if !products.valid
        _invalidate_llowfs!(products)
    end
    return products
end

"""Run-immutable calibration contract for C-RED 2 SCC reconstruction."""
struct SpidersSCCPlan{T,R<:Matrix{T},M<:Matrix{T}}
    reference_difference::R
    image_mask::BitMatrix
    correction_mask::BitMatrix
    image_to_actuators::M
end

"""Replaceable scratch owned by one C-RED 2 SCC execution context."""
struct SpidersSCCWorkspace{T,V<:Vector{T}}
    selected_residual::V
end

"""Caller-visible products from one C-RED 2 chopped difference."""
mutable struct SpidersSCCProducts{T,V<:Vector{T}}
    actuator_coefficients::V
    valid::Bool
end

function prepare_scc(
    reference_difference::AbstractMatrix{T},
    image_mask::AbstractMatrix{Bool},
    correction_mask::AbstractMatrix{Bool},
    image_to_actuators::AbstractMatrix{T},
) where {T<:AbstractFloat}
    image_shape = size(reference_difference)
    size(image_mask) == image_shape || throw(DimensionMismatch(
        "C-RED 2 image mask must match the reference difference",
    ))
    size(correction_mask) == image_shape || throw(DimensionMismatch(
        "C-RED 2 correction mask must match the reference difference",
    ))
    valid_pixels = count(image_mask)
    valid_pixels > 0 || throw(ArgumentError(
        "C-RED 2 image mask must select at least one pixel",
    ))
    count(correction_mask) > 0 || throw(ArgumentError(
        "C-RED 2 correction mask must select at least one pixel",
    ))
    for index in eachindex(image_mask, correction_mask)
        if correction_mask[index] && !image_mask[index]
            throw(ArgumentError(
                "C-RED 2 correction mask must be a subset of the image mask",
            ))
        end
    end
    size(image_to_actuators, 1) == valid_pixels || throw(DimensionMismatch(
        "image_to_actuators must have one row per selected C-RED 2 pixel",
    ))
    all(isfinite, reference_difference) || throw(ArgumentError(
        "C-RED 2 reference difference must contain only finite values",
    ))
    all(isfinite, image_to_actuators) || throw(ArgumentError(
        "C-RED 2 image-to-actuators calibration must contain only finite values",
    ))
    plan = SpidersSCCPlan(
        Matrix{T}(reference_difference),
        BitMatrix(image_mask),
        BitMatrix(correction_mask),
        Matrix{T}(image_to_actuators),
    )
    workspace = SpidersSCCWorkspace(zeros(T, valid_pixels))
    products = SpidersSCCProducts(
        zeros(T, size(image_to_actuators, 2)),
        false,
    )
    return plan, workspace, products
end

function process_scc!(
    products::SpidersSCCProducts{T},
    workspace::SpidersSCCWorkspace{T},
    plan::SpidersSCCPlan{T},
    chopped_difference::AbstractMatrix{T},
) where {T}
    Base.require_one_based_indexing(chopped_difference)
    size(chopped_difference) == size(plan.reference_difference) ||
        throw(DimensionMismatch(
            "C-RED 2 SCC difference shape changed after preparation",
        ))
    selected = 0
    @inbounds for index in eachindex(
        plan.image_mask,
        plan.correction_mask,
        chopped_difference,
        plan.reference_difference,
    )
        plan.image_mask[index] || continue
        selected += 1
        workspace.selected_residual[selected] = if plan.correction_mask[index]
            chopped_difference[index] - plan.reference_difference[index]
        else
            zero(T)
        end
    end
    mul!(
        products.actuator_coefficients,
        transpose(plan.image_to_actuators),
        workspace.selected_residual,
    )
    products.valid = all(isfinite, products.actuator_coefficients)
    if !products.valid
        fill!(products.actuator_coefficients, zero(T))
    end
    return products
end
