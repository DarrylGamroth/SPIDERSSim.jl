"""Run-immutable C-RED 2 frame-pairing contract."""
struct SpidersSCCPairPlan
    frame_shape::Tuple{Int,Int}
end

"""Persistent most-recent fringed and unfringed C-RED 2 frames."""
mutable struct SpidersSCCPairState{T,A<:Matrix{T}}
    fringed::A
    unfringed::A
    fringed_sequence::UInt64
    unfringed_sequence::UInt64
    last_sequence::UInt64
    have_fringed::Bool
    have_unfringed::Bool
end

"""Caller-visible rolling fringed-minus-unfringed SCC product."""
mutable struct SpidersSCCPairProducts{T,A<:Matrix{T}}
    difference::A
    sequence::UInt64
    fringed_sequence::UInt64
    unfringed_sequence::UInt64
    valid::Bool
end

function prepare_scc_pairing(
    ::Type{T}=Float32;
    frame_shape::Tuple{Int,Int}=(416, 380),
) where {T<:AbstractFloat}
    all(>(0), frame_shape) || throw(ArgumentError(
        "C-RED 2 SCC frame dimensions must be positive",
    ))
    plan = SpidersSCCPairPlan(frame_shape)
    state = SpidersSCCPairState(
        zeros(T, frame_shape),
        zeros(T, frame_shape),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        false,
        false,
    )
    products = SpidersSCCPairProducts(
        zeros(T, frame_shape),
        UInt64(0),
        UInt64(0),
        UInt64(0),
        false,
    )
    return plan, state, products
end

function reset_scc_pairing!(
    state::SpidersSCCPairState{T},
    products::SpidersSCCPairProducts{T},
) where {T}
    fill!(state.fringed, zero(T))
    fill!(state.unfringed, zero(T))
    state.fringed_sequence = UInt64(0)
    state.unfringed_sequence = UInt64(0)
    state.last_sequence = UInt64(0)
    state.have_fringed = false
    state.have_unfringed = false
    fill!(products.difference, zero(T))
    products.sequence = UInt64(0)
    products.fringed_sequence = UInt64(0)
    products.unfringed_sequence = UInt64(0)
    products.valid = false
    return products
end

function accept_scc_frame!(
    products::SpidersSCCPairProducts{T},
    state::SpidersSCCPairState{T},
    plan::SpidersSCCPairPlan,
    frame::AbstractMatrix{T},
    sequence::UInt64,
    phase::SpidersChopperPhase,
) where {T}
    Base.require_one_based_indexing(frame)
    size(frame) == plan.frame_shape || throw(DimensionMismatch(
        "C-RED 2 SCC frame shape changed after pairing preparation",
    ))
    sequence > state.last_sequence || throw(ArgumentError(
        "C-RED 2 SCC frame sequence must increase monotonically",
    ))

    if state.last_sequence != UInt64(0) && sequence - state.last_sequence != UInt64(1)
        state.have_fringed = false
        state.have_unfringed = false
    end
    if phase === SpidersFringed
        copyto!(state.fringed, frame)
        state.fringed_sequence = sequence
        state.have_fringed = true
    else
        copyto!(state.unfringed, frame)
        state.unfringed_sequence = sequence
        state.have_unfringed = true
    end
    state.last_sequence = sequence

    adjacent = state.have_fringed && state.have_unfringed &&
        abs(Int128(state.fringed_sequence) - Int128(state.unfringed_sequence)) == 1
    if adjacent
        @. products.difference = state.fringed - state.unfringed
        products.sequence = sequence
        products.fringed_sequence = state.fringed_sequence
        products.unfringed_sequence = state.unfringed_sequence
        products.valid = all(isfinite, products.difference)
    else
        fill!(products.difference, zero(T))
        products.sequence = sequence
        products.fringed_sequence = state.fringed_sequence
        products.unfringed_sequence = state.unfringed_sequence
        products.valid = false
    end
    if !products.valid
        fill!(products.difference, zero(T))
    end
    return products
end
