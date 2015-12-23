# Copyright 2015 The Australian National University
#
# Performs a signed Euclidean distance transform on an array of arbitrary
# dimension using the Hirata/Meijster algorithm.
#
# Olaf Delgado-Friedrichs dec 15


function init{T <: Real}(a::Array{T, 1})
    n = size(a)[1]

    out = zeros(T, n)
    a[1] > 0.0 && (out[1] = n)

    for i in 2:n
        a[i] > 0.0 && (out[i] = 1 + out[i - 1])
    end

    for i in n-1:-1:1
        out[i] > out[i + 1] && (out[i] = 1 + out[i + 1])
    end

    return out
end


function f{T}(a::Array{T, 1}, u, i)
    return (u - i)^2 + a[i]^2
end


function propagate{T <: Real}(a::Array{T, 1})
    m = size(a)[1]
    m <= 1 && return a

    q = 1
    s = zeros(Int64, m)
    t = zeros(T, m)

    s[q] = 1
    t[q] = 0

    for u in 2:m
        while q > 0 && f(a, t[q], s[q]) > f(a, t[q], u)
            q -= 1
        end

        if q < 1
            q = 1
            s[q] = u
            t[q] = 0
        else
            i = s[q]
            w = (u^2 - i^2 + a[u]^2 - a[i]^2) / 2(u - i)
            if w < m
                q += 1
                s[q] = u
                t[q] = w
            end
        end
    end

    result = zeros(T, m)
    for u in m:-1:1
        while q > 1 && u < t[q]
            q -= 1
        end
        result[u] = sqrt(f(a, u, s[q]))
    end

    return result
end


immutable Slivers
    dims::Array{Int64,1}
    k::Int64
end


function Base.start(it::Slivers)
    state = ones(Int64, length(it.dims))
    state[it.k] = 0
    return state
end


function Base.next(it::Slivers, state::Array{Int64, 1})
    n = length(it.dims)

    i = n
    while i > 0 && (i == it.k || state[i] == it.dims[i])
        i -= 1
    end

    if i == 0
        newState = Int64[]
    else
        newState = copy(state)
        newState[i] += 1

        for j in i+1:n
            if j != it.k
                newState[j] = 1
            end
        end
    end

    return (map(x -> x == 0 ? Colon() : x, state), newState)
end


function Base.done(it::Slivers, state::Array{Int64, 1})
    return length(state) == 0
end


slivers(a::Array, k::Int) = Slivers(collect(size(a)), k)


pos(x) = max( x, zero(typeof(x)))
neg(x) = max(-x, zero(typeof(x)))

go(f::Function, a::Array) = f(map(pos, a)) - f(map(neg, a))


function compute{T, n}(a::Array{T, n}, inside = x -> x > zero(T))
    out = map(x -> inside(x) ? 1.0 : -1.0, a)

    for idcs in slivers(out, 1)
        out[idcs...] = go(init, vec(out[idcs...]))
    end

    for k in 2:n
        for idcs in slivers(out, k)
            out[idcs...] = go(propagate, vec(out[idcs...]))
        end
    end

    for i in eachindex(out)
        if out[i] > 0
            out[i] -= 0.5
        elseif out[i] < 0
            out[i] += 0.5
        end
    end

    return out
end
