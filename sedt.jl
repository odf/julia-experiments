# Copyright 2015 The Australian National University
#
# Performs a signed Euclidean distance transform on an array of arbitrary
# dimension using the Hirata/Meijster algorithm.
#
# Olaf Delgado-Friedrichs dec 15


import Iterators


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


function propagate{T <: Real}(a::Array{T, 1})
    m = size(a)[1]
    m <= 1 && return a

    f(u, i::Int) = (u - i)^2 + a[i]^2

    q = 1
    s = zeros(Int64, m)
    t = zeros(T, m)

    s[q] = 1
    t[q] = 0

    for u in 2:m
        while q > 0 && f(t[q], s[q]) > f(t[q], u)
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
        result[u] = sqrt(f(u, s[q]))
    end

    return result
end


function slivers(a::Array, k::Int)
    dims = size(a)
    idcs = map(i -> i == k ? 0 : 1:dims[i], 1:length(dims))
    iter = Iterators.product(idcs...)
    fix  = k -> k == 0 ? Colon() : k

    return Iterators.imap(i -> map(fix, i), iter)
end


function go{T}(f::Function, a::Array{T})
    z = zero(T)
    return f(map(x -> max(x, z), a)) - f(map(x -> max(-x, z), a))
end


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
