function init(a::Array{Bool, 1})
    n = size(a)[1]

    out = zeros(Float64, n)
    a[1] && (out[1] = n)

    for i in 2:n
        a[i] && (out[i] = 1 + out[i - 1])
    end

    for i in n-1:-1:1
        out[i] > out[i + 1] && (out[i] = 1 + out[i + 1])
    end

    return out
end


function propagate{T <: Real}(a::Array{T, 1})
    m = size(a)[1]
    m <= 1 && return a

    g(i::Int) = a[i + 1]
    f(u, i::Int) = (u - i)^2 + g(i)^2

    q = 1
    s = zeros(Int64, m)
    t = zeros(T, m)

    s[q] = 0
    t[q] = 0

    for u in 1:m-1
        while q > 0 && f(t[q], s[q]) > f(t[q], u)
            q -= 1
        end

        if q < 1
            q = 1
            s[q] = u
            t[q] = 0
        else
            i = s[q]
            w = (u^2 - i^2 + g(u)^2 - g(i)^2) / 2(u - i)
            if w < m
                q += 1
                s[q] = u
                t[q] = w
            end
        end
    end

    result = zeros(T, m)
    for u in m:-1:1
        while q > 1 && u - 1 < t[q]
            q -= 1
        end
        result[u] = sqrt(f(u - 1, s[q]))
    end

    return result
end


function adjust{T <: Real, n}(a::AbstractArray{T, n})
    d = size(a)[1]
    out = zeros(Float64, size(a)...)

    for i in 1:d
        idcs = map(k -> k == 1 ? i : Colon(), 1:n)
        out[idcs...] = adjust(slice(a, idcs...))
    end

    return out
end


function adjust{T <: Real}(a::AbstractArray{T, 1})
    return (propagate(map(x -> max( x, zero(T)), a)) -
            propagate(map(x -> max(-x, zero(T)), a)))
end


function compute{T, n}(a::Array{T, n}, inside = x -> x > zero(T))
    d = size(a)[end]
    out = zeros(Float64, size(a)...)

    for i in 1:d
        idcs = map(k -> k == n ? i : Colon(), 1:n)
        out[idcs...] = compute(a[idcs...], inside)
    end

    return adjust(out)
end


function compute{T}(a::Array{T, 1}, inside = x -> x > zero(T))
    return init(map(inside, a)) - init(map(x -> !inside(x), a))
end
