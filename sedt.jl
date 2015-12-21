positives{T <: Number}(a::Array{T}) = map(x -> max(x, zero(T)), a)
negatives{T <: Number}(a::Array{T}) = map(x -> max(-x, zero(T)), a)


combine{T <: Number}(neg::Array{T}, pos::Array{T}) = pos - neg


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


function compute{T}(input::Array{T}, inside = x -> x > zero(T))
    output = Array(Float64, size(input)...)
    (xd, yd, zd) = size(input)

    for y in 1:yd
        for z in 1:zd
            a = input[:,y,z]
            output[:,y,z] = combine(init(map(x -> !inside(x), a)),
                                    init(map(inside, a)))
        end
    end

    for x in 1:xd
        for z in 1:zd
            a = vec(output[x,:,z])
            output[x,:,z] = combine(propagate(negatives(a)),
                                    propagate(positives(a)))
        end
    end

    for x in 1:xd
        for y in 1:yd
            a = vec(output[x,y,:])
            output[x,y,:] = combine(propagate(negatives(a)),
                                    propagate(positives(a)))
        end
    end

    return output
end
