module Donquicote

using PrettyTables, Statistics, StatsBase

function newton(f, ∂f, x0; xatol = 1e-7)

    diff = 1
    while diff >= xatol
        x1 = x0 - f(x0) / ∂f(x0);
        diff = abs(x1 - x0);
        x0 = x1;
    end
    return x0
end

function haley(f, ∂f, ∂2f, x0; xatol = 1e-7)
    diff = 1
    while diff >= xatol
        x1 = x0 - (2*f(x0)*∂f(x0)) / (2*∂f(x0)^2 - f(x0)*∂2f(x0));
        diff = abs(x1 - x0);
        x0 = x1;
    end
    return x0

end

function schroeder(f, ∂f, ∂2f, x0; xatol = 1e-7)
    diff = 1
    while diff >= xatol
        x1 = x0 - f(x0) / ∂f(x0)  - (∂2f(x0)*f(x0)^2)/(2*∂f(x0)^3);
        diff = abs(x1 - x0);
        x0 = x1;
    end
    return x0

end

function weightedmean(x, w)
    wxsum = wsum = 0
    for (x,w) in zip(x,w)
        wx = w*x
        if !ismissing(wx) | isfinite(wx)
            wxsum += wx
            wsum += w
        end
    end
    return wxsum / wsum
end

function weightedstd(x, w)
    μ = weightedmean(x, w)
    wxsum = wsum = 0
    for (x,w) in zip(x,w)
        wx = w*(x - μ)^2
        if !ismissing(wx) | isfinite(wx)
            wxsum += wx
            wsum += w
        end
    end
    σ = sqrt((1/wsum)*wxsum)
    return σ
end

function gini(x::Vector{Float64})
    n = length(x)
    sarray = sort(x)
    return 2*(sum(collect(1:n).*sarray))/(n*sum(sarray))-1
end

logm(x) = x > 0.0 ? log(x) : missing

lvar(x) = var(logm.(x))

nanmean(x) = mean(filter(!isnan,x))

format(x) = round(x, digits = 3)

function savehtml(filename, data)
    open("$filename.html", "w") do f
        pretty_table(f, data, backend = :html)
    end
end

export newton, haley, schroeder
export gini, nanmean, weightedmean, weightedstd, format, logm, lvar
export ones, zeros
export savehtml

end