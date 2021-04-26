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