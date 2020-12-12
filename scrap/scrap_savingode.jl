# Does saving ODE solution take time?

using BenchmarkTools
using OrdinaryDiffEq

function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end

function testsave()
    u0 = [1.0;0.0;0.0]
    tspan = (0.0,100.0)
    prob = ODEProblem(lorenz!,u0,tspan)
    sol = solve(prob, Tsit5())
end
@benchmark testsave()

# vs

function testnosave()
    u0 = [1.0;0.0;0.0]
    tspan = (0.0,100.0)
    prob = ODEProblem(lorenz!,u0,tspan)
    sol = solve(prob, Tsit5(), save_on=false, save_start=false, save_end=true)
end
@benchmark testnosave()

# yes, nosave is faster by approx 2x
