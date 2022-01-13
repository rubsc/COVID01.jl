# Basics of the mathematica SIR Model


using DifferentialEquations
using Plots

# This is a test for the Lorenz ODE system
function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end


u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
sol = solve(prob)

plot(sol,vars=(1,2,3))

# Now the basic SIR ODE system

beta = 0.25
gamma = 0.07
N = 7000
function SIR_basic!(du,u,p,t)
    du[1] = -beta* u[1]/N*u[2]
    du[2] = beta* u[1]/N*u[2] - gamma*u[2]
    du[3] = gamma*u[2]
end
   
   
   u0 = [N-15.0;15.0;0.0]
   tspan = (0.0,100.0)
   prob = ODEProblem(SIR_basic!,u0,tspan)
   sol = solve(prob)
   
   plot(sol,vars=(1,2,3))
   