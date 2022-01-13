# SEIR Model

using DifferentialEquations


beta = 0.215; gamma = 0.07; N = 7000;
alpha = 0.25;

"""SEIR Model
Extension of the basic SIR model

"""
function SEIR!(du,u,p,t)
    du[1] = -beta* u[1]/N*u[3]
    du[2] = beta*u[1]/N*u[3] - alpha*u[2]
    du[3] = alpha* u[2] - gamma*u[3]
    du[4] = gamma*u[3]
end
   
   
   u0 = [N-15.0;15.0;0.0;0.0]
   tspan = (0.0,100.0)
   prob = ODEProblem(SEIR!,u0,tspan)
   sol = solve(prob)
   
   plot(sol)
   