#SIR-X Model


using DifferentialEquations


beta = 0.215; gamma = 0.07; N = 7000;
eta0 = 0.003; eta = 0.003

"""SIR-X Model
Extension of the basic SIR model

"""
function SIR_X!(du,u,p,t)
    du[1] = -beta* u[1]/N*u[2] - eta0*u[1]                          #S
    du[2] = beta*u[1]/N*u[2] - gamma*u[2] - eta0*u[2] - eta*u[2]    #I
    du[3] = gamma*u[2] + eta0*u[1]                                  #R
    du[4] = (eta0 + eta)*u[2]                                       #X
end
   
   
   u0 = [N-15.0;15.0;0.0;0.0]
   tspan = (0.0,100.0)
   prob = ODEProblem(SIR_X!,u0,tspan)
   sol = solve(prob)
   
   plot(sol)
   