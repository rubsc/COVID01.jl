# A stochastic extension to the SIR and SEIR model


"""
We consider the 2 dimensional system of SDEs
 
dX_t = F(X_t) dt + G(X_t) dW_t

where X_t = (S_t; I_t) and F:R^2 -> R^2 and G(⋅) ∈ R^4.

The component R for the standard SIR model need not be considered.
"""

using DifferentialEquations
###########################################
#First some testing
f(du,u,p,t) = du .= 1.01u
function g(du,u,p,t)
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[2,1] = 1.2u[2]
  du[2,2] = 0.2u[2]
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,2))
sol = solve(prob,LambaEM())
plot(sol)



############################################################
beta = 0.25; gamma = 0.07; N = 7000;
function f(du,u,p,t)
    du[1] = -beta * u[2]*u[1]/N;
    du[2] = beta*u[1]*u[2]/N - gamma*u[2]
end

function g(du,u,p,t)
    du[1,1] = -sqrt(beta*u[1]*u[2]/N)
    du[1,2] = 0
    du[2,1] = sqrt(beta*u[1]*u[2]/N)
    du[2,2] = -sqrt(gamma*u[2])
  end

prob = SDEProblem(f,g,[N-15;15],(0.0,200.0),noise_rate_prototype=zeros(2,2))
sol = solve(prob,LambaEM())
   
plot(sol)
   