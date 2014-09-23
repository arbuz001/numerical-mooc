function g1 = euler_step(u,dt,g,a1,ve,ms,rho,A,Cd)
g1 = u + dt*f(u,g,a1,ve,ms,rho,A,Cd);
end