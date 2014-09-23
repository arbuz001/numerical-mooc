function g2 = f(u,g,a1,ve,ms,rho,A,Cd)
h = u(1);
v = u(2);
mp = u(3);

	if(mp >= 0)
		g2 = [v, 1/(ms+mp)*(-g*(ms+mp) + a1*ve - 1/2*rho*v*abs(v)*A*Cd), -a1];
	else
		g2 = [v, 1/(ms+mp)*(-g*(ms+mp) + 0.0*ve - 1/2*rho*v*abs(v)*A*Cd), 0.0];
	end
end