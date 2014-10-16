function get_option_price = f(v,S,ds,Smin,t,dt)

	[dim1,dim2] = size(v);

	if((dim1 > 1) & (dim2 > 1))
		i = floor((S - Smin)/ds)
		i_greater = ceil((S - Smin)/ds)
		
		j_less = floor( t / dt)
		j_greater = ceil( t / dt)
		
		get_option_price = v (j_greater,i_greater);
	else
		i_less = floor((S - Smin)/ds)
		i_greater = ceil((S - Smin)/ds)
		
		% get_option_price = v (i_less)
		get_option_price = v (i_less) + (v(i_greater) - v(i_less))*(S - i_less*ds)/ds;
	end

end