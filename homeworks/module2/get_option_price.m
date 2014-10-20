function yy = get_option_price(v,x,dx,xmin)
	% recoursively loop over dimentions
	Ndim = ndims(v);
	dim = size(v);
	
	if((Ndim == 2) && (dim(1) == 1)) 
		yy = spline(xmin + dx*(1:length(v)),v,x);
	else
		k = 1;
		
		dim_k = dim(k);
		x_k = x(2:length(x));
		dx_k = dx(2:length(x));
		xmin_k= xmin(2:length(x));
		
		idx = floor((x(k) - xmin(k))/dx(k));
		if(idx == ceil((x(k) - xmin(k))/dx(k)))
			% make sure index is not less than 1
			idx = max(idx,1);
			
			% make sure index is not greater than max index in that dimention
			idx = min(idx,dim(k)-1);
			
			v_k = v(idx,:);
			
			yy  = get_option_price(v_k,x_k,dx_k,xmin_k);
		else
			% make sure index is not less than 1
			idx = max(idx,1);
			
			% make sure index is not greater than max index in that dimention
			idx = min(idx,dim(k)-1);
			
			v1 = v(idx,:);
			v2 = v(idx+1,:);
			
			v_k  = v1 + (v2-v1)*x1;
			
			yy  = get_option_price(v_k,x_k,dx_k,xmin_k);
		end
	end
end