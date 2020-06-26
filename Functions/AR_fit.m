function mld_list = AR_fit(cx_select)
	%{
	Function to determine the AR model for surrogate data
	INPUT:
		cx_select: empirical time series ; T * N
	OUTPUT:
		mld_list: AR model 
	%}

	thres = -2; % threshold for rejecting higher BIC
	q_pool = 1:50; % the set of AR orders to test
	q_list = {}; mld_list = {};
	for i = 1:size(cx_select,2)
		x = cx_select(:,i);
		BIC = []; 
		for k = q_pool
		    mdl = ar(x,k);
		    BIC(k) = mdl.Report.Fit.BIC;
		end
		idx = find(diff(BIC)<thres == 0);
		q = q_pool(idx(1));
		q_list{i} = q;
		mld_list{i} = ar(x,q);
		i
	end
end