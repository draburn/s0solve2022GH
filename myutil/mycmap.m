function cMap = mycmap(n=1000)
	cMap = 0.6 + (0.4*jet(n));
	cMap(1,:) *= 0.4;
	cMap(end,:) *= 0.4;
	cMap(end,:) += 0.6;
return;
end
