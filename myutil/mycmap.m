function cMap = mycmap(n=1000)
	%%%cMap = 0.7 + (0.3*jet(n));
	cMap = 0.7 + (0.3*colormap("default"));
	cMap(1,:) *= 0.5;
	cMap(end,:) *= 0.5;
	cMap(end,:) += 0.5;
return;
end
