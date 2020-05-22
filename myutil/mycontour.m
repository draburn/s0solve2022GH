function mycontour( ...
  gridX, ...
  gridY, ...
  gridZ, ...
  zLo = min(min(gridZ)), ...
  zHi = max(max(gridZ)), ...
  cMapFull = mycmap, ...
  funchViz = @(z)(z), ...
  numCL = 20, ...
  zFlagVals = [], ...
  prm = [] )
	%
	vLo = funchViz(zLo);
	vHi = funchViz(zHi);
	gridV = cap(funchViz(gridZ),vLo,vHi);
	v0 = min(min(gridV));
	v1 = max(max(gridV));
	%
	numCFull = size(cMapFull,1);
	i0 = 1 + floor( (numCFull-1.0)*(v0-vLo)/(vHi-vLo) )
	i1 = 1 + ceil( (numCFull-1.0)*(v1-vLo)/(vHi-vLo) )
	cMap = cMapFull(i0:i1,:);
	%
	contourf( gridX, gridY, gridV, numCL );
	colormap(cMap);
	grid on;
return;
end
