% Note: contourf() accesses colorMap() based on scaled values of the gridZ;
%  simple attempts to use arbitrary contour levels with desired color values
%  doesn't work.

function vizPt_curve2d_contour( ...
  gridX, ...
  gridY, ...
  gridOmega, ...
  numCL, ...
  omegaVizMin, ...
  omegaLo, ...
  omega0, ...
  omegaVizMax, ...
  funchZ = @(omega_dummy)( sqrt(omega_dummy) ), ...
  strContourFunc = "contourf", ...
  prm = []  )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	commondefs;
	thisFile = "vizPt_curve2d_contour";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	gridZ = funchZ(cap(gridOmega,omegaVizMin,omegaVizMax));
	zMin = min(min(gridZ));
	zLo = funchZ(omegaLo);
	z0 = funchZ(omega0);
	zMax = max(max(gridZ));
	%
	assert( numCL >= 3 );
	numCLBelow = round( (numCL-3.0)*(z0-zMin)/(zMax-zMin) );
	numCLAbove = numCL - 3 - numCLBelow;
	clBelow = linspace( zMin, z0, numCLBelow+2 );
	clAbove = linspace( z0, zMax, numCLAbove+2 );
	contourLevels = [ clBelow(1:end-1), z0, clAbove(2:end) ];
	contourLevels(1) -= (zMax-zMin)*0.001/numCL;
	contourLevels(end) += (zMax-zMin)*0.001/numCL;
	assert(isrealarray(contourLevels,[1,numCL]));
	%echo__contourLevels = contourLevels
	%
	if ( strcmp(strContourFunc,"contour") )
		cMap = 0.3 * jet(numCL-1);
		cMap(1:numCLBelow+1,:) += 0.3;
		cMap(numCLBelow+2,:) = [1.0,0.0,0.0];
		cMap(numCLBelow+3:end,:) += 0.7;
		%echo__cMap = cMap
		%
		contour(gridX,gridY,gridZ,contourLevels);
		colormap(cMap);
	elseif ( strcmp(strContourFunc,"contourf") )
		cMap = 0.3 * jet(numCL-1);
		cMap(1:numCLBelow+1,:) += 0.3;
		cMap(numCLBelow+2:end,:) += 0.7;
		%echo__cMap = cMap
		%
		contourf(gridX,gridY,gridZ,contourLevels);
		colormap(cMap);
	else
		error("Unsupported value of strContourFunc.");
	end
	grid on;
end
