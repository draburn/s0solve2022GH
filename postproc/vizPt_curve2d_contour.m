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
  funchZ = @(omega_dummy)( omega_dummy ), ...
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
	contourLevels = linspace( zMin, zMax, numCL );
	assert(isrealarray(contourLevels,[1,numCL]));
	%
	numCol = numCL; % Must be true for contourf+colormap to work well.
	cMap = jet(numCol);
	%zEps = 0.01*(zMax-zMin)/numCol;
	%
	for n=1:numCol
		zCLo = contourLevels(n);
		if ( n<numCol )
			zCHi = contourLevels(n+1);
		else
			zCHi = (2*contourLevels(n)) - contourLevels(n-1);
		end
		%
		if ( zCLo < (zMax-zMin)/numCol )
			% Contains zero.
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.35;
		elseif ( zCLo > z0 )
			% Over z0.
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.4;
		elseif ( zCHi > z0 )
			% Contains z0. (Could also contain zLo, but, meh.)
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.75;
		elseif ( zCLo > zLo )
			% Between zLo and z0 without containing either.
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.6;
		elseif ( zCHi > zLo )
			% Contains zLo. (If also contains z0, see above.)
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.45;
		else
			% Under zLo.
			cMap(n,:) *= 0.2;
			cMap(n,:) += 0.8;
		end
	end
	%
	contourf(gridX,gridY,gridZ,contourLevels);
	colormap(cMap);
	%echo__contourLevels = contourLevels
	%echo__cMap = cMap
	%
	grid on;
end
