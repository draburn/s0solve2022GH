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
	if (0)
		echo__omeaVizMin = omegaVizMin
		echo__omegaLo = omegaLo
		echo__omega0 = omega0
		echo__omegaVizMax = omegaVizMax
		echo__omegaGridMin = min(min(gridOmega))
		echo__omegaGridMax = max(max(gridOmega))
	end
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
	numCol0 = 1000;
	numCol1 = 1000;
	numCol2 = 1000;
	cMap0 = jet(numCol0);
	cMap1 = jet(numCol1);
	cMap2 = jet(numCol1);
	cMap = jet(numCol);
	%zEps = 0.01*(zMax-zMin)/numCol;
	%
	for n=1:numCol
		zCLo = contourLevels(n);
		if ( n<numCol )
			zCMid = (contourLevels(n)+contourLevels(n+1))/2.0;
			zCHi = contourLevels(n+1);
		else
			zCMid = contourLevels(n);
			zCHi = (2*contourLevels(n)) - contourLevels(n-1);
		end
		%
		if ( zCLo < 0.5*(zMax-zMin)/numCol )
			% Contains zero.
			cMap(n,:) = 0.35 + 0.3*cMap0(1,:);
		elseif ( zCLo > z0 )
			% Over z0.
			m = 1 + round( (numCol2-1)*(zCMid-z0)/(zMax-z0) );
			cMap(n,:) = 0.65 + 0.2*cMap2(m,:);
		elseif ( zCHi > z0 )
			% Contains z0. (Could also contain zLo, but, meh.)
			cMap(n,:) = 0.15 + 0.8*(2.0-cMap1(end,:)-cMap2(1,:))/2.0;
		elseif ( zCLo > zLo )
			% Between zLo and z0 without containing either.
			m = 1 + round( (numCol1-1)*(zCMid-zLo)/(z0-zLo) );
			cMap(n,:) = 0.8 + 0.15*cMap1(m,:);
		elseif ( zCHi > zLo )
			% Contains zLo. (If also contains z0, see above.)
			cMap(n,:) = 0.25 + 0.5*(1.0-cMap1(1,:));
		else
			% Under zLo.
			m = 1 + round( (numCol1-1)*(zCMid-0.0)/(zLo-0.0) );
			cMap(n,:) = 0.9 + 0.1*cMap0(m,:);
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
