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
		if ( zCLo < 0.5*(zMax-zMin)/numCol ...
		  || ((n==1) && zCLo < (zMax-zMin)/numCol) )
			% Contains zero. (Maybe also contains zLo, z0.)
			cMap(n,:) = 0.35 + 0.3*cMap0(1,:);
		elseif ( zCHi < zLo )
			% Between zLo and z0 without containing either.
			% Note that color scale runs from 0 to z0, not 0 to zLo.
			m = 1 + round( (numCol1-1)*(zCMid-0.0)/(z0-0.0) );
			cMap(n,:) = 0.9 + 0.1*cMap1(m,:);
		elseif ( zCLo < zLo ...
		      || ((n==1) && zCLo + (zMax-zMin)/numCol ) )
			% Contains zLo. (Maybe also contains z0.)
			% Take color similar to (between zLo and z0) range,
			% to provide sense of how low value is.
			m = 1 + round( (numCol1-1)*(zLo-0.0)/(z0-0.0) );
			cMap(n,:) = 0.6 + 0.15*cMap1(m,:);
		elseif ( zCHi < z0 )
			% Between zLo and z0 without containing either.
			m = 1 + round( (numCol1-1)*(zCMid-zLo)/(z0-zLo) );
			cMap(n,:) = 0.8 + 0.15*cMap1(m,:);
		elseif ( zCLo < z0 )
			% Contains z0.
			cMap(n,:) = 0.15 + 0.8*(2.0-cMap1(end,:)-cMap2(1,:))/2.0;
		else
			% Over z0.
			m = 1 + round( (numCol2-1)*(zCMid-z0)/(zMax-z0) );
			cMap(n,:) = 0.65 + 0.2*cMap2(m,:);
		end
	end
	%
	contourf(gridX,gridY,gridZ,contourLevels);
	colormap(cMap);
	%
	grid on;
end
