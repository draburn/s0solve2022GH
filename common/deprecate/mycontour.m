function mycontour( ...
  gridX, ...
  gridY, ...
  gridF, ...
  numCuts, ...
  colorMap_base, ...
  fMin_base, ...
  fMax_base, ...
  color_underMin = [0.0,0.0,0.0], ...
  color_overMax = [1.0,1.0,1.0], ...
  fVals_flagged = [], ...
  colorMap_flagged = [], ...
  funchZ = @(f)( f ), ...
  strContourFunc = "contourf"  )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	commondefs;
	thisFile = "mycontour";
	%
	%
	numColorsBase = size(colorMap_base,1);
	numFlaggedVals = max(size(fVals_flagged));
	%
	if ( strcmp("colorbar",strContourFunc) || strcmp("colorbarf",strContourFunc) )
		% Display original f values, don't apply funchZ.
		fLo = fMin_base - 0.3*(fMax_base-fMin_base);
		fHi = fMax_base + 0.3*(fMax_base-fMin_base);
		xVals = [0,1];
		yVals = linspace( fLo, fHi, numColorsBase+1 );
		[ gridX, gridY ] = ndgrid( xVals, yVals );
		gridZ = gridY;
		zMin_base = fMin_base;
		zMax_base = fMax_base;
		zVals_flagged = fVals_flagged;
	else
		gridZ = funchZ( gridF );
		zMin_base = funchZ( fMin_base );
		zMax_base = funchZ( fMax_base );
		zVals_flagged = funchZ( fVals_flagged );
	end
	zMin = min(min(gridZ));
	zMax = max(max(gridZ));
	zLevels = linspace( min(min(gridZ)), max(max(gridZ)), numCuts+2 );
	zLevels(1) += sqrt(eps)*(zMax-zMin)/numCuts;
	zLevels(end) -= sqrt(eps)*(zMax-zMin)/numCuts;
	%echo__zLevels = zLevels
	%
	if ( strcmp("contour",strContourFunc) || strcmp("colorbar",strContourFunc) )
		for n=1:numCuts+2
			zTemp = zLevels(n);
			lFloat = (zTemp-zMin_base)/(zMax_base-zMin_base);
			if (lFloat<0)
				colorMap_tailored(n,:) = color_underMin;
			elseif (lFloat>1)
				colorMap_tailored(n,:) = color_overMax;
			else
				m = 1+round((numColorsBase-1.0)*lFloat+0.5);
				colorMap_tailored(n,:) = colorMap_base(m,:);
			end
		end
		%
		% Modification for flagged values might be wrong.
		for n=1:numFlaggedVals
			zTemp = zVals_flagged(n);
			lFloat = (zTemp-zMin)/(zMax-zMin);
			m = round( (numCuts+2)*lFloat );
			if ( m >= 1 && m <= numCuts+2 )
				colorMap_tailored(m,:) = colorMap_flagged(n,:);
			end
		end
	elseif ( strcmp("contourf",strContourFunc) || strcmp("colorbarf",strContourFunc) )
		for n=1:numCuts+1
			zTemp = (zLevels(n)+zLevels(n+1))/2.0;
			lFloat = (zTemp-zMin_base)/(zMax_base-zMin_base);
			if (lFloat<0)
				colorMap_tailored(n,:) = color_underMin;
			elseif (lFloat>1)
				colorMap_tailored(n,:) = color_overMax;
			else
				m = 1+round((numColorsBase-1.0)*lFloat);
				colorMap_tailored(n,:) = colorMap_base(m,:);
			end
		end
		% Seems last element won't be shown, but needs to be included
		% to prevent color map getting re-scaled.
		colorMap_tailored(numCuts+2,:) = [1.0,1.0,0.0];
		%
		% Modification for flagged values might be wrong.
		for n=1:numFlaggedVals
			zTemp = zVals_flagged(n);
			lFloat = (zTemp-zMin)/(zMax-zMin);
			m = round( 0.5 + (numCuts+1.0)*lFloat );
			if ( m >= 1 && m <= numCuts+2 )
				colorMap_tailored(m,:) = colorMap_flagged(n,:);
			end
		end
	end
	%echo__colorMap_tailored = colorMap_tailored
	%
	if (strcmp("contour",strContourFunc))
		contour( gridX, gridY, gridZ, zLevels );
	elseif (strcmp("contourf",strContourFunc))
		contourf( gridX, gridY, gridZ, zLevels );
	elseif (strcmp("colorbar",strContourFunc))
		contour( gridX, gridY, gridZ, zLevels );
		set(get(gcf,"children"),"xticklabel",[]);
	elseif (strcmp("colorbarf",strContourFunc))
		contourf( gridX, gridY, gridZ, zLevels );
		set(get(gcf,"children"),"xticklabel",[]);
	else
		error("Unsupported value of strContourFunc.");
	end
	colormap( colorMap_tailored );
	grid on;
end
