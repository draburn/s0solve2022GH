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
	gridZ = funchZ( gridF );
	zMin = min(min(gridZ));
	zMax = max(max(gridZ));
	fLevels = linspace( min(min(gridZ)), max(max(gridZ)), numCuts+2 );
	fLevels(1) += sqrt(eps)*(zMax-zMin)/numCuts;
	fLevels(end) -= sqrt(eps)*(zMax-zMin)/numCuts;
	%echo__fLevels = fLevels
	%
	if (strcmp("contour",strContourFunc))
		for n=1:numCuts+2
			fTemp = fLevels(n);
			lFloat = (fTemp-fMin_base)/(fMax_base-fMin_base);
			if (lFloat<0)
				colorMap_tailored(n,:) = color_underMin;
			elseif (lFloat>1)
				colorMap_tailored(n,:) = color_overMax;
			else
				m = 1+round((numColorsBase-1.0)*lFloat+0.5);
				colorMap_tailored(n,:) = colorMap_base(m,:);
			end
		end
	elseif (strcmp("contourf",strContourFunc))
		for n=1:numCuts+1
			fTemp = (fLevels(n)+fLevels(n+1))/2.0;
			lFloat = (fTemp-fMin_base)/(fMax_base-fMin_base);
			if (lFloat<0)
				colorMap_tailored(n,:) = color_underMin;
			elseif (lFloat>1)
				colorMap_tailored(n,:) = color_overMax;
			else
				m = 1+round((numColorsBase-1.0)*lFloat);
				colorMap_tailored(n,:) = colorMap_base(m,:);
			end
		end
		% Last element doesn't matter, but, needs to be included.
		colorMap_tailored(numCuts+2,:) = colorMap_tailored(numCuts+1,:);
	end
	%echo__colorMap_tailored = colorMap_tailored
	%
	if (strcmp("contour",strContourFunc))
		contour( gridX, gridY, gridZ, fLevels );
	elseif (strcmp("contourf",strContourFunc))
		contourf( gridX, gridY, gridZ, fLevels );
	else
		error("Unsupported value of strContourFunc.");
	end
	colormap( colorMap_tailored );
	grid on;
end
