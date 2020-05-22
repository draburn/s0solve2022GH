%function [ zMin, zMax, zAvg, zVar ] = contourfunch( ...
%  funchZ = @(x,y)( x.^2 + y.^2), ...
%  ax = [ -5.0, 5.0, -5.0, 5.0 ], ...
%  numXVals = 51, ...
%  numYVals = 51, ...
%  numContours = 20, ...
%  funchSupportsMultiArg = true, ...
%  useContourF = true )
%
% 2020.05.21: This has been largely deprecated by gridfunch().

function [ zMin, zMax, zAvg, zVar ] = contourfunch( ...
  funchZ = @(x,y)( x.^2 + y.^2), ...
  ax = [ -5.0, 5.0, -5.0, 5.0 ], ...
  numXVals = 51, ...
  numYVals = 51, ...
  numContours = 20, ...
  funchSupportsMultiArg = true, ...
  useContourF = true )
	%
	xVals = linspace(ax(1),ax(2),numXVals);
	yVals = linspace(ax(3),ax(4),numYVals);
	[ gridX, gridY ] = ndgrid( xVals, yVals );
	%
	%gridZ = funchZ(gridX,gridY); might work, but, not worth adding case.
	rvecX = reshape(gridX,[1,numXVals*numYVals]);
	rvecY = reshape(gridY,[1,numXVals*numYVals]);
	if (funchSupportsMultiArg)
		rvecZ = funchZ(rvecX,rvecY);
	else
		rvecZ = zeros(1,numXVals*numYVals);
		for n=1:numXVals*numYVals
			rvecZ(n) = funchZ(rvecX(n),rvecY(n));
		end
	end
	gridZ = reshape(rvecZ,[numXVals,numYVals]);
	%
	zMin = min(rvecZ);
	zMax = max(rvecZ);
	zAvg = sum(rvecZ)/(numXVals*numYVals);
	zSqAvg = sum(rvecZ.^2)/(numXVals*numYVals);
	zVar = sqrt(zSqAvg-(zAvg^2));
	%
	if (useContourF)
		contourf(gridX,gridY,gridZ,numContours);
	else
		contour(gridX,gridY,gridZ,numContours);
	end
	colormap(mycmap);
	grid on;
return;
end
