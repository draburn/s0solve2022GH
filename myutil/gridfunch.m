%function [ gridX, gridY, gridZ ] = gridfunch( ...
%  funchZ = @(x,y)( x.^2 + y.^2), ...
%  multiArgLevel = 1, ...
%  ax = [ -5.0, 5.0, -5.0, 5.0 ], ...
%  numXVals = 51, ...
%  numYVals = 51 )

function [ gridX, gridY, gridZ ] = gridfunch( ...
  funchZ = @(x,y)( x.^2 + y.^2), ...
  multiArgLevel = 1, ...
  ax = [ -5.0, 5.0, -5.0, 5.0 ], ...
  numXVals = 51, ...
  numYVals = 51 )
	%
	xVals = linspace(ax(1),ax(2),numXVals);
	yVals = linspace(ax(3),ax(4),numYVals);
	[ gridX, gridY ] = ndgrid( xVals, yVals );
	%
	switch (multiArgLevel)
	case {0}
		rvecX = reshape(gridX,[1,numXVals*numYVals]);
		rvecY = reshape(gridY,[1,numXVals*numYVals]);
		rvecZ = zeros(1,numXVals*numYVals);
		for n=1:numXVals*numYVals
			rvecZ(n) = funchZ(rvecX(n),rvecY(n));
		end
		gridZ = reshape(rvecZ,[numXVals,numYVals]);
	case {1}
		rvecX = reshape(gridX,[1,numXVals*numYVals]);
		rvecY = reshape(gridY,[1,numXVals*numYVals]);
		rvecZ = funchZ(rvecX,rvecY);
		gridZ = reshape(rvecZ,[numXVals,numYVals]);
	case {2}
		gridZ = funchZ(gridX,gridY);
	otherwise
		error(["Invalid multiArgLevel ( " num2str(multiArgLevel) " )."]);
	end
	%
return;
end
