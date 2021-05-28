function [ xOfCand, meritOfCand,  datOut ] = analyzeExt( xVals, fVals, prm=[], datIn=[] );
	%
	analyzeExt__init; thisFile = "analyzeExt";
	msg( thisFile, __LINE__, "THIS CODE IS UNWORTHY." );
	%
	x0 = xVals(indexOfGMin);
	x1 = xVals(indexOfGMin+1) - xVals(indexOfGMin-1);
	yVals = ( xVals - x0 ) / x1;
	diffYVals = diff(yVals);
	%
	meritOfCand = 0.5; % Worry about this later.
	datOut = [];
	%
	zVals = [];
	hVals = [];
	stencilSize=3;
	for n=1:numPts+1-stencilSize
		vecG = gVals(n:n+stencilSize-1)';
		vecY = yVals(n:n+stencilSize-1)';
		matY = [ ones(size(vecY)), vecY, vecY.^2 ];
		vecC = matY \ vecG;
		zValsForH(n) = (sum(vecY)-(vecY(1)+vecY(end))/2.0)/(stencilSize-1.0);
		%if (0==mod(stencilSize,2))
		%	zValsForH(n) = ( ...
		%	   vecY(floor((stencilSize+1)/2.0)) ...
		%     + vecY(ceil( (stencilSize+1)/2.0)) )/2.0;
		%else
		%	zValsForH(n) = vecY(round((stencilSize+1)/2.0));
		%end
		%zValsForH(n) = (vecY(1)+vecY(end))/2.0;
		hVals(n) = zValsForH(n) + 0.5*vecC(2)/vecC(3);
	end
	hVals *= x1;
	xValsForH = x0 + x1*zValsForH;
	%
	datOut.xValsForH = xValsForH;
	datOut.hVals = hVals;
	datOut.xValsLoForH = xVals(1:end+1-stencilSize);
	datOut.xValsHiForH = xVals(stencilSize:end);
	%
	if ( 3 == numPts )
		vecG = gVals(1:3)';
		vecY = yVals(1:3)';
		matY = [ ones(size(vecY)), vecY, vecY.^2 ];
		vecC = matY \ vecG;
		yOfCand = -0.5*vecC(2)/vecC(3);
		xOfCand = x0 + x1*yOfCand;
		return;
	end
	%
	% PLACEHOLDER: FORCE QUAD.
	vecG = gVals(indexOfGMin-1:indexOfGMin+1)';
	vecY = yVals(indexOfGMin-1:indexOfGMin+1)';
	matY = [ ones(size(vecY)), vecY, vecY.^2 ];
	vecC = matY \ vecG;
	yOfCand = -0.5*vecC(2)/vecC(3);
	xOfCand = x0 + x1*yOfCand;
	%
return;
end
