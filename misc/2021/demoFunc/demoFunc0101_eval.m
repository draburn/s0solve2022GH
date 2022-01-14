function vecF = demoFunc0101_eval( vecX, funcPrm )
	if (1==size(vecX,2))
		vecY = vecX - funcPrm.x0;
	else
		vecY = vecX - repmat(funcPrm.x0,1,size(vecX,2));
		% "vecX" and "vecY" aren't actually vectors at this point.
	end
	vecF = (funcPrm.m1 * vecY) ...
	  + ( (funcPrm.m2a*vecY) .* (funcPrm.m2b*vecY) ) ...
	  + ( (funcPrm.m3a*vecY) .* (funcPrm.m3b*vecY) .* (funcPrm.m3c*vecY) );
return;
end

%!test
%!	sizeX = 100;
%!	vecX0 = zeros(sizeX,1);
%!	seedPrm = demoFunc0101_genSeedPrm("easy");
%!	seedPrm.sizeX = sizeX;
%!	seedPrm.sizeF = sizeX;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%!	vecF0 = demoFunc0101_eval( vecX0, funcPrm );
%!	vecX4 = funcPrm.x0+(1e-4);
%!	vecX8 = funcPrm.x0+(1e-8);
%!	vecF4 = demoFunc0101_eval( vecX4, funcPrm );
%!	vecF8 = demoFunc0101_eval( vecX8, funcPrm );
%!	vecFSolu = demoFunc0101_eval( funcPrm.x0, funcPrm );
%!	res0 = sqrt(sum(vecF0.^2))
%!	res4 = sqrt(sum(vecF4.^2))
%!	res8 = sqrt(sum(vecF8.^2))
%!	resSolu = sqrt(sum(vecFSolu.^2))

%!test
%!	sizeX = 2;
%!	vecX0 = zeros(sizeX,1);
%!	seedPrm = demoFunc0101_genSeedPrm("xhard");
%!	seedPrm.sizeX = sizeX;
%!	seedPrm.sizeF = sizeX;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%!	x0 = 0.0;
%!	y0 = 0.0;
%!	x1 = 1.0;
%!	y1 = 1.0;
%!	funcPrm.x0(1) = x1;
%!	funcPrm.x0(2) = y1;
%!	%
%!	numXVals = 200;
%!	numYVals = 201;
%!	xSpanCoeff = 500.0;
%!	ySpanCoeff = 500.0;
%!	%
%!	xMid = 0.5 * (x0+x1);
%!	yMid = 0.5 * (y0+y1);
%!	xDel = 0.5 * abs(x1 - x0) * xSpanCoeff;
%!	yDel = 0.5 * abs(y1 - y0) * ySpanCoeff;
%!	xLo = xMid - xDel;
%!	xHi = xMid + xDel;
%!	yLo = yMid - yDel;
%!	yHi = yMid + yDel;
%!	xVals = xLo + ( (xHi-xLo)*(0:numXVals-1)/(numXVals-1.0) );
%!	yVals = yLo + ( (yHi-yLo)*(0:numYVals-1)/(numYVals-1.0) );
%!	[ meshX, meshY ] = meshgrid( xVals, yVals );
%!	xGridVals = reshape( meshX, 1, [] );
%!	yGridVals = reshape( meshY, 1, [] );
%!	vecRVals = [ xGridVals; yGridVals ];
%!	tic();
%!	%vecFVals = [ yGridVals; yGridVals ];
%!	vecFVals = demoFunc0101_eval( vecRVals, funcPrm );
%!	toc();
%!	omegaVals = 0.5*sqrt( (vecFVals(1,:).^2) + (vecFVals(2,:).^2) );
%!	meshOmega = reshape( omegaVals, numYVals, numXVals );
%!	contour( meshX, meshY, sqrt(2.0*meshOmega), 20 );
%!	axis equal;
%!	xlabel( "x" );
%!	ylabel( "y" );
%!	grid on;
