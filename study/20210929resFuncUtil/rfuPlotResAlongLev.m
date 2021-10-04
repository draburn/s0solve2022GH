function rfuPlotResAlongLev( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	%commondefs;
	%thisFile = "rfuPlotResAlongLev";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%assert( isrealscalar(verbLev) );
	% DRaburn 2021.09.29:
	%  We don't need this "common" stuff in the current code.
	%
	numFigs0 = mygetfield( prm, "numFigs0", 0 );
	assert( isposintscalar(numFigs0+1) );
	numFigs = numFigs0;
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	sizeX = size(vecX0,1);
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = mygetfield( prm, "vecF0", [] );
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	prm_calcNablaF = mygetfield( prm, "prm_calcNablaF", prm );
	prm_calcNablaSqF = mygetfield( prm, "prm_calcNablaSqF", prm );
	numPts = mygetfield( prm, "numPts", 20 );
	assert( isposintscalar(numPts) );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALCULATE STUFF
	%
	nablaF = rfuCalcNablaF( funchF, vecX0, prm_calcNablaF );
	vecDeltaNewt = -nablaF \ vecF0
	deltaNewtNorm = sqrt(sum(vecDeltaNewt.^2));
	deltaNewtNorm = 1.0;
	%
	vecG = nablaF' * vecF0;
	matH = nablaF' * nablaF;
	matI = eye(sizeX,sizeX);
	vecMu = (eps^-0.50).^linspace(-1.0,1.0,numPts)';
	%
	matDelta = NaN + zeros(sizeX,numPts);
	vecLambda = NaN + zeros(numPts);
	matX = NaN + zeros(sizeX,numPts);
	matF = NaN + zeros(sizeF,numPts);
	matFModel1 = NaN + zeros(sizeF,numPts);
	for n=1:numPts
		mu_temp = vecMu(n);
		vecDelta_temp = -( matH + mu_temp*matI ) \ vecG;
		lambda_temp = sqrt(sum(vecDelta_temp.^2)) / deltaNewtNorm;
		vecX_temp = vecX0 + vecDelta_temp;
		vecF_temp = funchF( vecX_temp );
		assert( isrealarray(vecF_temp,[sizeF,1]) );
		vecFModel1_temp = vecF0 + nablaF * vecDelta_temp;
		%
		matDelta(:,n) = vecDelta_temp;
		vecLambda(n) = lambda_temp;
		matX(:,n) = vecX_temp;
		matF(:,n) = vecF_temp;
		matFModel1(:,n) = vecFModel1_temp;
		%
		clear vecFModel1_temp;
		clear vecF_temp;
		clear vecX_temp;
		clear vecDelta_temp;
	end
	clear n;
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAKE PLOT(S).
	%
	matF = abs(matF);
	matFModel1 = abs(matFModel1);
	%
	numFigs++; figure(numFigs);
	n = 1;
	plot( ...
	  vecLambda, matFModel1(1,:), 'x-', 'markersize', 10, ...
	  vecLambda, matF(n,:), '+-', 'markersize', 10 );
	grid on;
	hold on;
	for n=2:sizeF
		plot( ...
		  vecLambda, matFModel1(n,:), 'x-', 'markersize', 10, ...
		  vecLambda, matF(n,:), '+-', 'markersize', 10 );
	end
	plot( vecLambda, 0*vecLambda, 'k-' );
	hold off
	clear n;
	%
	%
	%
return;
end


%!test
%!	commondefs;
%!	thisFile = "rfuPlotResAlongLev test";
%!	setprngstates(0);
%!	sizeX = 2
%!	sizeF = 3
%!	matA = randn(sizeF,sizeX);
%!	vecX0 = randn(sizeX,1);
%!	%funchF = @(x)( matA*(x-vecX0) );
%!	funchF = @(x)( (matA*(x-vecX0).^2).^2 );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	prm = [];
%!	tic();
%!	rfuPlotResAlongLev( funchF, vecX0, prm )
%!	toc();
