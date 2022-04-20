	vecX = datOut_800.dat_pp20220418.vecX;
	vecF = datOut_800.dat_pp20220418.vecF;
	matV = datOut_800.dat_pp20220418.matV;
	matW = datOut_800.dat_pp20220418.matW;
	%
	%
	%
	sizeV = size(matV,2)
	assert( isrealarray(matV,[sizeX,sizeV]) );
	assert( isrealarray(matW,[sizeX,sizeV]) );
	%
	sizeV = size(matW,2);
	matIV = ones(sizeV,sizeV);
	matWTW = matW'*matW;
	rcond(matWTW)
	%return;
	wSqScale = max(diag(matWTW));
	%
	numPVals = 10001;
	pVals = linspace( 0.0, 0.99, numPVals );
	pVals = ( 1.0 - (1.0-(pVals.^4)).^4);
	vecMG = -(matW'*vecF);
	vecDeltaVals = zeros(sizeX,numPVals);
	vecFModelVals = zeros(sizeF,numPVals);
	vecFVals = zeros(sizeF,numPVals);
	for n=1:numPVals
		p = pVals(n);
		matH_temp = p*matWTW + (1.0-p)*wSqScale*matIV;
		vecY_temp = matH_temp \ (p*vecMG);
		vecDelta_temp = matV * vecY_temp;
		vecDeltaVals(:,n) = vecDelta_temp;
		vecFModelVals(:,n) = vecF + matW*vecY_temp;
		vecFVals(:,n) = funchF( vecX + matV*vecY_temp );
	endfor
	%return
	% ALSO CHECK: RESET PRECONDITIONER.
	%
	%
	%
	figure(100);
	plot( ...
	  pVals, sqrt(sumsq(vecDeltaVals,1)), 's-' );
	%semilogy( ...
	%  sqrt(sumsq(vecDeltaVals,1)), sqrt(sumsq(vecFModelVals,1)), 'o-', ...
	%  sqrt(sumsq(vecDeltaVals,1)), sqrt(sumsq(vecFVals,1)), 'x-' );
	grid on;
	%
	figure(101);
	loglog( ...
	  sqrt(sumsq(vecDeltaVals,1)), sqrt(sumsq(vecFModelVals,1)), 'o-', ...
	  sqrt(sumsq(vecDeltaVals,1)), sqrt(sumsq(vecFVals,1)), 'x-' );
	grid on;
%
