function [ omega, vecNablaOmega, matNabla2Omega ] =  testfunc2021_funcOmega( vecX, testFuncPrm )
	assert( isrealarray(vecX,[testFuncPrm.sizeX,1]) );
	%
	vecY = vecX - testFuncPrm.vecXE;
	vecF = testFuncPrm.vecFE + testFuncPrm.matJ*vecY;
	for n=1:testFuncPrm.sizeF
		vecF(n) += 0.5*vecY'*(testFuncPrm.ary3K(:,:,n))*vecY;
	end
	omega = (vecF'*vecF)/2.0;
	if ( 1 == nargout )
		return;
	end
	%
	matJ = testFuncPrm.matJ;
	for n=1:testFuncPrm.sizeF
		matJ(n,:) += ( testFuncPrm.ary3K(:,:,n)*vecY )';
	end
	vecNablaOmega = matJ'*vecF;
	if ( 2 == nargout )
		return;
	end
	%
	matNabla2Omega = matJ'*matJ;
	for n=1:testFuncPrm.sizeF
		matNabla2Omega += vecF(n)*testFuncPrm.ary3K(:,:,n);
	end
return;
end


%!test
%!	numFigs0 = 0;
%!	numFigs = numFigs0;
%!	setAxisEqual = true;
%!	%
%!	sizeX = 2;
%!	testFuncPrm = testfunc2021_genPrm();
%!	funchOmega = @(dummyX) testfunc2021_funcOmega( dummyX, testFuncPrm );
%!	%
%!	numVals = 10 + round(5.0*abs(randn()));
%!	vecXVals = randn(sizeX,numVals);
%!	%
%!	%
%!	omegaVals = zeros(1,numVals);
%!	vecNablaOmegaVals = zeros(sizeX,numVals);
%!	matNabla2OmegaVals = zeros(sizeX,sizeX,numVals);
%!	for n=1:numVals
%!		[ omega, vecNablaOmega, matNabla2Omega ] = funchOmega( vecXVals(:,n) );
%!		assert( isrealscalar(omega) );
%!		assert( isrealarray(vecNablaOmega,[sizeX,1]) );
%!		assert( isrealarray(matNabla2Omega,[sizeX,sizeX]) );
%!		assert( issymmetric(matNabla2Omega) );
%!		omegaVals(n) = omega;
%!		vecNablaOmegaVals(:,n) = vecNablaOmega;
%!		matNabla2OmegaVals(:,:,n) = matNabla2Omega;
%!	end
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	epsNabla2Omega = sqrt(eps*sumsq(reshape(matNabla2OmegaVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check vecNablaOmega and matNabla2Omega.
%!	for n=1:numVals
%!		epsX = 1e-6;
%!		vecNablaOmega_fd = zeros(sizeX,1);
%!		matNabla2Omega_fd = zeros(sizeX,sizeX);
%!		for m=1:sizeX
%!			vecXP = vecXVals(:,n); vecXP(m) += epsX;
%!			vecXM = vecXVals(:,n); vecXM(m) -= epsX;
%!			[ omegaP, vecNablaOmegaP ] = funchOmega( vecXP );
%!			[ omegaM, vecNablaOmegaM ] = funchOmega( vecXM );
%!			vecNablaOmega_fd(m) = ( omegaP - omegaM ) / (2.0*epsX);
%!			matNabla2Omega_fd(m,:) = ( vecNablaOmegaP - vecNablaOmegaM ) / (2.0*epsX);
%!		end
%!		assert( reldiff(vecNablaOmega_fd,vecNablaOmegaVals(:,n),epsNablaOmega) < 100.0*epsX );
%!		assert( reldiff(matNabla2Omega_fd,matNabla2OmegaVals(:,:,n),epsNabla2Omega) < 100.0*epsX );
%!	end
%!	%
%!	isVectorized = false;
%!	ax = [ -5.0, 5.0, -5.0, 5.0 ];
%!	numXVals = [ 51, 55 ];
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = ...
%!	  genVizGrids( funchOmega, isVectorized, ax, numXVals );
%!	%
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = sqrt(sqrt(gridF)); strZ = "sqrt(sqrt(omega))";
%!	%gridZ = log( 1.0 + gridF.^2 ); strZ = "log( 1 + omega^2 )";
%!	contourf( gridX1, gridX2, gridZ );
%!	colormap( 0.3 + 0.7*colormap("default") );
%!	if (setAxisEqual)
%!		axis equal;
%!		axis equal; % Needed twice b/c of bug in Octave?
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( __FILE__, __LINE__, sprintf( "Please check figure(s) %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
%!	%
%!	return;
