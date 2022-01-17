% Function...
%   vecNHat is the surface normal.
%   vecUHat is the direction along which vecS(vecX) does not change for vecX near the surface;
%   (For a circle, vecNHat and vecUHat are the same.)

function [ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = ...
  demoFunc_surfEllip( vecXVals, vecXCent, matA=[], bigR=1.0, debugMode=false )
	if ( 2 > nargin || 5 < nargin )
		msg( __FILE__, __LINE__, "Bad nargin." );
		print_usage();
		return; % Superfluous?
	elseif ( 4 < nargout )
		msg( __FILE__, __LINE__, "Bad nargout." );
		print_usage();
		return; % Superfluous?
	end
	%
	sizeX = size(vecXVals,1);
	numVals = size(vecXVals,2);
	%
	assert( isscalar(debugMode) );
	if (debugMode)
		assert( isrealarray(vecXVals,[sizeX,numVals]) );
		assert( isrealarray(vecXCent,[sizeX,1]) );
		if (~isempty(matA))
			assert( isrealarray(matA,[sizeX,sizeX]) );
		end
		assert( isrealscalar(bigR) );
		assert( 0.0 < bigR );
	end
	%
	%
	vecDVals = vecXVals - vecXCent;
	dVals = sqrt(sumsq(vecDVals,1));
	numValsAtCent = sum(0.0==dVals);
	if (0<numValsAtCent)
		% As necessary, modify vecDVals to not be at center and update dVals.
		vecDVals(1,(0.0==dVals)) = 1.0;
		dVals(0.0==dVals) = 1.0;
	end
	%
	if (isempty(matA))
		vecUHatVals = vecDVals./dVals;
		vecSVals = vecXCent + (bigR*vecUHatVals);
		if ( 3 <= nargout )
			vecNHatVals = vecUHatVals;
			if ( 4 <= nargout )
			if ( 1 == numVals )
				matNablaSTVals = (bigR/dVals)*( eye(sizeX,sizeX) - (vecNHatVals*(vecNHatVals')) );
			else
				% There may be a non-loop way to evaluate this.
				% But, optimization here is not important.
				for n=1:numVals
					matNablaSTVals(:,:,n) = (bigR/dVals(n))*( ...
					  eye(sizeX,sizeX) - (vecNHatVals(:,n)*(vecNHatVals(:,n)')) );
				end
			end
			end
		end
	else
		% matA is not empty.
		vecBVals = matA * vecDVals;
		bVals = sqrt(sumsq(vecBVals,1));
		if ( 0 < sum(0.0==bVals) )
			error( "matA is singular." );
		end
		vecSVals = vecXCent + (bigR*(vecDVals./bVals));
		if ( 3 <= nargout )
			vecUHatVals = vecDVals./dVals;
			vecPVals = matA' * vecBVals;
			pVals = sqrt(sumsq(vecPVals,1));
			if ( 0 < sum(0.0==pVals) )
				error( "matA is singular." );
			end
			vecNHatVals = vecPVals./pVals;
			if ( 4 <= nargout )
			if ( 1 == numVals )
				matNablaSTVals = (bigR/bVals)*eye(sizeX,sizeX) ...
				  -  ((bigR/(bVals^3))*vecPVals) * (vecDVals');
			else
				% There may be a non-loop way to evaluate this.
				% But, optimization here is not important.
				for n=1:numVals
					matNablaSTVals(:,:,n) = (bigR/bVals(n))*eye(sizeX,sizeX) ...
					  -  ((bigR/(bVals(n)^3))*vecPVals(:,n)) * (vecDVals(:,n)');
				end
			end
			end
		end
	end
return;
end


%!test
%!	msg( __FILE__, __LINE__, "Performing basic execution test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	%
%!	% Test with minimal input.
%!	for n=1:numVals
%!		[ vecS, vecNHat, vecUHat, matNablaST ] = demoFunc_surfEllip( vecXVals(:,n), vecXCent );
%!		assert( isrealarray(vecS,[sizeX,1]) );
%!		assert( isrealarray(vecNHat,[sizeX,1]) );
%!		assert( isrealarray(vecUHat,[sizeX,1]) );
%!		assert( isrealarray(matNablaST,[sizeX,sizeX]) );
%!	end
%!	%
%!	% Test with maximal input.
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	bigR = 0.01 + abs(randn);
%!	debugMode = true;
%!	[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = demoFunc_surfEllip( vecXVals, vecXCent, matA, bigR, debugMode );
%!	assert( isrealarray(vecSVals,[sizeX,numVals]) );
%!	assert( isrealarray(vecNHatVals,[sizeX,numVals]) );
%!	assert( isrealarray(vecUHatVals,[sizeX,numVals]) );
%!	assert( isrealarray(matNablaSTVals,[sizeX,sizeX,numVals]) );
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Performing unscaled finite-differencing test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	matA = [];
%!	bigR = 0.01 + abs(randn());
%!	funchSurf = @(dummyX) demoFunc_surfEllip( dummyX, vecXCent, matA, bigR );
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
%!	epsS = sqrt(eps*sumsq(reshape(vecSVals,[],1)))/numVals;
%!	epsNHat = sqrt(eps);
%!	epsUHat = sqrt(eps);
%!	epsNablaST = sqrt(eps*sumsq(reshape(matNablaSTVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check vectorized evaluation.
%!	for n=1:numVals
%!		[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecXVals(:,n) );
%!		assert( reldiff(vecS,vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHat,vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHat,vecUHatVals(:,n),epsUHat) < sqrt(eps) );
%!		assert( reldiff(matNablaST,matNablaSTVals(:,:,n),epsNablaST) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	% Check derivative...
%!	matNablaSTVals_fd = zeros(sizeX,sizeX,numVals);
%!	epsX = 1e-6;
%!	for n=1:sizeX
%!		vecXPVals = vecXVals;
%!		vecXPVals(n,:) += epsX;
%!		vecSPVals = funchSurf( vecXPVals );
%!		vecXMVals = vecXVals;
%!		vecXMVals(n,:) -= epsX;
%!		vecSMVals = funchSurf( vecXMVals );
%!		matNablaSTVals_fd(n,:,:) = ( vecSPVals - vecSMVals ) / (2.0*epsX);
%!	end
%!	for n=1:numVals
%!		assert( reldiff(matNablaSTVals_fd(:,:,n),matNablaSTVals(:,:,n),epsNablaST) < 100.0*epsX );
%!	end
%!	msg( __FILE__, __LINE__, "BYE!" ); return;
%!	%
%!	%
%!	% Check A = I does nothing...
%!	[ vecSVals_eye, vecNHatVals_eye, vecUHatVals_eye, matNablaSTVals_eye ] = funchSurf( vecXVals, vecXCent, matA, bigR );
%!	for n=1:numVals
%!		assert( reldiff(vecSVals_eye(:,n),vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHatVals_eye(:,n),vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHatVals_eye(:,n),vecUHatVals(:,n),epsUHat) < sqrt(eps) );
%!		assert( reldiff(matNablaSTVals_eye(:,:,n),matNablaSTVals(:,:,n),epsNablaST) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Performing scaled finite-differencing test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	sizeF = round(1.0 + abs(sizeX*randn));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + matA0'*matA0;
%!	bigR = 0.01 + abs(randn());
%!	funchSurf = @(dummyX) demoFunc_surfEllip( dummyX, vecXCent, matA, bigR );
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
%!	epsS = sqrt(eps*sumsq(reshape(vecSVals,[],1)))/numVals;
%!	epsNHat = sqrt(eps);
%!	epsUHat = sqrt(eps);
%!	epsNablaST = sqrt(eps*sumsq(reshape(matNablaSTVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check vectorized evaluation.
%!	for n=1:numVals
%!		[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecXVals(:,n) );
%!		assert( reldiff(vecS,vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHat,vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHat,vecUHatVals(:,n),epsUHat) < sqrt(eps) );
%!		assert( reldiff(matNablaST,matNablaSTVals(:,:,n),epsNablaST) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	% Check derivative...
%!	matNablaSTVals_fd = zeros(sizeX,sizeX,numVals);
%!	epsX = 1e-6;
%!	for n=1:sizeX
%!		vecXPVals = vecXVals;
%!		vecXPVals(n,:) += epsX;
%!		vecSPVals = funchSurf( vecXPVals );
%!		vecXMVals = vecXVals;
%!		vecXMVals(n,:) -= epsX;
%!		vecSMVals = funchSurf( vecXMVals );
%!		matNablaSTVals_fd(n,:,:) = ( vecSPVals - vecSMVals ) / (2.0*epsX);
%!	end
%!	for n=1:numVals
%!		assert( reldiff(matNablaSTVals_fd(:,:,n),matNablaSTVals(:,:,n),epsNablaST) < 100.0*epsX );
%!	end
%!	msg( __FILE__, __LINE__, "BYE!" ); return;
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );
