% Function...
%  [ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = 
%    funcSurfEllip( vecXVals, vecXCent, bigR=1.0, matA=[], debugMode=false )
% Calculates "surface quantities" for an elliptical surface: ||A*(x-xCent)|| = R.
% Input...
%  vecXVals: collection of position vectors; size() is definitionally [ sizeX, numVals ].
%  vecXCent: position vector of center of min (or max) of omega; size() must be [ sizeX, 1 ].
%  bigR: scalar value of ellipse size; corresponds to radius if A=I; default is 1.0.
%  matA: scaling matrix; size() must be [ sizeX, sizeX ]; default is I.
%  debugMode: boolean flag to force debug checks; default is FALSE.
% Output...
%  vecSVals: vecXVals pulled to the surface; size() is [ sizeX, numVals ]. 
%  vecNHatVals: if requested, outward surface normal at each vecS; size() is [ sizeX, numVals ].
%  vecUHatVals: if requested, for each vecS, the local direction along which vecS does not change,
%   pointing outwards: this would be the outward surface normal if vecS were always the closest
%   point on the surface, but that problem appears to be non-trivial; size is [ sizeX, numVals ].
%  matNablaSTVals: if requested, the quantity nabla s';
%    size() is [ sizeX, sizeX ] if numVals is 1 and [ sizeX, sizeX, numVals ] if numVals >= 2.


function [ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = ...
  funcSurfEllip( vecXVals, vecXCent, bigR=1.0, matA=[], debugMode=false )
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
	if (debugMode)
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecXVals,[sizeX,numVals]) );
		assert( isrealarray(vecXCent,[sizeX,1]) );
		assert( isrealscalar(bigR) );
		assert( 0.0 < bigR );
		if (~isempty(matA))
			assert( isrealarray(matA,[sizeX,sizeX]) );
		end
	end
	%
	%
	vecDVals = vecXVals - vecXCent;
	dVals = sqrt(sumsq(vecDVals,1));
	numValsAtCent = sum(0.0==dVals);
	if (0<numValsAtCent)
		% As necessary, modify vecD to not be at center and update d.
		vecDVals(1,(0.0==dVals)) = 1.0;
		dVals(0.0==dVals) = 1.0;
	end
	%
	if (isempty(matA))
		vecUHatVals = vecDVals./dVals;
		vecSVals = vecXCent + (bigR*vecUHatVals);
		if ( 2 <= nargout )
			vecNHatVals = vecUHatVals;
			if ( 4 <= nargout )
			if ( 1 == numVals )
				matNablaSTVals = (bigR/dVals)*( eye(sizeX,sizeX) - (vecNHatVals*(vecNHatVals')) );
			else
				% There may be a non-loop way to evaluate this.
				% But, optimization here is not important.
				parfor n=1:numVals
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
		if ( 2 <= nargout )
			vecPVals = matA' * vecBVals;
			pVals = sqrt(sumsq(vecPVals,1));
			if ( 0 < sum(0.0==pVals) )
				error( "matA is singular." );
			end
			vecNHatVals = vecPVals./pVals;
			vecUHatVals = vecDVals./dVals;
			if ( 4 <= nargout )
			if ( 1 == numVals )
				matNablaSTVals = (bigR/bVals)*eye(sizeX,sizeX) - ((bigR/(bVals^3))*vecPVals) * (vecDVals');
			else
				% There may be a non-loop way to evaluate this.
				% But, optimization here is not important.
				parfor n=1:numVals
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
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 10 + round(10.0*abs(randn()));
%!	vecXCent = randn(sizeX,1);
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	%
%!	% Test with minimal input.
%!	for n=1:numVals
%!		[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurfEllip( vecXVals(:,n), vecXCent );
%!		assert( isrealarray(vecS,[sizeX,1]) );
%!		assert( isrealarray(vecNHat,[sizeX,1]) );
%!		assert( isrealarray(vecUHat,[sizeX,1]) );
%!		assert( isrealarray(matNablaST,[sizeX,sizeX]) );
%!	end
%!	%
%!	% Test with maximal input.
%!	bigR = 0.01 + abs(randn);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	debugMode = true;
%!	[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funcSurfEllip( vecXVals, vecXCent, bigR, matA, debugMode );
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
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 10 + round(10.0*abs(randn()));
%!	vecXCent = randn(sizeX,1);
%!	bigR = 0.01 + abs(randn());
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent, bigR );
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
%!	% Check surface self-consistency.
%!	[ vecSVals_s, vecNHatVals_s, vecUHatVals_s ] = funchSurf( vecSVals );
%!	for n=1:numVals
%!		assert( reldiff(vecSVals_s(:,n),vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHatVals_s(:,n),vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHatVals_s(:,n),vecUHatVals(:,n),epsUHat) < sqrt(eps) );
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
%!	%
%!	%
%!	% Check nHat and uHat wrt nablaST...
%!	for n=1:numVals
%!		assert( norm(matNablaSTVals(:,:,n)*vecNHatVals(:,n)) < sqrt(epsUHat*epsNablaST) );
%!		assert( norm(vecUHatVals(:,n)'*matNablaSTVals(:,:,n)) < sqrt(epsUHat*epsNablaST) );
%!	end
%!	%
%!	%
%!	% Check A = I does nothing...
%!	matA = eye(sizeX,sizeX);
%!	[ vecSVals_eye, vecNHatVals_eye, vecUHatVals_eye, matNablaSTVals_eye ] = funchSurf( vecXVals, vecXCent, bigR, matA );
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
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 10 + round(10.0*abs(randn()));
%!	vecXCent = randn(sizeX,1);
%!	bigR = 0.01 + abs(randn());
%!	sizeF = 1 + round(sizeX*abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + matA0'*matA0;
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent, bigR, matA );
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
%!	% Check surface self-consistency.
%!	[ vecSVals_s, vecNHatVals_s, vecUHatVals_s ] = funchSurf( vecSVals );
%!	for n=1:numVals
%!		assert( reldiff(vecSVals_s(:,n),vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHatVals_s(:,n),vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHatVals_s(:,n),vecUHatVals(:,n),epsUHat) < sqrt(eps) );
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
%!	%
%!	%
%!	% Check nHat and uHat wrt nablaST...
%!	for n=1:numVals
%!		assert( norm(matNablaSTVals(:,:,n)*vecNHatVals(:,n)) < sqrt(epsUHat*epsNablaST) );
%!		assert( norm(vecUHatVals(:,n)'*matNablaSTVals(:,:,n)) < sqrt(epsUHat*epsNablaST) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Testing conceptual understanding of nHat and vHat." );
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 10 + round(10.0*abs(randn()));
%!	vecXCent = randn(sizeX,1);
%!	bigR = 0.01 + abs(randn());
%!	sizeF = 1 + round(sizeX*abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + matA0'*matA0;
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent, bigR, matA );
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
%!	epsS = sqrt(eps*sumsq(reshape(vecSVals,[],1)))/numVals;
%!	epsNHat = sqrt(eps);
%!	epsUHat = sqrt(eps);
%!	epsNablaST = sqrt(eps*sumsq(reshape(matNablaSTVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check nHat...
%!	epsX = 1e-4;
%!	for n=1:numVals
%!		vecX = vecXVals(:,n);
%!		vecS = vecSVals(:,n);
%!		vecNHat = vecNHatVals(:,n);
%!		%
%!		vecY = randn(sizeX,1);
%!		vecXPY = vecX + epsX*vecY/(norm(vecY)+eps);
%!		[ vecSPY, vecNHatPY ] = funchSurf( vecXPY );
%!		assert( abs( (vecSPY-vecS)'*(vecNHatPY+vecNHatPY) ) < 10.0*epsX );
%!	end
%!	%
%!	%
%!	% Check uHat...
%!	epsX = 1e-4;
%!	vecXVals_p = vecSVals + epsX*vecUHatVals;
%!	[ vecSVals_p, vecNHatVals_p, vecUHatVals_p ] = funchSurf( vecXVals_p );
%!	for n=1:numVals
%!		assert( reldiff(vecSVals_p(:,n),vecSVals(:,n),epsS) < sqrt(eps) );
%!		assert( reldiff(vecNHatVals_p(:,n),vecNHatVals(:,n),epsNHat) < sqrt(eps) );
%!		assert( reldiff(vecUHatVals_p(:,n),vecUHatVals(:,n),epsUHat) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );
%!


%!test
%!	msg( __FILE__, __LINE__, "Generating visualization." );
%!	%setprngstates(75614704); % aixs() makes second graph dissappear.
%!	setprngstates();
%!	numFigs0 = 0;
%!	%
%!	numFigs = numFigs0;
%!	setAxisEqual = true;
%!	%
%!	sizeX = 2;
%!	vecXCent = randn(sizeX,1);
%!	bigR = 0.5 + abs(randn());
%!	sizeF = 1 + round(sizeX*abs(randn));;
%!	matA0 = randn(sizeF,sizeX);
%!	matA = 0.5*eye(sizeX,sizeX) + matA0'*matA0;
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent, bigR, matA );
%!	%
%!	% Trace surface with "pts".
%!	numPts = 1001;
%!	thetaPts = linspace(0.0,2.0*pi,numPts);
%!	vecXPts = vecXCent + bigR*[ cos(thetaPts); sin(thetaPts) ];
%!	vecSPts = funchSurf( vecXPts );
%!	distScale = max([ max(abs( vecSPts(1,:) - vecXCent(1) )), max(abs( vecSPts(2,:) - vecXCent(2) )) ]);
%!	%
%!	% Generate several "vals".
%!	numVals = 50 + round(10*abs(randn()));
%!	vecXVals = vecXCent + randn(sizeX,numVals)*bigR/norm(matA);
%!	vecSVals = funchSurf( vecXVals );
%!	%
%!	%
%!	numFigs++; figure(numFigs);
%!	n = 1;
%!	plot( ...
%!	  vecSPts(1,:), vecSPts(2,:), 'o-', 'markersize', 2, 'color', [0.5,0.5,0.5], ...
%!	  [ vecXVals(1,:); vecSVals(1,:) ], [ vecXVals(2,:); vecSVals(2,:) ], 'x-', 'markersize', 10, ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko', 'linewidth', 2, 'markersize', 7 );
%!	if (setAxisEqual)
%!		axis equal;
%!		axis equal; % Needed twice b/c of bug in Octave?
%!	end
%!	%axis( axis() + distScale*0.3*[-1,1,-1,1] );
%!	%ax1 = axis()
%!	%axis(ax1);
%!	grid on;
%!	%
%!	n = min([ 1, ceil( numPts * rand() ) ]);
%!	vecX = vecXPts(:,n);
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecX );
%!	%
%!	% Look at surface tangent.
%!	vecZ = 1e-4*(2.0*rand(sizeX,1)-1.0);
%!	[ vecSAtPlusZ ] = funchSurf( vecX + vecZ );
%!	[ vecSAtMinusZ ] = funchSurf( vecX - vecZ );
%!	vecTM = (vecSAtPlusZ+vecSAtMinusZ)/2.0 - distScale*(vecSAtPlusZ-vecSAtMinusZ)/norm(vecSAtPlusZ-vecSAtMinusZ);
%!	vecTP = (vecSAtPlusZ+vecSAtMinusZ)/2.0 + distScale*(vecSAtPlusZ-vecSAtMinusZ)/norm(vecSAtPlusZ-vecSAtMinusZ);
%!	%
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  vecSPts(1,:), vecSPts(2,:), 'ko-', 'markersize', 2, ...
%!	  vecXCent(1), vecXCent(2), 'k+', 'linewidth', 5, 'markersize', 10, ...
%!	  [ vecX(1); vecS(1) ], [ vecX(2); vecS(2) ], 's-', 'linewidth', 8, 'markersize', 10, ...
%!	  [ vecX(1); vecXCent(1) ], [ vecX(2); vecXCent(2) ], 'x-', 'linewidth', 4, 'markersize', 10, ...
%!	  [ vecS(1); vecS(1)+distScale*vecUHat(1) ], [ vecS(2); vecS(2)+distScale*vecUHat(2) ], '^-', 'linewidth', 3, 'markersize', 8, ...
%!	  [ vecS(1); vecS(1)+distScale*vecNHat(1) ], [ vecS(2); vecS(2)+distScale*vecNHat(2) ], 'v-', 'linewidth', 3, 'markersize', 8, ...
%!	  [ vecTM(1); vecTP(1) ], [ vecTM(2); vecTP(2) ], '-', 'linewidth', 2 );
%!	if (1)
%!	legend( ...
%!	  'surface', ...
%!	  'center', 'sample pt pulled to surf', 'sample pt pulled to center', ...
%!	  's-no-change direction', ...
%!	  'surface normal direction', 'surface tangent', ...
%!	  'location', 'NorthWestOutside' );
%!	end
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	%axis(ax1); %Can make graph disappear, b/c of Octave glitch?
%!	grid on;
%!	%
%!	msg( __FILE__, __LINE__, sprintf("Please check figures %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
