function [ vecXVals, datOut ] = vizFOCQLevCurve( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm=[] )
	%
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	doChecks = mygetfield( prm, "doChecks", true );
	%
	figNum = mygetfield( prm, "figNum", 10 );
	vecXFig = zeros(sizeX,1); vecXFig(1) = 1.0;
	vecYFig = zeros(sizeX,1); vecYFig(2) = 1.0;
	vecXFig = mygetfield( prm, "vecXFig", vecXFig );
	vecYFig = mygetfield( prm, "vecYFig", vecYFig );
	figure(figNum);
	clf;
	plot( vecXFig'*vecX0, vecYFig'*vecX0, 'ks', 'linewidth', 3, 'markersize', 25 );
	grid on;
	hold on;
	%
	curveSelector = mygetfield( prm, "curveSelector", 0 );
	assert( isscalar(doChecks) )
	if ( doChecks );
		assert( isbool(doChecks) );
		assert( 2 <= sizeX );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		assert( isrealarray(vecPhiHat,[sizeX,1]) );
		assert( abs(norm(vecPhiHat)-1.0) < (eps^0.75)*sizeX );
		assert( isrealarray(vecEta,[sizeF,1]) );
		%
		assert( isrealscalar(curveSelector) );
	end
	%
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeX-1,sizeX-1) );
	%
	vecXVals = [];
	datOut = [];
	%
	%
	% DRaburn 2022.01.29:
	%  Contrary to notes, take matPsi to span all of sizeX except vecPhiHat.
	matIX = eye(sizeX,sizeX);
	matPsi = orth( matIX - (vecPhiHat*(vecPhiHat')) );
	if ( doChecks )
		assert( isrealarray(matPsi,[sizeX,sizeX-1]) );
		matShy = [ vecPhiHat, matPsi ];
		%echo__shyShyT = matShy*(matShy')
		%echo__shyTShy = (matShy')*matShy
		assert( sum(sum(abs(matShy*(matShy')-matIX))) < (eps^0.50)*sizeX );
		assert( sum(sum(abs((matShy')*matShy-matIX))) < (eps^0.50)*sizeX );
		clear matShy;
	end
	vecLambda = matJ0 * vecPhiHat;
	matW = matJ0 * matPsi;
	szsq = sz*sz;
	matD = matSY'*matSY;
	matWTW = matW'*matW;
	%matIY = eye(sizeX-1,sizeX-1);
	%echo__vecEta = vecEta
	%
	%
	numVals = 101;
	vecXVals = zeros(sizeX,numVals);
	vecXVals(:,1) = vecX0;
	for n=2:numVals
		p = (n-1.0)/(numVals-1.0);
		matM = ((1.0-p)*matD) + (p*matWTW);
		matA = matIX - (p*(matW*(matM\(matW'))));
		%
		c0 = p * ( vecF0' * matA * vecLambda );
		c1 = (1.0-p)*szsq + p*( vecLambda' * matA * vecLambda ) + 2.0*p*( vecF0' * matA * vecEta )
		c2 = 3.0*p*( vecLambda' * matA * vecEta );
		c3 = 2.0*p*( vecEta' * matA * vecEta );
		%
		zPts = calcCubicRoots( c0, c1, c2, c3 );
		numPts = size(zPts,2);
		assert( isrealarray(zPts,[1,numPts]) );
		vecZetaPts = vecF0 + (vecLambda*zPts) + (vecEta*(zPts.^2));
		vecYPts = -p*( matM\(matW'*vecZetaPts) );
		vecXPts = vecX0 + (vecPhiHat*zPts) + (matPsi*vecYPts);
		switch (numPts)
		case 1
			vecXVals(:,n) = vecXPts;
			plot( vecXFig'*vecXPts, vecYFig'*vecXPts, 'ko', 'linewidth', 1, 'markersize', 10 );
		case 3
			if ( curveSelector < -0.5 )
				vecXVals(:,n) = vecXPts(:,1);
			elseif ( curveSelector < 0.5 )
				vecXVals(:,n) = vecXPts(:,2);
			else
				vecXVals(:,n) = vecXPts(:,3);
			end
			plot( ...
			  vecXFig'*vecXPts(:,1), vecYFig'*vecXPts(:,1), 'rv', 'linewidth', 1, 'markersize', 10, ...
			  vecXFig'*vecXPts(:,2), vecYFig'*vecXPts(:,2), 'bx', 'linewidth', 1, 'markersize', 10, ...
			  vecXFig'*vecXPts(:,3), vecYFig'*vecXPts(:,3), 'g^', 'linewidth', 1, 'markersize', 10 );
		otherwise
			error("Invalid case.");
		end
	end
	%
	switch (numPts)
	case 1
		vecXVals(:,n) = vecXPts;
		plot( vecXFig'*vecXPts, vecYFig'*vecXPts, 'ko', 'linewidth', 2, 'markersize', 20 );
	case 3
		if ( curveSelector < -0.5 )
			vecXVals(:,n) = vecXPts(:,1);
		elseif ( curveSelector < 0.5 )
			vecXVals(:,n) = vecXPts(:,2);
		else
			vecXVals(:,n) = vecXPts(:,3);
		end
		plot( ...
		  vecXFig'*vecXPts(:,1), vecYFig'*vecXPts(:,1), 'rv', 'linewidth', 2, 'markersize', 20, ...
		  vecXFig'*vecXPts(:,2), vecYFig'*vecXPts(:,2), 'bx', 'linewidth', 2, 'markersize', 20, ...
		  vecXFig'*vecXPts(:,3), vecYFig'*vecXPts(:,3), 'g^', 'linewidth', 2, 'markersize', 20 );
	otherwise
		error("Invalid case.");
	end
	%
	% CONSIDER EXTRAPOLATION HERE,
	%  AS PER OTHER CURVES!
return;
end


%!test
%!	setprngstates(0);
%!	%
%!	sizeF = 2 + round(abs(randn()));
%!	sizeX = 2 + round((sizeF-2.0)*rand());
%!	%
%!	vecF0 = randn(sizeF,1);
%!	matJ0 = randn(sizeF,sizeX);
%!	vecX0 = randn(sizeX,1);
%!	%
%!	matH00 = matJ0'*matJ0;
%!	[ matPsi00, matLambda00 ] = eig(matH00);
%!	[ lambda00AbsMin00, nOfAbsMin00 ] = min(abs(diag(matLambda00)));
%!	vecPhiHat = matPsi00(:,nOfAbsMin00);
%!	vecEta = randn(sizeF,1);
%!	%
%!	prm = [];
%!	calcFOCQLevCurve( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
