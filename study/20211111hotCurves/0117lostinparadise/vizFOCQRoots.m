function [ datOut ] = vizFOCQRoots( vecF0, vecLambda, vecEta, matW, prm=[] )
	%
	% model...
	%  vecF = vecF0 + (vecLambda*z) + (vecEta*z^2) + matW*vecY;
	%  xi = p*0.5*vecF'*vecF + (1-p)*0.5*( z^2 + vecY'*vecY );
	% look at critical points of xi(p) over z, vecY.
	%
	sizeF = size(matW,1);
	sizeY = size(matW,2);
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( isrealarray(vecLambda,[sizeF,1]) );
	assert( isrealarray(vecEta,[sizeF,1]) );
	assert( isrealarray(matW,[sizeF,sizeY]) );
	%
	figNum = mygetfield( prm, "figNum", 20 );
	%pVals = mygetfield( prm, "pVals", linspace(0.0,1.0,6) );
	pVals = mygetfield( prm, "pVals", linspace(0.1,1.0,4) );
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeY,sizeY) );
	zVals = mygetfield( prm, "zVals", linspace(-2.0,4.0,1001) );
	%
	matIF = eye(sizeF,sizeF);
	numPVals = length(pVals);
	numZVals = length(zVals);
	%
	figure(figNum);
	clf;
	hold on;
	datXi = zeros(numZVals,numPVals);
	for ip=1:numPVals
		p = pVals(ip);
		cellAry_legend(1,ip) = { sprintf( "p = %0.3f", p ) };
		%
		matM = (1.0-p)*(matSY'*matSY) + p * (matW'*matW);
		matA = matIF - p*(matW*(matM\(matW')));
		c0 = p*(vecF0'*matA*vecLambda);
		c1 = (1.0-p)*(sz^2) + p*(vecLambda'*matA*vecLambda) + 2.0*p*(vecF0'*matA*vecEta);
		c2 = 3.0*p*(vecLambda'*matA*vecEta);
		c3 = 2.0*p*(vecEta'*matA*vecEta);
		funchXi = @(dummyZ)( c0 + c1*dummyZ + c2*(dummyZ.^2) + c3*(dummyZ.^3) );
		%
		datXi(:,ip) = funchXi(zVals);
		plot( zVals, datXi(:,ip), '-', 'linewidth', 3 );
	end
	grid on;
	legend( cellAry_legend, 'location', 'southeast' );
	%
	ax = axis();
	plot( ...
	  [0.0,0.0], [ax(3),ax(4)], 'k-', ...
	  [ax(1),ax(2)], [0.0,0.0], 'k-' );
	%
	for ip=1:numPVals
		p = pVals(ip);
		%
		matM = (1.0-p)*(matSY'*matSY) + p * (matW'*matW);
		matA = matIF - p*(matW*(matM\(matW')));
		c0 = p*(vecF0'*matA*vecLambda);
		c1 = (1.0-p)*(sz^2) + p*(vecLambda'*matA*vecLambda) + 2.0*p*(vecF0'*matA*vecEta);
		c2 = 3.0*p*(vecLambda'*matA*vecEta);
		c3 = 2.0*p*(vecEta'*matA*vecEta);
		funchXi = @(dummyZ)( c0 + c1*dummyZ + c2*(dummyZ.^2) + c3*(dummyZ.^3) );
		%
		zRoots = calcCubicRoots( c0, c1, c2, c3 );
		numRoots = length(zRoots);
		switch (numRoots)
		case 1
			plot( ...
			  zRoots(1), 0.0, 'ko', 'markersize', 10, ...
			  zRoots(1), 0.0, 'k*', 'markersize', 10 );
		case 2
			plot( ...
			  zRoots(1), 0.0, 'co', 'markersize', 10, ...
			  zRoots(1), 0.0, 'c*', 'markersize', 10, ...
			  zRoots(2), 0.0, 'mo', 'markersize', 10, ...
			  zRoots(2), 0.0, 'm*', 'markersize', 10 );
		case 3
			plot( ...
			  zRoots(1), 0.0, 'ro', 'markersize', 10, ...
			  zRoots(1), 0.0, 'r*', 'markersize', 10, ...
			  zRoots(2), 0.0, 'go', 'markersize', 10, ...
			  zRoots(2), 0.0, 'g*', 'markersize', 10, ...
			  zRoots(3), 0.0, 'bo', 'markersize', 10, ...
			  zRoots(3), 0.0, 'b*', 'markersize', 10 );
		otherwise
			error( "?" );
		end
	end
	%
	hold off;
return;
end


%!test
%!	setprngstates(0);
%!	%
%!	sizeX = 2;
%!	sizeF = sizeX;
%!	sizeY = sizeX-1;
%!	matW = randn(sizeF,sizeY);
%!	vecF0 = randn(sizeF,1);
%!	vecLambda = randn(sizeF,1);
%!	vecEta = randn(sizeF,1);
%!	%
%!	matQ = orth( matW );
%!	makeQTLambdaZero = true;
%!	if (makeQTLambdaZero)
%!		vecLambda = vecLambda - matQ*(matQ'*vecLambda);
%!	end
%!	%
%!	prm = [];
%!	vizFOCQRoots( vecF0, vecLambda, vecEta, matW, prm )
