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
	%pVals = mygetfield( prm, "pVals", linspace(0.0,1.0,11) );
	%pVals = mygetfield( prm, "pVals", linspace(0.1,1.0,4) );
	%pVals = mygetfield( prm, "pVals", linspace(0.99,1.0,11) );
	%pVals = mygetfield( prm, "pVals", [0.0,0.05,0.1,0.5,0.9,0.95,1.0] );
	%pVals = mygetfield( prm, "pVals", [0.0,0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99,1.0] );
	%pVals = mygetfield( prm, "pVals", [0.0381] );
	%%%pVals = mygetfield( prm, "pVals", [0.0,0.01,0.038,0.0381,0.05,0.1,0.2,0.5,0.56,0.8,0.9,0.95,0.98,0.99,1.0] );
	pVals = mygetfield( prm, "pVals", [0.0,0.01,0.038,0.0381,0.05,0.1,0.2,0.56,0.8,0.9,1.0] );
	%pVals = mygetfield( prm, "pVals", [0.56] );
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeY,sizeY) );
	%
	useCnstAHack = mygetfield( prm, "useCnstAHack", false );
	%
	matIF = eye(sizeF,sizeF);
	numPVals = length(pVals);
	%
	%
	% Find min and max roots.
	zLo = -0.1;
	zHi = +0.1;
	for ip=1:numPVals
		p = pVals(ip);
		%
		matM = (1.0-p)*(matSY'*matSY) + p * (matW'*matW);
		matA = matIF - p*(matW*(matM\(matW')));
		if (useCnstAHack)
			matA = matIF;
		end
		c0 = p*(vecF0'*matA*vecLambda);
		c1 = (1.0-p)*(sz^2) + p*(vecLambda'*matA*vecLambda) + 2.0*p*(vecF0'*matA*vecEta);
		c2 = 3.0*p*(vecLambda'*matA*vecEta);
		c3 = 2.0*p*(vecEta'*matA*vecEta);
		funchXi = @(dummyZ)( c0 + c1*dummyZ + c2*(dummyZ.^2) + c3*(dummyZ.^3) );
		%
		zRoots = calcCubicRoots( c0, c1, c2, c3 );
		zLo = min([zLo,min(zRoots)]);
		zHi = max([zHi,max(zRoots)]);
	end
	zVals = mygetfield( prm, "zVals", linspace(zLo-0.2*(zHi-zLo),zHi+0.2*(zHi-zLo),1001) );
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
		if (useCnstAHack)
			matA = matIF;
		end
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
	%legend( cellAry_legend, 'location', 'southeast' );
	%%%legend( cellAry_legend, 'location', 'northwest' );
	%
	ax = axis();
	%ax(3) = -1.0;
	%ax(4) = 1.0;
	%axis(ax);
	plot( ...
	  [0.0,0.0], [ax(3),ax(4)], 'k-', ...
	  [ax(1),ax(2)], [0.0,0.0], 'k-' );
	%
	xiMin = 0.0;
	xiMax = 0.0;
	for ip=1:numPVals
		p = pVals(ip);
		%
		matM = (1.0-p)*(matSY'*matSY) + p * (matW'*matW);
		matA = matIF - p*(matW*(matM\(matW')));
		if (useCnstAHack)
			matA = matIF;
		end
		c0 = p*(vecF0'*matA*vecLambda);
		c1 = (1.0-p)*(sz^2) + p*(vecLambda'*matA*vecLambda) + 2.0*p*(vecF0'*matA*vecEta);
		c2 = 3.0*p*(vecLambda'*matA*vecEta);
		c3 = 2.0*p*(vecEta'*matA*vecEta);
		funchXi = @(dummyZ)( c0 + c1*dummyZ + c2*(dummyZ.^2) + c3*(dummyZ.^3) );
		%
		zRoots = calcCubicRoots( c0, c1, c2, c3 );
		numRoots = length(zRoots);
		switch (numRoots)
		case 0
			warning( "No root." );
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
			error( "Impossible." );
		end
	end
	%
	hold off;
return;
end


%!test
%!	%setprngstates();
%!	%setprngstates(30577088);
%!	%setprngstates(70188288);
%!	%setprngstates(47280832);
%!	%setprngstates(32188304);
%!	setprngstates();
%!	%
%!	prm = [];
%!	caseNum = 20
%!	switch(caseNum)
%!	case 0
%!		sizeX = 2;
%!		sizeF = sizeX;
%!		sizeY = sizeX-1;
%!		matW = randn(sizeF,sizeY);
%!		vecF0 = randn(sizeF,1);
%!		vecLambda = randn(sizeF,1);
%!		vecEta = randn(sizeF,1);
%!		makeQTLambdaZero = true;
%!	case 10
%!		matW = [ 1; 0.1 ]
%!		vecF0 = [ 0; 1.1 ]
%!		vecLambda = [ 0; -1 ]
%!		vecEta = [ 0; 0.2 ]
%!		makeQTLambdaZero = true;
%!	case 20
%!		% Hack, but, achieves "1 -> 3 -> 1 reconnect".
%!		matW = [ 1; 0 ]
%!		vecF0 = [ 100.0; 0.0 ]
%!		vecLambda = [ -1.0; -10.0 ]
%!		vecEta = [ 0; 1.0 ]
%!		makeQTLambdaZero = false;
%!		prm.useCnstAHack = true;
%!	case 21
%!		matW = [ 1; 0 ]
%!		vecF0 = [ 100.0; 0.0 ]
%!		vecLambda = [ -1.0; -10.0 ]
%!		vecEta = [ 0; 1.0 ]
%!		makeQTLambdaZero = false;
%!		prm.useCnstAHack = false;
%!	case 22
%!		matW = [ 1; 0.01 ]
%!		vecF0 = [ 100.0; 0.0 ]
%!		vecLambda = [ -1.0; -10.0 ]
%!		vecEta = [ 0; 1.0 ]
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = true;
%!	otherwise
%!		error( "Ivalid caseNum." );
%!	end
%!	%
%!	matQ = orth( matW );
%!	if (makeQTLambdaZero)
%!		vecLambda = vecLambda - matQ*(matQ'*vecLambda)
%!	end
%!	%
%!	vizFOCQRoots( vecF0, vecLambda, vecEta, matW, prm )
