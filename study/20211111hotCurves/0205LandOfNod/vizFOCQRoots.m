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
	%pVals = mygetfield( prm, "pVals", linspace(0.0,1.0,11) );
	pVals = mygetfield( prm, "pVals", [0.0,0.01,0.02,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.98,0.99,1.00] );
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeY,sizeY) );
	%
	useCnstAHack = mygetfield( prm, "useCnstAHack", false );
	useDivPHack = mygetfield( prm, "useDivPHack", false );
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
		if (useDivPHack)
			datXi(:,ip)/=(eps+p);
		end
		plot( zVals, datXi(:,ip), '-', 'linewidth', 3 );
	end
	grid on;
	%legend( cellAry_legend, 'location', 'southeast' );
	legend( cellAry_legend, 'location', 'northwest' );
	%
	ax = axis();
	%ax(3) = -100.0;
	%ax(4) = 100.0;
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
	%
	%
	%
	pVals = linspace(0,1,101);
	figure(figNum+1);
	plot( 0, 0, 'ko', 'markersize', 10 );
	hold on;
	plot( 0, 0, 'kx', 'markersize', 10 );
	for ip=1:length(pVals)
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
			  zRoots(1), p, 'ko', 'markersize', 10, ...
			  zRoots(1), p, 'k*', 'markersize', 10 );
		case 2
			plot( ...
			  zRoots(1), p, 'co', 'markersize', 10, ...
			  zRoots(1), p, 'c*', 'markersize', 10, ...
			  zRoots(2), p, 'mo', 'markersize', 10, ...
			  zRoots(2), p, 'm*', 'markersize', 10 );
		case 3
			plot( ...
			  zRoots(1), p, 'ro', 'markersize', 10, ...
			  zRoots(1), p, 'r*', 'markersize', 10, ...
			  zRoots(2), p, 'go', 'markersize', 10, ...
			  zRoots(2), p, 'g*', 'markersize', 10, ...
			  zRoots(3), p, 'bo', 'markersize', 10, ...
			  zRoots(3), p, 'b*', 'markersize', 10 );
		otherwise
			error( "Impossible." );
		end
	end
	grid on;
	hold off;
	%
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
%!	caseNum = 53
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
%!	case 30
%!		% 1->3->1 reconnect in 3D without only sizeY hack!
%!		matW = [ 0; 0; 1 ] % NOT 3x2???
%!		vecF0 = [ 100.0; 0.0; 0.0 ]
%!		vecLambda = [ -1.0; -10.0; 0.0 ]
%!		vecEta = [ 0; 1.0; 0.0 ]
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!	case 31
%!		% 1->3->1 reconnect in 3D with JTJ having nullspace dim >= 1..
%!		% No good.
%!		matW = [ 0, 0.01; 0, 0; 1, 2 ]
%!		vecF0 = [ 100.0; 0.0; 0.0 ]
%!		vecLambda = [ -1.0; -10.0; 0.0 ]
%!		vecEta = [ 0; 1.0; 0.0 ]
%!		makeQTLambdaZero = false;
%!		prm.useCnstAHack = false;
%!	case 40
%!		% Try with d > 0...
%!		matW = [ 0, 0.01; 0.01, 0; 1, 2 ]
%!		vecF0 = [ 100.0; 0.0; 0.0 ]
%!		vecLambda = [ 1.0; -10.0; 0.0 ]
%!		vecEta = [ 0; -1.0; 0.0 ]
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!	case 50
%!		% Try 2D 1->3->1 reconn, with reporting values.
%!		% Pops up at end!
%!		matW = [ 1; 0 ]
%!		vecF0 = [ 3; -1 ]
%!		vecLambda = [ 0; 3 ]
%!		vecEta = [ -1; -2]
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!	case 51
%!		% Try 2D 1->3->1 reconn, with reporting values.
%!		% YES! YES! YESSSS!!!!!
%!		matW = [ 1; 0 ]
%!		theta = 0.1
%!		eta = 1
%!		lambda = -3
%!		g = -30
%!		f = 1
%!		vecF0 = [ g; f ];
%!		vecLambda = [ 0; lambda ];
%!		vecEta = [ theta; eta ];
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!		prm.useDivPHack = false;
%!	case 52
%!		% Look at marginal connection.
%!		matW = [ 1; 0 ]
%!		c = 0.79;
%!		theta = 0.1*c
%!		eta = 1*c
%!		lambda = -3*c
%!		g = -30*c
%!		f = 1*c
%!		vecF0 = [ g; f ];
%!		vecLambda = [ 0; lambda ];
%!		vecEta = [ theta; eta ];
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!		prm.useDivPHack = true;
%!		discrim = 12*eta^2 - 48*f*eta^3
%!		ze = -( 1 + sqrt(3-12*f*eta)*[+1,-1]/3 )/(2*eta)
%!	case 53
%!		% Get case with smallest eigenvalue?
%!		matW = [ 3.1; 0 ]
%!		theta = 0.1
%!		eta = 1
%!		lambda = -3
%!		g = -100
%!		f = 1
%!		vecF0 = [ g; f ];
%!		vecLambda = [ 0; lambda ];
%!		vecEta = [ theta; eta ];
%!		makeQTLambdaZero = true;
%!		prm.useCnstAHack = false;
%!		prm.useDivPHack = false;
%!		discrim = 12*eta^2 - 48*f*eta^3
%!		ze = -( 1 + sqrt(3-12*f*eta)*[+1,-1]/3 )/(2*eta)
%!	otherwise
%!		error( "Ivalid caseNum." );
%!	end
%!	%prm.pVals = mygetfield( prm, "pVals", linspace(0.1,1.0,10) );
%!	prm.pVals = mygetfield( prm, "pVals", linspace(0,1,11) );
%!	%prm.pVals = mygetfield( prm, "pVals", [0.8,0.9,1.0] );
%!	sizeF = size(matW,1);
%!	%
%!	matQ = orth( matW );
%!	if (makeQTLambdaZero)
%!		vecLambda = vecLambda - matQ*(matQ'*vecLambda)
%!	end
%!	%
%!	reportVals = true;
%!	if (reportVals)
%!		echo__matW = matW
%!		echo__vecF0 = vecF0
%!		echo__vecLambda = vecLambda
%!		echo__vecEta = vecEta
%!		%
%!		matA1 = eye(sizeF,sizeF) - matQ*(matQ');
%!		echo_a = 2.0*vecEta'*[ vecEta, matA1*vecEta ];
%!		echo_b = 3.0*vecEta'*[ vecLambda, matA1*vecLambda ];
%!		echo_c = vecLambda'*[ vecLambda, matA1*vecLambda ] + 2.0*vecF0'*[ vecEta, matA1*vecEta ];
%!		echo_cPlusMu = [ +Inf, echo_c(2) ];
%!		echo_d = vecF0'*[ vecLambda, matA1*vecLambda ];
%!		echo_coeff = [ echo_a; echo_b; echo_c; echo_cPlusMu; echo_d ]
%!	end
%!	%
%!	vizFOCQRoots( vecF0, vecLambda, vecEta, matW, prm )
