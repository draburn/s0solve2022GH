	clear;
	commondefs;
	thisFile = "study_curveModels";
	setprngstates(0);
	numFigs = 0;
	tic();
	%
	sizeX = 2000;
	sizeF = sizeX;
	matJ = randn(sizeF,sizeX);
	vecF0 = randn(sizeF,1);
	matJAnti = randn(sizeF,sizeX);
	%
	matH = matJ'*matJ;
	matH = matJ'*matJ - 0.0001*(matJAnti'*matJAnti);
	vecG = matJ'*vecF0;
	%omega0 = 0.5*(vecF0'*vecF0);
	
	if (0)
		sizeX = 2
		sizeF = 0
		vecG = [ 2; 1 ]
		matH = [ 1, 2; 2, 1 ]
	end
	
	hFrobNorm = sqrt(sum(sum(matH.^2)));
	tic();
	[ matPsi, matLam ] = eig(matH); vecLam = diag(matLam);
	toc();
	msg( thisFile, __LINE__, sprintf( "Eigenvalue range: %g ~ %g; hFrobNorm = %g.", ...
	  min(vecLam), max(vecLam), hFrobNorm ) );
	%
	matI = eye(sizeX,sizeX);
	%
	return
	numVals = 100;
	muViz0 = 20.0;
	muViz1 = hFrobNorm
	%muViz1 = sqrt(hFrobNorm)
	%muViz1 = 25.0;
	muVals = muViz0 + (muViz1-muViz0)*linspace(0.0,1.0,numVals).^2;
	kappaVals = zeros(1,numVals);
	kappaPVals = zeros(1,numVals);
	kappaPPVals = zeros(1,numVals);
	iotaVals = zeros(1,numVals);
	iotaPVals = zeros(1,numVals);
	xVals = zeros(1,numVals);
	aVals = zeros(1,numVals);
	bVals = zeros(1,numVals);
	cVals = zeros(1,numVals);
	dVals = zeros(1,numVals);
	for n=1:numVals
		mu = muVals(n);
		matM = matH + (mu*matI);
		[ matR, cholFlag ] = chol(matM);
		%assert(0==cholFlag)
		if (0==cholFlag)
			vecDelta = -(matR \ (matR'\vecG) );
			vecDeltaP = -(matR \ (matR'\vecDelta) );
			kappa = 0.5*(vecDelta'*vecDelta);
			kappaP = vecDelta'*vecDeltaP;
			iota = -( (vecG'*vecDelta) + (0.5*(vecDelta'*matH*vecDelta)) );
			iotaP = -( (vecG'*vecDeltaP) + (vecDelta'*matH*vecDeltaP) );
			%
			kappaVals(n) = kappa;
			kappaPVals(n) = kappaP;
			iotaVals(n) = iota;
			iotaPVals(n) = iotaP;
			%
			% Model: kappa = A/(mu+B)^2, iota = C/(mu+B) + D/(mu+B)^2.
			% X = mu + B.
			cnstX = -2.0*kappa/kappaP;
			cnstA = kappa * (cnstX^2);
			cnstB = cnstX - mu;
			cnstC = cnstX*( (2.0*iota) + (cnstX*iotaP) );
			cnstD = -(cnstX^2)*( iota + (cnstX*iotaP) );
			%
			xVals(n) = cnstX;
			aVals(n) = cnstA;
			bVals(n) = cnstB;
			cVals(n) = cnstC;
			dVals(n) = cnstD;
			%
			% Model: kappa as above, iota = Q/(mu+B)^P.
			cnstP = -cnstX*iotaP/iota;
			cnstQ = iota*(cnstX^cnstP);
			pVals(n) = cnstP;
			qVals(n) = cnstQ;
			%
			%
			vecDeltaPP = -2.0*(matR \ (matR'\vecDeltaP) );
			kappaPP = (vecDelta'*vecDeltaPP) + (vecDeltaP'*vecDeltaP);
			kappaPPVals(n) = kappaPP;
			%
			% Model: kappa = R / (mu+S)^T, iota = U / (mu+S)^V.
			cnstT = (kappaP^2)/( (kappa*kappaPP) - (kappaP^2) );
			cnstY = -cnstT*kappa/kappaP;
			cnstS = cnstY - mu;
			if (1)
			% Mod 2021.10.27.
			if (cnstS > 0.0)
				cnstS = 0.0;
				cnstY = mu;
				cnstT = -cnstY*kappaP/kappa;
			end
			end
			cnstR = kappa * (cnstY^cnstT);
			tVals(n) = cnstT;
			sVals(n) = cnstS;
			rVals(n) = cnstR;
			%
			cnstV = -cnstY*iotaP/iota;
			cnstU = iota*(cnstY^cnstV);
			vVals(n) = cnstV;
			uVals(n) = cnstU;
			%
			% Model kappa = kappa1 * (mu1/mu)^K, iota = iota1 * (mu1/mu)^L.
			cnstK = -mu*kappaP/kappa;
			cnstL = -mu*iotaP/iota;
			kVals(n) = cnstK;
			lVals(n) = cnstL;
		end
	end
	%
	%
	%
	%n1A = 1+floor(0.1*(numVals-1))
	%n1A = 1+floor(0.5*(numVals-1))
	%n1A = 1+floor(0.6*(numVals-1))
	%n1A = 1+floor(0.8*(numVals-1))
	n1A = numVals
	mu1A = muVals(n1A)
	cnstA1A = aVals(n1A)
	cnstB1A = bVals(n1A)
	cnstC1A = cVals(n1A)
	cnstD1A = dVals(n1A)
	kappa1AVals = cnstA1A./((muVals+cnstB1A).^2);
	iota1AVals = cnstC1A./(muVals+cnstB1A) + cnstD1A./((muVals+cnstB1A).^2);
	%
	cnstP1A = pVals(n1A)
	cnstQ1A = qVals(n1A)
	iotaPQ1AVals = cnstQ1A./((muVals+cnstB1A).^cnstP1A);
	%
	cnstR1A = rVals(n1A)
	cnstS1A = sVals(n1A)
	cnstT1A = tVals(n1A)
	kappaRST1AVals = cnstR1A./((muVals+cnstS1A).^cnstT1A);
	cnstU1A = uVals(n1A)
	cnstV1A = vVals(n1A)
	iotaUV1AVals = cnstU1A./((muVals+cnstS1A).^cnstV1A);
	%
	cnstL1A = lVals(n1A)
	cnstK1A = kVals(n1A)
	kappaK1AVals = kappaVals(n1A)*((mu1A./muVals).^cnstK1A);
	iotaL1AVals = iotaVals(n1A)*((mu1A./muVals).^cnstL1A);
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  muVals, iotaVals, 'o-', ...
	  muVals, iota1AVals, 'x-', ...
	  mu1A, iotaVals(n1A), 's', 'markersize', 20, ...
	  muVals, iotaPQ1AVals, '^-', ...
	  muVals, iotaUV1AVals, '+-', ...
	  muVals, iotaL1AVals, '*-' );
	grid on;
	xlabel( 'mu' );
	ylabel( 'iota' );
	title( 'iota vs mu' );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  muVals, kappaVals, 'o-', ...
	  muVals, kappa1AVals, 'x-', ...
	  mu1A, kappaVals(n1A), 's', 'markersize', 20, ...
	  muVals, kappa1AVals, '-', ...
	  muVals, kappaRST1AVals, '+-', ...
	  muVals, kappaK1AVals, '*-' );
	grid on;
	xlabel( 'mu' );
	ylabel( 'kappa' );
	title( 'kappa vs mu' );
	
	
	return
	
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  muVals, 1.0./iotaVals, 'o-', ...
	  muVals, 1.0./iota1AVals, 'x-', ...
	  mu1A, 1.0./iotaVals(n1A), 's', 'markersize', 20, ...
	  muVals, 1.0./iotaPQ1AVals, '^-', ...
	  muVals, 1.0./iotaUV1AVals, '+-', ...
	  muVals, 1.0./iotaL1AVals, '*-' );
	grid on;
	xlabel( 'mu' );
	ylabel( '1/iota' );
	title( '1/iota vs mu' );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  muVals, 1.0./kappaVals, 'o-', ...
	  muVals, 1.0./kappa1AVals, 'x-', ...
	  mu1A, 1.0./kappaVals(n1A), 's', 'markersize', 20, ...
	  muVals, 1.0./kappa1AVals, '-', ...
	  muVals, 1.0./kappaRST1AVals, '+-', ...
	  muVals, 1.0./kappaK1AVals, '*-' );
	grid on;
	xlabel( 'mu' );
	ylabel( '1/kappa' );
	title( '1/kappa vs mu' );
	%
	%
	%
	numFigs++; figure(numFigs);
	plot( ...
	  kappaVals, iotaVals, 'o-', ...
	  kappa1AVals, iota1AVals, 'x-', ...
	  kappaVals(n1A), iotaVals(n1A), 's', 'markersize', 20, ...
	  kappa1AVals, iotaPQ1AVals, '-', ...
	  kappaRST1AVals, iotaUV1AVals, '+-' );
	grid on;
	xlabel( 'kappa' );
	ylabel( 'iota' );
	title( 'iota vs kappa' );

