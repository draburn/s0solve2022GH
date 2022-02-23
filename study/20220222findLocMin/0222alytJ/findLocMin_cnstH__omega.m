	error( "WORK IN PROGRESS." );
	
	muHi = 10.0*( abs(mu) + abs(hNorm) + abs(gNorm/deltaNormMin) );
	funchFZero = @(dummyMu)( norm( (matH+dummyMu*matI)\vecG ) - deltaNormTrgt );
	mu = fzero( funchFZero, [muLo, muHi] );
	matM = matH + mu*matI;
	matR = chol( matM );
	vecDelta = -( matR \ ( matR' \ vecG ) );
	%
	clear muHi;
	clear funchFZero;
