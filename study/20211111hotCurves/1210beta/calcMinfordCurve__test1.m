clear;
thisFile = "calcMinfordCurve__test1";
commondefs;
setprngstates(0);
numFigs = 0;
startTime = time();
msg( thisFile, __LINE__, "Please run __test3 for more thorough testing." );
return;
%
%
switch (20)
case 1
	sizeF = 2;
	sizeX = 2;
	funchOmega = @(x)(0.5*(x'*x));
	funchG = @(x)(x);
	funchH = @(x)(eye(2,2));
case 2
	sizeF = 2;
	sizeX = 2;
	omega0 = 0.0;
	vecG0 = [0.0;0.0];
	matH0 = [1.0,0.0;0.0,4.0];
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
case 10
	sizeF = 2;
	sizeX = 2;
	vecX0 = randn(sizeX,1);
	omega0 = 0.0;
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'*matHGen;
	funchOmega0 = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG0 = @(x)( vecG0 + matH0*x );
	funchH0 = @(x)( matH0 );
	funchOmega = @(x)(funchOmega0(x-vecX0));
	funchG = @(x)(funchG0(x-vecX0));
	funchH = @(x)(funchH0(x-vecX0));
case 11
	sizeF = 3;
	sizeX = 3;
	omega0 = abs(randn());
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'*matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
case 20
	sizeF = 10;
	sizeX = 10;
	omega0 = abs(randn());
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'+matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
otherwise
	error( "Invalid case." );
end
%
%
%
numTrials = 10;
msg( thisFile, __LINE__, sprintf( "Performing %d trials...", numTrials ) );
time0 = time();
for t=1:numTrials
	vecXC = randn(sizeX,1);
	bigR = abs(randn);
	vecX = randn(sizeX,1);
	matS = abs(randn(sizeX,sizeX));
	%
	if (0)
		assert( 2==sizeX )
		vecXC = [0.0;0.0];
		bigR = 1.0;
		vecX = [1.0;0.0];
		matS = [1.0,0.0;0.0,1.0];
	elseif (0)
		assert( 2==sizeX )
		matS = [1.0,0.0;0.0,1.0];
	end
	%
	if ( 1==t )
		msg( thisFile, __LINE__, "Sample calculation..." );
		echo__vecXC = vecXC
		echo__bigR = bigR
		echo__vecX = vecX
		echo__matS = matS
	end
	%
	%
	%
	vecXSurf = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
	vecXSurf2 = calcMinfordCurve__evalXSurf( vecXC, bigR, vecXSurf );
	if ( 1==t )
		echo__vecXSurf = vecXSurf
		echo__vecXSurf2 = vecXSurf2
	end
	assert( norm(vecXSurf-vecXSurf2)<=1e-4*(norm(vecXSurf)+norm(vecXSurf2)) );
	%
	%
	%
	vecXSurfScale = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS );
	vecXSurfScale2 = calcMinfordCurve__evalXSurf( vecXC, bigR, vecXSurfScale, matS );
	if ( 1==t )
		echo__vecXSurfScale = vecXSurfScale
		echo__vecXSurfScale2 = vecXSurfScale2
	end
	assert( norm(vecXSurfScale-vecXSurfScale2)<=1e-4*(norm(vecXSurfScale)+norm(vecXSurfScale2)) );
	%
	%
	%
	omega_unSurf = funchOmega(vecX);
	if ( 1==t )
		echo__omega_unSurf = omega_unSurf
	end
	%
	%
	%
	omegaSurf = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX );
	altOmegaSurf = funchOmega( calcMinfordCurve__evalXSurf( vecXC, bigR, vecX ) );
	if ( 1==t )
		echo__omegaSurf = omegaSurf
		echo__altOmegaSurf = altOmegaSurf
	end
	assert( norm(omegaSurf-altOmegaSurf)<=1e-4*(norm(omegaSurf)+norm(altOmegaSurf)) );
	%
	%
	%
	omegaSurfScale = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX, matS );
	altOmegaSurfScale = funchOmega( calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS ) );
	if ( 1==t )
		echo__omegaSurfScale = omegaSurfScale
		echo__altOmegaSurfScale = altOmegaSurfScale
	end
	assert( norm(omegaSurfScale-altOmegaSurfScale)<=1e-4*(norm(omegaSurfScale)+norm(altOmegaSurfScale)) );
	%
	%
	%
	vecGSurf = calcMinfordCurve__evalGSurf( funchG, vecXC, bigR, vecX );
	vecGSurf_fd = zeros(sizeX,1);
	for n=1:sizeX
		epsFD = 1e-4;
		vecXP = vecX;
		vecXM = vecX;
		vecXP(n) += epsFD;
		vecXM(n) -= epsFD;
		omegaP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXP );
		omegaM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXM );
		vecGSurf_fd(n) = ( omegaP - omegaM ) / (2.0*epsFD);
	end
	if ( 1==t )
		echo__vecGSurf = vecGSurf
		echo__vecGSurf_fd = vecGSurf_fd
	end
	assert( norm(vecGSurf-vecGSurf_fd)<=1e-4*(norm(vecGSurf)+norm(vecGSurf_fd)) );
	%
	%
	%
	vecGSurfScaled = calcMinfordCurve__evalGSurf( funchG, vecXC, bigR, vecX, matS );
	vecGSurfScaled_fd = zeros(sizeX,1);
	for n=1:sizeX
		epsFD = 1e-4;
		vecXP = vecX;
		vecXM = vecX;
		vecXP(n) += epsFD;
		vecXM(n) -= epsFD;
		omegaP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXP, matS );
		omegaM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXM, matS );
		vecGSurfScaled_fd(n) = ( omegaP - omegaM ) / (2.0*epsFD);
	end
	if ( 1==t )
		echo__vecGSurfScaled = vecGSurfScaled
		echo__vecGSurfScaled_fd = vecGSurfScaled_fd
	end
	assert( norm(vecGSurfScaled-vecGSurfScaled_fd)<=1e-4*(norm(vecGSurfScaled)+norm(vecGSurfScaled_fd)) );
	%
	%
	%
	matHSurf = calcMinfordCurve__evalHSurf( funchG, funchH, vecXC, bigR, vecX );
	omega00 = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX );
	matHSurf_fd2 = zeros(sizeX,sizeX);
	for m=1:sizeX
	for n=1:sizeX
		epsFD = 1e-4;
		vecXPP = vecX;
		vecXMP = vecX;
		vecXPM = vecX;
		vecXMM = vecX;
		vecXPP(m) += epsFD; vecXPP(n) += epsFD;
		vecXMP(m) -= epsFD; vecXMP(n) += epsFD;
		vecXPM(m) += epsFD; vecXPM(n) -= epsFD;
		vecXMM(m) -= epsFD; vecXMM(n) -= epsFD;
		omegaPP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPP );
		omegaMP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMP );
		omegaPM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPM );
		omegaMM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMM );
		matHSurf_fd2(m,n) = ( omegaPP + omegaMM - omegaMP - omegaPM ) / (4.0*(epsFD)^2);
	end
	end
	for m=1:sizeX
		epsFD = 1e-4;
		vecXP = vecX;
		vecXM = vecX;
		vecXP(m) += epsFD;
		vecXM(m) -= epsFD;
		omegaP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXP );
		omegaM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXM );
		matHSurf_fd2(m,m) = ( omegaP + omegaM - 2.0*omega00 ) / (epsFD^2);
	end
	if (1==t)
		echo__matHSurf = matHSurf
		echo__matHSurf_fd2 = matHSurf_fd2
	end
	assert( sqrt(sum(sum((matHSurf_fd2-matHSurf).^2))) <= ...
	  1e-3*(  sqrt(sum(sum((matHSurf_fd2).^2)))  + sqrt(sum(sum((matHSurf).^2)))  ) );
	%
	%
	%
	matHSurfScaled = calcMinfordCurve__evalHSurf( funchG, funchH, vecXC, bigR, vecX, matS );
	omega00Scaled = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX, matS );
	matHSurfScaled_fd2 = zeros(sizeX,sizeX);
	for m=1:sizeX
	for n=1:sizeX
		epsFD = 1e-4;
		vecXPP = vecX;
		vecXMP = vecX;
		vecXPM = vecX;
		vecXMM = vecX;
		vecXPP(m) += epsFD; vecXPP(n) += epsFD;
		vecXMP(m) -= epsFD; vecXMP(n) += epsFD;
		vecXPM(m) += epsFD; vecXPM(n) -= epsFD;
		vecXMM(m) -= epsFD; vecXMM(n) -= epsFD;
		omegaPP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPP, matS );
		omegaMP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMP, matS );
		omegaPM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPM, matS );
		omegaMM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMM, matS );
		matHSurfScaled_fd2(m,n) = ( omegaPP + omegaMM - omegaMP - omegaPM ) / (4.0*(epsFD)^2);
	end
	end
	for m=1:sizeX
		epsFD = 1e-4;
		vecXP = vecX;
		vecXM = vecX;
		vecXP(m) += epsFD;
		vecXM(m) -= epsFD;
		omegaP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXP, matS );
		omegaM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXM, matS );
		matHSurfScaled_fd2(m,m) = ( omegaP + omegaM - 2.0*omega00Scaled ) / (epsFD^2);
	end
	if (1==t)
		echo__matHSurfScaled = matHSurfScaled
		echo__matHSurfScaled_fd2 = matHSurfScaled_fd2
	end
	assert( sqrt(sum(sum((matHSurfScaled_fd2-matHSurfScaled).^2))) <= ...
	  1e-3*(  sqrt(sum(sum((matHSurfScaled_fd2).^2))) + sqrt(sum(sum((matHSurfScaled).^2)))  ) );
end
msg( thisFile, __LINE__, sprintf( "Passed %d trials in %0.3fs.", numTrials, time()-time0 ) );
