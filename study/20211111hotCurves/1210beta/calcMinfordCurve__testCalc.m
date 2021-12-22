function calcMinfordCurve__testCalc( sizeX, funchOmega, funchG, funchH, prm=[] )
	commondefs;
	thisFile = "calcMinfordCurve__testCalc";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%
	vecXC = randn(sizeX,1);
	bigR = abs(randn);
	vecX = randn(sizeX,1);
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Test values..." );
		echo__vecXC = vecXC
		echo__bigR = bigR
		echo__vecX = vecX
	end
	if ( abs(vecX-vecXC) < eps*bigR )
		msg( thisFile, __LINE__, "vecX is very nearly vecXC." );
		msg( thisFile, __LINE__, "This could happen by chance but is unlikely." );
		msg( thisFile, __LINE__, "Aborting test." );
		error( "Generated vecX is on top of vecXC." );
	end
	%
	doSurfDirTests = true;
	%
	%
	% Test that xSurf moves the point, but only once.
	vecXSurf = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
	vecXSurf2 = calcMinfordCurve__evalXSurf( vecXC, bigR, vecXSurf );
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing xSurf..." );
		echo__vecXSurf = vecXSurf
		echo__vecXSurf2 = vecXSurf2
	end
	if ( abs(vecXSurf-vecX) < eps*bigR )
		msg( thisFile, __LINE__, "__evalXSurf did not meaningfully move the point." );
		msg( thisFile, __LINE__, "This could happen by chance but is unlikely." );
		msg( thisFile, __LINE__, "Aborting test." );
		error( "Unscaled vecXSurf is on top of vecX." );
	end
	assert( norm(vecXSurf-vecXSurf2)<=1e-4*(norm(vecXSurf)+norm(vecXSurf2)) );
	clear vecXSurf;
	clear vecXSurf2;
	%
	%
	% Test omega.
	omegaSurf = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX );
	omegaSurf_alt = funchOmega( calcMinfordCurve__evalXSurf( vecXC, bigR, vecX ) );
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing omegaSurf..." );
		echo__omega = funchOmega( vecX );
		echo__omegaSurf = omegaSurf
		echo__omegaSurf_alt = omegaSurf_alt
	end
	assert( norm(omegaSurf-omegaSurf_alt)<=1e-4*(norm(omegaSurf)+norm(omegaSurf_alt)) );
	clear omegaSurf;
	clear omegaSurf_alt;
	%
	%
	% Test gradient; note that result is inexact.
	vecGSurf = calcMinfordCurve__evalGSurf( funchG, vecXC, bigR, vecX );
	if (doSurfDirTests)
		vecXS = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
		vecGFullSpace = funchG( vecXS );
		vecD = vecXS - vecXC;
		assert( abs(vecD'*vecGSurf) <= 1e-8*norm(vecD)*norm(vecGFullSpace) );
		clear vecXS;
		clear vecGFullSpace;
		clear vecD;
	end
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
		clear epsFD;
		clear vecXP;
		clear vecXM;
		clear omegaP;
		clear omegaM;
	end
	clear n;
	gNorm = norm(vecGSurf);
	gNorm_fd = norm(vecGSurf_fd);
	gRes = norm(vecGSurf-vecGSurf_fd);
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing gSurf..." );
		echo__vecGSurf = vecGSurf
		echo__vecGSurf_fd = vecGSurf_fd
		echo__gNorm = gNorm
		echo__gNorm_fd = gNorm_fd
		echo__gRes = gRes
	end
	assert( gRes <= 1e-4 * ( gNorm + gNorm_fd ) );
	clear vecGSurf;
	clear vecGSurf_fd;
	clear gNorm;
	clear gNorm_fd;
	clear gRes;
	%
	%
	% Test Hessian; results are inexact.
	% We could also compare Hessian via 1st deriv of grad, but, mathematically,
	% that should be redundant.
	matHSurf = calcMinfordCurve__evalHSurf( funchG, funchH, vecXC, bigR, vecX );
	if (doSurfDirTests)
		vecXS = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
		matHFullSpace = funchH( vecXS );
		vecD = vecXS - vecXC;
		assert( abs(vecD'*matHSurf*vecD) <= 1e-8*norm(vecD)^2*sqrt(sum(sum(matHFullSpace.^2))) );
		clear vecXS;
		clear matHFullSpace;
		clear vecD;
	end
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
		clear epsFD;
		clear vecXPP;
		clear vecXMP;
		clear vecXPM;
		clear vecXMM;
		clear omegaPP;
		clear omegaMP;
		clear omegaPM;
		clear omegaMM;
	end
	end
	clear m;
	clear n;
	hessNorm = sqrt(sum(sum((matHSurf).^2)));
	hessNorm_fd2 = sqrt(sum(sum((matHSurf_fd2).^2)));
	hessRes = sqrt(sum(sum((matHSurf_fd2-matHSurf).^2)));
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing hSurf..." );
		echo__matHSurf = matHSurf
		echo__matHSurf_fd2 = matHSurf_fd2
		echo__hessNorm = hessNorm
		echo__hessNorm_fd2 = hessNorm_fd2
		echo__hessRes = hessRes
	end
	assert( hessRes <= 1e-3*( hessNorm + hessNorm_fd2 ) );
	clear matHSurf;
	clear matHSurf_fd2;
	clear hessNorm;
	clear hessNorm_fd2;
	clear hessRes;
	%
	%
	% Now, with scaling...
	matS = abs(randn(sizeX,sizeX));
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing with scaling..." );
		echo__matS = matS
	end
	%
	%
	% We can also test R/matS...
	% First, let's repeat the unscaled tests, but with scaling...
	%
	% Test that xSurf moves the point, but only once.
	vecXSurf = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS );
	vecXSurf2 = calcMinfordCurve__evalXSurf( vecXC, bigR, vecXSurf, matS );
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing xSurf..." );
		echo__vecXSurf = vecXSurf
		echo__vecXSurf2 = vecXSurf2
	end
	if ( abs(vecXSurf-vecX) < eps*bigR )
		msg( thisFile, __LINE__, "__evalXSurf with scaling did not meaningfully move the point." );
		msg( thisFile, __LINE__, "This could happen by chance but is unlikely(?)." );
		msg( thisFile, __LINE__, "Aborting test." );
		error( "Scaled vecXSurf is on top of vecX." );
	end
	assert( norm(vecXSurf-vecXSurf2)<=1e-4*(norm(vecXSurf)+norm(vecXSurf2)) );
	clear vecXSurf;
	clear vecXSurf2;
	%
	%
	% Test omega.
	omegaSurf = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX, matS );
	omegaSurf_alt = funchOmega( calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS ) );
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing scaled omegaSurf..." );
		echo__omega = funchOmega( vecX );
		echo__omegaSurf = omegaSurf
		echo__omegaSurf_alt = omegaSurf_alt
	end
	assert( norm(omegaSurf-omegaSurf_alt)<=1e-4*(norm(omegaSurf)+norm(omegaSurf_alt)) );
	clear omegaSurf;
	clear omegaSurf_alt;
	%
	%
	% Test gradient; note that result is inexact.
	vecGSurf = calcMinfordCurve__evalGSurf( funchG, vecXC, bigR, vecX, matS );
	if (doSurfDirTests)
		vecXS = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
		vecGFullSpace = funchG( vecXS );
		vecD = vecXS - vecXC;
		assert( abs(vecD'*vecGSurf) <= 1e-8*norm(vecD)*norm(vecGFullSpace) );
		clear vecXS;
		clear vecGFullSpace;
		clear vecD;
	end
	vecGSurf_fd = zeros(sizeX,1);
	for n=1:sizeX
		epsFD = 1e-4;
		vecXP = vecX;
		vecXM = vecX;
		vecXP(n) += epsFD;
		vecXM(n) -= epsFD;
		omegaP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXP, matS );
		omegaM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXM, matS );
		vecGSurf_fd(n) = ( omegaP - omegaM ) / (2.0*epsFD);
		clear epsFD;
		clear vecXP;
		clear vecXM;
		clear omegaP;
		clear omegaM;
	end
	clear n;
	gNorm = norm(vecGSurf);
	gNorm_fd = norm(vecGSurf_fd);
	gRes = norm(vecGSurf-vecGSurf_fd);
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing scaled gSurf..." );
		echo__vecGSurf = vecGSurf
		echo__vecGSurf_fd = vecGSurf_fd
		echo__gNorm = gNorm
		echo__gNorm_fd = gNorm_fd
		echo__gRes = gRes
	end
	assert( gRes <= 1e-4 * ( gNorm + gNorm_fd ) );
	clear vecGSurf;
	clear vecGSurf_fd;
	clear gNorm;
	clear gNorm_fd;
	clear gRes;
	%
	%
	%
	% Test Hessian; results are inexact.
	% We could also compare Hessian via 1st deriv of grad, but, mathematically,
	% that should be redundant.
	matHSurf = calcMinfordCurve__evalHSurf( funchG, funchH, vecXC, bigR, vecX, matS );
	if (doSurfDirTests)
		vecXS = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
		matHFullSpace = funchH( vecXS );
		vecD = vecXS - vecXC;
		assert( abs(vecD'*matHSurf*vecD) <= 1e-8*norm(vecD)^2*sqrt(sum(sum(matHFullSpace.^2))) );
		clear vecXS;
		clear matHFullSpace;
		clear vecD;
	end
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
		omegaPP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPP, matS );
		omegaMP = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMP, matS );
		omegaPM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXPM, matS );
		omegaMM = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecXMM, matS );
		matHSurf_fd2(m,n) = ( omegaPP + omegaMM - omegaMP - omegaPM ) / (4.0*(epsFD)^2);
		clear epsFD;
		clear vecXPP;
		clear vecXMP;
		clear vecXPM;
		clear vecXMM;
		clear omegaPP;
		clear omegaMP;
		clear omegaPM;
		clear omegaMM;
	end
	end
	clear m;
	clear n;
	hessNorm = sqrt(sum(sum((matHSurf).^2)));
	hessNorm_fd2 = sqrt(sum(sum((matHSurf_fd2).^2)));
	hessRes = sqrt(sum(sum((matHSurf_fd2-matHSurf).^2)));
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing scaled hSurf..." );
		echo__matHSurf = matHSurf
		echo__matHSurf_fd2 = matHSurf_fd2
		echo__hessNorm = hessNorm
		echo__hessNorm_fd2 = hessNorm_fd2
		echo__hessRes = hessRes
	end
	assert( hessRes <= 1e-3*( hessNorm + hessNorm_fd2 ) );
	clear matHSurf;
	clear matHSurf_fd2;
	clear hessNorm;
	clear hessNorm_fd2;
	clear hessRes;	
	%
	%
	%
	% Do some conceptual checks...
	%
	% Test that using matS = I produces same result as no matS.
	% We could do test of this form on omega, g, and h, but,
	%  mathematically, that's (probably(?)) unnecessary.
	% Also, test that xSurf moves the point, and differently than when S = I, but only once.
	vecXSurf_noS = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX );
	vecXSurf_eye = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, eye(sizeX,sizeX) );
	vecXSurf = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS );
	modFactor = abs(randn());
	vecXSurf_mod = calcMinfordCurve__evalXSurf( vecXC, bigR*modFactor, vecX, matS*modFactor );
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Testing xSurf with and without scaling..." );
		echo__vecX = vecX
		echo__vecXSurf_noS = vecXSurf_noS
		echo__vecXSurf_eye = vecXSurf_eye
		echo__vecXSurf = vecXSurf
		echo__modFactor = modFactor
		echo__vecXSurf_mod = vecXSurf_mod
	end
	if ( abs(vecXSurf-vecXSurf_noS) < eps*bigR )
		msg( thisFile, __LINE__, "Scaling did not meaningfully change where the point moves." );
		msg( thisFile, __LINE__, "This could happen by chance but is unlikely." );
		msg( thisFile, __LINE__, "Aborting test." );
		error( "Scaled vecXSurf is on top of unscaled vecXSurf." );
	end
	assert( norm(vecXSurf_noS-vecXSurf_eye)<=1e-4*(norm(vecXSurf_noS)+norm(vecXSurf_eye)) );
	assert( norm(vecXSurf-vecXSurf_mod)<=1e-4*(norm(vecXSurf)+norm(vecXSurf_mod)) );
	clear vecXSurf_noS;
	clear vecXSurf_eye;
	clear vecXSurf;
	clear modFactor;
	clear vecXSurf_mod;
	%
	%
	%whos
	if ( verbLev >= VERBLEV__PROGRESS )
		msg( thisFile, __LINE__, "Passed all test calculations for this case." );
	end
return;
end
