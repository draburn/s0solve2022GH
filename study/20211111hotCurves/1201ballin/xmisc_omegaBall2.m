% Run x1 first, then run xmisc_omegaBall.
%
thisFile = "xmisc_omegaBall2";
funchOmega = @(x)( 0.5*sum(funchF(x).^2,1) );
funchH = @(x)( testFunc_evalH( x, testFuncPrm ) );
%
msg( thisFile, __LINE__, "Let's do a comparison..." );
index1 = 5;
index2 = 5;
%vecX = [ x1Mesh(index1,index2); x2Mesh(index1,index2) ]
normR = norm(vecX-vecXC)
%assert( normR > bigR )
%
omegaBall_meshCalcu = omegaBallMesh(index1,index2)
omegaBall_newModule = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecX )
%
%
vecGBall_dOmega = zeros(sizeX,1);
for n=1:sizeX
	epsFD = 1e-4;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	omegaP = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXP );
	omegaM = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXM );
	vecGBall_dOmega(n) = ( omegaP - omegaM ) / (2.0*epsFD);
end
vecGBall_dOmega = vecGBall_dOmega
vecGBall_newMod = omegaBall_evalGrad( funchG, vecXC, bigR, bigC, vecX )
%
%
matHBall_newModule = omegaBall_evalHess( funchG, funchH, vecXC, bigR, bigC, vecX )
%
matHBall_dg = zeros(sizeX,sizeX);
for n=1:sizeX
	epsFD = 1e-4;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	gradP = omegaBall_evalGrad( funchG, vecXC, bigR, bigC, vecXP );
	gradM = omegaBall_evalGrad( funchG, vecXC, bigR, bigC, vecXM );
	matHBall_dg(:,n) = ( gradP - gradM ) / (2.0*epsFD);
end
matHBall_dg = matHBall_dg
%
omega0 = omegaBall_newModule;
matHBall_d2omega = zeros(sizeX,sizeX);
for n=1:sizeX
	epsFD = 1e-3;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	omegaP = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXP );
	omegaM = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXM );
	matHBall_d2omega(n,n) = ( omegaP + omegaM - 2.0*omega0 ) / (epsFD^2);
end
for n=1:sizeX
for m=1:n-1
	epsFD = 1e-3;
	vecXPP = vecX;
	vecXPM = vecX;
	vecXMP = vecX;
	vecXMM = vecX;
	vecXPP(n) += epsFD; vecXPP(m) += epsFD;
	vecXPM(n) += epsFD; vecXPM(m) -= epsFD;
	vecXMP(n) -= epsFD; vecXMP(m) += epsFD;
	vecXMM(n) -= epsFD; vecXMM(m) -= epsFD;
	omegaPP = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXPP );
	omegaPM = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXPM );
	omegaMP = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXMP );
	omegaMM = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecXMM );
	matHBall_d2omega(n,m) = ( omegaPP + omegaMM - omegaPM - omegaMP ) / (4.0*(epsFD^2));
end
end
for n=1:sizeX
for m=n+1:sizeX
	matHBall_d2omega(n,m) = matHBall_d2omega(m,n);
end
end
matHBall_d2omega = matHBall_d2omega

%vecXStart = vecXC
%vecXStart = vecXC-[0.1;0.1]
%vecXStart = [ 2.5; 0.2 ]
%vecXStart = [ 2.4; 0.1 ]
%vecXStart = [ 2.0; 6.5 ]
vecXStart = [ 1.5; 7.0 ]
vecRStart = vecXStart - vecXC;
vecXStart = vecXC + bigR*vecRStart/norm(vecRStart);
msg( thisFile, __LINE__, "Starting alpha..." );
vecXAlpha = omegaBall_findMin_alpha( funchOmega, funchG, vecXC, bigR, bigC, vecXStart )
msg( thisFile, __LINE__, "Finished alpha." );
msg( thisFile, __LINE__, "Starting beta..." );
vecXBeta = omegaBall_findMin_beta( funchOmega, funchG, vecXC, bigR, vecXStart )
msg( thisFile, __LINE__, "Finished beta." );
myPrm.forceOnSurf = true;
msg( thisFile, __LINE__, "Starting beta force..." );
vecXBetaForce = omegaBall_findMin_beta( funchOmega, funchG, vecXC, bigR, vecXStart, myPrm )
msg( thisFile, __LINE__, "Finished beta force." );
%
bigR = bigR
dBetaForce = norm(vecXBetaForce-vecXC)
foores = dBetaForce - bigR
assert( abs(foores) < 0.01*bigR )
