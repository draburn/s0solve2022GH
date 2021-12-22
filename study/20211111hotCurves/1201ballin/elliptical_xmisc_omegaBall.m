thisFile = "elliptical_xmisc_omegaBall";
numFigs = 10;
%
vecXC = [ 0.0; -0.5 ]
%matS = [ 1.0, -1.0; -1.0, 3.0 ];
setprngstates();
matS = diag(1.0+abs(randn(2,1)));
ax = [ vecXC(1), vecXC(1), vecXC(2), vecXC(2) ] ...
  + 1.5*sqrt(sum(sum(matS.^2)))*[-1,1,-1,1];
xpost; thisFile = "elliptical_xmisc_omegaBall";
%
r1Mesh = x1Mesh - vecXC(1);
r2Mesh = x2Mesh - vecXC(2);
rNormMesh = sqrt( r1Mesh.^2 + r2Mesh.^2 );
rHat1Mesh = r1Mesh./(eps*max(max(rNormMesh))+rNormMesh);
rHat2Mesh = r2Mesh./(eps*max(max(rNormMesh))+rNormMesh);
%
yE1Mesh = vecXC(1) + matS(1,1)*rHat1Mesh + matS(1,2)*rHat2Mesh;
yE2Mesh = vecXC(2) + matS(2,1)*rHat1Mesh + matS(2,2)*rHat2Mesh;
matYE = [ reshape(yE1Mesh,1,[]); reshape(yE2Mesh,1,[]) ];
matFE = funchF(matYE,testFuncPrm);
fE1Mesh = reshape(matFE(1,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
fE2Mesh = reshape(matFE(2,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
omegaEMesh = 0.5*( fE1Mesh.^2 + fE2Mesh.^2 );
%
numEllipPts = 1003;
vecTheta = 2.0*pi*linspace(0,1.0,numEllipPts);
matEllipPts = zeros(2,numEllipPts);
matEllipPts(1,:) = vecXC(1) + matS(1,1)*cos(vecTheta) + matS(1,2)*sin(vecTheta);
matEllipPts(2,:) = vecXC(2) + matS(2,1)*cos(vecTheta) + matS(2,2)*sin(vecTheta);
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh), contourPlot_numLevels );
axis equal;
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( vecXC(1), vecXC(2), 'ko', 'linewidth', 3, 'markersize', 10 );
plot( vecXC(1), vecXC(2), 'k+', 'linewidth', 3, 'markersize', 10 );
plot( matEllipPts(1,:), matEllipPts(2,:), 'k-', 'linewidth', 1, 'markersize', 1 );
grid on;
title( "elliptical: omega vs x1, x2" );
xlabel( "x1" );
ylabel( "x2" );
hold off;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaEMesh), contourPlot_numLevels );
axis equal;
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( vecXC(1), vecXC(2), 'ko', 'linewidth', 3, 'markersize', 10 );
plot( vecXC(1), vecXC(2), 'k+', 'linewidth', 3, 'markersize', 10 );
plot( matEllipPts(1,:), matEllipPts(2,:), 'k-', 'linewidth', 1, 'markersize', 1 );
grid on;
title( "elliptical: omegaE vs x1, x2" );
xlabel( "x1" );
ylabel( "x2" );
hold off;
%
%
%
msg( thisFile, __LINE__, "Let's do a comparison..." );
index1 = 11
index2 = 20
vecX = [ x1Mesh(index1,index2); x2Mesh(index1,index2) ]
%
mio_vecR = vecX - vecXC;
mio_rNorm = norm(mio_vecR);
mio_vecXMod = vecXC + matS*mio_vecR/mio_rNorm
%
omegaE_meshCalc = omegaEMesh(index1,index2)
omegaE_analytic = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecX )
%
%
%
vecGE_dOmega = zeros(sizeX,1);
for n=1:sizeX
	epsFD = 1e-4;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	omegaP = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXP );
	omegaM = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXM );
	vecGE_dOmega(n) = ( omegaP - omegaM ) / (2.0*epsFD);
end
vecGE_dOmega = vecGE_dOmega
vecGE_analyt = elliptical_omegaBall_evalGrad( funchG, vecXC, matS, vecX )
%
%
funchH = @(x)( testFunc_evalH( x, testFuncPrm ) );
matHE_analytic = elliptical_omegaBall_evalHess( funchG, funchH, vecXC, matS, vecX )
%
matHE_dg = zeros(sizeX,sizeX);
for n=1:sizeX
	epsFD = 1e-4;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	gradP = elliptical_omegaBall_evalGrad( funchG, vecXC, matS, vecXP );
	gradM = elliptical_omegaBall_evalGrad( funchG, vecXC, matS, vecXM );
	matHE_dg(:,n) = ( gradP - gradM ) / (2.0*epsFD);
end
matHE_dg = matHE_dg
%
omega0 = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecX );
matHE_d2omega = zeros(sizeX,sizeX);
for n=1:sizeX
	epsFD = 1e-3;
	vecXP = vecX;
	vecXM = vecX;
	vecXP(n) += epsFD;
	vecXM(n) -= epsFD;
	omegaP = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXP );
	omegaM = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXM );
	matHE_d2omega(n,n) = ( omegaP + omegaM - 2.0*omega0 ) / (epsFD^2);
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
	omegaPP = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXPP );
	omegaPM = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXPM );
	omegaMP = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXMP );
	omegaMM = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecXMM );
	matHE_d2omega(n,m) = ( omegaPP + omegaMM - omegaPM - omegaMP ) / (4.0*(epsFD^2));
end
end
for n=1:sizeX
for m=n+1:sizeX
	matHE_d2omega(n,m) = matHE_d2omega(m,n);
end
end
matHE_d2omega = matHE_d2omega