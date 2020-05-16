myclear;
tic();
matA = eye(2,2);
matB = eye(2,2);
vecX0 = zeros(2,1);
indexNNonLin = 1;
c1 = 1.0;
cx = 0.1;
c0 = (c1^2)*(1.0+cx)/4.0;
funchF = @(x) blm0429_func( x, matA, matB, vecX0, indexNNonLin, c1, c0 );
funchFNorm = @(x)sqrt(sum(funchF([x]).^2,1));
%
figIndex++; figure(figIndex);
funchVizXY = @(x,y)asinh(funchFNorm([x;y]));
contourfunch( funchVizXY, [-0.75,0.15,-0.05,0.05] );
colormap(0.6+0.4*jet);
grid on;
%
figIndex++; figure(figIndex);
%rvecX1 = linspace(-0.75,0.15,1001);
rvecX1 = linspace(-0.474,-0.472,1001);
[ minVal, indexOfMin ] = min(funchFNorm([1;0]*rvecX1))
%plot( rvecX1, funchFNorm([1;0]*rvecX1), 'o-' );
%grid on;
%
vecX = [rvecX1(indexOfMin);0];
%
fdjaco_prm.epsFD = 1E-6;
fdjaco_prm.fdOrder = 2;
matJ = fdjaco( funchF, vecX, fdjaco_prm );
vecF = funchF(vecX);
%
[matU,matS,matV] = svd(matJ);
%
matUNN = matU(:,1:end-1);
matSNN = matS(1:end-1,1:end-1);
matVNN = matV(:,1:end-1);
matJNN = matUNN*matSNN*(matVNN');
%
vecVNull = matV(:,end);
vecUNull = matU(:,end);
sNull = matS(end,end);
%
d = 100;
vecXP = vecX + (d*vecVNull);
vecXM = vecX - (d*vecVNull);
vecFP = funchF(vecXP);
vecFM = funchF(vecXM);
vecA = ( vecFP + vecFM - (2.0*vecF) ) / (2*d*d);
vecB = ( vecFP - vecFM - (2.0*d*sNull*vecUNull) ) / (2*d*d*d);
funchFModel = @(x)( vecF + matJ*(x-vecX) + ...
  (vecA*( ((vecVNull') * (x-vecX)).^2) ) + ...
  (vecB*( ((vecVNull') * (x-vecX)).^3) ) );
assert( norm(funchFModel(vecX)-vecF) <= eps^0.5 * norm(vecF) );
assert( norm(funchFModel(vecXP)-vecFP) <= eps^0.5 * norm(vecFP) );
assert( norm(funchFModel(vecXM)-vecFM) <= eps^0.5 * norm(vecFM) );
%
numVals = 1001;
rvecLambda = linspace(-1E0,1E0,numVals);
%rvecLambda = linspace(-0.4728,-0.47276,numVals);
ff = vecF'*vecF
fu = vecF'*vecUNull*sNull
fa = vecF'*vecA
fb = vecF'*vecB
for n=1:numVals
	% Take vecDelta = vecVNull * lambda + matVNN * vecR,
	% find vecR which minimizes ||funchFModel||.
	lambda = rvecLambda(n);
	vecResNL = vecF + (sNull*vecUNull*lambda) + (vecA*(lambda^2)) + (vecB*(lambda^3));
	vecR = - (matUNN*matSNN) \ vecResNL;
	matX(:,n) = vecX + (vecVNull * lambda) + (matVNN*vecR);
	matFModel(:,n) = funchFModel(matX(:,n));
	matF(:,n) = funchF(matX(:,n));
	matResNL(:,n) = vecResNL;
end
toc
%
%figIndex++; figure(figIndex);
%rvecLambda -= 5.19548109;
semilogy( ...
  rvecLambda, sqrt(sum(matFModel.^2,1)), 'o-', ...
  rvecLambda, sqrt(sum(matF.^2,1)), 'x-', ...
  rvecLambda, sqrt(sum(matResNL.^2,1)), '^-', ...
  rvecLambda, abs(vecF'*matResNL), 'v-' );
grid on;
%
[ minVal, indexOfMin ] = min(sqrt(sum(matFModel.^2,1)))

return;