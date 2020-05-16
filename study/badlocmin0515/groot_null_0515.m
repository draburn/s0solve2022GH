if (0)
myclear;
test_groot;
end
tic();
figIndex = 0;

matJ = grootDatOut.iterDat(end).matJ;
vecX = grootDatOut.iterDat(end).vecX;
vecF = funchF(vecX);
vecXRoot = grootDatOut.vecX0;
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
numVals = 101;
rvecLambda = linspace(-1E1,1E1,numVals);
%rvecLambda = linspace(4.8,5.4,numVals);
%rvecLambda = linspace(5.1954810895,5.1954810905,numVals);
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