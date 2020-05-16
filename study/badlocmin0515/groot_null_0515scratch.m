myclear;
setprngstates();
%
%matJ = [1,1,1;1,0,0;0,1,0];
matJ = randn(3,3);
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
vecX = randn(3,1);
%funchF = @(x)( [1;1;1] + matJ*x + [1;-1;0]*((vecVNull'*x).^2) + [0;1;1]*((vecVNull'*x).^3) );
funchF = @(x)( [1;1;1] + matJ*x + [1;-1;0]*((vecVNull'*(x-vecX)).^2) ...
  + [0;1;1]*((vecVNull'*(x-vecX)).^3) );
%
vecF = funchF(vecX);
%
d = 1.0;
vecXP = vecX + (d*vecVNull);
vecXM = vecX - (d*vecVNull);
vecFP = funchF(vecXP);
vecFM = funchF(vecXM);
vecA = ( vecFP + vecFM - (2.0*vecF) ) / (2*d*d)
vecB = ( vecFP - vecFM - (2.0*d*sNull*vecUNull) ) / (2*d*d*d)
funchFModel = @(x)( vecF + matJ*(x-vecX) + ...
  (vecA*( ((vecVNull') * (x-vecX)).^2) ) + ...
  (vecB*( ((vecVNull') * (x-vecX)).^3) ) );
assert( norm(funchFModel(vecX)-vecF) <= eps^0.5 * norm(vecF) );
assert( norm(funchFModel(vecXP)-vecFP) <= eps^0.5 * norm(vecFP) );
assert( norm(funchFModel(vecXM)-vecFM) <= eps^0.5 * norm(vecFM) );
