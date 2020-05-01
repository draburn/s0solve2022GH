myclear;
%setprngstates(0);
thisFile = "blm0430";
%
sizeN = 5;
%
if (0)
	matA = eye(sizeN,sizeN);
	matB = eye(sizeN,sizeN);
	vecX0 = zeros(sizeN,1);
	indexNNonLin = 1;
	vecX1 = zeros(sizeN,1);
	vecX2 = zeros(sizeN,1);
	vecX3 = zeros(sizeN,1);
	vecX1(1) = 1.0;
	vecX2(2) = 1.0;
	vecX3(1) = 2.0;
	vecX1 = randn(sizeN,1);
	vecX2 = randn(sizeN,1);
	vecX3 = randn(sizeN,1);
else
	matA = randn(sizeN,sizeN);
	matB = randn(sizeN,sizeN);
	vecX0 = randn(sizeN,1);
	indexNNonLin = ceil((10*eps) + ((sizeN-(20*eps))*rand) );
	vecX1 = randn(sizeN,1);
	vecX2 = randn(sizeN,1);
	vecX3 = randn(sizeN,1);
end
%
c1 = 1.0;
cx = 0.01;
c0 = (c1^2)*(1.0+cx)/4.0;
funchF = @(dummy_matX)( blm0429_func_alt( dummy_matX, matA, matB, vecX0, indexNNonLin, c1, c0 ) );
%
vecF1 = funchF(vecX1);
vecF2 = funchF(vecX2);
vecF3 = funchF(vecX3);
fdjaco_prm.fdOrder = 2;
fdjaco_prm.epsFD = eps^0.3;
matJ1 = fdjaco( funchF, vecX1, fdjaco_prm );
matJ2 = fdjaco( funchF, vecX2, fdjaco_prm );
toc();
%
% Our model...
%   vecF = vecF0 + (matJ0 * vecX) + (vecAHat * ( p*(z^2) + q*(z^3) ));
%   z = vecBHat' * vecX;
%   matJ = matJ0 + (vecAHat * (vecBHat') * ( 2.0*p*z + 3.0*q*(z^2) ));
% Hence...
matDJ = matJ2 - matJ1
if (0)
	vecATemp = matDJ*randn(sizeN,1);vecAHat = vecATemp/sqrt(sum(vecATemp.^2))
	vecATemp = matDJ*randn(sizeN,1);vecAHat = vecATemp/sqrt(sum(vecATemp.^2))
	vecATemp = matDJ*randn(sizeN,1);vecAHat = vecATemp/sqrt(sum(vecATemp.^2))
	vecATemp = matDJ*randn(sizeN,1);vecAHat = vecATemp/sqrt(sum(vecATemp.^2))
	vecATemp = matDJ*randn(sizeN,1);vecAHat = vecATemp/sqrt(sum(vecATemp.^2))
	return;
end
vecATest = vecX2 - vecX1; % Could be any random vector; I like this one.
vecATemp = matDJ * vecATest;
assert( norm(vecATemp) > 0.0 );
vecBTest = vecATemp; % Could be any random vector; I like this one.
vecBTemp = (matDJ')*vecBTest;
assert( norm(vecBTemp) > 0.0 );
vecAHat = vecATemp / norm(vecATemp)
vecBHat = vecBTemp / norm(vecBTemp)
if (0)
	% Let's check...
	vecATrue = matA(:,indexNNonLin);
	vecBTrue = matB(indexNNonLin,:)';
	trueSign = sign(vecATrue'*vecAHat);
	vecResAHat = vecAHat - trueSign * (vecATrue/norm(vecATrue))
	vecResBHat = vecBHat - trueSign * (vecBTrue/norm(vecBTrue))
	return
end
%
z1 = vecBHat' * vecX1;
z2 = vecBHat' * vecX2;
z3 = vecBHat' * vecX3;
%
adF2 = (vecAHat')*( vecF2 - (vecF1+( matJ1 * (vecX2-vecX1) )) );
adF3 = (vecAHat')*( vecF3 - (vecF1+( matJ1 * (vecX3-vecX1) )) );
matZ = [ ...
  (z2-z1)*(z2-z1), (z2-z1)*( (z2^2) + (z2*z1) - (2*(z1^2)) ); ...
  (z3-z1)*(z3-z1), (z3-z1)*( (z3^2) + (z3*z1) - (2*(z1^2)) ) ];
vecCoeff = matZ \ [ adF2; adF3 ];
p = vecCoeff(1)
q = vecCoeff(2)
%
matJ0 = matJ1 - (vecAHat*(vecBHat')*( (2.0*p*z1) + (3.0*q*(z1^2)) ))
vecF0 = vecF1 - (matJ1*vecX1) + (vecAHat*( (p*(z1^2)) + (2.0*q*(z1^3)) ))
vecF0 = vecF1 - (matJ0*vecX1) - (vecAHat*( (p*(z1^2)) + (q*(z1^3)) ))
%%%%
if (0)
	% Let's check...
	vecResF1 = vecF0 + (matJ0*vecX1) + vecAHat*( p*z1^2 + q*z1^3) - vecF1
	vecResF2 = vecF0 + (matJ0*vecX2) + vecAHat*( p*z2^2 + q*z2^3) - vecF2
	vecResF3 = vecF0 + (matJ0*vecX3) + vecAHat*( p*z3^2 + q*z3^3) - vecF3
	matResJ1 = matJ0 + vecAHat*(vecBHat')*( 2.0*p*z1 + 3.0*q*z1^2 ) - matJ1
	matResJ2 = matJ0 + vecAHat*(vecBHat')*( 2.0*p*z2 + 3.0*q*z2^2 ) - matJ2
	resF1 = sqrt(sum((vecResF1).^2))
	resF2 = sqrt(sum((vecResF2).^2))
	resF3 = sqrt(sum((vecResF3).^2))
	resJ1 = sqrt(sum(sum((matResJ1).^2)))
	resJ2 = sqrt(sum(sum((matResJ2).^2)))
	relresF1 = resF1 / sqrt(sum(vecF1.^2))
	relresF2 = resF2 / sqrt(sum(vecF2.^2))
	relresF3 = resF3 / sqrt(sum(vecF3.^2))
	relresJ1 = resJ1 / sqrt(sum(sum(matJ1.^2)))
	relresJ2 = resJ2 / sqrt(sum(sum(matJ2.^2)))
	return
end
%%%%
%
vecV = matJ0 \ vecAHat;
vecW = matJ0 \ vecF0;
%
btvq = (vecBHat')*vecV*q;
btvp = (vecBHat')*vecV*p;
btw = (vecBHat')*vecW;
rvecZRoots = roots([ btvq, btvp, 1.0, btw ])
rvecIsRealIsh = abs(imag(rvecZRoots)) < sqrt(eps);
assert(1==sum(rvecIsRealIsh));
%
zRoot = real(rvecZRoots(rvecIsRealIsh))
vecXGuess = -vecW - (vecV*( (p*(zRoot^2)) + (q*(zRoot^3)) ))
echo__vecX0 = vecX0
%
rvecResX = vecXGuess - vecX0
resX = sqrt(sum(rvecResX.^2))
