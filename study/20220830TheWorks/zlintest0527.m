clear;
sizeX = 60;
sizeF = 50;
sizeV = 11;
sizeVLocal = 5;
%
matV = utorthdrop(randn(sizeX,sizeV));
matW0 = randn(sizeF,sizeV);
matA0 = matW0'*matW0;
matVLocal = utorthdrop(matV*randn(sizeV,sizeVLocal));
matWLocal = randn(sizeF,sizeVLocal);
%
matW = matW0;
matA = matA0;
for n=1:sizeVLocal
	vecV = matVLocal(:,n);
	vecW = matWLocal(:,n);
	vecU = matV'*vecV;
	matEU = eye(sizeV,sizeV) - vecU*(vecU');
	matW = matW + ( vecW - matW*vecU) * (vecU');
	matA = matEU'*matA*matEU;
	matA = (matA'+matA)/2.0;
endfor
%
matU = matV'*matVLocal;
matEU = eye(sizeV,sizeV) - matU*(matU');
matW_alt = matW0 + ( matWLocal - matW0*matU ) * (matU');
matA_alt = matEU' * matA0 * matEU;
%
assert( reldiff(matW,matW_alt) < sqrt(eps) );
assert( reldiff(matA,matA_alt) < sqrt(eps) );
