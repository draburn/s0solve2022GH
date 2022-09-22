clear;
mydefs;
setprngstates();
sizeX = 10+round(50*rand)
sizeF = 10+round(50*rand)
sizeK = 1+round(min([sizeX,sizeF])*rand)
sizeL = 1+round(min([sizeX,sizeF])*rand)
%
matV = randn(sizeX,sizeK);
matJSecret = zeros(sizeF,sizeX);
for m=1:sizeF
for l=1:sizeL
	n = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	nzeListSecret(m) = n;
	matJSecret(m,n) += randn();
endfor
endfor
%
matV = utorthdrop( matV );
matW = matJSecret*matV;
%
%matJSecret
%
%matOMP1 = zeros(sizeF,sizeX);
matEta1 = zeros(sizeF,sizeX);
matChi1 = zeros(sizeF,sizeX);
for m=1:sizeF
for n=1:sizeX
	sumWV = 0.0;
	sumWW = 0.0;
	sumVV = 0.0;
	for k=1:sizeK
		w = matW(m,k);
		v = matV(n,k);
		sumWV += w*v;
		sumVV += v*v;
		sumWW += w*w;
	endfor
	%matOMP1(m,n) = sumWV^2;
	matEta1(m,n) = sumWV /(eps+sumVV);
	matChi1(m,n) =  (sumWV^2)/((eps+sumVV)*(eps+sumWW));
endfor
endfor
%matOMP1
%matEta1
%matChi1
%
matEta2 = (matW * (matV')) * diag(1.0./(eps+sumsq(matV,2)));
matChi2 = diag(1.0./(eps+sumsq(matW,2))) * ((matW * (matV')).^ 2) * diag(1.0./(eps+sumsq(matV,2)));
%
rdEta = reldiff( matEta1, matEta2 )
rdChi = reldiff( matChi1, matChi2 )
