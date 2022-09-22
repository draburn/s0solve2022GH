clear;
mydefs;
numFigs = 0;
if (1)
	setprngstates(0);
	sizeX = 10
	sizeF = 8
	sizeK = 6
	sizeL = 4
	
elseif (1)
	setprngstates(0);
	sizeX = 5
	sizeF = 4
	sizeK = 3
	sizeL = 2
else
	setprngstates();
	sizeX = 10+round(50*rand)
	sizeF = 10+round(50*rand)
	sizeK = 1+round(min([sizeX,sizeF])*rand)
	sizeL = 1+round(min([sizeX,sizeF])*rand)
endif
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
%matV(:,1) = 0.0;
matV(8,1) = 10.0;
%matV(:,1) = 1.0;
matV = utorthdrop( matV );
matW = matJSecret*matV;
%
%matJSecret
%matV
%sumsq(matV,1)
%matSInv = sumsq(matV,1) - matV.^2 % Autobroadcast.
matS = 1./(sqrt(eps) + abs(sumsq(matV,1) - (matV.^2))); % % Autobroadcast. abs() should be superfluous.
%
%matOMP1 = zeros(sizeF,sizeX);
matEta1 = zeros(sizeF,sizeX);
matChi1 = zeros(sizeF,sizeX);
matEtaS1 = zeros(sizeF,sizeX);
matChiS1 = zeros(sizeF,sizeX);
%matEtaSX1 = zeros(sizeF,sizeX);
%matChiSX1 = zeros(sizeF,sizeX); Moot now that I knwo how to calc matChiS.
for m=1:sizeF
for n=1:sizeX
	sumWV = 0.0;
	sumWW = 0.0;
	sumVV = 0.0;
	sumWVS = 0.0;
	sumWWS = 0.0;
	sumVVS = 0.0;
	sumVVSS = 0.0;
	for k=1:sizeK
		w = matW(m,k);
		v = matV(n,k);
		s = matS(n,k);
		sumWV += w*v;
		sumVV += v*v;
		sumWW += w*w;
		sumWVS += w*v*s;
		sumVVS += v*v*s;
		sumWWS += w*w*s;
		sumVVSS += v*v*s*s;
	endfor
	%matOMP1(m,n) = sumWV^2;
	matEta1(m,n) = sumWV /(eps+sumVV);
	matChi1(m,n) =  (sumWV^2)/((eps+sumVV)*(eps+sumWW));
	matEtaS1(m,n) = sumWVS /(eps+sumVVS);
	matChiS1(m,n) =  (sumWVS^2)/((eps+sumVVS)*(eps+sumWWS));
	%matChiSX1(m,n) = (sumWVS^2)/((eps+sumVVSS)*(eps+sumWW));
endfor
endfor
%matOMP1
matEta1Sq = matEta1.^2
matChi1
%matEtaS1
matChiS1
%matChiSX1
%
matEta2 = (matW * (matV')) * diag(1.0./(eps+sumsq(matV,2)));
matChi2 = diag(1.0./(eps+sumsq(matW,2))) * ((matW * (matV')).^ 2) * diag(1.0./(eps+sumsq(matV,2)));
%%%matChiS2 =  ( ((matW * ((matV.*matS)')).^ 2) ./ ( eps + ((matW.^2)*(matS')) );

matFoo_wvs = matW * (( matV.*matS )');
matFoo_vvsInv = diag( 1.0 ./ (eps+sum( (matV.^2).*matS, 2 )) );
matFoo_wws = (matW.^2) * (matS');
matChiS2 = ((matFoo_wvs.^2)*matFoo_vvsInv) ./ (eps+matFoo_wws);

%matChiSX2 = diag(1.0./(eps+sumsq(matW,2))) * ((matW * ((matV.*matS)')).^ 2) * diag(1.0./(eps+sumsq(matV.*matS,2)));
%
rdEta = reldiff( matEta1, matEta2 )
rdChi = reldiff( matChi1, matChi2 )
rdChiS = reldiff( matChiS1, matChiS2 )
%rdChiSX = reldiff( matChiSX1, matChiSX2 )
%rdSX = reldiff( matChi1, matChiSX2 )
%rdX = reldiff( matChiS1, matChiSX2 )

numFigs++; figure(numFigs);
matEta1Scl = diag(1./max(matEta1'.^2))*(matEta1.^2);
image(64*matEta1Scl);
colormap(hot(64));

numFigs++; figure(numFigs);
image(64*matChi1);
colormap(hot(64));
%
numFigs++; figure(numFigs);
image(64*matChiS1);
colormap(hot(64));
%
%numFigs++; figure(numFigs);
%imagesc(matChiSX1);

%matJSSq = matJSecret.^2;
%matJSSqScl = diag(1./max(matJSSq'))*matJSSq
%
numFigs++; figure(numFigs);
%image(64*(matJSSqScl.^0.25));
image(64*double(matJSecret~=0.0));
colormap(hot(64));
