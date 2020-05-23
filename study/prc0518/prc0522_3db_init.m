myclear;
setprngstates(51618464);
sizeX = 3;
sizeF = 3;
%
funcPrm.vecF0 = randn(sizeF,1);
funcPrm.matJ0 = randn(sizeF,sizeX);
funcPrm.polyOrder = 5;
funcPrm.c2 = randn(sizeF,sizeX,sizeX)/20;
funcPrm.c3 = randn(sizeF,sizeX,sizeX,sizeX)/100;
funcPrm.c4 = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/500;
funcPrm.c5 = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/2000;
%
funcPrm.expPrm.numTerms = round(exp(2+randn));
for n=1:funcPrm.expPrm.numTerms
	funcPrm.expPrm.matC(:,n) = randn(sizeF,1);
	funcPrm.expPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.expPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.expPrm.rvecXScale(1,n) = abs(randn);
end
%
funcPrm.vecF0Denom = ones(sizeF,1);
funcPrm.denomExpPrm.numTerms = round(exp(2+randn));
for n=1:funcPrm.denomExpPrm.numTerms
	funcPrm.denomExpPrm.matC(:,n) = abs(randn(sizeF,1));
	funcPrm.denomExpPrm.matX0(:,n) = randn(sizeX,1);
	funcPrm.denomExpPrm.rvecPowA(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecPowB(1,n) = 1.0+abs(randn);
	funcPrm.denomExpPrm.rvecXScale(1,n) = abs(randn);
end
%
vecXRoot = randn(sizeX,1);
vecFPreRoot = demofunc0518_eval(vecXRoot,funcPrm);
%
funchF = @(x)( demofunc0518_eval(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
funchJ = @(x)( fdjaco(funchF,x) );
vecXGroot0 = zeros(sizeX,1);
if (0)
	grootPrm.funchJ = funchJ;
	[ vecXM, retCode, grootDat ] = groot( funchF, vecXGroot0, grootPrm );
else
	% "Resolve"...
	groot_vecFMHat = [
	  0.520266905658494
	  0.147619556424842
	  0.841148508550938 ];
	groot_matW = eye(sizeF,sizeF) + 10*groot_vecFMHat*(groot_vecFMHat');
	groot_funchF = @(x)( groot_matW * funchF(x) );
	groot_funchJ = @(x)( fdjaco(groot_funchF,x) );
	grootPrm.funchJ = groot_funchJ;
	[ vecXM, retCode, grootDat ] = groot( groot_funchF, vecXGroot0, grootPrm );
end
