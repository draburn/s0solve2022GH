clear;
thisFile = "test_fminsurf";
commondefs;
setprngstates(0);
numFigs = 0;
startTime = time();
%
sizeX = 2
%
function [ bigL, vecDL ] = funchBigL( vecX, matA, vecXRoot )
	%matA = [ 1, 2; 3, 4 ];
	%vecXRoot = [ 2; 2 ];
	vecF = matA*(vecX-vecXRoot);
	vecDL = matA'*vecF;
	bigL = 0.5*vecF'*vecF;
end
%
function [ vecS, matDST ] = funchVecS( vecX, bigR )
	bigR = 1.0;
	magX = norm(vecX);
	if (0.0==magX)
		vecX(1) = bigR;
		magX = bigR;
	end
	sizeX = size(vecX,1);
	vecS = bigR*vecX/magX;
	matDST = ((bigR/magX)*eye(sizeX,sizeX)) - (bigR/(magX^3))*(vecX*vecX');
end
%
function [ bigP, vecDP ] = funchBigP( vecD )
	bigP = 0.5*(vecD'*vecD);
	vecDP = vecD;
end
%
sizeF = sizeX;
matA = randn(sizeF,sizeX);
vecXRoot = randn(sizeX,1);
bigR = 1.0;
%
vecX0 = randn(sizeX,1);
%
tic();
%vecXF = fminsurf( @funchBigL, @funchVecS, @funchBigP, vecX0 )
funchBigL_at = @(x)( funchBigL( x, matA, vecXRoot ) );
funchVecS_at = @(x)( funchVecS( x, bigR ) );
funchBigP_at = @(x)( funchBigP( x ) );
vecXF = fminsurf( funchBigL_at, funchVecS_at, funchBigP_at, vecX0 );
toc();
%
echo__vecX0 = vecX0
[ vecS, matDST ] = funchVecS_at(vecX0)
[ bigL, vecDL ] = funchBigL_at(vecS)
[ bigP, vecDP ] = funchBigP_at(vecX0-vecS)
%
echo__vecXF = vecXF
[ vecS, matDST ] = funchVecS_at(vecXF)
[ bigL, vecDL ] = funchBigL_at(vecS)
[ bigP, vecDP ] = funchBigP_at(vecXF-vecS)
