blm0509init;
thisFile = "blm0509calc";
%
fdjaco_prm.epsFD = 1E-6;
fdjaco_prm.fdOrder = 2;
funchJ = @(x)( fdjaco( funchF, x, fdjaco_prm ) );
vecXStart = [-0.49;0.01];
%
vecX = vecXStart;
for n=1:50
	vecF = funchF(vecX);
	matJ = funchJ(vecX);
	if ( rcond(matJ'*matJ) < 1e-16 )
		break;
	end
	[matU,matS,matV] = svd(matJ);
	vecDeltaTemp = -(matJ'*matJ)\(matJ'*vecF);
	m = 0;
	while (1)
		vecXTemp = vecX+vecDeltaTemp;
		vecFTemp = funchF(vecXTemp);
		if (norm(vecFTemp)<norm(vecF))
			vecDelta = vecDeltaTemp;
			break;
		end
		vecDeltaTemp*=0.1;
		m++;
		if ( m>1000 )
			break;
		end
	end
	if ( m>1000 )
		break;
	end
	vecX += vecDelta;
	msg( thisFile, __LINE__, sprintf( ...
	 "  %2d;   %10.3e,   %10.3e;   %10.3e,   %10.3e", ...
	 n, norm(vecF), norm(vecDelta),
	 matS(1,1), matS(2,2) ));
end
%
vecXM = vecX
vecFM = vecF
matJM = matJ
[matUM,matSM,matVM] = svd(matJM);
%%[matUM,matSM,matVM] = svd(matJM'*matJM);

matUM(:,1)' * vecFM
% THIS WAS SUPPOSED TO BE NEAR ZERO!
