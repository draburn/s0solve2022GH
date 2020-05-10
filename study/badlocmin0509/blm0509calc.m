blm0509init;
thisFile = "blm0509calc";
%
fdjaco_prm.epsFD = 1E-6;
fdjaco_prm.fdOrder = 2;
funchJ = @(x)( fdjaco( funchF, x, fdjaco_prm ) );
%vecXStart = [-0.49;0.01];
%vecXStart = [-0.48;-0.45];
vecXStart = [-0.5008417558849074;-0.0424983458406762];
%
vecX = vecXStart;
for n=1:50
	vecF = funchF(vecX);
	matJ = funchJ(vecX);
	vecG = matJ'*vecF;
	%if ( rcond(matJ'*matJ) < 1e-16 )
	%	break;
	%end
	%[matU,matS,matV] = svd(matJ);
	%vecDeltaTemp = -(matJ'*matJ)\(matJ'*vecF);
	vecDeltaTemp = -vecG;
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
	 "  %2d;   %10.3e,   %10.3e,   %10.3e", ...
	 n, norm(vecF), norm(vecG), norm(vecDelta) ));
end
%
vecXM = vecX
vecFM = vecF
vecGM = vecG
matJM = matJ
[matUM,matSM,matVM] = svd(matJM);
%%[matUM,matSM,matVM] = svd(matJM'*matJM);

norm(matUM(:,1)' * vecFM)/norm(vecFM)
% THIS SHOULD BE VERY SMALL.

if (0)
	rvecLambda = linspace(-10.0,10.0,1001);
	%rvecLambda = linspace(-1.0,1.0,10001);
	%rvecLambda = linspace(-0.475,-0.473,100000);
	matX = matVM(:,2)*rvecLambda + repmat(vecXM,size(rvecLambda));
	matF = funchF(matX);
	semilogy(rvecLambda,sqrt(sum(matF.^2,1)),'o-');
	%semilogy(rvecLambda,sqrt(sum((matUM(:,2)' * matF).^2,1)),'o-');
	grid on;
	% THIS WILL ONLY SHOW A STEEP DIP IF THE
	% "ROOT-ISH CURVE" IS A LINE.
end
