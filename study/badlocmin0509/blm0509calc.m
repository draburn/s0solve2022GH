blm0509init;
thisFile = "blm0509calc";
%
fdjaco_prm.epsFD = 1E-6;
fdjaco_prm.fdOrder = 2;
funchJ = @(x)( fdjaco( funchF, x, fdjaco_prm ) );
%
%vecXStart = [-0.49;0.01];
vecXStart = [ -0.48684527346537659; 0.00973760207499515 ];
%
%vecXStart = [-0.48;-0.45];
%vecXStart = [-0.5008417558849074;-0.0424983458406762];
%
vecX = vecXStart;
for n=1:100
	vecF = funchF(vecX);
	matJ = funchJ(vecX);
	vecG = matJ'*vecF;
	vecG *= norm(vecG)/(norm(matJ*vecG));
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
		vecDeltaTemp*=0.5;
		m++;
		if ( m>100 )
			break;
		end
	end
	if ( m>100 )
		break;
	end
	vecX += vecDelta;
	%if (0==mod(n,100))
	msg( thisFile, __LINE__, sprintf( ...
	 "  %2d;   %10.3e,   %10.3e,   %10.3e", ...
	 n, norm(vecF), norm(vecG), norm(vecDelta) ));
	%end
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

matCM = matUM(:,1); %Columnspace of Jm.
matVMNonNull = matVM(:,1); %Non-null space of Jm.
vecVMHatNull = matVM(:,2); %Null space of Jm, 1D is assumed.
yPrev = 0.0;
for k=1:300
	z = 0.001*k;
	funchXCurve = @(y,z)( repmat(vecXM,[1,size(y,2)]) ...
	 + (matVMNonNull*y) + (vecVMHatNull*z) );
	funchFSub = @(y)( (matCM') * funchF( funchXCurve(y,z) ));
	y = yPrev;
	n = 0;
	while (1)
		epsP = 1e-6;
		epsM = 1e-6;
		f = funchFSub(y);
		if (norm(f)<1e-14)
			break;
		end
		fp = funchFSub(y+epsP);
		fm = funchFSub(y-epsM);
		dfdy = (fp-fm)/(eps+epsP+epsM);
		dy = -(dfdy)\f;
		m = 0;
		while (1)
			yTemp = y+dy;
			fTemp = funchFSub(yTemp);
			if (norm(fTemp)<norm(f))
				break;
			end
			m++;
			if (m>50)
				%msg(thisFile,__LINE__,"Failed to converge.");
				break;
			end
			dy *= 0.1;
		end
		if (m>50)
			break;
		end
		n++;
		if (n>100)
			break;
		end
		y+=dy;
	end
	rvecZ(k) = z;
	rvecY(k) = y;
	matX(:,k) = funchXCurve(y,z);
	matF(:,k) = funchF(matX(:,k));
	matFSub(:,k) = funchFSub(y);
end
