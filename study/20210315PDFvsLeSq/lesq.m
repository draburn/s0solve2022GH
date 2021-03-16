% DRabun 2021.03.15
% A simple least-squares fit.

function dat = lesq( pts, numCoeff )
	numPts = size(pts,1);
	numCol = size(pts,2);
	%
	matX = zeros(numPts,numCoeff);
	matX(:,1) = 1.0;
	for m=2:numCoeff
		matX(:,m) = pts(:,1).^(m-1);
	end
	vecY = pts(:,2);
	if ( 3 <= numCol )
		vecW = pts(:,3);
		matW = diag( sqrt(vecW) );
	else
		vecW = ones(numPts,1);
		matW = eye(numPts,numPts);
	end
	%
	matWX = matW * matX;
	vecWY = matW * vecY;
	%
	matWXTWX = matWX' * matWX;
	vecWXTWY = matWX' * vecWY;
	%
	vecC = matWXTWX \ vecWXTWY;
	%
	vecR = (matX*vecC) - vecY;
	vecWR = matW * vecR;
	rho = norm(vecWR);
	%
	%
	vecC_flip = zeros(numCoeff,1);
	for m=1:numCoeff
		vecC_flip(m) = vecC(numCoeff+1-m);
	end
	rootsCmplx = roots(vecC_flip)
	%
	numRoots = 0;
	vecRoots = [];
	for m=1:size(rootsCmplx,1)
	if (isreal(rootsCmplx(m)))
		numRoots++;
		vecRoots(numRoots) = rootsCmplx(m);
	end
	end
	%
	dat.pts = pts;
	dat.numCoeff = numCoeff;
	dat.fitType = "lesq";
	dat.vecC = vecC;
	dat.vecRoots = vecRoots;
	%
return;
end

%!test
%!	setprngstates(0);
%!	pts = abs(randn(3,2))
%!	pts(:,3) = 1.0;
%!	dat1 = lesq( pts, 1 );
%!	dat2 = lesq( pts, 2 );
%!	dat3 = lesq( pts, 3 );
%!	xMin = min(pts(:,1));
%!	xMin = min([ xMin, min(dat1.vecRoots) ]);
%!	xMin = min([ xMin, min(dat2.vecRoots) ]);
%!	xMin = min([ xMin, min(dat3.vecRoots) ]);
%!	xMax = max(pts(:,1));
%!	xMax = max([ xMax, max(dat1.vecRoots) ]);
%!	xMax = max([ xMax, max(dat2.vecRoots) ]);
%!	xMax = max([ xMax, max(dat3.vecRoots) ]);
%!	xVals = linspace( xMin-0.2*(xMax-xMin), xMax+0.2*(xMax-xMin), 1000 );
%!	yVals1 = dat1.vecC(1) + (0.0*xVals);
%!	yVals2 = dat2.vecC(1) + (dat2.vecC(2)*xVals);
%!	yVals3 = dat3.vecC(1) + (dat3.vecC(2)*xVals) + (dat3.vecC(3)*(xVals.^2));
%!	plot( ...
%!	  pts(:,1), pts(:,2), 'ro', ...
%!	  xVals, yVals2, 'b-', ...
%!	  dat2.vecRoots, 0*dat2.vecRoots, 'bx', ...
%!	  xVals, yVals3, 'g-', ...
%!	  dat3.vecRoots, 0*dat3.vecRoots, 'gx', ...
%!	  xVals, 0*xVals, 'k-' );
%!	grid on;
