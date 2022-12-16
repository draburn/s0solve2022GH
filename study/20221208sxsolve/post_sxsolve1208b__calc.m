function datOut = post_sxsolve1208b__calc( vecX, matV, f, vecG, matH, matB=[], prm=[] )
	datOut = [];
	datOut.lev = __calc_lev( vecX, matV, f, vecG, matH, [], prm );
	datOut.levScl = __calc_lev( vecX, matV, f, vecG, matH, matB, prm );
return;
endfunction
function datOut = __calc_lev( vecX, matV, f, vecG, matH, matB, prm )
	datOut = [];
	sizeX = size(vecX,1);
	if (isempty(matB))
		matB = eye(size(matH));
	endif
	matBInv = inv(matB);
	%
	vecGScl = matBInv * vecG;
	matHScl = matBInv * matH  * matBInv;
	clear vecG;
	clear matH;
	hScl = sqrt(sum(sum(matHScl.^2)));
	matI = eye(size(matHScl));
	matHI = hScl*matI;
	matA = matHScl - matHI;
	%
	rvecS_src = linspace( 0.0, 1.0, 101 );
	rvecS_src = 1.0 - (1.0-rvecS_src.^1).^3;
	%
	datOut.rvecS = 0.0;
	datOut.matZ = zeros(sizeX,1);
	datOut.rvecFModel = f;
	%
	for n = 2 : length(rvecS_src);
		s = rvecS_src(n);
		[ matR, cholFlag ] = chol( matHI + s * matA );
		if ( 0 ~= cholFlag )
			break;
		endif
		vecZ = matR \ ( matR' \ (-s*vecGScl) );
		fModel = f + vecGScl'*vecZ + (vecZ'*matHScl*vecZ)/2.0;
		datOut.rvecS = [ datOut.rvecS, s ];
		datOut.matZ = [ datOut.matZ, vecZ ];
		datOut.rvecFModel = [ datOut.rvecFModel, fModel ];
	endfor
	%
	datOut.matDelta = matV * matBInv * datOut.matZ;
	datOut.rvecDeltaNorm = sqrt(sum( datOut.matDelta.^2, 1 ));
	datOut.rvecBDeltaNorm = sqrt(sum( (matB*datOut.matDelta).^2 , 1 ));
	datOut.matX = vecX + datOut.matDelta;
	[ datOut.rvecF, datOut.matG ] = prm.funchFG( datOut.matX );
return;
endfunction
