function groot1d_extFitViz( xVals_sorted, fVals_sorted, prm )
	thisFile = "groot1d_extFitViz";
	%
	gVals_sorted = abs(fVals_sorted);
	[ gMin, indexOfGMin ] = min( gVals_sorted );
	numPts = size(xVals_sorted,2);
	%
	if ( gMin < 1e-6 )
		msg( thisFile, __LINE__, sprintf( "numPts = %d", numPts ) );
		msg( thisFile, __LINE__, sprintf( "indexOfGMin = %d", indexOfGMin ) );
		for stencilSize = 5:5
			%for nLo = indexOfGMin+1-stencilSize : indexOfGMin
			%for nLo = indexOfGMin+2-stencilSize : indexOfGMin-1
			for nLo = indexOfGMin+3-stencilSize : indexOfGMin-2
				nHi = nLo + stencilSize - 1;
				msg( thisFile, __LINE__, sprintf( "Considering stencil %d, %d ~ %d.", ...
				  stencilSize, nLo, nHi ) );
				if (  (1 <= nLo)  &&  (numPts >= nHi) )
					prm_extFit_viz.numFigs = 10000 + 1000*numPts + 100*stencilSize + nLo-(indexOfGMin+1-stencilSize);
					msg( thisFile, __LINE__, sprintf( "  Making figure %d.", prm_extFit_viz.numFigs+1 ) );
					xLo = 2.14;%-0.892; %xVals_sorted(nLo);
					xHi = 2.15;%-0.890; %xVals_sorted(nHi);
					pLo = 1.5;
					pHi = 5.0;
					extFit_viz( ...
					  xVals_sorted(nLo:nHi), ...
					  gVals_sorted(nLo:nHi), ...
					  xLo, ...
					  xHi, ...
					  pLo, ...
					  pHi, ...
					  prm_extFit_viz );
					clear xHi;
					clear xLo;
				else
					msg( thisFile, __LINE__, "  Nope." );
				end
			end % nLo loop.
		end % Stencil loop
	end % gMin check.
	%
return;
