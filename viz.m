figIndex = 0;
%
figIndex++; figure(figIndex);
vizPt_vs( p, "deltaNorm", "omegaLin" );
%
figIndex++; figure(figIndex);
vizPt_vs( p, "deltaNorm", "omega" );
%.
prm.figIndex0 = figIndex;
vizPt_curve( p, [8,4], prm );
