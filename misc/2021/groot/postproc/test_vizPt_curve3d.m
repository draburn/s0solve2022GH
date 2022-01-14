clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt_curve3d";
figIndex = 0;
%
tic();
%studyPtDat_mod = studyPtDat;
%studyPtDat_mod.curveDat = studyPtDat.curveDat([1,2,4,8]);
%vizPt_curve2d(studyPtDat_mod, 1, 2, "omegalin");
n1 = 1;
n2 = 3;
n3 = 8;
numPts3 = size( studyPtDat.curveDat(n3).matDelta, 2 );
%
figIndex++; figure(figIndex);
vizPt_curve2d( studyPtDat, n1, n2, "omega" );
figIndex++; figure(figIndex);
vizPt_curve2d( studyPtDat, n1, n2, "resS1" );
%
figIndex++; figure(figIndex);
m = round(numPts3*0.5);
prm.vecXP = studyPtDat.vecX0 + studyPtDat.curveDat(n3).matDelta(:,m);
prm.strPName = sprintf( "(%s plane %d/%d)", studyPtDat.curveDat(n3).curveName, m, numPts3 );
vizPt_curve2d( studyPtDat, n1, n2, "omega", prm );
figIndex++; figure(figIndex);
vizPt_curve2d( studyPtDat, n1, n2, "resS1", prm );
%
figIndex++; figure(figIndex);
m = numPts3;
prm.vecXP = studyPtDat.vecX0 + studyPtDat.curveDat(n3).matDelta(:,end);
prm.strPName = sprintf( "(%s plane %d/%d)", studyPtDat.curveDat(n3).curveName, m, numPts3 );
vizPt_curve2d( studyPtDat, n1, n2, "omega", prm );
figIndex++; figure(figIndex);
vizPt_curve2d( studyPtDat, n1, n2, "resS1", prm );
toc;
