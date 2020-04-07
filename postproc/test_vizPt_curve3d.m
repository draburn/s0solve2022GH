clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt_curve3d";
%
tic();
%studyPtDat_mod = studyPtDat;
%studyPtDat_mod.curveDat = studyPtDat.curveDat([1,2,4,8]);
%vizPt_curve2d(studyPtDat_mod, 1, 2, "omegalin");
figure(1);
vizPt_curve2d( studyPtDat, 1, 2, "omega" );
figure(2);
prm.vecXP = studyPtDat.vecX0 + studyPtDat.curveDat(8).matDelta(:,end)*0.5;
vizPt_curve2d( studyPtDat, 1, 2, "omega", prm );
figure(3);
prm.vecXP = studyPtDat.vecX0 + studyPtDat.curveDat(8).matDelta(:,end)*1.0;
vizPt_curve2d( studyPtDat, 1, 2, "omega", prm );
toc;
