clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt_contour";
%
tic();
%studyPtDat_mod = studyPtDat;
%studyPtDat_mod.curveDat = studyPtDat.curveDat([1,2,4,8]);
%vizPt_curve2d(studyPtDat_mod, 1, 2, "omegalin");
figure(1);
vizPt_curve2d(studyPtDat, 1, 8, "omegalin");
figure(2);
vizPt_curve2d(studyPtDat, 1, 8, "omega");
toc;
