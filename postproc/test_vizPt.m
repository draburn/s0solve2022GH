clear;
%
test_studyPt;
%close("all");
commondefs;
thisFile = "test_vizPt";
%
tic();
figure(1);
vizPt_vs(studyPtDat,"deltaNorm","omegalin",prm);
figure(2);
vizPt_vs(studyPtDat,"deltaNorm","omega",prm);
figure(3);
vizPt_curve1d(studyPtDat,[4,2],prm);
prm.figIndex0 = 3;
vizPt_curve(studyPtDat,[4,2],prm);
toc;
