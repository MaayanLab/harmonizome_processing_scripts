function outcm = cmcolmerge(incm, method)

outcm = cmtranspose(cmrowmerge(cmtranspose(incm), method));


