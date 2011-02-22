setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib

if (${?DYLD_LIBRARY_PATH}) then
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}../NtuplePackage/lib
endif

setenv PATH ${PATH}:${THISDIR}/bin
setenv PATH ${PATH}:${THISDIR}/../NtuplePackage/bin

setenv NTUPLEPKGINCLUDE ${THISDIR}/../NtuplePackage/interface
setenv NTUPLEPKGLIB ${THISDIR}/../NtuplePackage/lib

