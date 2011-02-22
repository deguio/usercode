setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv DYLD_LIBRARY_PATH ${THISDIR}/lib
setenv PATH ${PATH}:${THISDIR}/bin


setenv NTUPLEPKGDIR ${THISDIR}/../NtuplePackage/

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${NTUPLEPKGDIR}/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${NTUPLEPKGDIR}/lib
setenv PATH ${PATH}:${NTUPLEPKGDIR}/bin
setenv PATH ${PATH}:${NTUPLEPKGDIR}/scripts

setenv NTUPLEPKGINCLUDE ${NTUPLEPKGDIR}/interface
setenv NTUPLEPKGLIB ${NTUPLEPKGDIR}/lib