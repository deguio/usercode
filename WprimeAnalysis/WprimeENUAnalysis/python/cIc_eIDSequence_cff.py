import FWCore.ParameterSet.Config as cms


# add electron id 

from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *


eIDClassesLoose = eidCutBasedClassesExt.clone()
eIDClassesLoose.electronQuality = 'loose'

eIDClassesMedium = eidCutBasedClassesExt.clone()
eIDClassesMedium.electronQuality = 'medium'

eIDClassesTight = eidCutBasedClassesExt.clone()
eIDClassesTight.electronQuality = 'tight'
cIc_eIDSequence = cms.Sequence(
                     eIDClassesLoose+ 
                     eIDClassesMedium+ 
                     eIDClassesTight )



