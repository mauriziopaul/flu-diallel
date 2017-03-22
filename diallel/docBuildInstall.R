#!/usr/local/bin/RScript

require("devtools")

docBuildInstall  <- function(dir, ...){
	document("../treatmentResponseDiallel", roclets=c('rd', 'namespace'))
	build("../treatmentResponseDiallel")
	install("../treatmentResponseDiallel", quick=TRUE, upgrade_dependencies=FALSE)
}
 