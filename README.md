# SEMA 
Streaming EM approximation algorithm 

At this moment the SEMA algorithm can fit mixed effects/multilevel models with fixed effects and random effects.

The function is written such that you only have to add the new incoming data and the current state of the parameters. At the first call, so when the current state is not yet available, the function creates a list with all required objects. Also when a new individual enters, the function automatically creates a list with all needed objects.

The SEMA package can be installed and used as follows: 

library(devtools)	
install_github(“L-Ippel/SEMA”)
library(SEMA)
?SEMA
