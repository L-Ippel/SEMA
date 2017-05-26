# SEMA 
Streaming EM approximation 

This project is still work in progress and we are currently working on making the SEMA function more user friendly. 
At this moment the SEMA algorithm can fit mixed effects/multilevel models with fixed effects and random effects.

The function is written such that you only have to add the new incoming data and the current state of the parameters. At the first call, so when the current state is not yet available, the function creates a list with all required objects. Also when a new individual enters, the function automatically creates a list with all needed objects.

Using the new formulation of the SEMA algorithm, random intercepts models can also be fitted (without using the previous entry)
