# K10571

### sim_linear.R simlulation study when g functions are linear. 

### sim_nonlinear.R simulation study when g functions are non linear. 

Ben: The noise level seem to be too large (sigma=1) for the nonlinear example (true function ranges between -1 and 1)? Did you change it for the plots you produced? I have added your non linear function to my notebook as well. - Luyao 

Luyao: I totally agree. -Ben

Ten figures are added into the folder called figures. Each figure plots the g function fitted using the real data through linear interpolations. ***I noticed that the scales of those 10 g functions very a lot.*** Therefore, the sigmas we use to evaluate the likelihood for each g should be different. Otherwise, one or two of all ten g will dominate the likelihood function. 
