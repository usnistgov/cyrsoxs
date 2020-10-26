Frequently Asked Questions
========================================
[I see the errror: `cudaMalloc(...) returned out of memory`. What should I do?](#cudaMalloc)

## I see the errror: `cudaMalloc(...) returned out of memory`. What should I do?

The reason that `cudaMalloc` error is appearing is that the memory on your current GPU is not 
sufficient to run the simulation. Current version assumes that all the data needed for simulation can 
resides in the GPU memory. 

  