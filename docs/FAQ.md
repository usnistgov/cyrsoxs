Frequently Asked Questions
========================================
[Why do I see the errror: `cudaMalloc(...) returned out of memory`?](#cudaMalloc)

## Why do I see the errror: `cudaMalloc(...) returned out of memory`?

The reason that `cudaMalloc` error is appearing is that the memory on your current GPU is not 
sufficient to run the simulation. Current version assumes that all the data needed for simulation can 
resides in the GPU memory. 

  