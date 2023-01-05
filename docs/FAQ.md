# Frequently Asked Questions

## Why do I see the error: `cudaMalloc(...) returned out of memory`?

The reason that `cudaMalloc` error is appearing is that the memory on your current GPU is not sufficient to run the simulation. The current version requires all of the data needed for a simulation to fit in GPU memory.
