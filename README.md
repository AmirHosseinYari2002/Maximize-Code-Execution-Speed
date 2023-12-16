<h1 align='center'> Optimize Taylor Series Expansion </h1>

<h2 align='center'> Sharif University of Technology </h2>

<h3 align='center'> Electrical Engineering Department </h3>

In this project, we want to calculate the Taylor series expansion of the exponential function.
To increase the speed of code execution, I have written 3 different functions that optimize the speed respectively.
- fx1: In this function, I took the common denominator between the two terms, which will cause the calculations to be halved. I have also used multithreading and SIMD.
- fx2: In this function, I added loop unroll to the code.
- fx3: In this function, I used SIMD and multithreading in another way.

By performing the above actions, I was able to increase the code execution speed up to 80x.

### SIMD
SIMD stands for "Single Instruction, Multiple Data". It is a type of parallel processing technique used in computer architecture and programming to accelerate data processing tasks by performing multiple operations on multiple data elements simultaneously.<br>
In SIMD, a single instruction operates on multiple data elements in parallel, allowing for efficient vectorized computation. This contrasts with traditional scalar processing, where a single instruction operates on a single data element at a time.

### Multithreading
Multithreading is a programming concept that allows multiple threads to execute concurrently within a single process. A thread is the smallest unit of execution within a program, and multithreading enables a program to perform multiple tasks or operations concurrently, taking advantage of modern multi-core processors and improving overall performance and responsiveness. <br>
In a multithreaded application, each thread operates independently and can execute different parts of the program's code simultaneously. This allows the program to perform tasks concurrently, handle multiple user interactions, and utilize available system resources more efficiently. 

### Loop Unrolling
Loop unrolling is an optimization technique used in computer programming and compiler design to improve the performance of loops. The goal of loop unrolling is to reduce the overhead of loop control and increase instruction-level parallelism, which can lead to better utilization of hardware resources and faster execution.
