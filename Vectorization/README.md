# Architecture of high performance processors
## Lab №1 - Vectorization
### Description:
Create manual *vectorization* in effective way.
### Task:
Multiply 2 matrices in following ways:
1. Result - C1 matrix
    * with **enabled vectorization**
    * with **disabled vectorization**
2. Result - C2 matrix
    * with **manual** vectorization using *assembly language* or *intrinsics*

* Size of **external** matrix *(chosen by student)* is **1 000 x 1 000**
* Size of **internal** matrix *(chosen by teacher)* is **4 x 4**
* Type of data *(chosen by teacher)* is **double**

### Requirements:
1. Minimum 2 functions in one project
    * **deafult** vectorization
    * **manual** vectorization
2. C1 matrix must equals C2 matrix
    * do not show matrices
3. SSE2 manual vectorization must not be slowly than deafult vectorization 
### Restrictions:
1. Can't transpose matrix
2. Check time using functions *high_resolution_clock* or *rdtsc, __rdtsc()*
3. Can't use *clock()* or *time()*
4. Execute only release version
### Remind!
* We should proof efficiency of code using MSDN disassembler or Inter Vtune(preffered)
* Speed up for float - 3-3.5, for double - 1.2-1.5

### Bonuses:
* Manual vectorization with any sizes
* Proof using Intel Vtune
* Find bottlenecks using vTune
* Uniq solution
* Optimization for L3 cache (also for bonuses L1 and L2) 

### Timeline of matrix multiplication exponent:
| Year | O(n) complexity    | Name              |
|------|-----------|----------------------------|
| 1969 | 2.8074    | Strassen                   |
| 1978 | 2.796     | Pan                        |
| 1979 | 2.780     | Bini, Capovani, Romani     |
| 1981 | 2.522     | Schönhage                  |
| 1981 | 2.517     | Romani                     |
| 1981 | 2.496     | Coppersmith, Winograd      |
| 1986 | 2.479     | Strassen                   |
| 1990 | 2.3755    | Coppersmith, Winograd      |
| 2010 | 2.3737    | Stothers                   |
| 2013 | 2.3729    | Williams                   |
| 2014 | 2.3728639 | Le Gall                    |
| 2020 | 2.3728596 | Alman, Williams            |

### Some notes of developing process:
1. During testing different theories, we understood,
that the most effective way to store matricies to improve
performance is use common array, which represents 4D matrix,
so the basic formula for access to matrix element with h, l, i, j indices
is **mtx[l][h][i][j] = arr[h*n*x*y + x*y*l + y*i + j]**, where
l, h, x and y is sizes of matix
2. Using functions decreasing performance
and I wrote this, because I thought that it's not true. I was wrong :(
3. The best algorithms for solve problem are 'divide and conquer' algorithms
4. I will use Strassens alogrithm because Algorithms (such as Coppersmith-Winograd, which works only for NxN matricies)
with a preferred asymptotic running time 
over the Strassen calculation are rarely utilized practically, 
on the grounds that the huge steady factors in their running occasions make them inappropriate. 
5. When I create Strassens' algorithm, I understand, that it's actually sucks
with huge amount of data 
6. This is reason, why I burned down: https://community.intel.com/t5/Intel-ISA-Extensions/Need-help-Why-my-avx-code-is-slower-than-SSE-code/td-p/1034874 
