# SplitCache

## SplitCache paper Implemented 
And also the default solid cache implemented.

> Please read the paper for better understanding what is its porpuse.

This code basicaly tries to count the number of cache misses when multipling matrix A and matrix B to generate matrix C.  
It also compares 6 different algorithms when doing this job. ( Inner product, Outer product, Gustavson and ...).
## Referencing:
*If you are using this code or any other ideas from the paper/code. I would be appreciated if you cite to this paper:*
  > SplitCache: ...

### Before starting:
- The Cache by defualt uses *LRU replacement policy*. If you want to use other replacement policies, you should edit the code.
- `k0 = 0xF0000000;` is the default starting address for saving matrices inside the RAM, you could change it to any address you want.
- The cache sizes are indicated by their rows inside the **Sizes.txt** file and also by the terminal arguments that indicates **Block Size** and **Number of Ways inside each row** (which will be discussed in the following)
- `string filename00 = "_result";` used for naming the *.csv* result files. 
- You need to adjust **Sizes.txt** and **terminal codes.txt** for your own input matrix dimentions and cache sizes. (which is described in the following)
- You can add multiple lines to **Sizes.txt** and **terminal codes.txt**. But they should be in the same struction as described.

### Compile
- First compile the `cpp` files using:
  ```
  g++ -std=c++11 -pthread <input cpp file> <output file name>
  ```
- If you want to use back-and-forth optimization, edit the line: 
```
#define my_optimize1 false
```
and change `false` to `true`

- If you are compiling this code on Windows OS make sure to uncomment the line: 
```
//typedef unsigned long long int ulong;	
```

### How to use
The file **Sizes.txt** in each folder contains the number of rows for the cache. 
- in *Solid Cache* for example, it is like this:
  ```
  256 0 0 256 0 0 256 0 0
  ```
which means that it contains 256 rows for the cache. We used 9 arguments here, you will see in the following why :)
- in *Splitted Cache* for example, it is like this:
```
128 32 32 128 32 32 128 32 32
```
which means: `<row of the Cache for Matrix **A** for *Inner and Rev. Inner product*> <row of the Cache for Matrix **B** for *Inner and Rev. Inner product*> <row of the Cache for Matrix **C** for *Inner and Rev. Inner product*> <row of the Cache for Matrix **A** for *Outer and Rev. Outer product*> <row of the Cache for Matrix **B** for *Outer and Rev. Outer product*> <row of the Cache for Matrix **C** for *Outer and Rev. Outer product*> <row of the Cache for Matrix **A** for *Gustavson and Rev. Gustavson product*> <row of the Cache for Matrix **B** for *Gustavson and Rev. Gustavson product*> <row of the Cache for Matrix **C** for *Gustavson and Rev. Gustavson product*>`

Multiple sizes can be used in each line for both *solid cache* and *SplitCache*.  
So for *solid cache*, size of matrix B and C are 0. Because they use the cache space of matrix A jointly.  
The file **terminal Codes.txt** contains the arguments that passed to the executable file through the terminal. for example:
```
./solidCache.exe 128 4096 2048 128 4
```
Here, it contains:
`./<name of the Exe. file> <M argument for Matrix> <K argument for Matrix> <N argument for Matrix> <Cache Block size in Bytes> <Cache ways per row>`

