# SplitCache
SplitCache paper Implemented  and also default solid cache


### Compile
- First compile the `cpp` files useing:
  ```
  g++ -std=c++11 -pthread <input cpp> <output file name.exe>
  ```
- If you want to use back-and-forth optimization edit line: 
```
#define my_optimize1 false
```
and change `false` to `true`

- If you are compiling this code on Windows OS make sure to uncomment line: 
```
//typedef unsigned long long int ulong;	
```

The Cache by defualt uses LRU replacement policy 

