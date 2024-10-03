## Building
First clone HElib with git.
```shell
git clone git@github.com:homenc/HElib.git
cd HElib
git checkout [3e337a6]
```
Then, copy `github_20241003.patch` into `HElib/` and apply the changes by
```shell
git apply github_20241003.patch
```
You need to replace the path `path_of_fatboot.cpp` in `HElib/src/extractDigits.cpp` by the directory path of `fatboot.cpp`.
A directory called `saved_ZZX` needs to be created in this directory to store the computed digit extraction polynomials.

Then build and install HElib following the instructions in [INSTALL.md](https://github.com/homenc/HElib/blob/master/INSTALL.md), some dependencies like NTL and gmp may be needed. Remember to build in either `Release` or `RelWithDebInfo` mode.

You need to change the `helib_DIR` variable in `CMakeLists.txt` to the correct installation path of HElib.

Finally, switch to the directory of `fatboot.cpp` and build the test program by
```shell
cmake -S . -B build
cmake --build build
```

## Running
If the building is successful, it will produce a binary file `build/fatboot`.
Running `./build/fatboot -h` will display the following help message
```shell
Usage: ./build/fatboot [i=<arg>] [h=<arg>] [t=<arg>] [newbts=<arg>] [newks=<arg>] [thick=<arg>] [repeat=<arg>] [rad2=<arg>] [baseline=<arg>] [cache=<arg>] [s2cFirst=<arg>]
  i        index of the chosen parameter set [ default=0 ]
  h        hwt of encapsulated key [ default=0 ]
  t        parameter t used in new bts [ default=0 ]
  newbts   if new bts is used [ default=0 ]
  newks    if new ks is used [ default=1 ]
  thick    if thick bts is used [ default=0 ]
  repeat   number of tests [ default=5 ]
  rad2     force to use radix-2 [ default=0 ]
  baseline to test baseline [ default=0 ]
  cache    use cache or not [ default=1 ]
  s2cFirst if SlotToCoeff is the first step in general bootstrapping [ default=0 ]
```
+ Three parameter sets are available (given in Table 2 of the paper with ID I, II, III). By setting the argument `i` to an integer in 0 to 2, the corresponding parameter set is chosen.
+ The Hamming weight of the encapsulated secret key needs to be assigned through the argument `h`. `h` is 26 for parameter sets I and II, and is 24 for parameter set II.
+ Keep argument `t` as the default value 0.
+ Always set `newbts=1` so that the optimization of Ma et al. is used.
+ Set `thick=1` to test general bootstrapping. Set `thick=0` to test thin bootstrapping.
+ Argument `repeat` indicates the number of tests to run. The performance data will be averaged over these tests. 
+ Argument `rad2` controls the decomposition style when $p\equiv 3\bmod 4$. `rad2=0` means Bruun style while `rad2=1` means Radix-2 style.
+ Set `baseline=1` to test the baselines for general or thin bootstrapping.
+ Set `cache=1` to precompute $\kappa(i)$ in linear transformations as double-CRTs. Refer to Section 5.2 for its discussion.
+ Set `s2cFirst=1` to evaluate SlotToCoeff as the first step in general bootstrapping. Refer to the comment below Table 5 for details.

