# 1 Summary

This repository holds supplementary code of the paper "Analysis and Improvements of Deep Learning-based Key Recovery Attack". It contains:

* three lookup tables (12_7_nd7_table, 14_11_5_4_nd7_table and 14_9_nd6_table) mentioned in this paper. Each table is a 64 MB binary file containing 2\*\*24 floating point numbers representing 2\*\*24 entries of neutral distinguisher, and the size of each floating point number is 4 bytes,
* a C++ implementation of Speck32/64 (speck.cpp and speck.h),
* the code written in C++ of the improved attack on 11-round Speck32/64 (improved_11_round_attack.cpp),
* the code written in C++ of the improved attack on 12-round Speck32/64 (improved_12_round_attack.cpp).

# 2 Compiling the attack code

There are two ways to compile the attack code:

Method 1: Compile the attack code using Makefile:

1. G++ that supports C++2011 is required to compile these C++ programs.
2. Run `make` or `make all` in this repository to build the C++ source files, and you will get two executable files:
   * improved_11_round_attack for improved_11_round_attack.cpp to attack 11-round Speck32/64,
   * improved_12_round_attack for improved_12_round_attack.cpp to attack 12-round Speck32/64.
3. If g++ is missing, for Windows OS, you can download and install MinGW-W64 through [MinGW-w64 - for 32 and 64 bit Windows - Browse /mingw-w64 at SourceForge.net](https://sourceforge.net/projects/mingw-w64/files/mingw-w64/).
4. If you have installed MinGW-W64 but can't find the command `make`, try to run `mingw32-make` or `mingw32-make all` instead.

Method 2: Compile the attack code using your own IDE:

1. Note that there's a main function both in improved_11_round_attack.cpp and improved_12_round_attack.cpp. So use improved_11_round_attack.cpp, speck.h and speck.cpp  to test 11-round attack and use improved_12_round_attack.cpp, speck.h and speck.cpp to test 12-round attack.
7. Please test the attacks in Release mode instead of Debug mode.
