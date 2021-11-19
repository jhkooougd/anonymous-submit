# 1 Summary

This repository holds supplementary code of the paper "Improved Deep Learning-based Attack with Application to Round Reduced Speck32/64". It contains:

* three lookup tables (12_7_nd7_table, 14_11_5_4_nd7_table and 14_9_nd6_table) mentioned in this paper,
* a C++ implementation of Speck32/64 (speck.cpp and speck.h),
* the code written in C++ of the improved attack on 11-round Speck32/64 (improved_11_round_attack.cpp),
* the code written in C++ of the second version of improved attack on 12-round Speck32/64 (improved_12_round_attack.cpp).

# 2 Compiling the attack code

G++ that supports C++2011 is required to compile these C++ programs.

Run `make` or `make all` in this repository to build the C++ source files, and you will get two executable files:

* improved_11_round_attack for improved_11_round_attack.cpp to attack 11-round Speck32/64,
* improved_12_round_attack for improved_12_round_attack.cpp to attack 12-round Speck32/64.

