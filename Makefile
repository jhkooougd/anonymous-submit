all : improved_11_round_attack naive_12_round_attack

improved_11_round_attack : speck.cpp improved_11_round_attack.cpp
	g++ -O2 -std=c++11 -o $@ $^

naive_12_round_attack : speck.cpp naive_12_round_attack.cpp
	g++ -O2 -std=c++11 -o $@ $^