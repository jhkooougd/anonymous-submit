all : improved_11_round_attack improved_12_round_attack

improved_11_round_attack : speck.cpp improved_11_round_attack.cpp
	g++ -O2 -std=c++11 -o $@ $^

improved_12_round_attack : speck.cpp improved_12_round_attack.cpp
	g++ -O2 -std=c++11 -o $@ $^