#pragma once
#include <stdint.h>
#include <utility>
#include <stdlib.h>
using namespace std;

#define WORD_SIZE 16
#define ALPHA 7
#define BETA 2
#define MASK_VAL 0xffff
#define M 4
#define MAX_NR 50
#define RAND_BYTE (rand() & 0xff)
#define RAND_WORD ((RAND_BYTE << 8) | RAND_BYTE)

typedef uint16_t word;
typedef pair<word, word> block;

struct neutral_bit {
    uint32_t bit_size;
    uint32_t bit_pos[3];
};

void enc_one_round(const block& p, const word& k, block& c);
void dec_one_round(const block& c, const word& k, block& p);
void expand_key(const word mk[], word keys[], const uint32_t& nr);
void encrypt(const block& p, const word keys[], const uint32_t& nr, block& c);
void decrypt(const block& c, const word keys[], const uint32_t& nr, block& p);
bool check_testvector();