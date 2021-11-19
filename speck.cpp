#include <stdio.h>
#include <time.h>
#include <random>
#include "speck.h"

word rol(const word& a, const uint32_t& b) {
    return (a << b) | (a >> (WORD_SIZE - b));
}

word ror(const word& a, const uint32_t& b) {
    return (a >> b) | (a << (WORD_SIZE - b));
}

void enc_one_round(const block& p, const word& k, block& c) {
    word& cl = c.first = p.first;
    word& cr = c.second = p.second;
    cl = ror(cl, ALPHA);
    cl += cr;
    cl ^= k;
    cr = rol(cr, BETA);
    cr ^= cl;
}

void dec_one_round(const block& c, const word& k, block& p) {
    word& pl = p.first = c.first;
    word& pr = p.second = c.second;
    pr ^= pl;
    pr = ror(pr, BETA);
    pl ^= k;
    pl -= pr;
    pl = rol(pl, ALPHA);
}

void expand_key(const word mk[], word keys[], const uint32_t& nr) {
    keys[0] = mk[M - 1];
    word l[MAX_NR];
    block inner_state;
    for (int i = 0, j = M - 2; i < M - 1; i++, j--) {
        l[i] = mk[j];
    }
    for (int i = 0, j = M - 1, k = 1; i < nr - 1; i++, j++, k++) {
        enc_one_round(block(l[i], keys[i]), i, inner_state);
        l[j] = inner_state.first;
        keys[k] = inner_state.second;
    }
}

void encrypt(const block& p, const word keys[], const uint32_t& nr, block& c) {
    c = p;
    for (int i = 0; i < nr; i++) {
        enc_one_round(c, keys[i], c);
    }
}

void decrypt(const block& c, const word keys[], const uint32_t& nr, block& p) {
    p = c;
    for (int i = nr - 1; i >= 0; i--) {
        dec_one_round(p, keys[i], p);
    }
}

bool check_testvector() {
    word mk[4] = {0x1918, 0x1110, 0x0908, 0x0100};
    block p = {0x6574, 0x694c};
    word keys[22];
    expand_key(mk, keys, 22);
    block c;
    encrypt(p, keys, 22, c);
    if (c == block(0xa868, 0x42f2)) {
        printf("Testvector verified.\n");
        return true;
    } else {
        printf("Testvector not verified.\n");
        return false;
    }
}