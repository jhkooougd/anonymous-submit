#include "speck.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <string>
#include <algorithm>

static vector<word> hd_7, hd_6;
static ofstream fout;

void generate_hd() {
    // hd_7
    for (uint32_t i = 0; i < (1 << 16); i++) {
        uint32_t hd = 0;
        uint32_t tmp = i;
        for (int j = 0; j < 16; j++) {
            hd += tmp & 1;
            tmp >>= 1;
        }
        if (hd <= 2) hd_7.push_back(i);
    }

    // hd_6
    for (uint32_t i = 0; i < (1 << 6); i++) {
        uint32_t hd = 0;
        uint32_t tmp = i;
        for (int j = 0; j < 6; j++) {
            hd += tmp & 1;
            tmp >>= 1;
        }
        if (hd <= 1) hd_6.push_back(i);
    }
}

void extrac_bits_to_uint(const uint32_t& n, const block x0[], const block x1[], uint32_t res[], const vector<uint32_t>& bits) {
    memset(res, 0, n * sizeof(uint32_t));
    for (int i = 0; i < n; i++) {
        for (auto x : bits) {
            res[i] = (res[i] << 1) | ((x0[i].first >> x) & 1);
        }
        for (auto x : bits) {
            res[i] = (res[i] << 1) | ((x0[i].second >> x) & 1);
        }
        for (auto x : bits) {
            res[i] = (res[i] << 1) | ((x1[i].first >> x) & 1);
        }
        for (auto x : bits) {
            res[i] = (res[i] << 1) | ((x1[i].second >> x) & 1);
        }
    }
}

void verify_search(double bs, const block c0[], const block c1[], const uint32_t& structure_size, const word& kg11, const word& kg10, const float nd[], const vector<uint32_t>& bits, word& bk10, word& bk11) {
    block *d0 = new block[structure_size], *d1 = new block[structure_size];
    block *e0 = new block[structure_size], *e1 = new block[structure_size];
    uint32_t *X = new uint32_t[structure_size];
    word sur_kg11, sur_kg10;
    double score;
    bk10 = kg10;
    bk11 = kg11;
    for (word x : hd_7) {
        sur_kg11 = kg11 ^ x;
        for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], sur_kg11, d0[i]), dec_one_round(c1[i], sur_kg11, d1[i]);
        for (word y : hd_6) {
            sur_kg10 = kg10 ^ y;
            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], sur_kg10, e0[i]), dec_one_round(d1[i], sur_kg10, e1[i]);
            extrac_bits_to_uint(structure_size, e0, e1, X, bits);
            score = 0;
            for (int i = 0; i < structure_size; i++) {
                score += nd[X[i]];
            }
            if (score > bs) {
                bs = score;
                bk10 = sur_kg10;
                bk11 = sur_kg11;
            }
        }
    }
    delete[] d0, d1, e0, e1, X;
}

bool load_table_from_file(const char *path, float table[], const uint32_t& table_size) {
    ifstream fin(path, ios::in | ios::binary);
    if (sizeof(float) != 4) {
        printf("the size of float is not 4 bytes!\n");
        return false;
    }
    if (!fin.is_open()) {
        printf("fail to open file!\n");
        return false;
    }
    fin.read((char*)table, 4 * table_size);
    fin.close();
    
    // test machine endianness
    uint32_t test_int = 0x12345678;
    unsigned char *p = (unsigned char*)(&test_int);
    if (*p != 0x78) {
        // big endian
        char *tmp = (char*)table;
        for (int i = 0; i < 4 * table_size; i += 4) {
            swap(tmp[i], tmp[i + 3]);
            swap(tmp[i + 1], tmp[i + 2]);
        }
    }
    return true;
}

void make_homogeneous_set(const block& diff, const vector<neutral_bit>& neutral_bits, block p0[], block p1[]) {
    p0[0].first = RAND_WORD, p0[0].second = RAND_WORD;
    uint32_t structure_size = 1;
    for (auto x : neutral_bits) {
        word dl = 0, dr = 0;
        for (int i = 0; i < x.bit_size; i++) {
            if (x.bit_pos[i] < WORD_SIZE) {
                dr |= 1 << x.bit_pos[i];
            } else {
                dl |= 1 << (x.bit_pos[i] - WORD_SIZE);
            }
        }
        for (int i = 0, j = structure_size; i < structure_size; i++, j++) {
            p0[j].first = p0[i].first ^ dl;
            p0[j].second = p0[i].second ^ dr;
        }
        structure_size <<= 1;
    }
    for (int i = 0; i < structure_size; i++) {
        p1[i].first = p0[i].first ^ diff.first;
        p1[i].second = p0[i].second ^ diff.second;
    }
}

void attack_with_three_NDs(const uint32_t& nr, const block& diff, const double c[], const float* nd_table[], const vector<uint32_t>* bits[], const vector<neutral_bit>& NBs,
                            word& dk_10, word& dk_11, double& time_cost, uint32_t& structure_consumption) {
    word key[M];
    word ks[MAX_NR];
    word sk11, sk10;
    double score;
    const uint32_t structure_size = 1 << NBs.size();
    block *p0 = new block[structure_size], *p1 = new block[structure_size];
    block *c0 = new block[structure_size], *c1 = new block[structure_size];
    block *d0 = new block[structure_size], *d1 = new block[structure_size];
    block *e0 = new block[structure_size], *e1 = new block[structure_size];
    uint32_t *X = new uint32_t[structure_size];
    uint32_t num = 0;
    clock_t start = clock();
    for (int i = 0; i < M; i++) key[i] = RAND_WORD;
    expand_key(key, ks, nr);
    sk11 = ks[nr - 1], sk10 = ks[nr - 2];
    printf("true subkey is (0x%x, 0x%x)\n", sk10 & 0x3f, sk11);
    while (true) {
        // generate a plaintext structure
        make_homogeneous_set(diff, NBs, p0, p1);
        for (int i = 0; i < structure_size; i++) dec_one_round(p0[i], 0, p0[i]), dec_one_round(p1[i], 0, p1[i]);
        for (int i = 0; i < structure_size; i++) encrypt(p0[i], ks, nr, c0[i]), encrypt(p1[i], ks, nr, c1[i]);
        num += 1;
        // stage 1: guess kg_11[5~0]
        for (word kg_11_L = 0; kg_11_L < (1 << 6); kg_11_L++) {
            for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_11_L, d0[i]), dec_one_round(c1[i], kg_11_L, d1[i]);
            extrac_bits_to_uint(structure_size, d0, d1, X, *bits[0]);
            score = 0;
            // access lookup table 1 and calculate key guess score
            for (int i = 0; i < structure_size; i++) {
                score += nd_table[0][X[i]];
            }
            if (score > c[0]) {
                // stage 2: guess kg_11[14~6]
                for (word kg_11_H = 0; kg_11_H < (1 << 9); kg_11_H++) {
                    word kg_11 = (kg_11_H << 6) | kg_11_L;
                    for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_11, d0[i]), dec_one_round(c1[i], kg_11, d1[i]);
                    extrac_bits_to_uint(structure_size, d0, d1, X, *bits[1]);
                    score = 0;
                    // access lookup table 2 and calculate key guess score
                    for (int i = 0; i < structure_size; i++) {
                        score += nd_table[1][X[i]];
                    }
                    if (score > c[1]) {
                        // stage 3: guess (kg_10[5~0], kg_11[15])
                        // kg_11[15] == 0: guess kg_10[5~0]
                        for (word kg_10 = 0; kg_10 < (1 << 6); kg_10++) {
                            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], kg_10, e0[i]), dec_one_round(d1[i], kg_10, e1[i]);
                            extrac_bits_to_uint(structure_size, e0, e1, X, *bits[2]);
                            score = 0;
                            // access lookup table 3 and calculate key guess score
                            for (int i = 0; i < structure_size; i++) {
                                score += nd_table[2][X[i]];
                            }
                            if (score > c[2]) {
                                word bk_10, bk_11;
                                verify_search(score, c0, c1, structure_size, kg_11, kg_10, nd_table[2], *bits[2], bk_10, bk_11);
                                clock_t end = clock();
                                dk_10 = bk_10 ^ (sk10 & 0x3f);
                                dk_11 = bk_11 ^ sk11;
                                time_cost = (end - start + 0.0) / CLOCKS_PER_SEC;
                                structure_consumption = num;
                                delete[] p0, p1, c0, c1, d0, d1, e0, e1, X;
                                return;
                            }
                        }
                        // kg_11[15] == 1: guess kg_10[5~0]
                        word kg_11_n = kg_11 | (1 << 15);
                        for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_11_n, d0[i]), dec_one_round(c1[i], kg_11_n, d1[i]);
                        for (word kg_10 = 0; kg_10 < (1 << 6); kg_10++) {
                            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], kg_10, e0[i]), dec_one_round(d1[i], kg_10, e1[i]);
                            extrac_bits_to_uint(structure_size, e0, e1, X, *bits[2]);
                            score = 0;
                            // access lookup table 3 and calculate key guess score
                            for (int i = 0; i < structure_size; i++) {
                                score += nd_table[2][X[i]];
                            }
                            if (score > c[2]) {
                                word bk_10, bk_11;
                                verify_search(score, c0, c1, structure_size, kg_11_n, kg_10, nd_table[2], *bits[2], bk_10, bk_11);
                                clock_t end = clock();
                                dk_10 = bk_10 ^ (sk10 & 0x3f);
                                dk_11 = bk_11 ^ sk11;
                                time_cost = (end - start + 0.0) / CLOCKS_PER_SEC;
                                structure_consumption = num;
                                delete[] p0, p1, c0, c1, d0, d1, e0, e1, X;
                                return;
                            }
                        }
                    }
                }
            }
        }
    }
}

void test(const uint32_t& t, const uint32_t& nr, const block& diff, const double c[], const float* nd_table[], const vector<uint32_t>* bits[], const vector<neutral_bit>& NBs) {
    double structure_sum = 0, time_sum = 0, acc = 0;
    for (int i = 0; i < t; i++) {
        double time_cost;
        uint32_t structure_consumption = 0;
        word dk_10, dk_11;
        printf("cur t is %d\n", i);
        attack_with_three_NDs(nr, diff, c, nd_table, bits, NBs, dk_10, dk_11, time_cost, structure_consumption);
        printf("the number of tested structures is %d\n", structure_consumption);
        printf("time consumption of current attack is %f\n", time_cost);
        printf("differences between kg and sk are (0x%x, 0x%x)\n", dk_10, dk_11);
        uint32_t d = 0;
        word tmp = dk_10;
        for (int j = 0; j < 6; j++) {
            d += tmp & 1;
            tmp >>= 1;
        }
        if (dk_11 == 0 && d <= 1) acc += 1;
        time_sum += time_cost;
        structure_sum += structure_consumption;
        fout << hex << dk_10 << ' ' << dk_11 << ' ' << time_cost << ' ' << dec << structure_consumption << endl;
    }
    printf("the average time consumption is %f\n", time_sum / t);
    printf("the average number of structure consumption is %f\n", structure_sum / t);
    printf("the accuracy is %f\n", acc / t);
}

int main() {
    check_testvector();
    string res_file = "./key_recovery_res/11_round_attack/attack_res";
    fout.open(res_file);
    srand(time(0));
    generate_hd();
    vector<uint32_t> bits1_for_ND7 = {12, 11, 10, 9, 8, 7};
    vector<uint32_t> bits2_for_ND7 = {14, 13, 12, 11, 5, 4};
    vector<uint32_t> bits3_for_ND6 = {12, 11, 10, 9, 8, 7};
    vector<uint32_t>* bits[3] = {&bits1_for_ND7, &bits2_for_ND7, &bits3_for_ND6};
    vector<neutral_bit> neutral_bits({{1,{20}}, {1,{21}}, {1,{22}}, {2,{9,16}}, {3,{2,11,25}}, {1,{14}}, {1,{15}}, {2,{6,29}}, {1,{23}}, {1,{30}}});
    float *table_1 = new float[1 << 24], *table_2 = new float[1 << 24], *table_3 = new float[1 << 24];
    if (!load_table_from_file("./12_7_nd7_table", table_1, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    if (!load_table_from_file("./14_11_5_4_nd7_table", table_2, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    if (!load_table_from_file("./12_7_nd6_table", table_3, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    float* tables[3] = {table_1, table_2, table_3};
    double c[3] = {30, 40, 80};
    test(1000, 11, {0x211, 0xa04}, c, (const float**)tables, (const vector<uint32_t>**)bits, neutral_bits);

    // only use the first eight neutral bits
    // vector<neutral_bit> neutral_bits_2({{1,{20}}, {1,{21}}, {1,{22}}, {2,{9,16}}, {3,{2,11,25}}, {1,{14}}, {1,{15}}, {2,{6,29}}});
    // double c_2[3] = {7.5, 10, 20};
    // test(1000, 11, {0x211, 0xa04}, c_2, (const float**)tables, (const vector<uint32_t>**)bits, neutral_bits_2);

    delete[] table_1, table_2, table_3;
    return 0;
}