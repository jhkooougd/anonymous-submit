#include "speck.h"
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <string.h>
#include <string>
#include <iostream>
#include <stdio.h>

vector<word> hd_7, hd_6;
ofstream fout;

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

void verify_search(double bs, const block c0[], const block c1[], const uint32_t& structure_size, const word& kg12, const word& kg11, const float nd[], const vector<uint32_t>& bits, word& bk11, word& bk12) {
    block *d0 = new block[structure_size], *d1 = new block[structure_size];
    block *e0 = new block[structure_size], *e1 = new block[structure_size];
    uint32_t *X = new uint32_t[structure_size];
    word sur_kg12, sur_kg11;
    double score;
    bk11 = kg11;
    bk12 = kg12;
    for (word x : hd_7) {
        sur_kg12 = kg12 ^ x;
        for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], sur_kg12, d0[i]), dec_one_round(c1[i], sur_kg12, d1[i]);
        for (word y : hd_6) {
            sur_kg11 = kg11 ^ y;
            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], sur_kg11, e0[i]), dec_one_round(d1[i], sur_kg11, e1[i]);
            extrac_bits_to_uint(structure_size, e0, e1, X, bits);
            score = 0;
            for (int i = 0; i < structure_size; i++) {
                score += nd[X[i]];
            }
            if (score > bs) {
                bs = score;
                bk11 = sur_kg11;
                bk12 = sur_kg12;
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
                            word& dk_11, word& dk_12, double& time_cost, uint32_t& structure_consumption) {
    word key[M];
    word ks[MAX_NR];
    word sk12, sk11;
    double score;
    const uint32_t structure_size = 1 << NBs.size();
    block *p0 = new block[structure_size], *p1 = new block[structure_size];
    block *c0 = new block[structure_size], *c1 = new block[structure_size];
    block *d0 = new block[structure_size], *d1 = new block[structure_size];
    block *e0 = new block[structure_size], *e1 = new block[structure_size];
    uint32_t *X = new uint32_t[structure_size];
    uint32_t num = 0;
    // ensure weak key setting
    do {
        for (int i = 0; i < M; i++) key[i] = RAND_WORD;
        expand_key(key, ks, nr);
    } while(nr == 12 && (((ks[2] >> 12) & 1) == ((ks[2] >> 11) & 1)));
    sk12 = ks[nr - 1], sk11 = ks[nr - 2];
    clock_t start = clock();
    while (true) {
        // generate a plaintext structure
        make_homogeneous_set(diff, NBs, p0, p1);
        for (int i = 0; i < structure_size; i++) dec_one_round(p0[i], 0, p0[i]), dec_one_round(p1[i], 0, p1[i]);
        for (int i = 0; i < structure_size; i++) encrypt(p0[i], ks, nr, c0[i]), encrypt(p1[i], ks, nr, c1[i]);
        num += 1;
        // stage 1: guess kg_12[5~0]
        for (word kg_12_L = 0; kg_12_L < (1 << 6); kg_12_L++) {
            for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_12_L, d0[i]), dec_one_round(c1[i], kg_12_L, d1[i]);
            extrac_bits_to_uint(structure_size, d0, d1, X, *bits[0]);
            score = 0;
            // access lookup table 1 and calculate key guess score
            for (int i = 0; i < structure_size; i++) {
                score += nd_table[0][X[i]];
            }
            if (score > c[0]) {
                // stage 2: guess kg_12[14~6]
                for (word kg_12_H = 0; kg_12_H < (1 << 9); kg_12_H++) {
                    word kg_12 = (kg_12_H << 6) | kg_12_L;
                    for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_12, d0[i]), dec_one_round(c1[i], kg_12, d1[i]);
                    extrac_bits_to_uint(structure_size, d0, d1, X, *bits[1]);
                    score = 0;
                    // access lookup table 2 and calculate key guess score
                    for (int i = 0; i < structure_size; i++) {
                        score += nd_table[1][X[i]];
                    }
                    if (score > c[1]) {
                        // stage 3: guess (kg_11[5~0], kg_12[15])
                        // kg_12[15] == 0: guess kg_11[5~0]
                        for (word kg_11 = 0; kg_11 < (1 << 6); kg_11++) {
                            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], kg_11, e0[i]), dec_one_round(d1[i], kg_11, e1[i]);
                            extrac_bits_to_uint(structure_size, e0, e1, X, *bits[2]);
                            score = 0;
                            // access lookup table 3 and calculate key guess score
                            for (int i = 0; i < structure_size; i++) {
                                score += nd_table[2][X[i]];
                            }
                            if (score > c[2]) {
                                word bk_11, bk_12;
                                verify_search(score, c0, c1, structure_size, kg_12, kg_11, nd_table[2], *bits[2], bk_11, bk_12);
                                clock_t end = clock();
                                dk_11 = bk_11 ^ (sk11 & 0x3f);
                                dk_12 = bk_12 ^ sk12;
                                time_cost = (end - start + 0.0) / CLOCKS_PER_SEC;
                                structure_consumption = num;
                                delete[] p0, p1, c0, c1, d0, d1, e0, e1, X;
                                return;
                            }
                        }
                        // kg_12[15] == 1: guess kg_11[5~0]
                        word kg_12_n = kg_12 | (1 << 15);
                        for (int i = 0; i < structure_size; i++) dec_one_round(c0[i], kg_12_n, d0[i]), dec_one_round(c1[i], kg_12_n, d1[i]);
                        for (word kg_11 = 0; kg_11 < (1 << 6); kg_11++) {
                            for (int i = 0; i < structure_size; i++) dec_one_round(d0[i], kg_11, e0[i]), dec_one_round(d1[i], kg_11, e1[i]);
                            extrac_bits_to_uint(structure_size, e0, e1, X, *bits[2]);
                            score = 0;
                            // access lookup table 3 and calculate key guess score
                            for (int i = 0; i < structure_size; i++) {
                                score += nd_table[2][X[i]];
                            }
                            if (score > c[2]) {
                                word bk_11, bk_12;
                                verify_search(score, c0, c1, structure_size, kg_12_n, kg_11, nd_table[2], *bits[2], bk_11, bk_12);
                                clock_t end = clock();
                                dk_11 = bk_11 ^ (sk11 & 0x3f);
                                dk_12 = bk_12 ^ sk12;
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
        word dk_11, dk_12;
        printf("cur t is %d\n", i);
        attack_with_three_NDs(nr, diff, c, nd_table, bits, NBs, dk_11, dk_12, time_cost, structure_consumption);
        printf("the number of tested structures is %d\n", structure_consumption);
        printf("time consumption of current attack is %f\n", time_cost);
        printf("differences between kg and sk are (0x%x, 0x%x)\n", dk_11, dk_12);
        uint32_t d = 0;
        word tmp = dk_11;
        for (int j = 0; j < 6; j++) {
            d += tmp & 1;
            tmp >>= 1;
        }
        if (dk_12 == 0 && d <= 1) acc += 1;
        time_sum += time_cost;
        structure_sum += structure_consumption;
        fout << hex << dk_11 << ' ' << dk_12 << ' ' << time_cost << ' ' << dec << structure_consumption << endl;
    }
    printf("the average time consumption is %f\n", time_sum / t);
    printf("the average number of structure consumption is %f\n", structure_sum / t);
    printf("the accuracy is %f\n", acc / t);
}

int main() {
    check_testvector();
    string res_file = "./key_recovery_res/12_round_attack/attack_res";
    fout.open(res_file);
    srand(time(0));
    generate_hd();
    vector<uint32_t> bits1_for_ND7 = {12, 11, 10, 9, 8, 7};
    vector<uint32_t> bits2_for_ND7 = {14, 13, 12, 11, 5, 4};
    vector<uint32_t> bits3_for_ND6 = {14, 13, 12, 11, 10, 9};
    vector<uint32_t>* bits[3] = {&bits1_for_ND7, &bits2_for_ND7, &bits3_for_ND6};
    vector<neutral_bit> neutral_bits({{1,{22}}, {1,{20}}, {1,{13}}, {2,{12,19}}, {2,{14,21}}, {2,{6,29}}, {1,{30}}, {3,{0,8,31}}});
    float *table_1 = new float[1 << 24], *table_2 = new float[1 << 24], *table_3 = new float[1 << 24];
    if (!load_table_from_file("./12_7_nd7_table", table_1, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    if (!load_table_from_file("./14_11_5_4_nd7_table", table_2, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    if (!load_table_from_file("./14_9_nd6_table", table_3, 1 << 24)) {
        printf("loading table from file failed!\n");
        return 0;
    }
    float* tables[3] = {table_1, table_2, table_3};
    double c[3] = {10, 8, 30};
    test(100, 12, {0x8020, 0x4101}, c, (const float**)tables, (const vector<uint32_t>**)bits, neutral_bits);
    delete[] table_1, table_2, table_3;
    fout.close();
    return 0;
}