import numpy as np
import speck as sp
from os import urandom
from os.path import isfile
import time
import generate_look_up_table as glut

word_size = sp.WORD_SIZE()


def extract_sensitive_bits(raw_x, bits=None):
    # get new-x according to sensitive bits
    id0 = [word_size - 1 - v for v in bits]
    id1 = [v + word_size * i for i in range(4) for v in id0]
    new_x = raw_x[:, id1]
    return new_x


# make a plaintext structure from a random plaintext pair
# diff: difference of the plaintext pair
# neutral_bits is used to form the plaintext structure
def make_plaintext_structure(diff=(0x211, 0xa04), neutral_bits=None):
    p0l = np.frombuffer(urandom(2), dtype=np.uint16)
    p0r = np.frombuffer(urandom(2), dtype=np.uint16)
    for i in neutral_bits:
        if isinstance(i, int):
            i = [i]
        d0 = 0
        d1 = 0
        for j in i:
            d = 1 << j
            d0 |= d >> 16
            d1 |= d & 0xffff
        p0l = np.concatenate([p0l, p0l ^ d0])
        p0r = np.concatenate([p0r, p0r ^ d1])
    p1l = p0l ^ diff[0]
    p1r = p0r ^ diff[1]
    return p0l, p0r, p1l, p1r


# x0l shape: [n]
# output shape: [n]
def extrac_bits_to_uint(x0l, x0r, x1l, x1r, bits=None):
    n = len(x0l)
    m = len(bits)
    r0l = np.zeros(n, dtype=np.uint32)
    r0r = np.zeros(n, dtype=np.uint32)
    r1l = np.zeros(n, dtype=np.uint32)
    r1r = np.zeros(n, dtype=np.uint32)
    for i in range(m):
        index = bits[i]
        offset = m - 1 - i
        r0l = r0l + (((x0l >> index) & 1) << offset)
        r0r = r0r + (((x0r >> index) & 1) << offset)
        r1l = r1l + (((x1l >> index) & 1) << offset)
        r1r = r1r + (((x1r >> index) & 1) << offset)
    res = (r0l << (3*m)) + (r0r << (2*m)) + (r1l << m) + r1r
    return res


def hw(x, hd):
    num = len(x)
    tp = np.zeros(num, dtype=np.uint8)
    for i in range(16):
        tp = tp + ((x >> i) & 1)
    res = np.array([v for v in range(num) if tp[v] <= hd], dtype=np.uint16)
    return res


hd_7 = hw(x=np.array(range(2**16), dtype=np.uint16), hd=2)
hd_6 = hw(x=np.array(range(2**6), dtype=np.uint16), hd=1)


def verify_search(bs, c0l, c0r, c1l, c1r, kg12, kg11, nd, bits):
    sur_kg12 = kg12 ^ hd_7
    sur_kg11 = kg11 ^ hd_6
    bk = [kg11, kg12]
    for k12 in sur_kg12:
        d0l, d0r = sp.dec_one_round((c0l, c0r), k12)
        d1l, d1r = sp.dec_one_round((c1l, c1r), k12)
        for k11 in sur_kg11:
            e0l, e0r = sp.dec_one_round((d0l, d0r), k11)
            e1l, e1r = sp.dec_one_round((d1l, d1r), k11)
            x = extrac_bits_to_uint(e0l, e0r, e1l, e1r, bits=bits)
            z = nd[x]

            s = np.sum(z)
            if s > bs:
                bs = s
                bk = [k11, k12]
    return bk


def attack_with_three_NDs(nr=12, diff=(0x8020, 0x4101), c=None, nd_path=None, bits=None, NBs=None):
    nd = []
    for file in nd_path:
        nd.append(np.load(file))
    if nr == 12:
        while 1:
            key = np.frombuffer(urandom(8), dtype=np.uint16).reshape(4, -1)
            ks = sp.expand_key(key, nr)
            if (ks[2][0] >> 12) & 1 != (ks[2][0] >> 11) & 1:
                break
    else:
        key = np.frombuffer(urandom(8), dtype=np.uint16).reshape(4, -1)
        ks = sp.expand_key(key, nr)
    sk12, sk11 = ks[nr - 1][0], ks[nr - 2][0]
    print('true subkey is ', (hex(sk11 & 0x3f), hex(sk12)))
    num = 0
    start = time.time()
    while 1:
        p0l, p0r, p1l, p1r = make_plaintext_structure(diff=diff, neutral_bits=NBs)
        p0l, p0r = sp.dec_one_round((p0l, p0r), 0)
        p1l, p1r = sp.dec_one_round((p1l, p1r), 0)
        c0l, c0r = sp.encrypt((p0l, p0r), ks)
        c1l, c1r = sp.encrypt((p1l, p1r), ks)

        num = num + 1
        print('\r num is {}'.format(num), end=' ')

        # partial with guess of sk12[5~0]
        for kg_12_L in range(2 ** 6):
            d0l, d0r = sp.dec_one_round((c0l, c0r), kg_12_L)
            d1l, d1r = sp.dec_one_round((c1l, c1r), kg_12_L)
            x1 = extrac_bits_to_uint(d0l, d0r, d1l, d1r, bits=bits[0])
            z1 = nd[0][x1]

            s1 = np.sum(z1)
            if s1 > c[0]:
                for kg_12_H in range(2**9):
                    kg_12 = (kg_12_H << 6) + kg_12_L
                    e0l, e0r = sp.dec_one_round((c0l, c0r), kg_12)
                    e1l, e1r = sp.dec_one_round((c1l, c1r), kg_12)
                    x2 = extrac_bits_to_uint(e0l, e0r, e1l, e1r, bits=bits[1])
                    z2 = nd[1][x2]

                    s2 = np.sum(z2)
                    if s2 > c[1]:
                        # the 15-th bit of kg_11 is 0
                        for kg_11 in range(2 ** 6):
                            f0l, f0r = sp.dec_one_round((e0l, e0r), kg_11)
                            f1l, f1r = sp.dec_one_round((e1l, e1r), kg_11)
                            x3 = extrac_bits_to_uint(f0l, f0r, f1l, f1r, bits=bits[2])
                            z3 = nd[2][x3]

                            s3 = np.sum(z3)
                            if s3 > c[2]:
                                kg11, kg12 = verify_search(s3, c0l, c0r, c1l, c1r, kg12=kg_12, kg11=kg_11, nd=nd[2], bits=bits[2])
                                end = time.time()
                                return [kg11 ^ (sk11 & 0x3f), kg12 ^ sk12], end - start, num
                        # the 15-th bit of kg_11 is 1
                        kg_12_n = kg_12 + (1 << 15)
                        g0l, g0r = sp.dec_one_round((c0l, c0r), kg_12_n)
                        g1l, g1r = sp.dec_one_round((c1l, c1r), kg_12_n)
                        for kg_11 in range(2 ** 6):
                            h0l, h0r = sp.dec_one_round((g0l, g0r), kg_11)
                            h1l, h1r = sp.dec_one_round((g1l, g1r), kg_11)
                            x3 = extrac_bits_to_uint(h0l, h0r, h1l, h1r, bits=bits[2])
                            z3 = nd[2][x3]

                            s3 = np.sum(z3)
                            if s3 > c[2]:
                                kg11, kg12 = verify_search(s3, c0l, c0r, c1l, c1r, kg12=kg_12_n, kg11=kg_11, nd=nd[2], bits=bits[2])
                                end = time.time()
                                return [kg11 ^ (sk11 & 0x3f), kg12 ^ sk12], end - start, num


def bit_count(x):
    count = 0
    while x > 0:
        count += 1
        x = x & (x - 1)
    return count


def test(t=100, nr=12, diff=(0x8020, 0x4101), c=None, nd_path=None, bits=None, NBs=None):
    structure_sum, time_sum, acc = 0, 0, 0
    res = []
    for i in range(t):
        print('cur t is ', i)
        kg_t, time_t, num_t = attack_with_three_NDs(nr=nr, diff=diff, c=c, nd_path=nd_path, bits=bits, NBs=NBs)
        print(' ')
        print('the number of tested sturctures is ', num_t)
        print('time consumption of current attack is ', time_t)
        print('differences between kg and sk are ', (hex(kg_t[0]), hex(kg_t[1])))
        if bit_count(kg_t[0]) <= 1 and bit_count(kg_t[1]) == 0:
            acc = acc + 1
        time_sum = time_sum + time_t
        structure_sum = structure_sum + num_t
        res.append(kg_t)
    print('the average time consumption is ', time_sum / t)
    print('the average number of structures is ', structure_sum / t)
    print('the accuracy is ', acc / t)
    np.save('./key_recovery_record/naive_attack_res.npy', res)


if __name__ == '__main__':
    # replace neural distinguisher with lookup tables
    if not isfile('./saved_model/0.5854_14_12_6_4_nd6.npy'):
        print('generating lookup tables.')
        glut.generate_required_tables()

    bits1_for_ND7 = [12, 11, 10, 9, 8, 7]
    bits2_for_ND7 = [14, 13, 12, 11, 5, 4]
    bits3_for_ND6 = [12, 11, 10, 9, 8, 7]
    bits = [bits1_for_ND7, bits2_for_ND7, bits3_for_ND6]
    table_1 = './saved_model/0.5416_12_7_nd7.npy'
    table_2 = './saved_model/0.5722_14_11_5_4_nd7.npy'
    table_3 = './saved_model/0.6388_12_7_nd6.npy'
    tables = [table_1, table_2, table_3]
    c = [10, 8, 30]
    neutral_bits = [22, 20, 13, [12, 19], [14, 21], [6, 29], 30, [0, 8, 31]]
    test(t=1000, nr=12, diff=(0x8020, 0x4101), c=c, nd_path=tables, bits=bits, NBs=neutral_bits)