import numpy as np
from keras.models import load_model


def uint_to_array(x, bit_len=16):
    y = np.zeros((len(x), bit_len), dtype=np.uint8)
    for j in range(bit_len):
        y[:, j] = (x >> (bit_len - 1 - j)) & 1
    return y


def generate_look_up_table(net_path=None, block_size=16):
    # load neural distinguisher
    nd = load_model(net_path)
    # get inference table
    x = np.array(range(2**block_size), dtype=np.uint32)
    new_x = uint_to_array(x, bit_len=block_size)
    y = nd.predict(new_x, batch_size=10**4, verbose=0)
    y = np.squeeze(y)
    y = np.log2(y / (1 - y))
    return y


def generate_required_tables():
    saved_folder = './saved_model/'
    net_path = saved_folder + '0.5416_12_7_student_7_distinguisher.h5'
    y = generate_look_up_table(net_path=net_path, block_size=24)
    np.save(saved_folder+'0.5416_12_7_nd7.npy', y)

    saved_folder = './saved_model/'
    net_path = saved_folder + '0.5722_14_11_5_4_student_7_distinguisher.h5'
    y = generate_look_up_table(net_path=net_path, block_size=24)
    np.save(saved_folder+'0.5722_14_11_5_4_nd7.npy', y)

    saved_folder = './saved_model/'
    net_path = saved_folder + '0.6388_12_7_student_6_distinguisher.h5'
    y = generate_look_up_table(net_path=net_path, block_size=24)
    np.save(saved_folder+'0.6388_12_7_nd6.npy', y)

    saved_folder = './saved_model/'
    net_path = saved_folder + '0.5854_14_12_6_4_student_6_distinguisher.h5'
    y = generate_look_up_table(net_path=net_path, block_size=24)
    np.save(saved_folder+'0.5854_14_12_6_4_nd6.npy', y)

if __name__ == '__main__':
    generate_required_tables()