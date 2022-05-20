import ctypes
import numpy as np
import cv2

lib = np.ctypeslib.load_library(libname='libising', loader_path='.')

### lib API def

# void set_rand_seed(unsigned int seed)
lib.set_rand_seed.argtypes = [
    ctypes.c_int
]

# void init_system(int in_dim0, int in_dim1)
lib.init_system.argtypes = [
    ctypes.c_int,
    ctypes.c_int
]

# void set_param(double in_betaJ, double in_betamuB, double in_freq)
lib.set_param.argtypes = [
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double
]

# void set_spin(int * in_spin)
lib.set_spin.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags="C_CONTIGUOUS")
]

# void set_random_spin()

# void get_spin(int * out_spin)
lib.get_spin.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags="C_CONTIGUOUS")
]

# double next_frame(double t, double time_per_frame)
lib.next_frame.argtypes = [
    ctypes.c_double,
    ctypes.c_double
]
lib.next_frame.restype = ctypes.c_double

# void free_system()

### end lib API def

seed = 20220518202200
dim0 = 512
dim1 = 512
betaJ = 0.5
betamuB = 0.0
freq = 1.0  # adjust time scale
spin = np.zeros([dim0, dim1], dtype=np.int32)

lib.set_rand_seed(seed)
lib.init_system(dim0, dim1)
lib.set_param(betaJ, betamuB, freq)
lib.set_random_spin()
# lib.set_spin(spin)

fps = 25  # frame per second
size = (dim0, dim1)  # video size
video = cv2.VideoWriter('Video.mp4', cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), fps, size, True)
time_per_frame = 2  # system time interval per frame

t = 0.0
for frame in range(40*25):
    if frame % 10 == 0:
        avg_spin = 2.0 * spin.sum() / (dim0 * dim1) - 1.0
        print(f'frame: {frame}, spin: {avg_spin}')

    t = lib.next_frame(t, time_per_frame)
    lib.get_spin(spin)

    out_img = np.zeros([size[0], size[1], 3], dtype=np.uint8)
    # spin: + = (255, 0, 0); - = (0, 0, 255)
    out_img[:, :, 0] = spin[:, :] * 255
    out_img[:, :, 2] = (1 - spin[:, :]) * 255
    # smoooooth
    out_img = cv2.blur(out_img, (5, 5))

    video.write(out_img)

video.release()
cv2.destroyAllWindows()

lib.free_system()
