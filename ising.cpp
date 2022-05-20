/**
 * A Kinetic Ising Model
 * H = -J d d + B d
 * v3, complexity per step is O(1)
 */

#include <iostream>
#include <cmath>

int dim[2];
int N;

int * spin = nullptr;

const int n_event_types = 10;

double rate[n_event_types];
int n_events[n_event_types];
int * event_pos[n_event_types];
int * grid_event_type;
int * grid_event_idx;

/// beta * J
double betaJ;
double betamuB;

/// adjust time scale
double freq = 1.0;

inline int idx2_to_idx1(int i, int j) {
    return dim[1] * i + j;
}

inline void idx1_to_idx2(int idx, int & i, int & j) {
    i = idx / dim[1];
    j = idx % dim[1];
}

int get_event_type(int i, int j) {
    // sum of neighbor spins
    // 0, 1, 2, 3, 4
    int sum = 0;
    sum += spin[idx2_to_idx1(i+1 < dim[0] ? i+1 : 0, j)];
    sum += spin[idx2_to_idx1(i-1 >= 0 ? i-1 : dim[0]-1, j)];
    sum += spin[idx2_to_idx1(i, j+1 < dim[1] ? j+1 : 0)];
    sum += spin[idx2_to_idx1(i, j-1 >= 0 ? j-1 : dim[1]-1)];
    // spin on site
    // 0, 1
    int self = spin[idx2_to_idx1(i, j)];
    // event type:
    // sum + self * 5
    int event_type = sum + self * 5;
    return event_type;
}

void remove_event(int idx) {
    // remove old entry
    int old_event_type = grid_event_type[idx];
    int old_event_idx = grid_event_idx[idx];
    int old_array_end = n_events[old_event_type] - 1;
    
    // fix replacing entry
    int replaced_grid = event_pos[old_event_type][old_array_end];
    grid_event_idx[replaced_grid] = old_event_idx;
    
    event_pos[old_event_type][old_event_idx] = replaced_grid;
    n_events[old_event_type] -= 1;
}

void insert_event(int event_type, int idx) {
    // add new entry
    int array_end = n_events[event_type];
    event_pos[event_type][array_end] = idx;
    grid_event_type[idx] = event_type;
    grid_event_idx[idx] = array_end;
    n_events[event_type] += 1;
}

void init_event_rates() {
    // init event rates
    for (int i = 0; i < n_event_types; i++) {
        int sum = i % 5;
        int self = i / 5;
        double betaE = -betaJ * (2 - 4 * self) * (2 * sum - 4) + betamuB * (2 - 4 * self);
        rate[i] = 0.5 * freq * (1 - tanh(0.5 * betaE));  // kinetic model is from eq 2.18 in https://link.springer.com/chapter/10.1007/978-3-662-06758-1_2
    }
}

void init_events() {
    // set event array length to 0
    for (int i = 0; i < n_event_types; i++) {
        n_events[i] = 0;
    }
    // init events
    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            // get event type
            int event_type = get_event_type(i, j);
            // add event
            insert_event(event_type, idx2_to_idx1(i, j));
        }
    }
}

void free_events() {
    // free memory
    for (int i = 0; i < n_event_types; i++) {
        delete [] event_pos[i];
    }
    delete [] grid_event_type;
    delete [] grid_event_idx;
}

void fix_events(int i, int j) {
    int idx1 = idx2_to_idx1(i, j);
    // remove old entry
    remove_event(idx1);
    // add new entry
    int event_type = get_event_type(i, j);
    insert_event(event_type, idx1);
}

double kmc() {
    double total_rate = 0.0;
    for (int i = 0; i < n_event_types; i++) {
        total_rate += rate[i] * n_events[i];
    }
    double p = (double)rand() / (double)RAND_MAX * total_rate;
    double rate_asum = 0.0;
    int chosen_event_type = n_event_types - 1;
    for (int i = 0; i < n_event_types; i++) {
        rate_asum += rate[i] * n_events[i];
        if (p <= rate_asum) {
            chosen_event_type = i;
            break;
        }
    }
    int event_idx = rand() % n_events[chosen_event_type];
    int idx = event_pos[chosen_event_type][event_idx];
    int i, j;
    idx1_to_idx2(idx, i, j);
    // flip
    spin[idx] = 1 - spin[idx];
    // fix events
    fix_events(i, j);
    fix_events(i+1 < dim[0] ? i+1 : 0, j);
    fix_events(i-1 >= 0 ? i-1 : dim[0]-1, j);
    fix_events(i, j+1 < dim[1] ? j+1 : 0);
    fix_events(i, j-1 >= 0 ? j-1 : dim[1]-1);
    // time step
    double u = (double)rand() / (double)RAND_MAX;
    double t = -log(u) / total_rate;
    return t;
}

void read_spin() {
    for (int k = 0; k < N; k++) {
        char c;
        std::cin >> c;
        if (c == '+') {
            spin[k] = 1;
        } else if (c == '-'){
            spin[k] = 0;
        } else {
            throw std::exception();
        }
    }
}

void print_spin() {
    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            int t = idx2_to_idx1(i, j);
            if (spin[t] == 1) {
                std::cout << '+';
            } else if (spin[t] == 0){
                std::cout << '-';
            }
        }
        std::cout << '\n';
    }
}

void debug() {
    std::cout << "spin:\n";
    print_spin();
    std::cout << "n events:\n";
    for (int i = 0; i < n_event_types; i++) {
        std::cout << n_events[i] << " ";
    }
    std::cout << "\n";
    std::cout << "event lists:\n";
    for (int i = 0; i < n_event_types; i++) {
        std::cout << "list " << i << ":\n";
        for (int k = 0; k < n_events[i]; k++) {
            std::cout << event_pos[i][k] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "event type:\n";
    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            std::cout << (char)(grid_event_type[idx2_to_idx1(i, j)] + '0');
        }
        std::cout << "\n";
    }
    std::cout << "event idx:\n";
    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            std::cout << grid_event_idx[idx2_to_idx1(i, j)] << " ";
        }
        std::cout << "\n";
    }
}

/// dylib API

#ifdef __cplusplus
extern "C" {
#endif

void set_rand_seed(unsigned int seed) {
    srand(seed);
}

void init_system(int in_dim0, int in_dim1) {
    dim[0] = in_dim0;
    dim[1] = in_dim1;
    N = dim[0] * dim[1];
    spin = new int[N];
    for (int i = 0; i < n_event_types; i++) {
        event_pos[i] = new int[N];
    }
    grid_event_type = new int[N];
    grid_event_idx = new int[N];
}

void set_param(double in_betaJ, double in_betamuB, double in_freq) {
    betaJ = in_betaJ;
    betamuB = in_betamuB;
    freq = in_freq;
    init_event_rates();
}

void set_spin(int * in_spin) {
    memcpy(spin, in_spin, N * sizeof(int));
    init_events();
}

void set_random_spin() {
    for (int i = 0; i < N; i++) {
        spin[i] = rand() & 1;
    }
    init_events();
}

void get_spin(int * out_spin) {
    memcpy(out_spin, spin, N * sizeof(int));
}

double next_frame(double t, double time_per_frame) {
    while (t < time_per_frame) {
        t += kmc();
    }
    return t - time_per_frame;
}

void free_system() {
    delete [] spin;
    for (int i = 0; i < n_event_types; i++) {
        delete [] event_pos[i];
    }
    delete [] grid_event_type;
    delete [] grid_event_idx;
}

#ifdef __cplusplus
}
#endif

/// end dylib API
