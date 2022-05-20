# A Kinetic Ising Model (KMC simulation)

A kinetic Monte Carlo simulator for a kinetic Ising model

## Simulation System

Ising model on 2D square lattice.

Ising model is not kinetic, but you can assign an arbitrary transition rate as long as it conforms to the detailed balance condition.

Transition rate:

W = 0.5 \* lambda \* (1 - tanh(0.5 \* beta \* dE))

Transition rate is from eq 2.18 in https://link.springer.com/chapter/10.1007/978-3-662-06758-1_2 .

Complexity for each KMC step is O(1).

## Usage

Compile `ising.cpp` to a dynamic-link library.

Example shell commands using gcc for unix-like system:

```shell
g++ -O2 -shared ising.cpp -o libising.so
```

Run `render.py` to render video.

```shell
python3 render.py
```

`render.py` will load libising as a ctypes lib, so compile `ising.cpp` before running the script. opencv is required for video rendering.
