### Description

This code simulates the Fractional-Step (FS) algorithm described in [1] and [2] on the Ziff–Gulari–Barshad model. Although the algorithm is parallel, the code is not and it runs all blocks sequentially. It was used as a proof of concept.



This video is a simulation on a 200x200 lattice with 4x4 blocks.

![ZGB](./figures/zgb-fs.gif)

The next figure is a comparison of the FS algorithm with the Gillespie (or SSA) algorithm on a 10x10 lattice with 2x2 blocks and dt=0.01.

![comparison](./figures/comparison-10.jpg)

The next figure is a comparison of the FS algorithm with the Gillespie (or SSA) algorithm on a 100x100 lattice with 2x2 block and dt=0.1.

![comparison](./figures/comparison-100.jpg)
