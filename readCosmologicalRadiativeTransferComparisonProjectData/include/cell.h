#ifndef CELL_H
#define CELL_H

#define NCELLS 128
#define CELL_WIDTH 1.f

#define get_array_index(i, j, k)    \
     i + NCELLS * j + (NCELLS * NCELLS) * k  \


#endif
