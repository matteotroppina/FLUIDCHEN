#include "UtilsGPU.hpp"
#include "Grid.hpp"

#ifdef __CUDACC__
void generate_gpu_masks(Grid &_grid, bool* &fluid_mask, uint8_t* &boundary_type, uint8_t* &border_position, bool * &d_fluid_mask, uint8_t * &d_boundary_type, uint8_t * &d_border_position){
    int size_linear = (_grid.domain().size_x + 2) * (_grid.domain().size_y + 2);
    fluid_mask = new bool[size_linear];
    boundary_type = new uint8_t[size_linear];
    border_position = new uint8_t[size_linear];

    for (int i = 0; i <= _grid.domain().size_x + 1; i++) {
        for (int j = 0; j <= _grid.domain().size_y + 1; j++) {
            int idx = i + j * (_grid.domain().size_x + 2);
            Cell cell = _grid.cell(i, j);
            if (cell.type() == cell_type::FLUID) {
                fluid_mask[idx] = 1;
                boundary_type[idx] = 0;
                border_position[idx] = 8; // 8 is a flag for no boundary
            } else {
                fluid_mask[idx] = 0;
                boundary_type[idx] = static_cast<uint8_t>(cell.type());

                if (cell.is_border(border_position::TOP)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::TOP); // 0
                }
                if (cell.is_border(border_position::BOTTOM)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::BOTTOM); // 1
                }
                if (cell.is_border(border_position::RIGHT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::RIGHT); // 2
                }
                if (cell.is_border(border_position::LEFT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::LEFT); // 3
                }

                // B_NW cell
                if (cell.is_border(border_position::TOP) && cell.is_border(border_position::LEFT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::TOP) + 4; // 4
                }
                // B_SE cell
                if (cell.is_border(border_position::BOTTOM) && cell.is_border(border_position::RIGHT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::BOTTOM) + 4; // 5
                }
                // B_NE cell
                if (cell.is_border(border_position::TOP) && cell.is_border(border_position::RIGHT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::RIGHT) + 4; // 6
                }
                // B_SW cell
                if (cell.is_border(border_position::BOTTOM) && cell.is_border(border_position::LEFT)) {
                    border_position[idx] = static_cast<uint8_t>(border_position::LEFT) + 4; // 7
                }

                if (cell.type() == cell_type::INNER_OBSTACLE) {
                    border_position[idx] = 8;
                }
            }
        }
    }
    // Copy data to GPU (only once)
    cudaMalloc(&d_fluid_mask, size_linear * sizeof(bool));
    cudaMemcpy(d_fluid_mask, fluid_mask, size_linear * sizeof(bool), cudaMemcpyHostToDevice);

    cudaMalloc(&d_boundary_type, size_linear * sizeof(uint8_t));
    cudaMemcpy(d_boundary_type, boundary_type, size_linear * sizeof(uint8_t), cudaMemcpyHostToDevice);

    cudaMalloc(&d_border_position, size_linear * sizeof(uint8_t));
    cudaMemcpy(d_border_position, border_position, size_linear * sizeof(uint8_t), cudaMemcpyHostToDevice);

    // Free memory
    delete[] fluid_mask;
    delete[] boundary_type;
    delete[] border_position;
}
// Allocating and copying data to GPU
// GPU has a different memory space
void allocate_gpu_memory(Grid &_grid, double * &d_p_matrix_new, double * &d_p_matrix, double * &d_rs_matrix) {
    int size_linear = (_grid.domain().size_x + 2) * (_grid.domain().size_y + 2);

    // Allocate memory on GPU
    cudaMalloc(&d_p_matrix_new, size_linear * sizeof(double));
    cudaMalloc(&d_p_matrix, size_linear * sizeof(double));
    cudaMalloc(&d_rs_matrix, size_linear * sizeof(double));
}

void free_gpu_memory(double * &d_p_matrix_new, double * &d_p_matrix, double * &d_rs_matrix, bool * &d_fluid_mask, uint8_t * &d_boundary_type, uint8_t * &d_border_position) {
    cudaFree(d_p_matrix_new);
    cudaFree(d_p_matrix);
    cudaFree(d_rs_matrix);
    cudaFree(d_fluid_mask);
    cudaFree(d_boundary_type);
}
#endif



