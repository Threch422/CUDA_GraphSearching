#include <iostream>
#include "helpers.h"

/* you can define data structures and helper functions here */

// counting the number of NO2 group
int no2_counting(int* n_o, int n_o_size) {	
	// Find number of NO2
	int no2_size = 0;
	int no2_count = 1;
	for (int i = 0; i < n_o_size; i++) {
		for (int j = i + 1; j < n_o_size; j++) {
			if (n_o[i] == n_o[j])
				no2_count++;			
		}
		
		if (no2_count == 2) {
			no2_size++;
		}
		no2_count = 1;
		
	}
	return no2_size;

}

// Trimming the invalid C-N edges, change to NO2 array
void no2_making(int* source, int* result, int n_o_size, int no2_size){
	// Trimming empty cell
	int j = 0;
	for (int i = 0; i < n_o_size; i++){
		if (source[i] != -1) {
			result[j] = source[i];
			result[j + no2_size] = source[i + n_o_size];
			result[j + 2 * no2_size]= source[i + 2 * n_o_size];
			result[j + 3 * no2_size] = source[i + 3 * n_o_size];
			j++;
		}
	} 
}

// Counting the number of C6RING
int c6_ring_count(int* source, int no2_size) {
	int counter = 0;
	for (int i = 0; i < no2_size; i++) {
		int inner_count = 0;
		for (int j = 0; j < 6; j++) {
			if (source[i + j * no2_size] != -1)
				inner_count++;
			else
				break;
		}
		if (inner_count == 6)
			counter++;
	}
	return counter;
}

// Trimming non-completed "c6chain", and store the reverse order of c6ring
void c6_ring_tidying(int* source, int* result, int no2_size, int c6ring_num) {
	int j = 0;
	for (int i = 0; i < no2_size; i++) {
		if (source[i + 2 * no2_size] != -1 && source[i + 3 * no2_size] != -1 && source[i + 4 * no2_size] != -1) {
				result[j] = source[i];
				result[j + c6ring_num] = source[i + no2_size];
				result[j + 2 * c6ring_num] = source[i + 2 * no2_size];
				result[j + 3 * c6ring_num] = source[i + 3 * no2_size];
				result[j + 4 * c6ring_num] = source[i + 4 * no2_size];
				result[j + 5 * c6ring_num] = source[i + 5 * no2_size];
				j++;
				result[j] = source[i];
				result[j + c6ring_num] = source[i + 5 * no2_size];
				result[j + 2 * c6ring_num] = source[i + 4 * no2_size];
				result[j + 3 * c6ring_num] = source[i + 3 * no2_size];
				result[j + 4 * c6ring_num] = source[i + 2 * no2_size];
				result[j + 5 * c6ring_num] = source[i + no2_size];
				j++;
		}
	}
}

// Kernel function, searching no2 group from C-N edges and N-O edges
__global__ void search_no2(int* d_c_n, int* d_n_o, int* result, int c_n_size, int n_o_size, int no2_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < n_o_size; i+= num_threads) {
			int o1 = d_n_o[i + n_o_size];
			int o2 = -1;
			int nitrogen_count = 1;

		for (int j = i + 1; j < n_o_size; j++) {
			if (d_n_o[i] == d_n_o[j]) {
				nitrogen_count++;
				o2 = d_n_o[j + n_o_size];
			}
		}
		if (nitrogen_count == 2) {
			for (int k = 0; k < c_n_size; k++) {
				if (d_c_n[k + c_n_size] == d_n_o[i]) {
					result[i] = d_c_n[k];
					result[i + n_o_size] = d_n_o[i];
					result[i + 2 * n_o_size] = o1;
					result[i + 3 * n_o_size] = o2;
				}
			}
		}
	}


}

// Getting nearby atom bounded with C atom in C_NO2
__global__ void search_c6ring(int* d_c_c, int* d_c_no2, int* d_c_no2_temp, int c_c_size, int no2_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < no2_size; i += num_threads) {
		d_c_no2_temp[i] = d_c_no2[i];
		for (int j = 0; j < c_c_size; j++) {
			if (d_c_no2[i] == d_c_c[j]) {
				if (d_c_no2_temp[i + no2_size] == -1)
					d_c_no2_temp[i + no2_size] = d_c_c[j + c_c_size];
				else
					d_c_no2_temp[i + 2 * no2_size] = d_c_c[j + c_c_size];
			}
		}
		d_c_no2_temp[i + 3 * no2_size] = d_c_no2[i + no2_size];
		d_c_no2_temp[i + 4 * no2_size] = d_c_no2[i + 2 * no2_size];
		d_c_no2_temp[i + 5 * no2_size] = d_c_no2[i + 3 * no2_size];
	}

}

// Mapping between different groups of C_NO2, extracting C6ring which is belongs to TNT structure
__global__ void search_c6ring2(int* d_c_no2_temp, int* d_c_ring_detect, int no2_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < no2_size; i += num_threads) {
		d_c_ring_detect[i] = d_c_no2_temp[i];
		d_c_ring_detect[i + no2_size] = d_c_no2_temp[i + no2_size];
		d_c_ring_detect[i + 5 * no2_size] = d_c_no2_temp[i + 2 * no2_size];
		int c1 = -1; int c2 = -1; int c3 = -1; int c4 = -1;

		for (int j = 0; j < no2_size; j++) {
			if (j == i)
				continue;
			if (d_c_no2_temp[i + no2_size] == d_c_no2_temp[j + no2_size]) {
				c1 = d_c_no2_temp[j];
				c2 = d_c_no2_temp[j + 2 * no2_size];
			}
			else if (d_c_no2_temp[i + no2_size] == d_c_no2_temp[j + 2 * no2_size]) {
				c1 = d_c_no2_temp[j];
				c2 = d_c_no2_temp[j + no2_size];
			}

			if (d_c_no2_temp[i + 2 * no2_size] == d_c_no2_temp[j + no2_size]) {
				c3 = d_c_no2_temp[j];
				c4 = d_c_no2_temp[j + 2 * no2_size];
			}
			else if (d_c_no2_temp[i + 2 * no2_size] == d_c_no2_temp[j + 2 * no2_size]) {
				c3 = d_c_no2_temp[j];
				c4 = d_c_no2_temp[j + no2_size];
			}
		}
		if (c1 == -1 || c2 == -1 || c3 == -1 || c4 == -1 || c1 == d_c_ring_detect[i] || c2 == d_c_ring_detect[i] || c3 == d_c_ring_detect[i] || c4 == d_c_ring_detect[i])
			continue;

		if (c1 == c3) {
			d_c_ring_detect[i + 2 * no2_size] = c1;
			d_c_ring_detect[i + 3 * no2_size] = c2;
			d_c_ring_detect[i + 4 * no2_size] = c4;
		}
		if (c1 == c4) {
			d_c_ring_detect[i + 2 * no2_size] = c1;
			d_c_ring_detect[i + 3 * no2_size] = c2;
			d_c_ring_detect[i + 4 * no2_size] = c3;
		}
		if (c2 == c3) {
			d_c_ring_detect[i + 2 * no2_size] = c1;
			d_c_ring_detect[i + 3 * no2_size] = c2;
			d_c_ring_detect[i + 4 * no2_size] = c4;
		}
		if (c2 == c4) {
			d_c_ring_detect[i + 2 * no2_size] = c1;
			d_c_ring_detect[i + 3 * no2_size] = c2;
			d_c_ring_detect[i + 4 * no2_size] = c3;
		}
	}
}

// Copy the c6ring for 8 times for later mapping :D
__global__ void tnt_c6ring_copying(int* d_c6ring, int* d_tnt_c6ring, int c6ring_num) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < c6ring_num / 8; i += num_threads) {
		for (int j = 0; j < 8; j++) {
			d_tnt_c6ring[i * 8 + j] = d_c6ring[i];
			d_tnt_c6ring[i * 8 + j + c6ring_num] = d_c6ring[i + c6ring_num / 8];
			d_tnt_c6ring[i * 8 + j + 2 * c6ring_num] = d_c6ring[i + 2 * c6ring_num / 8];
			d_tnt_c6ring[i * 8 + j + 3 * c6ring_num] = d_c6ring[i + 3 * c6ring_num / 8];
			d_tnt_c6ring[i * 8 + j + 4 * c6ring_num] = d_c6ring[i + 4 * c6ring_num / 8];
			d_tnt_c6ring[i * 8 + j + 5 * c6ring_num] = d_c6ring[i + 5 * c6ring_num / 8];
		}
	}
}

// Here comes mapping of tnt
__global__ void tnt_mapping(int* d_tnt, int* d_c_no2, int* d_tnt_c6ring, int tnt_num, int no2_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;
	// Copy the tnt_c6ring to d_tnt first
	for (int i = tid; i < tnt_num; i += num_threads) {
		d_tnt[i] = d_tnt_c6ring[i];
		d_tnt[i + tnt_num] = d_tnt_c6ring[i + tnt_num];
		d_tnt[i + 2 * tnt_num] = d_tnt_c6ring[i + 2 * tnt_num];
		d_tnt[i + 3 * tnt_num] = d_tnt_c6ring[i + 3 * tnt_num];
		d_tnt[i + 4 * tnt_num] = d_tnt_c6ring[i + 4 * tnt_num];
		d_tnt[i + 5 * tnt_num] = d_tnt_c6ring[i + 5 * tnt_num];

		// mapping is only required on O atom as we have already done N mapping in constructing C6RING
		/* Assuming we have 6 O atom, which are 9 10 11 12 13 14
		*  map_id = 0 -->  9 10 11 12 13 14
		*  map_id = 1 -->  9 10 11 12 14 13
		*  map_id = 2 -->  9 10 12 11 13 14
		*  map_id = 3 -->  9 10 12 11 14 13
		*  map_id = 4 -->  10 9 11 12 13 14
		*  map_id = 5 -->  10 9 11 12 14 13
		*  map_id = 6 -->  10 9 12 11 13 14
		*  map_id = 7 -->  10 9 12 11 14 13
		*/
		int map_id = tid % 8;
		int o1, o2, o3, o4, o5, o6;
		// Mapping N atom to d_tnt
		for (int j = 0; j < no2_size; j++) {
			if (d_tnt[i] == d_c_no2[j]) {
				d_tnt[i + 6 * tnt_num] = d_c_no2[j + no2_size];
				o1 = d_c_no2[j + 2 * no2_size];
				o2 = d_c_no2[j + 3 * no2_size];
			}

			if (d_tnt[i + 2 * tnt_num] == d_c_no2[j]) {
				d_tnt[i + 7 * tnt_num] = d_c_no2[j + no2_size];
				o3 = d_c_no2[j + 2 * no2_size];
				o4 = d_c_no2[j + 3 * no2_size];
			}

			if (d_tnt[i + 4 * tnt_num] == d_c_no2[j]) {
				d_tnt[i + 8 * tnt_num] = d_c_no2[j + no2_size];
				o5 = d_c_no2[j + 2 * no2_size];
				o6 = d_c_no2[j + 3 * no2_size];
			}
		}
		switch(map_id) {
			case 0:
				d_tnt[i + 9 * tnt_num] = o1;
				d_tnt[i + 10 * tnt_num] = o2;
				d_tnt[i + 11 * tnt_num] = o3;
				d_tnt[i + 12 * tnt_num] = o4;
				d_tnt[i + 13 * tnt_num] = o5;
				d_tnt[i + 14 * tnt_num] = o6;
				break;
			case 1:
				d_tnt[i + 9 * tnt_num] = o1;
				d_tnt[i + 10 * tnt_num] = o2;
				d_tnt[i + 11 * tnt_num] = o3;
				d_tnt[i + 12 * tnt_num] = o4;
				d_tnt[i + 13 * tnt_num] = o6;
				d_tnt[i + 14 * tnt_num] = o5;
				break;
			case 2:
				d_tnt[i + 9 * tnt_num] = o1;
				d_tnt[i + 10 * tnt_num] = o2;
				d_tnt[i + 11 * tnt_num] = o4;
				d_tnt[i + 12 * tnt_num] = o3;
				d_tnt[i + 13 * tnt_num] = o5;
				d_tnt[i + 14 * tnt_num] = o6;
				break;
			case 3:
				d_tnt[i + 9 * tnt_num] = o1;
				d_tnt[i + 10 * tnt_num] = o2;
				d_tnt[i + 11 * tnt_num] = o4;
				d_tnt[i + 12 * tnt_num] = o3;
				d_tnt[i + 13 * tnt_num] = o6;
				d_tnt[i + 14 * tnt_num] = o5;
				break;
			case 4:
				d_tnt[i + 9 * tnt_num] = o2;
				d_tnt[i + 10 * tnt_num] = o1;
				d_tnt[i + 11 * tnt_num] = o3;
				d_tnt[i + 12 * tnt_num] = o4;
				d_tnt[i + 13 * tnt_num] = o5;
				d_tnt[i + 14 * tnt_num] = o6;
				break;
			case 5:
				d_tnt[i + 9 * tnt_num] = o2;
				d_tnt[i + 10 * tnt_num] = o1;
				d_tnt[i + 11 * tnt_num] = o3;
				d_tnt[i + 12 * tnt_num] = o4;
				d_tnt[i + 13 * tnt_num] = o6;
				d_tnt[i + 14 * tnt_num] = o5;
				break;
			case 6:
				d_tnt[i + 9 * tnt_num] = o2;
				d_tnt[i + 10 * tnt_num] = o1;
				d_tnt[i + 11 * tnt_num] = o4;
				d_tnt[i + 12 * tnt_num] = o3;
				d_tnt[i + 13 * tnt_num] = o5;
				d_tnt[i + 14 * tnt_num] = o6;
				break;
			case 7:
				d_tnt[i + 9 * tnt_num] = o2;
				d_tnt[i + 10 * tnt_num] = o1;
				d_tnt[i + 11 * tnt_num] = o4;
				d_tnt[i + 12 * tnt_num] = o3;
				d_tnt[i + 13 * tnt_num] = o6;
				d_tnt[i + 14 * tnt_num] = o5;
				break;								
		}
	}
}

/**
 * please remember to set final_results and final_result_size 
 * before return.
 */
void tnt_counting(int num_blocks_per_grid, int num_threads_per_block,
        int* c_c, int* c_n, int* c_h, int* n_o,
        int c_c_size, int c_n_size, int c_h_size, int n_o_size,
        int* &final_results, int &final_result_size) {
	
	// Memory allocation for GPU
	int* d_c_c, *d_c_n, *d_c_h, *d_n_o;
	cudaMalloc((void **) &d_c_c, 2 * c_c_size * sizeof(int));
	cudaMalloc((void **) &d_c_n, 2 * c_n_size * sizeof(int));
	cudaMalloc((void **) &d_c_h, 2 * c_h_size * sizeof(int));
	cudaMalloc((void **) &d_n_o, 2 * n_o_size * sizeof(int));

	// Copy the chem edges to GPU
	cudaMemcpy(d_c_c, c_c, 2 * c_c_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c_n, c_n, 2 * c_n_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c_h, c_h, 2 * c_h_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n_o, n_o, 2 * n_o_size * sizeof(int), cudaMemcpyHostToDevice);
	
	// Getting the number of NO2 group by using N-O edges
	int no2_size = no2_counting(n_o, n_o_size);

	// Prep for marking N atom position and get their O
	int* no2_arr_temp, *d_no2;
	no2_arr_temp = (int*)malloc(4 * n_o_size * sizeof(int));
	cudaMalloc((void **) &d_no2, 4 * n_o_size * sizeof(int));
	
	for (int i = 0; i < 4 * n_o_size; i++)
		no2_arr_temp[i] = -1;
	
	cudaMemcpy(d_no2, no2_arr_temp,  4 * n_o_size * sizeof(int), cudaMemcpyHostToDevice);

	search_no2<<<num_blocks_per_grid, num_threads_per_block>>>(d_c_n, d_n_o, d_no2, c_n_size, n_o_size, no2_size);

	cudaMemcpy(no2_arr_temp,  d_no2, 4 * n_o_size * sizeof(int), cudaMemcpyDeviceToHost);

	int* c_no2 = (int*)malloc(4 * no2_size * sizeof(int));

	no2_making(no2_arr_temp, c_no2, n_o_size, no2_size);
	// Finish making C_NO2

	// Putting neighbour atom of C atom in C_NO2 group
	int* c_no2_temp, *d_c_no2, *d_c_no2_temp;
	c_no2_temp = (int*)malloc(6 * no2_size * sizeof(int));
	for (int i = 0; i < 6 * no2_size; i++)
		c_no2_temp[i] = -1;
	cudaMalloc((void **) &d_c_no2, 4 * no2_size * sizeof(int));
	cudaMalloc((void **) &d_c_no2_temp, 6 * no2_size * sizeof(int));
	cudaMemcpy(d_c_no2, c_no2,  4 * no2_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c_no2_temp, c_no2_temp,  6 * no2_size * sizeof(int), cudaMemcpyHostToDevice);

	search_c6ring<<<num_blocks_per_grid, num_threads_per_block>>>(d_c_c, d_c_no2, d_c_no2_temp, c_c_size, no2_size);
	cudaMemcpy(c_no2_temp, d_c_no2_temp,  6 * no2_size * sizeof(int), cudaMemcpyDeviceToHost);

	// Mapping C atom in the connected C atom with those in C_NO2 group, finding valid tnt_C6ring
	int* c_ring_detect, *d_c_ring_detect;
	c_ring_detect = (int*)malloc(6 * no2_size * sizeof(int));
	for (int i = 0; i < 6 * no2_size; i++)
		c_ring_detect[i] = -1;
	cudaMalloc((void **) &d_c_ring_detect, 6 * no2_size * sizeof(int));
	cudaMemcpy(d_c_ring_detect, c_ring_detect, 6 * no2_size * sizeof(int), cudaMemcpyHostToDevice);

	search_c6ring2<<<num_blocks_per_grid, num_threads_per_block>>>(d_c_no2_temp, d_c_ring_detect, no2_size);
	cudaMemcpy(c_ring_detect, d_c_ring_detect, 6 * no2_size * sizeof(int), cudaMemcpyDeviceToHost);

	// Cleaning invalid C6-chain, extract valid c6ring only
	int c6ring_num = c6_ring_count(c_ring_detect, no2_size) * 2;
	int* c6ring = (int*)malloc(6 * c6ring_num * sizeof(int));
	c6_ring_tidying(c_ring_detect, c6ring, no2_size, c6ring_num);

	// TNT structre generation and mapping
	int* d_c6ring;
	cudaMalloc((void **) &d_c6ring, 6 * c6ring_num * sizeof(int));
	cudaMemcpy(d_c6ring, c6ring, 6 * c6ring_num * sizeof(int), cudaMemcpyHostToDevice);
	int tnt_num = c6ring_num * 8;
	int* tnt_c6ring = (int*)malloc(6 * tnt_num * sizeof(int));
	int* d_tnt_c6ring;
	cudaMalloc((void **) &d_tnt_c6ring, 6 * tnt_num * sizeof(int));
	tnt_c6ring_copying<<<num_blocks_per_grid, num_threads_per_block>>>(d_c6ring, d_tnt_c6ring, tnt_num);
	cudaMemcpy(tnt_c6ring, d_tnt_c6ring, 6 * tnt_num * sizeof(int), cudaMemcpyDeviceToHost);

	int* tnt = (int*)malloc(15 * tnt_num * sizeof(int));
	int* d_tnt;
	cudaMalloc((void**) &d_tnt, 15 * tnt_num * sizeof(int));
	tnt_mapping<<<num_blocks_per_grid, num_threads_per_block>>>(d_tnt, d_c_no2, d_tnt_c6ring, tnt_num, no2_size);
	cudaMemcpy(tnt, d_tnt, 15 * tnt_num * sizeof(int), cudaMemcpyDeviceToHost);





	// Result
	final_results = tnt;
	final_result_size = tnt_num;

	// Free the memory used to prevent memory leak
	free(no2_arr_temp);
	free(tnt_c6ring);
	free(c6ring);
	free(c_ring_detect);
	free(c_no2);
	free(c_no2_temp);

	cudaFree(d_c6ring);
	cudaFree(d_tnt);
	cudaFree(d_tnt_c6ring);
	cudaFree(d_c_ring_detect);
	cudaFree(d_c_no2_temp);
	cudaFree(d_no2);
	cudaFree(d_c_no2);
	cudaFree(d_c_c);
	cudaFree(d_c_h);
	cudaFree(d_c_n);
	cudaFree(d_n_o);     
	

}
