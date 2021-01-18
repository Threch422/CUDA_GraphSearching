#include <iostream>
#include <vector>
#include <algorithm>
#include "helpers.h"

/* you can define data structures and helper functions here */

__global__ void c6ring1(int* d_c_c, int* result, int c_c_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;


	for (int i = tid; i < c_c_size; i += num_threads) {
		int c1, c2, c3, c4, c5;
		bool found = false;
		c1 = d_c_c[i + c_c_size];
		for (int a = 0; a < c_c_size; a++) {
			if (d_c_c[a] == c1) {
				c2 = d_c_c[a + c_c_size];
				if (c2 != d_c_c[i]) {
					for (int b = 0; b < c_c_size; b++) {
						if (d_c_c[b] == c2) {
							c3 = d_c_c[b + c_c_size];
							if (c3 != c1 && c3 != d_c_c[i]) {
								for (int c = 0; c < c_c_size; c++) {
									if (d_c_c[c] == c3) {
										c4 = d_c_c[c + c_c_size];
										if (c4 != c2 && c4 != c1 && c4 != d_c_c[i]) {
											for (int d = 0; d < c_c_size; d++) {
												if (c4 == d_c_c[d]) {
													c5 = d_c_c[d + c_c_size];
													if (c5 != c3 && c5 != c2 && c5 != c1 && c5 != d_c_c[i]) {
														for (int e = 0; e < c_c_size; e++) {
															if (c5 == d_c_c[e] && d_c_c[e + c_c_size] == d_c_c[i]) {
																	result[i] = d_c_c[i];
																	result[i + c_c_size] = c1;
																	result[i + 2 * c_c_size] = c2;
																	result[i + 3 * c_c_size] = c3;
																	result[i + 4 * c_c_size] = c4;
																	result[i + 5 * c_c_size] = c5;
																	found = true;
																	break;
															}
														}
														if (found)
															break;
													}
												}
											}
											if (found)
												break;
										}
									}
								}
								if (found)
									break;
							}
						}
					}
					if (found)
						break;
				}
			}
		}
	}
}

__global__ void cn_elim(int* d_c6elim, int* d_c_n, int* d_n_o, int c6size, int c_n_size, int n_o_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < c6size; i += num_threads) {
		for (int a = 0; a < c_n_size; a++) {
			if (d_c6elim[i] == d_c_n[a]) {
				int count1 = 0;
				for (int b = 0; b < n_o_size; b++) {
					if (d_c_n[a + c_n_size] == d_n_o[b]) {
						count1++;
					}
				}
				if (count1 < 2) {
					d_c_n[a] = -1;
					d_c_n[a + c_n_size] = -1;
				}
			}

			if (d_c6elim[i + 2 * c6size] == d_c_n[a]) {
				int count1 = 0;
				for (int b = 0; b < n_o_size; b++) {
					if (d_c_n[a + c_n_size] == d_n_o[b]) {
						count1++;
					}
				}
				if (count1 < 2) {
					d_c_n[a] = -1;
					d_c_n[a + c_n_size] = -1;
				}
			}

			if (d_c6elim[i + 4 * c6size] == d_c_n[a]) {
				int count1 = 0;
				for (int b = 0; b < n_o_size; b++) {
					if (d_c_n[a + c_n_size] == d_n_o[b]) {
						count1++;
					}
				}
				if (count1 < 2) {
					d_c_n[a] = -1;
					d_c_n[a + c_n_size] = -1;
				}
			}			
		}
	}
} 

__global__ void c6ring2(int* d_c6ring_1, int* d_c_n, int* d_n_o, int c_n_size, int n_o_size, int c_c_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < c_c_size; i += num_threads) {
		if (d_c6ring_1[i] != -1) {
			int count0 = 0;
			int count1 = 0;
			int count2 = 0;			
				
			for (int a = 0; a < c_n_size; a++) {
				if (d_c6ring_1[i] == d_c_n[a]) {
					for (int b = 0; b < n_o_size; b++) {
						if (d_c_n[a + c_n_size] == d_n_o[b]) {
							count0++;
						}
					}					
				}							
				if (d_c6ring_1[i + 2 * c_c_size] == d_c_n[a]) {
					for (int b = 0; b < n_o_size; b++) {
						if (d_c_n[a + c_n_size] == d_n_o[b]) {
							count1++;
						}
					}					
				}
				if (d_c6ring_1[i + 4 * c_c_size] == d_c_n[a]) {
					for (int b = 0; b < n_o_size; b++) {
						if (d_c_n[a + c_n_size] == d_n_o[b]) {
							count2++;
						}
					}					
				}			
			}
			if (count0 < 2 || count1 < 2 || count2 < 2) {
				d_c6ring_1[i] = -1;
				d_c6ring_1[i + c_c_size] = -1;
				d_c6ring_1[i + 2 * c_c_size] = -1;
				d_c6ring_1[i + 3 * c_c_size] = -1;
				d_c6ring_1[i + 4 * c_c_size] = -1;
				d_c6ring_1[i + 5 * c_c_size] = -1;
			}				
		}
	}
}

int c6ring_count(int* c6ring, int c_c_size) {
	int count = 0;
	for (int i = 0; i < c_c_size; i++) {
		if (c6ring[i] != -1)
			count++;
	}
	return count;
}

void c6ring_elim(int* source, int* result, int c_c_size, int count) {
	int j = 0;
	for (int i = 0; i < c_c_size; i++) {
		if (source[i] != -1) {
			result[j] = source[i];
			result[j + count] = source[i + c_c_size];
			result[j + 2 * count] = source[i + 2 * c_c_size];
			result[j + 3 * count] = source[i + 3 * c_c_size];
			result[j + 4 * count] = source[i + 4 * c_c_size];
			result[j + 5 * count] = source[i + 5 * c_c_size];
			j++;
		}
	}
}

void prefix(int* source, int c6size, int* result) {
	result[0] = 0;
	for (int i = 1; i <= c6size; i++) {
		result[i] = source[i - 1] + result[i - 1];
	}
}

int total_mapping(int* source, int c6size) {
	int count = 0;
	for (int i = 0; i < c6size; i++) {
		count += source[i];
	}
	return count;
}

__global__ void count_mapping(int* source, int* d_c_n, int* d_n_o, int c_n_size, int n_o_size, int c6size, int* result) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < c6size; i += num_threads) {
		int count0 = 0;
		for (int a = 0; a < c_n_size; a++) {
			int count1 = 0; int count2 = 0; int count3 = 0;
			if (source[i] == d_c_n[a]) {
				for (int b = 0; b < n_o_size; b++) {
					if (d_c_n[a + c_n_size] == d_n_o[b]) {
						count1++;
					}
				}
				for (int b = 0; b < c_n_size; b++) {
					if (source[i + 2 * c6size] == d_c_n[b]) {
						for (int c = 0; c < n_o_size; c++) {
							if (d_c_n[b + c_n_size] == d_n_o[c]) {
								count2++;
							}
						}
						for (int c = 0; c < c_n_size; c++) {
							if (source[i + 4 * c6size] == d_c_n[c]) {
								for (int d = 0; d < n_o_size; d++) {
									if (d_c_n[c + c_n_size] == d_n_o[d]) {
										count3++;
									}
								}
								count0 = count0 + (count1 * (count1 - 1)) * (count2 * (count2 - 1)) * (count3 * (count3 - 1));
								count3 = 0;
							}
						}
						count2 = 0;
					}
				}
				count1 = 0;
			}
		}
		result[i] = count0;
	}
}
/*
__global__ void c_n_ext(int* source, int* d_n_o, int s_size, int n_o_size, int* result) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < s_size; i += num_threads) {
		for ()
	}

}
*/

void n_ext(int* c_n, std::vector<int> &source, int c_n_size, int s_size, std::vector<int> &result) {
	for (int i = 0; i < s_size; i++) {
		for (int j = 0; j < c_n_size; j++) {
			if (source[i] == c_n[j]) {
				if (std::find(result.begin(), result.end(), c_n[j + c_n_size]) != result.end()) {
					continue;
				}
				else {
					result.push_back(c_n[j + c_n_size]);
				}
			}

		}
	}
}

__global__ void n_prefix(int* n_ex, int* d_n_o, int* result, int n_size, int n_o_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	int count = 0;
	for (int i = tid; i < n_size; i += num_threads) {
		for (int j = 0; j < n_o_size; j++) {
			if (n_ex[i] == d_n_o[j]) {
				count++;
			}
		}
		count = count * (count - 1);
		result[i] = count;
	}
}

// Getting all NO2 permutation
__global__ void no2_permutation(int* d_n_arr, int* d_n_o, int* d_n_prefix, int n_arr_size, int n_o_size, int* result, int result_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < n_arr_size; i += num_threads) {
		int idx = d_n_prefix[i];
		int idx1 = -1;
		while(idx < d_n_prefix[i+1]) {
			int count = 0; int o1 = -1; int o2 = -1;
			for (int a = 0; a < n_o_size; a++) {
				if (d_n_arr[i] == d_n_o[a]) {
					if (count == 0 && o1 == -1 && a > idx1) {
						o1 = d_n_o[a + n_o_size];
						idx1 = a;
						count++;
					}
					else if (count == 1 && o2 == -1) {
						o2 = d_n_o[a + n_o_size];
						count++;
					}
					if (count == 2 && o1 != -1 && o2 != -1) {
						result[idx] = d_n_arr[i];
						result[idx + result_size] = o1;
						result[idx + 2 * result_size] = o2;
						idx++;
						result[idx] = d_n_arr[i];
						result[idx + result_size] = o2;
						result[idx + 2 * result_size] = o1;
						idx++;
						o2 = -1;
						count--;
					}
				}
			}
		}
	}
}

__global__ void tnt_mapping(int* d_no2, int* c6ring, int* d_c_n, int* d_map_prefix, int no2_size, int c6size, int c_n_size, int* result, int result_size) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;

	for (int i = tid; i < c6size; i += num_threads) {
		int idx = d_map_prefix[i];
		while (idx < d_map_prefix[i + 1]) {
			int o1 = -1; int o2 = -1; int o3 = -1; int o4 = -1; int o5 = -1; int o6 = -1;
			for (int a = 0; a < c_n_size; a++) {
				if (c6ring[i] == d_c_n[a]) {
					for (int b = 0; b < no2_size; b++) {
						if (d_c_n[a + c_n_size] == d_no2[b]) {
							o1 = d_no2[b + no2_size];
							o2 = d_no2[b + 2 * no2_size];
							for (int c = 0; c < c_n_size; c++) {
								if (c6ring[i + 2 * c6size] == d_c_n[c]) {
									for (int d = 0; d < no2_size; d++) {
										if (d_c_n[c + c_n_size] == d_no2[d]) {
											o3 = d_no2[d + no2_size];
											o4 = d_no2[d + 2 * no2_size];
											for (int e = 0; e < c_n_size; e++) {
												if (c6ring[i + 4 * c6size] == d_c_n[e]) {
													for (int f = 0; f < no2_size; f++) {
														if (d_c_n[e + c_n_size] == d_no2[f]) {
															o5 = d_no2[f + no2_size];
															o6 = d_no2[f + 2 * no2_size];

															result[idx] = c6ring[i];
															result[idx + result_size] = c6ring[i + c6size];
															result[idx + 2 * result_size] = c6ring[i + 2 * c6size];
															result[idx + 3 * result_size] = c6ring[i + 3 * c6size];
															result[idx + 4 * result_size] = c6ring[i + 4 * c6size];
															result[idx + 5 * result_size] = c6ring[i + 5 * c6size];	
															result[idx + 6 * result_size] = d_no2[b];
															result[idx + 7 * result_size] = d_no2[d];
															result[idx + 8 * result_size] = d_no2[f];
															result[idx + 9 * result_size] = o1;
															result[idx + 10 * result_size] = o2;
															result[idx + 11 * result_size] = o3;
															result[idx + 12 * result_size] = o4;																																																																																																							
															result[idx + 13 * result_size] = o5;
															result[idx + 14 * result_size] = o6;
															idx++;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
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
	
	// Finding c6ring valid in TNT structure form C-C edges
	int* c6ring_1 = (int*)malloc(6 * c_c_size * sizeof(int));
	for (int i = 0; i < 6 * c_c_size; i++) {
		c6ring_1[i] = -1;
	}
	int* d_c6ring_1;
	cudaMalloc((void **) &d_c6ring_1, 6 * c_c_size * sizeof(int));
	cudaMemcpy(d_c6ring_1, c6ring_1, 6 * c_c_size * sizeof(int), cudaMemcpyHostToDevice);
	c6ring1<<<num_blocks_per_grid, num_threads_per_block>>>(d_c_c, d_c6ring_1, c_c_size);
	c6ring2<<<num_blocks_per_grid, num_threads_per_block>>>(d_c6ring_1, d_c_n, d_n_o, c_n_size, n_o_size, c_c_size);
	cudaMemcpy(c6ring_1, d_c6ring_1, 6 * c_c_size * sizeof(int), cudaMemcpyDeviceToHost);
	int c6count = c6ring_count(c6ring_1, c_c_size);
	int* c6elim = (int*)malloc(6 * c6count * sizeof(int));
	c6ring_elim(c6ring_1, c6elim, c_c_size, c6count);

	// Elimiate those invalid C6-Chain
	int* d_c6elim;
	cudaMalloc((void **) &d_c6elim, 6 * c6count * sizeof(int));
	cudaMemcpy(d_c6elim, c6elim, 6 * c6count * sizeof(int), cudaMemcpyHostToDevice);
	cn_elim<<<num_blocks_per_grid, num_threads_per_block>>>(d_c6elim, d_c_n, d_n_o, c6count, c_n_size, n_o_size);
	cudaMemcpy(c_n, d_c_n, 2 * c_n_size * sizeof(int), cudaMemcpyDeviceToHost);

	// Do the mapping counting of each valid TNT C6ring
	int* mapping = (int*)malloc(c6count * sizeof(int));
	int* d_mapping;
	cudaMalloc((void **) &d_mapping, c6count * sizeof(int));
	//cudaMemcpy(d_mapping, mapping, c6count * sizeof(int), cudaMemcpyHostToDevice);
	count_mapping<<<num_blocks_per_grid, num_threads_per_block>>>(d_c6elim, d_c_n, d_n_o, c_n_size, n_o_size, c6count, d_mapping);
	cudaMemcpy(mapping, d_mapping, c6count * sizeof(int), cudaMemcpyDeviceToHost);
	int* map_prefix = (int*)malloc((c6count + 1) * sizeof(int));
	prefix(mapping, c6count, map_prefix);

	cn_elim<<<num_blocks_per_grid, num_threads_per_block>>>(d_c6elim, d_c_n, d_n_o, c6count, c_n_size, n_o_size);
	cudaMemcpy(c_n,  d_c_n, 2 * n_o_size * sizeof(int), cudaMemcpyDeviceToHost);

	// N extraction and O permutation
	std::vector<int> v;
	for (int i = 0; i < c6count; i++) {
		for (int j = 0; j < 3; j++) {
			if (std::find(v.begin(), v.end(), c6elim[i + j * 2 * c6count]) != v.end()) {
				continue;
			}
			else {
				v.push_back(c6elim[i + j * 2 * c6count]);
			}
		}
	}
	std::vector<int> v1;
	n_ext(c_n, v, c_n_size, v.size(), v1);											// Vector of N atom in TNT 
	int* n_arr = (int*)malloc(v1.size() * sizeof(int));
	copy(v1.begin(), v1.end(), n_arr);												// Copy vector to array

	int* d_n_arr;
	cudaMalloc((void **) & d_n_arr, v1.size() * sizeof(int));
	cudaMemcpy(d_n_arr, n_arr, v1.size() * sizeof(int), cudaMemcpyHostToDevice);
	int* n_pref = (int*)malloc(v1.size() * sizeof(int));							// Num of permutation of each N
	int* d_n_pref;
	cudaMalloc((void **) & d_n_pref, v1.size() * sizeof(int));
	n_prefix<<<num_blocks_per_grid, num_threads_per_block>>>(d_n_arr, d_n_o, d_n_pref, v1.size(), n_o_size);
	cudaMemcpy(n_pref, d_n_pref, v1.size()* sizeof(int), cudaMemcpyDeviceToHost);
	int* n_prefix = (int*)malloc((v1.size() + 1) * sizeof(int));					// Prefix sum of Permutation of N, used to set location
	prefix(n_pref, v1.size(), n_prefix);
	int* d_n_prefix;
	cudaMalloc((void **) & d_n_prefix, (v1.size() + 1) * sizeof(int));
	cudaMemcpy(d_n_prefix, n_prefix, (v1.size() + 1) * sizeof(int), cudaMemcpyHostToDevice);
	int* no2_perm = (int*)malloc(3 * n_prefix[v1.size()] * sizeof(int));
	int*d_no2_perm;
	cudaMalloc((void **) & d_no2_perm, 3 * n_prefix[v1.size()] * sizeof(int));
	no2_permutation<<<num_blocks_per_grid, num_threads_per_block>>>(d_n_arr, d_n_o, d_n_prefix, v1.size(), n_o_size, d_no2_perm, n_prefix[v1.size()]);
	cudaMemcpy(no2_perm, d_no2_perm, 3 * n_prefix[v1.size()] * sizeof(int), cudaMemcpyDeviceToHost);


	// Start Mapping
	int total_map_count = total_mapping(mapping, c6count);
	int* result = (int*)malloc(15 * total_map_count * sizeof(int));
	int* d_map_prefix, *d_result;
	cudaMalloc((void **) &d_map_prefix, (c6count + 1) * sizeof(int));
	cudaMemcpy(d_map_prefix, map_prefix, (c6count + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMalloc((void **) &d_result, 15 * total_map_count * sizeof(int));

	tnt_mapping<<<num_blocks_per_grid, num_threads_per_block>>>(d_no2_perm, d_c6elim, d_c_n, d_map_prefix, n_prefix[v1.size()], c6count, c_n_size, d_result, total_map_count);
	cudaMemcpy(result, d_result, 15 * total_map_count * sizeof(int), cudaMemcpyDeviceToHost);

	// Result
	final_results = result;
	final_result_size = total_map_count;

	// Release memory used
	free(c6ring_1);
	free(mapping);
	free(map_prefix);
	free(no2_perm);
	free(n_prefix);
	free(n_arr);
	free(n_pref);
	free(c6elim);

	cudaFree(d_c6ring_1);
	cudaFree(d_c6elim);
	cudaFree(d_mapping);
	cudaFree(d_map_prefix);
	cudaFree(d_result);
	cudaFree(d_n_arr);
	cudaFree(d_n_pref);
	cudaFree(d_no2_perm);
	cudaFree(d_result);

	cudaFree(d_c_h);
	cudaFree(d_c_c);
	cudaFree(d_c_n);	
	cudaFree(d_n_o);
	



}
