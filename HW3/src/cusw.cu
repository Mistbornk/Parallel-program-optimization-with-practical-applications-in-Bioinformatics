#include "sw.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#define DIR_NONE  0  // 來自 0，起點
#define DIR_DIAG  1  // 來自斜對角（配對/不配）
#define DIR_LEFT  2  // 來自左（gap in query）
#define DIR_UP    3  // 來自上（gap in ref）

extern "C" void cuda_warmup() {
    cudaFree(0); // 觸發初始化
}

struct MaxScore {
    int score;
    int i;
    int j;
};

__device__ void update_max_score(MaxScore* d_max, int h, int i, int j) {
    int old = atomicMax(&(d_max->score), h);
    if (h > old) {
        d_max->i = i;
        d_max->j = j;
    }
}

__device__ int dev_max4(int a, int b, int c, int d) {
    return max(max(a, b), max(c, d));
}

// kernel to compute diagonal wavefront, only using 3 row buffers
__global__ void cuda_sw_kernel(
    const char* ref, const char* query,
    int* H_prev2, int* H_prev1, int* H_curr,
    int* E_prev, int* E_curr,
    int* F_prev, int* F_curr,
    uint8_t* d_dir,
    int k, int M, int N,
    int match, int mismatch,
    int gap_open, int gap_extend,
    MaxScore* d_max
) {
    // i + j = k + 1
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int i_start = max(1, k - N + 1);  // i 起點
    int i_end   = min(M, k);          // i 終點
    int i = i_start + tid;
    int j = k - i + 1;

    if (i <=0 || i > i_end || j <= 0 || j > N) return; 

    int idx_up = i;
    int idx_left = i - 1;
    int idx_diag = i - 1;

    int sub_score = (ref[i-1] == query[j-1]) ? match : -mismatch;
    int e = max(H_prev1[idx_left] - gap_open, E_prev[idx_left] - gap_extend);
    int f = max(H_prev1[idx_up] - gap_open, F_prev[idx_up] - gap_extend);
    int h = dev_max4(0, H_prev2[idx_diag] + sub_score, e, f);

    H_curr[i] = h;
    E_curr[i] = e;
    F_curr[i] = f;

    // 決定方向
    uint8_t dir = DIR_NONE;
    int h_diag = H_prev2[idx_diag] + sub_score;
    if (h == e) dir = DIR_LEFT;
    else if (h == f) dir = DIR_UP;
    else if (h == h_diag) dir = DIR_DIAG;
    d_dir[j * (M + 1) + i] = dir;  // 對應於 direction[j][i]

    update_max_score(d_max, h, i, j);
}

SmithWaterman cuda_smith_waterman(std::string_view ref, std::string_view query,
    int match, int mismatch, int gap_open, int gap_extend
) {
    int M = ref.size();
    int N = query.size();
    int size = max(M, N);

    std::vector<int> H_prev2(size  + 1, 0); std::vector<int> H_prev1(size  + 1, 0); std::vector<int> H_curr(size  + 1, 0);
    std::vector<int> E_prev(size  + 1, 0); std::vector<int> E_curr(size  + 1, 0);
    std::vector<int> F_prev(size  + 1, 0); std::vector<int> F_curr(size  + 1, 0);
    std::vector<uint8_t> flat_dir((M + 1) * (N + 1));
    
    int *d_H_prev2; int *d_H_prev1; int *d_H_curr;
    int *d_E_prev; int *d_E_curr;
    int *d_F_prev; int *d_F_curr;    
    char *d_ref, *d_query;
    MaxScore* d_max;
    MaxScore h_max = {0, 0, 0};
    uint8_t* d_dir;
    
    // malloc
    cudaMalloc(&d_H_prev2, sizeof(int) * (size  + 1)); cudaMalloc(&d_H_prev1, sizeof(int) * (size  + 1)); cudaMalloc(&d_H_curr,  sizeof(int) * (size  + 1));
    cudaMalloc(&d_E_prev,  sizeof(int) * (size  + 1)); cudaMalloc(&d_E_curr,  sizeof(int) * (size  + 1));
    cudaMalloc(&d_F_prev,  sizeof(int) * (size  + 1)); cudaMalloc(&d_F_curr,  sizeof(int) * (size  + 1));
    cudaMalloc(&d_ref, sizeof(char) * M); cudaMalloc(&d_query, sizeof(char) * N); cudaMalloc(&d_max, sizeof(MaxScore));
    cudaMalloc(&d_dir, sizeof(uint8_t) * (M + 1) * (N + 1));

    // memcpy
    cudaMemcpy(d_H_prev2, H_prev2.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice); cudaMemcpy(d_H_prev1, H_prev1.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice); cudaMemcpy(d_H_curr, H_curr.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(d_E_prev, E_prev.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice); cudaMemcpy(d_E_curr, E_curr.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F_prev, F_prev.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice); cudaMemcpy(d_F_curr, F_curr.data(), sizeof(int) * (size  + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ref, ref.data(), M, cudaMemcpyHostToDevice); cudaMemcpy(d_query, query.data(), N, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_max, &h_max, sizeof(MaxScore), cudaMemcpyHostToDevice);

    // kernel implement
    for (int k = 1; k< M + N - 1; ++k) {
        int i_start = max(1, k - N + 1);
        int i_end = min(M, k);
        int thread_count = i_end - i_start + 1;
    
        int threadsPerBlock = 256;
        int numBlocks = (thread_count + threadsPerBlock - 1) / threadsPerBlock;
    
        cuda_sw_kernel<<<numBlocks, threadsPerBlock>>>(
            d_ref, d_query,
            d_H_prev2, d_H_prev1, d_H_curr,
            d_E_prev, d_E_curr,
            d_F_prev, d_F_curr,
            d_dir,
            k, M, N, match, mismatch, gap_open, gap_extend,
            d_max
        );

        std::swap(d_H_prev2, d_H_prev1);
        std::swap(d_H_prev1, d_H_curr);
        std::swap(d_E_prev, d_E_curr);
        std::swap(d_F_prev, d_F_curr);
    }
    cudaDeviceSynchronize();

    cudaMemcpy(&h_max, d_max, sizeof(MaxScore), cudaMemcpyDeviceToHost);
    cudaMemcpy(flat_dir.data(), d_dir, sizeof(uint8_t) * (M + 1) * (N + 1), cudaMemcpyDeviceToHost);
    
    // clean up
    cudaFree(d_H_prev2); cudaFree(d_H_prev1); cudaFree(d_H_curr);
    cudaFree(d_E_prev); cudaFree(d_E_curr);
    cudaFree(d_F_prev); cudaFree(d_F_curr);
    cudaFree(d_ref); cudaFree(d_query);
    cudaFree(d_max);
    cudaFree(d_dir);


    //direction[j][i] = flat_dir[j * (M + 1) + i];
    std::string aligned_ref, aligned_query, match_line;
    int i = h_max.i, j = h_max.j;
    while (i > 0 && j > 0) {
        uint8_t dir = flat_dir[j * (M + 1) + i];
        if (dir == DIR_LEFT) {
            aligned_ref += ref[i - 1];
            aligned_query += '-';
            match_line += ' ';
            --i;
        } else if (dir == DIR_UP) {
            aligned_ref += '-';
            aligned_query += query[j - 1];
            match_line += ' ';
            --j;
        } else if (dir == DIR_DIAG) {
            aligned_ref += ref[i - 1];
            aligned_query += query[j - 1];
            match_line += (ref[i - 1] == query[j - 1] ? '|' : '*');
            --i; --j;
        } else break;  // DIR_NONE or invalid
    }
    std::reverse(aligned_ref.begin(), aligned_ref.end());
    std::reverse(aligned_query.begin(), aligned_query.end());
    std::reverse(match_line.begin(), match_line.end());
    
    return SmithWaterman{
        .score = static_cast<int16_t>(h_max.score),
        .aligned_seq1 = std::move(aligned_ref),
        .aligned_seq2 = std::move(aligned_query),
        .match_line = std::move(match_line),
        .start1 = static_cast<size_t>(i),
        .end1 = static_cast<size_t>(h_max.i),
        .start2 = static_cast<size_t>(j),
        .end2 = static_cast<size_t>(h_max.j)
    };

    //return SmithWaterman{
    //    .score = h_max.score,
    //    .aligned_seq1 = "",
    //    .aligned_seq2 = "",
    //    .match_line = "",
    //    .start1 = 0,
    //    .end1 = static_cast<size_t> (h_max.i),
    //    .start2 = 0,
    //    .end2 = static_cast<size_t> (h_max.j)
    //};
}