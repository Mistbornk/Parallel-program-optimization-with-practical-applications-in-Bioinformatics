// optimized_cuda_sw.cu

#include "sw.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

struct MaxScore {
    int score;
    int i;
    int j;
};

__device__ int dev_max4(int a, int b, int c, int d) {
    return max(max(a, b), max(c, d));
}

__global__ void cuda_sw_kernel_expand(
    int* H, int* E, int* F, const char* ref, const char* query, int k,
    int M, int N, int match, int mismatch, int gap_open, int gap_extend, MaxScore* d_max
) {
    // i: ref, j:query, k:thread or 斜對角 idx (start on 1)
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int i = tid + 1;
    int j = k - i + 1;

    if (i > 0 && j > 0 && i <= M && j <= N) {
        int idx = j * (M+1) + i;
        int idx_up = (j-1) * (M+1) + i;
        int idx_left = j * (M+1) + (i-1);
        int idx_diag = (j-1) * (M+1) + (i-1);

        int sub_score = (ref[i-1] == query[j-1]) ? match : -mismatch;
        int e = max(H[idx_left] - gap_open, E[idx_left] - gap_extend);
        int f = max(H[idx_up] - gap_open, F[idx_up] - gap_extend);
        int h = dev_max4(0, H[idx_diag] + sub_score, e, f);

        H[idx] = h;
        E[idx] = e;
        F[idx] = f;

        // 儲存最大值（只在值變大時記錄座標）
        int old = atomicMax(&(d_max->score), h);
        if (h > old) {
            d_max->i = i;
            d_max->j = j;
        }
    }
}

__global__ void cuda_sw_kernel_shrink(
    int* H, int* E, int* F, const char* ref, const char* query, int k,
    int M, int N, int match, int mismatch, int gap_open, int gap_extend, MaxScore* d_max
) {
    int i_start = max(1, k - N + 1);
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int i = tid + i_start;
    int j = k + 1 - i;
    //printf("[k=%d] thread %d → i=%d, j=%d\n", k, threadIdx.x, i, j);
    
    if (i > 0 && j > 0 && i <= M && j <= N) {
        int idx = j * (M+1) + i;
        int idx_up = (j-1) * (M+1) + i;
        int idx_left = j * (M+1) + (i-1);
        int idx_diag = (j-1) * (M+1) + (i-1);

        int sub_score = (ref[i-1] == query[j-1]) ? match : -mismatch;
        int e = max(H[idx_left] - gap_open, E[idx_left] - gap_extend);
        int f = max(H[idx_up] - gap_open, F[idx_up] - gap_extend);
        int h = dev_max4(0, H[idx_diag] + sub_score, e, f);

        H[idx] = h;
        E[idx] = e;
        F[idx] = f;

        // 儲存最大值（只在值變大時記錄座標）
        int old = atomicMax(&(d_max->score), h);
        if (h > old) {
            d_max->i = i;
            d_max->j = j;
        }
    }
}


SmithWaterman cuda_smith_waterman(std::string_view ref, std::string_view query,
    int match, int mismatch, int gap_open, int gap_extend
) {
    int M = ref.size();
    int N = query.size();

    std::vector<int> H((M+1)*(N+1), 0);
    std::vector<int> E((M+1)*(N+1), 0);
    std::vector<int> F((M+1)*(N+1), 0);

    int max_score = 0;
    int max_i = 0, max_j = 0;

    int *dpH_dev, *dpE_dev, *dpF_dev;
    char *dev_ref, *dev_query;
    MaxScore* d_max;
    MaxScore h_max = {0, 0, 0};


    cudaMalloc(&dpH_dev, (M+1)*(N+1)*sizeof(int));
    cudaMalloc(&dpE_dev, (M+1)*(N+1)*sizeof(int));
    cudaMalloc(&dpF_dev, (M+1)*(N+1)*sizeof(int));
    cudaMalloc(&dev_ref, M * sizeof(char));
    cudaMalloc(&dev_query, N * sizeof(char));
    cudaMalloc(&d_max, sizeof(MaxScore));

    cudaMemcpy(dpH_dev, H.data(), (M+1)*(N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dpE_dev, E.data(), (M+1)*(N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dpF_dev, F.data(), (M+1)*(N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ref, ref.data(), M, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_query, query.data(), N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_max, &h_max, sizeof(MaxScore), cudaMemcpyHostToDevice);

	// kernel implement
    int threadsPerBlock = 1024;
    for (int k = 1; k <= min(M, N); ++k) {
        int totalThreads = k;
        int numBlocks = (totalThreads + threadsPerBlock - 1) / threadsPerBlock;
    
        cuda_sw_kernel_expand<<<numBlocks, threadsPerBlock>>>(
            dpH_dev, dpE_dev, dpF_dev, dev_ref, dev_query, 
            k, M, N, match, mismatch, gap_open, gap_extend, d_max
        );
        cudaDeviceSynchronize();
    }
    
    for (int k = min(M, N) + 1; k <= M + N - 1; ++k) {
        int i_start = max(1, k - N + 1);
        int i_end = min(M, k);
        int thread_count = i_end - i_start + 1;
        int numBlocks = (thread_count + threadsPerBlock - 1) / threadsPerBlock;
    
        cuda_sw_kernel_shrink<<<numBlocks, threadsPerBlock>>>(
            dpH_dev, dpE_dev, dpF_dev, dev_ref, dev_query,
            k, M, N, match, mismatch, gap_open, gap_extend, d_max
        );
        cudaDeviceSynchronize();
    }    

    cudaMemcpy(H.data(), dpH_dev, (M+1)*(N+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(E.data(), dpE_dev, (M+1)*(N+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(F.data(), dpF_dev, (M+1)*(N+1)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_max, d_max, sizeof(MaxScore), cudaMemcpyDeviceToHost);

    cudaFree(dpH_dev);
    cudaFree(dpE_dev);
    cudaFree(dpF_dev);
    cudaFree(dev_ref);
    cudaFree(dev_query);
    cudaFree(d_max);

    // Print H matrix
    //printf("\nScoring matrix H:\n");
    //for (int i = 0; i <= M; ++i) {
    //    for (int j = 0; j <= N; ++j) {
    //        printf("%3d ", H[j*(M+1)+i]);
    //    }
    //    printf("\n");
    //}

    std::string aligned_ref, aligned_query, match_line;
    int i = max_i = h_max.i, j = max_j = h_max.j;
    max_score = h_max.score;
    while (i > 0 && j > 0 && H[j * (M+1) + i] > 0) {
        int idx      = j * (M+1) + i;
        int idx_diag = (j-1) * (M+1) + (i-1);
        int idx_up   = (j-1) * (M+1) + i;
        int idx_left = j * (M+1) + (i-1);
    
        int score_diag = H[idx_diag] + ((ref[i-1] == query[j-1]) ? match : -mismatch);
        int score_E = std::max(H[idx_left] - gap_open, E[idx_left] - gap_extend);
        int score_F = std::max(H[idx_up] - gap_open, F[idx_up] - gap_extend);
    
        if (H[idx] == score_E) {
            aligned_ref += ref[i-1];
            aligned_query += '-';
            match_line += ' ';
            i--;
        } else if (H[idx] == score_F) {
            aligned_ref += '-';
            aligned_query += query[j-1];
            match_line += ' ';
            j--;
        } else {
            aligned_ref += ref[i-1];
            aligned_query += query[j-1];
            match_line += (ref[i-1] == query[j-1] ? '|' : '*');
            i--; j--;
        }
    }    

    std::reverse(aligned_ref.begin(), aligned_ref.end());
    std::reverse(aligned_query.begin(), aligned_query.end());
    std::reverse(match_line.begin(), match_line.end());

    return SmithWaterman{
        .score = static_cast<int16_t>(max_score),
        .aligned_seq1 = std::move(aligned_ref),
        .aligned_seq2 = std::move(aligned_query),
        .match_line = std::move(match_line),
        .start1 = static_cast<size_t>(i),
        .end1 = static_cast<size_t>(max_i),
        .start2 = static_cast<size_t>(j),
        .end2 = static_cast<size_t>(max_j)
    };
}