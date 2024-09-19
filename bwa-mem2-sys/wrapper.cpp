#include "ext/bwa-mem2/src/fastmap.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

uint64_t tprof[LIM_R][LIM_C];

worker_t* new_worker_t(int32_t nreads, int32_t nthreads)
{
    worker_t *w = (worker_t *) malloc(sizeof(worker_t));
    int32_t memSize = nreads;
    int32_t readLen = READ_LEN;

    /* Mem allocation section for core kernels */
    w->regs = NULL; w->chain_ar = NULL; w->seedBuf = NULL;

    w->regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
    w->chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w->seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

    assert(w->seedBuf  != NULL);
    assert(w->regs     != NULL);
    assert(w->chain_ar != NULL);

    w->seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

    /* SWA mem allocation */
    int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
    for(int l=0; l<nthreads; l++)
    {
        w->mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w->mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w->mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w->mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);

        w->mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
        w->mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;

        assert(w->mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
        assert(w->mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
        assert(w->mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
        assert(w->mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
    }

    for(int l=0; l<nthreads; l++) {
        w->mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w->mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w->mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w->mmc.wsize[l] = wsize;

        assert(w->mmc.seqPairArrayAux[l] != NULL);
        assert(w->mmc.seqPairArrayLeft128[l] != NULL);
        assert(w->mmc.seqPairArrayRight128[l] != NULL);
    }

    for (int l=0; l<nthreads; l++)
    {
        // BATCH_MUL accounts for smems per reads in smem kernel
        w->mmc.wsize_mem[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
        w->mmc.wsize_mem_s[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
        w->mmc.wsize_mem_r[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
        w->mmc.matchArray[l]    = (SMEM *) _mm_malloc(w->mmc.wsize_mem[l] * sizeof(SMEM), 64);
        w->mmc.min_intv_ar[l]   = (int32_t *) malloc(w->mmc.wsize_mem[l] * sizeof(int32_t));
        w->mmc.query_pos_ar[l]  = (int16_t *) malloc(w->mmc.wsize_mem[l] * sizeof(int16_t));
        w->mmc.enc_qdb[l]       = (uint8_t *) malloc(w->mmc.wsize_mem[l] * sizeof(uint8_t));
        w->mmc.rid[l]           = (int32_t *) malloc(w->mmc.wsize_mem[l] * sizeof(int32_t));
        w->mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64); // candidate not for reallocation, deferred for next round of changes.
    }

    return w;
}

void destroy_worker_t(worker_t *w, int32_t nthreads)
{
    free(w->chain_ar);
    free(w->regs);
    free(w->seedBuf);

    for(int l=0; l<nthreads; l++) {
        _mm_free(w->mmc.seqBufLeftRef[l*CACHE_LINE]);
        _mm_free(w->mmc.seqBufRightRef[l*CACHE_LINE]);
        _mm_free(w->mmc.seqBufLeftQer[l*CACHE_LINE]);
        _mm_free(w->mmc.seqBufRightQer[l*CACHE_LINE]);
    }

    for(int l=0; l<nthreads; l++) {
        free(w->mmc.seqPairArrayAux[l]);
        free(w->mmc.seqPairArrayLeft128[l]);
        free(w->mmc.seqPairArrayRight128[l]);
    }

    for(int l=0; l<nthreads; l++) {
        _mm_free(w->mmc.matchArray[l]);
        free(w->mmc.min_intv_ar[l]);
        free(w->mmc.query_pos_ar[l]);
        free(w->mmc.enc_qdb[l]);
        free(w->mmc.rid[l]);
        _mm_free(w->mmc.lim[l]);
    }

    free(w);
}

alloc_ref_string(char* binary_seq_file) {
    FILE *fr = fopen(binary_seq_file, "r");

    if (fr == NULL) {
        fprintf(stderr, "Error: can't open %s input file\n", binary_seq_file);
        exit(EXIT_FAILURE);
    }

    int64_t rlen = 0;
    fseek(fr, 0, SEEK_END);
    rlen = ftell(fr);
    ref_string = (uint8_t*) _mm_malloc(rlen, 64);

    rewind(fr);

    /* Reading ref. sequence */
    err_fread_noeof(ref_string, 1, rlen, fr);

    uint64_t timer  = __rdtsc();
    tprof[REF_IO][0] += timer - tim;

    fclose(fr);

    fclose(fr);
}