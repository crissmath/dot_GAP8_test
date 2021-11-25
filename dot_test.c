#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pulp.h"
#include "Gap8.h"
#include "rt/rt_api.h"
  
   #define NELEM 128
   #define INS 50000000

   #define use_simd
   //#define use_normal
  
   unsigned char vecA[NELEM];
   unsigned char vecB[NELEM];

  
  unsigned int GOLDEN_VALUE = 20743432;

unsigned int dotproduct(unsigned int acc, unsigned char* vA, unsigned char* vB, unsigned int N);
unsigned int dotproduct_simd(unsigned int acc, unsigned char* vA, unsigned char* vB, unsigned int N);
  
   void init_vec(unsigned char* V, int n, unsigned char off)
   {
       int i;
       for(i=0;i<n;i++)
           V[i] = (unsigned char)i + off;
   }
  
   int main()
   {
  
     unsigned int acc = 0;
  
     if(get_core_id() == 0) {
       init_vec(vecA,NELEM,0);
       init_vec(vecB,NELEM,1);
       #ifdef use_simd
       printf("SIMD\n");
       acc += dotproduct_simd(acc, vecA, vecB, NELEM);
       #else
       #ifdef use_normal
       printf("Normal\n");
       acc += dotproduct(acc, vecA, vecB, NELEM);
       #endif
       #endif
       if(acc != GOLDEN_VALUE)
         printf("dot product is is %d instead of %d\n",acc, GOLDEN_VALUE);
       else
         printf("Nice! Well done! 0 errors\n");
     }
  
   }



unsigned int dotproduct_simd(unsigned int acc, unsigned char* vA, unsigned char* vB, unsigned int N)
{
    uint i, base = 0;
    float t1, t2;
    unsigned char elemA, elemB;
    unsigned int instr, cycles, ldstall, jrstall, imstall;



    unsigned char acc_1,acc_2,acc_3,acc_4;
     v4u A = {1,2,3,4};
     v4u B = {2,3,4,5};

    rt_perf_t *perf;
    perf = rt_alloc(RT_ALLOC_L2_CL_DATA, sizeof(rt_perf_t));

    rt_perf_init(perf);
    rt_perf_conf(perf, (1<<RT_PERF_ACTIVE_CYCLES) | (1<<RT_PERF_INSTR) |
                       (1<< RT_PERF_LD_STALL )    | (1<< RT_PERF_JR_STALL ) |
                       (1 << RT_PERF_IMISS ) ) ;

    rt_perf_reset(perf);

    //start the monitoring
    rt_perf_start(perf);
    t1   = rt_time_get_us();
    for(i = 0; i<INS; i++){ 
        //A = *((v4u *) (&vA[base])); // load A
        //B = *((v4u *) (&vB[base])); // load B
        //printf(" i = %d\n", i);
        //printf("A[%d] = %u, A[%d] = %u, A[%d] = %u, A[%d] = %u\n", base, A[base], base+1, A[base+1], base+2, A[base+2], base+3, A[base+3]);
        //printf("B[%d] = %u, B[%d] = %u, B[%d] = %u, B[%d] = %u\n", base, B[base], base+1, B[base+1], base+2, B[base+2], base+3, B[base+3]);
        acc_1 = gap8_sumdotpu4(A, B, acc_1);
        acc_2 = gap8_sumdotpu4(A, B, acc_2);
        acc_3 = gap8_sumdotpu4(A, B, acc_3);
        acc_4 = gap8_sumdotpu4(A, B, acc_4);
        //printf("acc = %d \n", acc);
        //acc += vA[i]*vB[i];
    }
    t2   = rt_time_get_us();

    int instruciones = 4*(8*INS); 
    float time = (t2 - t1)/1e6;
    float Gflops = instruciones/(time*1e9);
    printf("intruc = %d, time = %f, Gflops = %f\n", instruciones, time, Gflops); 
    //stop the HW counter used for monitoring
    rt_perf_stop(perf);

    //get the total measurement
    rt_perf_save(perf);

    instr   = rt_perf_get(perf, RT_PERF_INSTR);
    cycles  = rt_perf_get(perf, RT_PERF_ACTIVE_CYCLES);
    ldstall = rt_perf_get(perf, RT_PERF_LD_STALL);
    jrstall = rt_perf_get(perf, RT_PERF_JR_STALL);
    imstall = rt_perf_get(perf, RT_PERF_IMISS);

    printf("Perf of dot product: \n \t cycles : %d \n \t instructions : %d \n \t load stalls : %d \n \t jump stalls : %d \n \t insrtructions stalls: %d \n\n", cycles, instr, ldstall, jrstall, imstall);

    return acc;
}



unsigned int dotproduct(unsigned int acc, unsigned char* vA, unsigned char* vB, unsigned int N)
{
    int i;
    unsigned char elemA, elemB;
    unsigned int instr, cycles, ldstall, jrstall, imstall;

    rt_perf_t *perf;
    perf = rt_alloc(RT_ALLOC_L2_CL_DATA, sizeof(rt_perf_t));

    rt_perf_init(perf);
    rt_perf_conf(perf, (1<<RT_PERF_ACTIVE_CYCLES) | (1<<RT_PERF_INSTR) |
                       (1<< RT_PERF_LD_STALL ) | (1<< RT_PERF_JR_STALL ) |
                       (1 << RT_PERF_IMISS ) ) ;

    rt_perf_reset(perf);

    //start the monitoring
    rt_perf_start(perf);

    for(i = 0; i<N; i++){
        acc += vA[i]*vB[i];
    }

    //stop the HW counter used for monitoring
    rt_perf_stop(perf);

    //get the total measurement
    rt_perf_save(perf);

    instr   = rt_perf_get(perf, RT_PERF_INSTR);
    cycles  = rt_perf_get(perf, RT_PERF_ACTIVE_CYCLES);
    ldstall = rt_perf_get(perf, RT_PERF_LD_STALL);
    jrstall = rt_perf_get(perf, RT_PERF_JR_STALL);
    imstall = rt_perf_get(perf, RT_PERF_IMISS);

    printf("Perf of dot product: \n \t cycles : %d \n \t instructions : %d \n \t load stalls : %d \n \t jump stalls : %d \n \t insrtructions stalls: %d \n\n", cycles, instr, ldstall, jrstall, imstall);

    return acc;
}