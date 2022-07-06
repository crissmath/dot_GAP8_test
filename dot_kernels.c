/***************************************
*          SIMD DOT 
****************************************/
void dotproduct_simd(signed int acc, signed char* vA, signed char* vB, signed int N)
{
    uint i, base = 0;
    float t1, t2;
    unsigned int instr, cycles, ldstall, jrstall, imstall;
    unsigned char acc_1,acc_2,acc_3,acc_4;
    v4u A;
    v4u B;

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
    for(i = 0; i<N; i=i+4){ 
        A = *((v4u *) (&vA[i])); // load A
        B = *((v4u *) (&vB[i])); // load B
        printf(" i = %d\n", i);
        printf("A[%d] = %u, A[%d] = %u, A[%d] = %u, A[%d] = %u\n", base, A[base], base+1, A[base+1], base+2, A[base+2], base+3, A[base+3]);
        printf("B[%d] = %u, B[%d] = %u, B[%d] = %u, B[%d] = %u\n", base, B[base], base+1, B[base+1], base+2, B[base+2], base+3, B[base+3]);
        acc += gap8_sumdotpu4(A, B, acc);
        //printf("acc = %d \n", acc);
        //acc += vA[i]*vB[i];
    }
    t2   = rt_time_get_us();

  printf(" acc = %d\n", acc);
  /*  
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
*/
    pi_cl_team_barrier();
  
}
/***************************************
*          DOT Product 
****************************************/
void dotproduct(signed int acc, signed char* vA, signed char* vB, signed int N)
{
    int i;
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

    printf(" acc = %d\n", acc);

    for(i = 0; i<N; i++){
        acc += vA[i]*vB[i];
        printf(" acc = %d\n", acc);
    }

    printf(" acc = %d\n", acc);

    //stop the HW counter used for monitoring
    rt_perf_stop(perf);

    //get the total measurement
    rt_perf_save(perf);

    instr   = rt_perf_get(perf, RT_PERF_INSTR);
    cycles  = rt_perf_get(perf, RT_PERF_ACTIVE_CYCLES);
    ldstall = rt_perf_get(perf, RT_PERF_LD_STALL);
    jrstall = rt_perf_get(perf, RT_PERF_JR_STALL);
    imstall = rt_perf_get(perf, RT_PERF_IMISS);

    printf("Perf of dot product: \n \t cycles : %d "
                                 "\n \t instructions : %d"
                                 "\n \t load stalls : %d"
                                 "\n \t jump stalls : %d"
                                 "\n \t insrtructions stalls: %d"
                                 "\n\n", cycles, instr, ldstall, jrstall, imstall);

    pi_cl_team_barrier();
}