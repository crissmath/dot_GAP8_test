#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Gap8.h"
#include "pulp.h"
#include "pmsis.h"
#include "rt/rt_api.h"
#include "bsp/ram.h"
#include "bsp/ram/hyperram.h"

#define VERBOSE
#define V_size 24             // vector size
#define BUFFER_SIZE_L2 (8)    // buffer size

// Com Buffers 
signed char *buff_comm;
signed char *buff_com_tmp;

/* init ram */
struct pi_device ram;           /* ram object */
struct pi_hyperram_conf conf;   /* config object */
//
signed char *L3_va;
signed char *L3_vb;
uint32_t    *L3_out;

// set voltage
#define pi_pmu_set_voltage(x, y)      ( rt_voltage_force(RT_VOLTAGE_DOMAIN_MAIN, x, NULL) )

// struct for dot product task
typedef struct   dot_arg
{
    signed char  *v_a;
    signed char  *v_b;
    uint32_t     *acc;
    uint32_t     dim;
}dot_arg_t;

dot_arg_t Arg;

// struct for dot product Cluster task
typedef struct dot_cl_arg{
    signed char * L2_va;
    signed char * L2_vb;
    signed char * L1_va;
    signed char * L1_vb;
    uint32_t    * L1_acc;
    uint32_t    * L2_acc;
    uint32_t    dim;
}dot_cl_arg_t;

dot_cl_arg_t cl_Arg;

/***************************************
*          SIMD DOT 
****************************************/
void dotproduct_simd(dot_arg_t *Arg)
{    
    signed char *vA    = Arg->v_a;
    signed char *vB    = Arg->v_b;
    uint32_t    *c     = Arg->acc;
    uint32_t     N     = Arg->dim;

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
    for(i = 0; i < N; i=i+4){ 
        A = *((v4u *) (&vA[i])); // load A
        B = *((v4u *) (&vB[i])); // load B
        /*printf(" i = %d\n", i);
        printf("A[%d] = %u, A[%d] = %u, A[%d] = %u, A[%d] = %u\n",
                 base, A[base], base+1, A[base+1], base+2, A[base+2], base+3, A[base+3]);
        printf("B[%d] = %u, B[%d] = %u, B[%d] = %u, B[%d] = %u\n",
                 base, B[base], base+1, B[base+1], base+2, B[base+2], base+3, B[base+3]);*/
        *c = gap8_sumdotpu4(A, B, *c);
        //printf("acc = %d \n", acc);
        //acc += vA[i]*vB[i];
    }
    t2   = rt_time_get_us();
  
    //int instruciones = 4*(8*INS); 
    //float time = (t2 - t1)/1e6;
    //float Gflops = instruciones/(time*1e9);
    //printf("intruc = %d, time = %f, Gflops = %f\n", instruciones, time, Gflops); 
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
/***************************************
*          DOT Product 
****************************************/
void dotproduct( dot_arg_t *Arg)
{
    signed char *vA    = Arg->v_a;
    signed char *vB    = Arg->v_b;
    uint32_t    *c     = Arg->acc;
    uint32_t     N     = Arg->dim;

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

    for(i = 0; i<N; i++){
        *c += vA[i]*vB[i];
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
    printf("Perf of dot product: \n \t cycles : %d "
                                 "\n \t instructions : %d"
                                 "\n \t load stalls : %d"
                                 "\n \t jump stalls : %d"
                                 "\n \t insrtructions stalls: %d"
                                 "\n\n", cycles, instr, ldstall, jrstall, imstall);
    pi_cl_team_barrier();
}

/***************************************
*          Zero init vectors
*****************************************/
   void zero_init(signed char* V, int n, signed char off)
   { 
       int i;
       for(i=0;i<n;i++)
           V[i] = 0;
   }
/***************************************
*          Init vectors
*****************************************/
   void init_vec(signed char* V, int n, signed char off)
   { 
       int i;
       for(i=0;i<n;i++)
           V[i] = (signed char)1 + off;
   }
/***************************************
*          Print vectors
*****************************************/
   void print_vec(signed char* V, int n)
   { 
       int i;
       for(i=0;i<n;i++){
           printf(" v[%d] = %u \n", i, V[i]);
       }
       printf("\n\n");
   }
/**************************************************
*               Dummy test
***************************************************/
void dummy( dot_arg_t *Arg){

    signed char *a    = Arg->v_a;
    signed char *b    = Arg->v_b;
    uint32_t    *c    = Arg->acc;
    uint32_t     N    = Arg->dim;

    print_vec(a, N);
    print_vec(b, N);
    printf(" c = %d\n", *c);
    printf(" size = %d\n", N);

    int i;
    for(i = 0; i<N; i++){
        *c += a[i] * b[i];
    }
    pi_cl_team_barrier();
}
/***************************************
*    Generate vector L3
*****************************************/
void generate_vector_L3( uint32_t n, uint8_t off, signed char *V,
                         struct pi_device *ram, signed char *buff_com) {
#ifdef VERBOSE
  printf("start generate...\n");
#endif
  static int factor = 1;
  uint32_t b, i, ncp, resto = 0, base = 0;
  // number of blocks
  ncp = (uint32_t)(n * sizeof(unsigned char)) / BUFFER_SIZE_L2;
  // last part of vector 
  resto = (uint32_t)(n * sizeof(unsigned char)) % BUFFER_SIZE_L2;
#ifdef VERBOSE
  printf("number of cp resto : %d %d\n", ncp, resto);
#endif
  for (b = 0; b < ncp; b++) {
    base = b * BUFFER_SIZE_L2;
    for (i = 0; i < BUFFER_SIZE_L2; i++) {
      //buff_com[i] = (signed char)((b * 2 + i % 4) - 1);
      buff_com[i] = (signed char)(1)+off;
#ifdef VERBOSE
  printf("buff_comm[%i]: %d \n", i, buff_com[i]);
#endif
    }
    // copy L2->L3
    pi_ram_write(ram, (uint32_t)V + base, buff_com, (uint32_t)BUFFER_SIZE_L2);
  }
  if (resto > 0) {
    for (i = 0; i < resto; i++) {
      buff_com[i] = (signed char)((i + 1) * factor);
    }
    pi_ram_copy(ram, (uint32_t)V + base, buff_com, (uint32_t)resto, 0);
  }
  factor++;
#ifdef VERBOSE
  printf("end generate ..\n");
#endif
}
/**********************************************************************
 *
 *                        Print Vector from L3
 *
 **********************************************************************/
/*! \brief This function printed the vector of size m*n allocated in the ram
 *
 * This function requires the size (n) and the buffer(buffer_com) where the
 * vector is going to be read.
 **********************************************************************/
void print_vector_L3(char *name, uint32_t n, unsigned char *V,
                     struct pi_device *ram, signed char *buff_com) {

  printf("init printmatrix !...\n");
  buff_com = (signed char *)pmsis_l2_malloc((uint32_t)BUFFER_SIZE_L2);
  if (buff_com == NULL) {
    printf("Error allocation buff_comm!\n");
    pmsis_exit(-1);
  }

  uint32_t b, i, ncp, resto = 0, base = 0;
  // nb blks
  ncp   = (uint32_t)( n * sizeof(signed char)) / BUFFER_SIZE_L2;
  // last part of blk
  resto = (uint32_t)( n * sizeof(signed char)) % BUFFER_SIZE_L2;
#ifdef VERBOSE
  printf("number of cp resto : %d %d\n", ncp, resto);
#endif
  for (b = 0; b < ncp; b++) {
    base = b * BUFFER_SIZE_L2;
    pi_ram_read(ram, (uint32_t)V + base, buff_com, (uint32_t)BUFFER_SIZE_L2);
    for (i = 0; i < (uint32_t)BUFFER_SIZE_L2; i++) {
      printf("%s[%d,%d] = %d;\n", name, b, i, buff_com[i]);
    }
  }
  if (resto > 0) {
    pi_ram_read(ram, (uint32_t)V + base, buff_com, (uint32_t)resto);
    for (i = 0; i < (uint32_t)resto; i++) {
      printf("%s[%d,%d] = %d;\n", name, b, i, buff_com[i]);
    }
  }
  pmsis_l2_malloc_free(buff_com, (uint32_t)BUFFER_SIZE_L2);
  printf("fin printmatrix...\n");
}
/************************************************************
*               Dot Cluster DMA  
*        core0(master) of the cluster
*************************************************************/
void cluster_dma(dot_cl_arg_t *cl_Arg)
{
    signed char *L2_va  = cl_Arg->L2_va;
    signed char *L2_vb  = cl_Arg->L2_vb;
    signed char *L1_va  = cl_Arg->L1_va;
    signed char *L1_vb  = cl_Arg->L1_vb;
    uint32_t    *L2_acc = cl_Arg->L2_acc;
    uint32_t    *L1_acc = cl_Arg->L1_acc;
    uint32_t      N     = cl_Arg->dim;

    uint32_t coreid = pi_core_id(), start = 0, end = 0;
    int acc_d = 0;
    L1_acc = &acc_d;
    //Core0(master) of cluster init DMA transfer L2->L1.
    if (!coreid)
    {
        printf("Core %d requesting va DMA transfer from l2_va to l1_va. size :%d \n", coreid, sizeof(N));
        pi_cl_dma_copy_t copy_va;
        copy_va.dir     = PI_CL_DMA_DIR_EXT2LOC;       // external to cl memory 
        copy_va.merge   = 0;                         
        copy_va.size    = N;   
        copy_va.id      = 0;                             
        copy_va.ext     = L2_va;            
        copy_va.loc     = L1_va;

        pi_cl_dma_memcpy(&copy_va);
        pi_cl_dma_wait(&copy_va);
        printf("Core %d : va Transfer done.\n", coreid);


        printf("Core %d requesting vb DMA transfer from l2_vb to l1_vb.\n", coreid);
        pi_cl_dma_copy_t copy_vb;
        copy_vb.dir   = PI_CL_DMA_DIR_EXT2LOC;       // external to cl memory 
        copy_vb.merge = 0;                         
        copy_vb.size  = N;   
        copy_vb.id    = 0;                             
        copy_vb.ext   = L2_vb;            
        copy_vb.loc   = L1_vb;

        pi_cl_dma_memcpy(&copy_vb);
        pi_cl_dma_wait(&copy_vb);
        printf("Core %d : vb Transfer done.\n", coreid);
    }
    //start = (coreid * ((uint32_t) BUFFER_SIZE_L2  / pi_cl_cluster_nb_cores()));
    //end   = (start - 1 + ((uint32_t) BUFFER_SIZE_L2 / pi_cl_cluster_nb_cores()));
    // sync
    pi_cl_team_barrier(0);

    // Each core computes on specific portion of buffer.
    //if( !coreid ){
    printf("Core %d: L1_acc = %d\n",coreid, *L1_acc);
    for(uint32_t i = 0 ; i < N; i++){
        printf("Core %d : computing...\n", coreid);
        printf("size: %d Core: %d %s[%d] = %d;\t%s[%d] = %d;\n",
                N, coreid, "L1_va", i, L1_va[i], "L1_vb", i, L1_vb[i]);
        
        *L1_acc += (L1_va[i] * L1_vb[i]);
        printf("L1_acc = %d\n", *L1_acc); 
      }
      printf("Core %d : L1_acc = %d, acc_d: %d\n",coreid, *L1_acc, acc_d);
    //}
    //sync
    pi_cl_team_barrier(0);

    // Core 0 of cluster initiates DMA transfer from L1 to L2.
    if (!coreid)
    {
        printf("Core %d requesting out DMA transfer from L1_acc to L2_acc.\n", coreid);
        pi_cl_dma_copy_t copy_acc;
        copy_acc.dir   = PI_CL_DMA_DIR_LOC2EXT;
        copy_acc.merge = 0;
        copy_acc.size  = (uint32_t) N;
        copy_acc.id    = 0;
        copy_acc.ext   = (uint32_t)L2_acc;
        copy_acc.loc   = (uint32_t)L1_acc;

        pi_cl_dma_memcpy(&copy_acc);
        pi_cl_dma_wait(&copy_acc);
        printf("Core %d : Transfer done.\n", coreid);
    }
     pi_cl_team_barrier(0);
}

/************************************************************
*               Delegate function
*  In this function we are in the core0(master) of the cluster
*  The cluster controller is core 0 of the cluster.
*  Core 1 to 7 are the slave cores.
*************************************************************/
void cluster_delegate(void *arg)
{
/*In this function we are in the core 0 of the cluster*/
    unsigned int time;
    printf("Cluster master core entry\n");
    /* vector allocation
    signed char *L1_va  = (signed char *) pmsis_l1_malloc((uint32_t)(V_size * sizeof(signed char)));
    signed char *L1_vb  = (signed char *) pmsis_l1_malloc((uint32_t)(V_size * sizeof(signed char)));
    uint32_t    *L1_out = (uint32_t *)    pmsis_l1_malloc((uint32_t)(sizeof(uint32_t)));
    if (L1_va_t == NULL) {
    printf("buff alloc failed v_A!\n");
    pmsis_exit(-1);
    }
    if (L1_vb_t == NULL) {
    printf("buff alloc failed v_B !\n");
    pmsis_exit(-1);
    }
    if (L1_out == NULL) {
    printf("buff alloc failed v_C !\n");
    pmsis_exit(-1);
    }*/

    // init vector
    //init_vec(L1_va_t, V_size, 0);
    //init_vec(L1_vb_t, V_size, 2);
    //L1_out = &acc;

    //printf("init vector\n");
    //printf("VA\n");
    //print_vec(va, V_size);
    //printf("VB\n");
    //print_vec(vb, V_size);
    //printf("c = %d\n\n\n", *(L1_out));

  /*
    Tasks to execute
  */
    printf("Run Task\n");
    //Arg.v_a = L1_va_t;
    //Arg.v_b = L1_vb_t;
    //Arg.acc = L1_out;
    //Arg.dim = V_size;
    /*pi_cl_team_fork( 1, (void *)dummy, (void *) &Arg);
    printf("acc=%d\n", acc);*/


    //uint nb_cores = pi_cl_cluster_nb_cores();
    uint nb_cores = 1;

    printf("Run Test DMA\n");
    pi_cl_team_fork( nb_cores, (void *)cluster_dma, arg);
    
    //printf("acc=%d\n", *L1_acc);

    //acc = 0;

    //printf("Run dotSIMD\n");
    //pi_cl_team_fork( nb_cores, (void *)dotproduct_simd, (void *) &Arg);
    //printf("acc=%d\n", *(L1_out));
}


/* FC Entry*/
int fc_main()
{
  printf("Entering Fabric Controller\n");
  
  uint32_t errors = 0;    
  uint32_t core_id = pi_core_id();        // in this point we only have core 0
  uint32_t cluster_id = pi_cluster_id();  // the cluster ID is 32 per default 
  printf("[cl_id = %d, core_id = %d] Hello World!\n", (int) cluster_id, (int) core_id);

//Set SoC Voltage to 1.2V
    uint32_t voltage_in_mV = 1200;
    if(pi_pmu_set_voltage(voltage_in_mV, 0)==-1){
        printf("Frequency set failed!\n");
        pmsis_exit(-1);
    }
// Set FC Freq at 250 MHz
    uint32_t fc_freq_in_Hz = 250 * 1000 * 1000;
    pi_freq_set(PI_FREQ_DOMAIN_FC, fc_freq_in_Hz);
    printf("Fabric Controller Frequency %d Hz\n", (int) pi_freq_get(PI_FREQ_DOMAIN_FC));

// ******* Alloc buffera L2 ********* 
// buff com
  buff_comm = (signed char *)pmsis_l2_malloc((uint32_t)BUFFER_SIZE_L2);
  if (buff_comm == NULL) {
  printf("buff alloc failed  buff_comm!\n");
  pmsis_exit(-1);
  }
// L2_va  
  signed char *L2_va = (signed char *)pmsis_l2_malloc((uint32_t)V_size);
  if (L2_va == NULL) {
  printf("buff alloc failed  L2_va!\n");
  pmsis_exit(-1);
  }
// L2_vb 
  signed char *L2_vb = (signed char *)pmsis_l2_malloc((uint32_t)V_size);
  if (L2_vb == NULL) {
  printf("buff alloc failed  L2_vb!\n");
  pmsis_exit(-1);
  }
// L2_acc  
  uint32_t *L2_acc = (uint32_t *)pmsis_l2_malloc(sizeof(uint32_t));
  if (L2_acc == NULL) {
  printf("buff alloc failed  buff_comm!\n");
  pmsis_exit(-1);
  }
#ifdef VERBOSE
    printf("L2:buff_comm allocated : %d %ld.\n", buff_comm,
           (uint32_t) V_size * sizeof(signed char));
#endif

//******  Alloc Buffers L3. *********
// ram config
  pi_hyperram_conf_init(&conf);
  pi_open_from_conf(&ram, &conf);
  if (pi_ram_open(&ram)) {
    printf("Error ram open !\n");
    pmsis_exit(-3);
  }
  
// L3_va
  if (pi_ram_alloc(&ram, &L3_va, (uint32_t) V_size * sizeof(signed char))) {
    printf("L3:Ram malloc failed va!\n");
    pmsis_exit(-4);
  }
  else {
#ifdef VERBOSE
    printf("L3:Ram allocated : %lx %ld.\n", L3_va,
           (uint32_t) V_size * sizeof(signed char));
#endif
  }
// L3_vb
  if (pi_ram_alloc(&ram, &L3_vb, (uint32_t) V_size * sizeof(signed char))) {
    printf("L3:Ram malloc failed v_b!\n");
    pmsis_exit(-4);
  }
  else {
#ifdef VERBOSE
    printf("L3:Ram allocated : %lx %ld.\n", L3_vb,
           (uint32_t) V_size * sizeof(signed char));
#endif
  }
// L3_out
  if (pi_ram_alloc(&ram, &L3_out, (uint32_t)sizeof(uint32_t))) {
    printf("L3:Ram malloc failed out!\n");
    pmsis_exit(-4);
  }  
  else {
#ifdef VERBOSE
    printf("L3:Ram allocated : %lx %ld.\n", L3_out,
           (uint32_t) V_size * sizeof(signed char));
#endif
  }

// generate vector L3
    generate_vector_L3( V_size, 1, L3_va, &ram, buff_comm);
    //print_vector_L3( "va", V_size, L3_va, &ram, buff_comm);

    generate_vector_L3( V_size, 1, L3_vb, &ram, buff_comm);
    //print_vector_L3( "vb", V_size, L3_vb, &ram, buff_comm);

// Copy L3->L2
printf("start: Copy L3->L2\n");
uint32_t k;
for (k = 0; k < V_size; k += BUFFER_SIZE_L2) {
      pi_ram_copy(&ram, (uint32_t)L3_va + k, L2_va + k, (uint32_t)BUFFER_SIZE_L2, 1);
      pi_ram_copy(&ram, (uint32_t)L3_vb + k, L2_vb + k, (uint32_t)BUFFER_SIZE_L2, 1);
}


#ifdef VERBOSE
/*
uint32_t r, i, b;
buff_com_tmp = (signed char *)pmsis_l2_malloc((uint32_t)BUFFER_SIZE_L2);
for (r = 0; r < (V_size/BUFFER_SIZE_L2); r ++) {
      k = r * BUFFER_SIZE_L2;
      pi_ram_read(&ram, (uint32_t)L2_va + k, buff_com_tmp, (uint32_t)BUFFER_SIZE_L2);
      for( i = 0; i < BUFFER_SIZE_L2; i++){
          printf("%s[%d,%d] = %d;\n", "va", k, i, buff_com_tmp[i]);
      }
  }
pmsis_l2_malloc_free(buff_com_tmp, (uint32_t)BUFFER_SIZE_L2);
*/
int i;
uint32_t *out;
uint32_t  dot = 0;
*out = &dot;

for( i = 0; i < V_size; i++){
          printf("%s[%d] = %d;\t%s[%d] = %d;\n",
                 "L2_va", i, L2_va[i], "L2_vb", i, L2_vb[i]);
      *out +=  L2_va[i] * L2_vb[i];
      }
printf("out=%d\n", *out);
printf("done: Copy L3->L2\n");
#endif


/****** Configure & open cluster *******/
struct pi_device cluster_dev;
struct pi_cluster_conf cl_conf;

// Init cluster configuration structure.
pi_cluster_conf_init(&cl_conf);
pi_open_from_conf(&cluster_dev, (void *)&cl_conf);
if (pi_cluster_open(&cluster_dev)){
    printf("Cluster open failed !\n");
    pmsis_exit(-1);
    }

// Set the max freq for the cluster @1.2V
uint32_t cl_freq_in_Hz = 175 * 1000 * 1000;
pi_freq_set(PI_FREQ_DOMAIN_CL, cl_freq_in_Hz);
printf("Cluster Frequency %d Hz\n", (int) pi_freq_get(PI_FREQ_DOMAIN_CL));


// --------------- alloc L1 buffers---------------------------- 
printf("start: L1 allocs\n");
// L1_va
signed char *L1_va = pmsis_l1_malloc((uint32_t)V_size);
if (L1_va == NULL)
{
  printf("L1_va alloc failed !\n");
  pi_cluster_close(&cluster_dev);
  pmsis_exit(-4);
}
// L1_va
signed char *L1_vb = pmsis_l1_malloc((uint32_t)V_size);
if (L1_vb == NULL)
{
  printf("L1_vb alloc failed !\n");
  pi_cluster_close(&cluster_dev);
  pmsis_exit(-4);
}
// L1_acc
uint32_t *L1_acc = pmsis_l1_malloc(sizeof(uint32_t));
if (L1_va == NULL)
{
  printf("L1_acc alloc failed !\n");
  pi_cluster_close(&cluster_dev);
  pmsis_exit(-4);
}
printf("done: L1 allocs\n");
//--------- init cl_Arg ---------------------
uint32_t acc = 0;
L2_acc = &acc;
printf("L2_acc = %d\n", *L2_acc);

cl_Arg.L2_va  = L2_va;
cl_Arg.L2_vb  = L2_vb;
cl_Arg.L1_va  = L1_va;
cl_Arg.L1_va  = L1_vb;
cl_Arg.L2_acc = L2_acc;
cl_Arg.L1_acc = L1_acc;
cl_Arg.dim    = V_size;
//------- Task setup --------------------------- 
struct pi_cluster_task cl_task = {0};
//memset(cl_task, 0, sizeof(struct pi_cluster_task));
cl_task.entry  = cluster_delegate;
cl_task.arg    = (void*)&cl_Arg;
//struct pi_cluster_task cl_task = {0};
//cl_task.entry = cluster_delegate;
//cl_task.arg = NULL;


// send task to cluster

pi_cluster_send_task_to_cl(&cluster_dev, &cl_task);


// Free L1 memory
//pi_l2_free(cl_task, sizeof(struct pi_cluster_task));
pi_cl_l1_free(&cluster_dev, L1_va , V_size);
pi_cl_l1_free(&cluster_dev, L1_vb , V_size);
pi_cl_l1_free(&cluster_dev, L1_acc, sizeof(uint32_t));

//printf("Close cluster after end of computation.\n");
//pi_cluster_close(&cluster_dev);

printf("L2_acc = %d\n", *L2_acc);

// Free L2 memory
//pi_l2_free(cl_task, sizeof(struct pi_cluster_task));
pi_l2_free(L2_va , V_size);
pi_l2_free(L2_vb , V_size);
pi_l2_free(L2_acc, sizeof(uint32_t));

// Terminate and exit the test
printf("Test done !\n");
pmsis_exit(errors);

  /*
  unsigned int acc = 0;
    if(get_core_id() == 0) {
       // choose the function use 
       #ifdef use_simd
       printf("SIMD\n");
       acc += dotproduct_simd(acc, v_A, v_B, V_size);
       #else

       #ifdef use_normal
       printf("Normal\n");
       acc += dotproduct(acc, v_A, v_B, V_size);
       #endif
       #endif

       if(acc != GOLDEN_VALUE)
         printf("dot product is is %d instead of %d\n",acc, GOLDEN_VALUE);
       else
         printf("Nice! Well done! 0 errors\n");
     }
    */
}// end main 

int main(void){
      printf("\n\n\t    *** Dot product test *** \n\n");
      return pmsis_kickoff((void *) fc_main);
}



