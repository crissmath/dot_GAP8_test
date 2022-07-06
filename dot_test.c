#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Gap8.h"
#include "pulp.h"
#include "pmsis.h"
#include "rt/rt_api.h"

#define V_size 4    
#define pi_pmu_set_voltage(x, y)      ( rt_voltage_force(RT_VOLTAGE_DOMAIN_MAIN, x, NULL) )
/*Select function*/
//#define use_simd
#define use_normal

/* Calculate value */  
unsigned int GOLDEN_VALUE = 20743432;

signed char *va;
signed char *vb;
uint32_t    *out;

typedef struct dot_arg{
    signed char *v_a;
    signed char *v_b;
    uint32_t    *acc;
    uint32_t    dim;
}dot_arg_t;

dot_arg_t Arg;

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
           V[i] = (signed char)i + off;
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
    uint32_t     N     = Arg->dim;

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

void cluster_delegate(void *arg)
{
    unsigned int time;
    uint32_t acc = 0;

    printf("Cluster master core entry\n");
    // vector allocation
    va  = (signed char *) pi_l2_malloc((uint32_t)(V_size * sizeof(signed char)));
    vb  = (signed char *) pi_l2_malloc((uint32_t)(V_size * sizeof(signed char)));
    out = (uint32_t *)    pi_l2_malloc((uint32_t)(sizeof(uint32_t)));
    if (va == NULL) {
    printf("buff alloc failed v_A!\n");
    pmsis_exit(-1);
    }
    if (vb == NULL) {
    printf("buff alloc failed v_B !\n");
    pmsis_exit(-1);
    }
    if (out == NULL) {
    printf("buff alloc failed v_C !\n");
    pmsis_exit(-1);
    }

    // init vector
    init_vec(va, V_size, 0);
    init_vec(vb, V_size, 10);
    out = &acc;

    printf("init vector\n");
    printf("VA\n");
    print_vec(va, V_size);
    printf("VB\n");
    print_vec(vb, V_size);
    printf("VC\n");
    printf("%d\n", *(out));

  /*
    Tasks to execute
  */
    printf("Run Task\n");
    Arg.v_a = va;
    Arg.v_b = vb;
    Arg.acc = out;
    Arg.dim = V_size;
    pi_cl_team_fork( 1, (void *)dummy, (void *) &Arg);

    printf("acc=%d\n", acc);

    /*
    time = pi_perf_read(PI_PERF_ACTIVE_CYCLES);
    printf("Computation done in %d cycles at %.2f operations per cycle....\n", 
            time, ((float) (V_size*V_size)/time));
    printf(" v_C = %d\n", v_C);
    */
 /*   printf("Run SIMD Dot product\n");
    pi_perf_conf(1 << PI_PERF_ACTIVE_CYCLES);
    pi_perf_reset(); 
    pi_perf_start();
    pi_cl_team_fork( 1, dotproduct_simd, (acc, v_A, v_B, V_size));
    pi_perf_stop();

    time = pi_perf_read(PI_PERF_ACTIVE_CYCLES);
    printf("Computation done in %d cycles at %.2f operations per cycle....\n", 
            time, ((float) (V_size*V_size)/time));
*/
}



int fc_main()
{
  printf("Entering main controller\n");
  uint32_t errors = 0;
  uint32_t core_id = pi_core_id();
  uint32_t cluster_id = pi_cluster_id();

  printf("[cl_id = %d, core_id = %d] Hello World!\n", (int) cluster_id, (int) core_id);

   /* 
        Set SoC Voltage to 1.2V
    */
    uint32_t voltage_in_mV = 1200;
    if(pi_pmu_set_voltage(voltage_in_mV, 0)==-1){
        printf("Frequency set failed!\n");
        pmsis_exit(-1);
    }
    /* 
        Set Cluster Freq at 250 MHz 
    */
    uint32_t fc_freq_in_Hz = 250 * 1000 * 1000;
    pi_freq_set(PI_FREQ_DOMAIN_FC, fc_freq_in_Hz);
    printf("Fabric Controller Frequency %d Hz\n", (int) pi_freq_get(PI_FREQ_DOMAIN_FC));
    /* 
        Configure & open cluster. 
    */
    struct pi_device cluster_dev   = {0};
    struct pi_cluster_conf cl_conf = {0};

    // Init cluster configuration structure.
    pi_cluster_conf_init(&cl_conf);
    cl_conf.id = 0;                             /* Set cluster ID. */
    pi_open_from_conf(&cluster_dev, &cl_conf);
    if (pi_cluster_open(&cluster_dev))
    {
        printf("Cluster open failed !\n");
        pmsis_exit(-1);
    }

    /* 
        Set the max freq for the cluster @1.2V
    */
    uint32_t cl_freq_in_Hz = 175 * 1000 * 1000;
    pi_freq_set(PI_FREQ_DOMAIN_CL, cl_freq_in_Hz);
    printf("Cluster Frequency %d Hz\n", (int) pi_freq_get(PI_FREQ_DOMAIN_CL));

    /* 
        Prepare cluster task and send it to cluster. 
    */
    struct pi_cluster_task cl_task = {0};
    cl_task.entry = cluster_delegate;
    cl_task.arg = NULL;
    // send task to cluster someone need odentifier 
    pi_cluster_send_task_to_cl(&cluster_dev, &cl_task);
    pi_cluster_close(&cluster_dev);

    // Terminate and exit the test
    printf("Test success !\n");
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



