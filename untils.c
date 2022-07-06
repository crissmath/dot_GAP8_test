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
/***************************************
*          Cluster delegate
*****************************************/