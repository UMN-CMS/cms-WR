/*********************************************************
  nn_1000_100.c
  --------------------------------------------------------
  generated at Tue Nov 30 10:25:11 2010
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/

#include <math.h>

#define Act_Logistic(sum, bias)  ( (sum+bias<10000.0) ? ( 1.0/(1.0 + exp(-sum-bias) ) ) : 0.0 )
#ifndef NULL
#define NULL (void *)0
#endif

typedef struct UT {
          float act;         /* Activation       */
          float Bias;        /* Bias of the Unit */
          int   NoOfSources; /* Number of predecessor units */
   struct UT   **sources; /* predecessor units */
          float *weights; /* weights from predecessor units */
        } UnitType, *pUnit;

  /* Forward Declaration for all unit types */
  static UnitType Units[12];
  /* Sources definition section */
  static pUnit Sources[] =  {
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 1, Units + 2, 
Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 

  };

  /* Weigths definition section */
  static float Weights[] =  {
3.994970, -1.334730, 
-4.047990, 1.401710, 
-4.582210, 2.679330, 
4.276280, -1.247670, 
4.600200, -1.809180, 
-2.638370, -3.166390, 
-3.356800, 7.479620, 
-3.490130, 76.019096, 
1.938330, -10.139990, -4.820500, 1.654340, 1.773820, 2.652800, -2.501840, -4.092090, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[12] = 
  {
    { 0.0, 0.0, 0, NULL , NULL },
    { /* unit 1 (Old: 1) */
      0.0, 0.116520, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 2 (Old: 2) */
      0.0, -0.902200, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 3 (Old: 3) */
      0.0, -0.678170, 2,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 4 (Old: 4) */
      0.0, 0.791190, 2,
       &Sources[2] , 
       &Weights[2] , 
      },
    { /* unit 5 (Old: 5) */
      0.0, 1.055550, 2,
       &Sources[4] , 
       &Weights[4] , 
      },
    { /* unit 6 (Old: 6) */
      0.0, -0.941540, 2,
       &Sources[6] , 
       &Weights[6] , 
      },
    { /* unit 7 (Old: 7) */
      0.0, -0.642600, 2,
       &Sources[8] , 
       &Weights[8] , 
      },
    { /* unit 8 (Old: 8) */
      0.0, 6.246410, 2,
       &Sources[10] , 
       &Weights[10] , 
      },
    { /* unit 9 (Old: 9) */
      0.0, 0.441250, 2,
       &Sources[12] , 
       &Weights[12] , 
      },
    { /* unit 10 (Old: 10) */
      0.0, -7.514060, 2,
       &Sources[14] , 
       &Weights[14] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, -0.494970, 8,
       &Sources[16] , 
       &Weights[16] , 
      }

  };



int heavyNuNN_1000_100(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};
  pUnit unit;


  /* layer definition section (names & member units) */

  static pUnit Input[2] = {Units + 1, Units + 2}; /* members */

  static pUnit Hidden1[8] = {Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10}; /* members */

  static pUnit Output1[1] = {Units + 11}; /* members */

  static int Output[1] = {11};

  for(member = 0; member < 2; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 8; member++) {
    unit = Hidden1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for (member = 0; member < 1; member++) {
    unit = Output1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for(member = 0; member < 1; member++) {
    out[member] = Units[Output[member]].act;
  }

  return(OK);
}
