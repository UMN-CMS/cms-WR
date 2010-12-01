/*********************************************************
  nn_1500_1000.c
  --------------------------------------------------------
  generated at Tue Nov 30 10:39:46 2010
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
-4.703830, -0.145050, 
-4.070430, -2.229690, 
-0.196760, -3.634840, 
-4.781370, -0.142890, 
-5.322780, 0.123080, 
-4.796770, -0.810860, 
1.505670, 3.635490, 
-11.034010, -1.668980, 
-5.312420, -1.660800, -6.375520, -1.602010, -11.711080, -4.005900, 2.424420, 15.228810, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[12] = 
  {
    { 0.0, 0.0, 0, NULL , NULL },
    { /* unit 1 (Old: 1) */
      0.0, 0.920510, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 2 (Old: 2) */
      0.0, -0.889600, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 3 (Old: 3) */
      0.0, 0.707680, 2,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 4 (Old: 4) */
      0.0, 1.561420, 2,
       &Sources[2] , 
       &Weights[2] , 
      },
    { /* unit 5 (Old: 5) */
      0.0, 1.434830, 2,
       &Sources[4] , 
       &Weights[4] , 
      },
    { /* unit 6 (Old: 6) */
      0.0, 0.626490, 2,
       &Sources[6] , 
       &Weights[6] , 
      },
    { /* unit 7 (Old: 7) */
      0.0, -0.622160, 2,
       &Sources[8] , 
       &Weights[8] , 
      },
    { /* unit 8 (Old: 8) */
      0.0, 2.727360, 2,
       &Sources[10] , 
       &Weights[10] , 
      },
    { /* unit 9 (Old: 9) */
      0.0, -0.671410, 2,
       &Sources[12] , 
       &Weights[12] , 
      },
    { /* unit 10 (Old: 10) */
      0.0, 1.984020, 2,
       &Sources[14] , 
       &Weights[14] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, 0.981720, 8,
       &Sources[16] , 
       &Weights[16] , 
      }

  };



int heavyNuNN_1500_1000(float *in, float *out, int init)
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
