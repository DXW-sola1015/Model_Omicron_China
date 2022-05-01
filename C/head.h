#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define N_AGE 14

typedef struct
{
  int S[N_AGE];
  int C[N_AGE];
  int E[N_AGE];
  double I[N_AGE];
  int Ia[N_AGE];
  int Is_Drug_preO[N_AGE];
  int Is_Drug_preR[N_AGE];
  int H_preO[N_AGE];
  int H_preR[N_AGE];
  int ICU_preO[N_AGE];
  int ICU_preR[N_AGE];
  int D[N_AGE];
  int R[N_AGE];
} STATE_SEIR;

typedef struct
{
  double E2I;
  double I2R;
  double Is2H;
  double H2D;
  double H2R;
  double ICU2D;
  double ICU2R;
} GAMMA;

typedef struct
{
  double Symp[N_AGE];
  double Drug[N_AGE];
  double Hosp[N_AGE];
  double ICU[N_AGE];
  double Death_ICU[N_AGE];
  double Death_HOS[N_AGE];

} PROB_PROGRESS;

typedef struct
{
  double Inf;
  double Symp;
  double Hosp;
  double Death;
  double OT;
  double rate4R2S;
  double rate4V;
} VE;

typedef struct
{
  int S[N_AGE];
  int I[N_AGE];
  int Is[N_AGE];
  int DRUG[N_AGE];
  int H[N_AGE];
  int ICU[N_AGE];
  int D[N_AGE];
} SAVE;

double FUNC_NGM(int n_col, double input_CM[][n_col]);

void FUNC_read_col_int(char *input_path, int n_row, int output_return[]);

void FUNC_read_col_double(char *input_path, int n_row, double output_return[]);

void FUNC_read_matrix_int(char *input_path, int n_row, int n_col,
                          int output_return[][n_col]);

void FUNC_read_matrix_double(char *input_path, int n_row, int n_col,
                             double output_return[][n_col]);

void FUNC_zero_SEIR(int input_i, STATE_SEIR *input_SEIR);

void FUNC_zero_SAVE(int input_i, SAVE *input_save);

void FUNC_read_VE(int n_col, int col_s, double input_VE[][n_col], VE *out_VE);

int FUNC_transmission(const gsl_rng *input_RNG, int input_i, int input_tstep, int input_n_tran,
                      STATE_SEIR *input_SEIR, VE input_VE, int whe_R, STATE_SEIR *target_SEIR,
                      GAMMA input_GAMMA, PROB_PROGRESS input_P, SAVE *output_STATE, SAVE *output_NEW,
                      double input_Kappa, double input_FOI[], double input_DRUGeffect, double input_rate4unwill,
                      int whe_output);

int FUNC_vac_allocate(int input_dose, int dim2, int input_index, int input_strategy[][dim2],
                      STATE_SEIR *target_pop, STATE_SEIR *vac_pop, int out_admin[]);

void FUNC_Delta(int input_t, int input_tstep, int input_i, int input_j, int input_num[], int input_denom,
                int input_Delta[][N_AGE], STATE_SEIR *input_SEIR, STATE_SEIR *Delta_SEIR);
