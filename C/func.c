#include "head.h"

double FUNC_NGM(int n_col, double input_CM[][n_col])
{
  double input_data[n_col * n_col];
  double tmp, output_return;
  int Index_cell = 0;

  for (int i = 0; i < n_col; i++)
  {
    for (int j = 0; j < n_col; j++)
    {
      input_data[Index_cell] = input_CM[i][j];
      Index_cell++;
    }
  }

  gsl_matrix_view m = gsl_matrix_view_array(input_data, n_col, n_col);
  gsl_vector_complex *eval = gsl_vector_complex_alloc(n_col);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n_col, n_col);
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n_col);
  gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free(w);
  gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  gsl_complex eval_value = gsl_vector_complex_get(eval, 0);
  output_return = GSL_REAL(eval_value);

  for (int i = 1; i < n_col; i++)
  {
    gsl_complex eval_value = gsl_vector_complex_get(eval, i);
    tmp = GSL_REAL(eval_value);
    if (output_return < tmp)
    {
      output_return = tmp;
    }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return output_return;
}

void FUNC_read_col_int(char *input_path, int n_row, int output_return[])
{
  FILE *fp;
  fp = fopen(input_path, "r");
  for (int i = 0; i < n_row; i++)
  {
    fscanf(fp, "%d", &output_return[i]);
  }
  fclose(fp);
}

void FUNC_read_col_double(char *input_path, int n_row, double output_return[])
{
  FILE *fp;
  fp = fopen(input_path, "r");
  for (int i = 0; i < n_row; i++)
  {
    fscanf(fp, "%lf", &output_return[i]);
  }
  fclose(fp);
}

void FUNC_read_matrix_int(char *input_path, int n_row, int n_col,
                          int output_return[][n_col])
{
  FILE *fp;
  fp = fopen(input_path, "r");
  for (int i = 0; i < n_row; i++)
  {
    for (int j = 0; j < n_col; j++)
    {
      fscanf(fp, "%d", &output_return[i][j]);
    }
  }
  fclose(fp);
}

void FUNC_read_matrix_double(char *input_path, int n_row, int n_col,
                             double output_return[][n_col])
{
  FILE *fp;
  fp = fopen(input_path, "r");
  for (int i = 0; i < n_row; i++)
  {
    for (int j = 0; j < n_col; j++)
    {
      fscanf(fp, "%lf", &output_return[i][j]);
    }
  }
  fclose(fp);
}

void FUNC_read_VE(int n_col, int col_s, double input_VE[][n_col], VE *out_VE)
{
  (*out_VE).Inf = input_VE[0][col_s];
  (*out_VE).OT = input_VE[1][col_s];
  (*out_VE).Symp = input_VE[2][col_s];
  (*out_VE).Hosp = input_VE[3][col_s];
  (*out_VE).Death = input_VE[4][col_s];
  (*out_VE).rate4V = input_VE[5][col_s];
  (*out_VE).rate4R2S = input_VE[6][col_s];
}

void FUNC_zero_SEIR(int input_i, STATE_SEIR *input_SEIR)
{
  (*input_SEIR).S[input_i] = 0;
  (*input_SEIR).C[input_i] = 0;
  (*input_SEIR).E[input_i] = 0;
  (*input_SEIR).I[input_i] = 0.0;
  (*input_SEIR).Ia[input_i] = 0;
  (*input_SEIR).Is_Drug_preO[input_i] = 0;
  (*input_SEIR).Is_Drug_preR[input_i] = 0;
  (*input_SEIR).H_preO[input_i] = 0;
  (*input_SEIR).H_preR[input_i] = 0;
  (*input_SEIR).ICU_preO[input_i] = 0;
  (*input_SEIR).ICU_preR[input_i] = 0;
  (*input_SEIR).D[input_i] = 0;
  (*input_SEIR).R[input_i] = 0;
}

void FUNC_zero_SAVE(int input_i, SAVE *input_save)
{
  (*input_save).S[input_i] = 0;
  (*input_save).I[input_i] = 0;
  (*input_save).Is[input_i] = 0;
  (*input_save).DRUG[input_i] = 0;
  (*input_save).H[input_i] = 0;
  (*input_save).ICU[input_i] = 0;
  (*input_save).D[input_i] = 0;
}

int FUNC_transmission(const gsl_rng *input_RNG, int input_i, int input_tstep, int input_n_tran,
                      STATE_SEIR *input_SEIR, VE input_VE, int whe_R, STATE_SEIR *target_SEIR,
                      GAMMA input_GAMMA, PROB_PROGRESS input_P, SAVE *output_STATE, SAVE *output_NEW,
                      double input_Kappa, double input_FOI[], double input_DRUGeffect, double input_rate4unwill,
                      int whe_output)
{
  int R2S, S2TRAN, S2E, C2E, E2I, E2Ia, E2Is, IsDRUG, IsnonDRUG;
  int IsDRUG2preO, IsDRUG2preR, IsnonDRUG2preO, IsnonDRUG2preR;
  int Ia2R, Is2R, Is2O, Is2ICU, Is2H, ICU2ICUpreD, ICU2ICUpreR, ICU2R, ICU2D;
  int H2HpreD, H2HpreR, H2R, H2D, UNWILL;
  int tmp_return;

  // R2S
  R2S = gsl_ran_binomial(input_RNG, 1.0 - exp(-input_VE.rate4R2S / (double)input_tstep), (*input_SEIR).R[input_i]);

  // S2E, C2E, AND TRANSITION
  S2TRAN = gsl_ran_binomial(input_RNG, 1.0 - exp(-input_VE.rate4V / (double)input_tstep), (*input_SEIR).S[input_i]);
  S2E = gsl_ran_binomial(input_RNG, input_FOI[input_i], (*input_SEIR).S[input_i]);
  if (S2TRAN != 0 && (S2TRAN + S2E) > (*input_SEIR).S[input_i])
  {
    S2TRAN = (int)((double)S2TRAN / (double)(S2TRAN + S2E));
    S2E = (*input_SEIR).S[input_i] - S2TRAN;
  }
  tmp_return = S2TRAN;
  C2E = gsl_ran_binomial(input_RNG, input_FOI[input_i], (*input_SEIR).C[input_i]);

  // E2Ia, E2IS, IsDRUG AND IsnonDRUG
  E2I = gsl_ran_binomial(input_RNG, input_GAMMA.E2I, (*input_SEIR).E[input_i]);
  E2Is = gsl_ran_binomial(input_RNG, input_P.Symp[input_i] * input_VE.Symp, E2I);
  E2Ia = E2I - E2Is;
  IsDRUG = gsl_ran_binomial(input_RNG, input_P.Drug[input_i], E2Is);
  IsnonDRUG = E2Is - IsDRUG;

  // Ia AND Is TO PREOURCOME
  IsDRUG2preO = gsl_ran_binomial(input_RNG, input_P.Hosp[input_i] * input_VE.Hosp * (1.0 - input_DRUGeffect), IsDRUG);
  IsDRUG2preR = IsDRUG - IsDRUG2preO;
  IsnonDRUG2preO = gsl_ran_binomial(input_RNG, input_P.Hosp[input_i] * input_VE.Hosp, IsnonDRUG);
  IsnonDRUG2preR = IsnonDRUG - IsnonDRUG2preO;

  // Ia AND Is TO OUTCOME
  Is2R = gsl_ran_binomial(input_RNG, input_GAMMA.I2R, (*input_SEIR).Is_Drug_preR[input_i]);
  Is2O = gsl_ran_binomial(input_RNG, input_GAMMA.Is2H, (*input_SEIR).Is_Drug_preO[input_i]);
  Is2ICU = gsl_ran_binomial(input_RNG, input_P.ICU[input_i], Is2O);
  Is2H = Is2O - Is2ICU;

  // OUTCOME TO D AND R
  Ia2R = gsl_ran_binomial(input_RNG, input_GAMMA.I2R, (*input_SEIR).Ia[input_i]);

  ICU2ICUpreD = gsl_ran_binomial(input_RNG, input_P.Death_ICU[input_i] * input_VE.Death, Is2ICU);
  ICU2ICUpreR = Is2ICU - ICU2ICUpreD;
  ICU2R = gsl_ran_binomial(input_RNG, input_GAMMA.ICU2R, (*input_SEIR).ICU_preR[input_i]);
  ICU2D = gsl_ran_binomial(input_RNG, input_GAMMA.ICU2D, (*input_SEIR).ICU_preO[input_i]);

  H2HpreD = gsl_ran_binomial(input_RNG, input_P.Death_HOS[input_i] * input_VE.Death, Is2H);
  H2HpreR = Is2H - H2HpreD;
  H2R = gsl_ran_binomial(input_RNG, input_GAMMA.H2R, (*input_SEIR).H_preR[input_i]);
  H2D = gsl_ran_binomial(input_RNG, input_GAMMA.H2D, (*input_SEIR).H_preO[input_i]);
  
  // if (whe_output == 1)
  // {
  //   printf("!");
  // }

  // PROGRESS
  UNWILL = (int)((1.0 - input_rate4unwill) * (double)input_n_tran + 0.5);
  (*input_SEIR).S[input_i] += (input_n_tran - UNWILL);
  (*input_SEIR).C[input_i] += UNWILL;

  (*input_SEIR).S[input_i] = (*input_SEIR).S[input_i] - S2TRAN - S2E + R2S;
  (*input_SEIR).C[input_i] = (*input_SEIR).C[input_i] - C2E;
  (*input_SEIR).E[input_i] = (*input_SEIR).E[input_i] + S2E + C2E - E2I;
  (*input_SEIR).Ia[input_i] = (*input_SEIR).Ia[input_i] + E2Ia - Ia2R;
  (*input_SEIR).Is_Drug_preO[input_i] = (*input_SEIR).Is_Drug_preO[input_i] + IsDRUG2preO + IsnonDRUG2preO - Is2O;
  (*input_SEIR).Is_Drug_preR[input_i] = (*input_SEIR).Is_Drug_preR[input_i] + IsDRUG2preR + IsnonDRUG2preR - Is2R;
  (*input_SEIR).ICU_preO[input_i] = (*input_SEIR).ICU_preO[input_i] + ICU2ICUpreD - ICU2D;
  (*input_SEIR).ICU_preR[input_i] = (*input_SEIR).ICU_preR[input_i] + ICU2ICUpreR - ICU2R;
  (*input_SEIR).H_preO[input_i] = (*input_SEIR).H_preO[input_i] + H2HpreD - H2D;
  (*input_SEIR).H_preR[input_i] = (*input_SEIR).H_preR[input_i] + H2HpreR - H2R;
  (*input_SEIR).D[input_i] = (*input_SEIR).D[input_i] + H2D + ICU2D;
  (*input_SEIR).R[input_i] = (*input_SEIR).R[input_i] + H2R + ICU2R + Is2R + Ia2R - R2S;

  if (whe_R == 1)
  {
    (*target_SEIR).S[input_i] = (*target_SEIR).S[input_i] + H2R + ICU2R + Is2R + Ia2R;
  }

  (*input_SEIR).I[input_i] = input_VE.OT * ((double)(*input_SEIR).Is_Drug_preR[input_i] +
                                            (double)(*input_SEIR).Is_Drug_preO[input_i] +
                                            input_Kappa * (double)(*input_SEIR).Ia[input_i]);

  // SAVE
  (*output_STATE).Is[input_i] = (*input_SEIR).Is_Drug_preO[input_i] + (*input_SEIR).Is_Drug_preR[input_i];
  (*output_STATE).H[input_i] = (*input_SEIR).H_preO[input_i] + (*input_SEIR).H_preR[input_i];
  (*output_STATE).ICU[input_i] = (*input_SEIR).ICU_preO[input_i] + (*input_SEIR).ICU_preR[input_i];
  (*output_STATE).D[input_i] = (*input_SEIR).D[input_i];

  (*output_NEW).I[input_i] = E2I;
  (*output_NEW).Is[input_i] = E2Is;
  (*output_NEW).H[input_i] = H2HpreR + H2HpreD;
  (*output_NEW).ICU[input_i] = ICU2ICUpreR + ICU2ICUpreD;
  (*output_NEW).D[input_i] = H2D + ICU2D;
  (*output_NEW).DRUG[input_i] = IsDRUG;
  return tmp_return;
}

int FUNC_vac_allocate(int input_dose, int dim2, int input_index, int input_strategy[][dim2],
                      STATE_SEIR *target_pop, STATE_SEIR *vac_pop, int out_admin[])
{
  int tmp_return = 0;
  int AVA[N_AGE], AVA_TOTAL = 0, tmp_admin;

  for (int i = 0; i < N_AGE; i++)
  {
    AVA[i] = (*target_pop).S[i] * input_strategy[i][input_index];
    AVA_TOTAL += AVA[i];
  }

  if (AVA_TOTAL > input_dose)
  {
    for (int i = 0; i < N_AGE; i++)
    {
      tmp_admin = (int)((double)input_dose * ((double)AVA[i] / (double)AVA_TOTAL) + 0.5);
      if (tmp_admin <= AVA[i])
      {
        out_admin[i] = tmp_admin;
        tmp_return += tmp_admin;
        (*target_pop).S[i] -= tmp_admin;
        (*vac_pop).S[i] += tmp_admin;
      }
      else
      {
        tmp_admin = AVA[i];
        out_admin[i] = tmp_admin;
        tmp_return += tmp_admin;
        (*target_pop).S[i] -= tmp_admin;
        (*vac_pop).S[i] += tmp_admin;
      }
    }
  }
  else
  {
    for (int i = 0; i < N_AGE; i++)
    {
      tmp_admin = AVA[i];
      out_admin[i] = tmp_admin;
      tmp_return += tmp_admin;
      (*target_pop).S[i] -= tmp_admin;
      (*vac_pop).S[i] += tmp_admin;
    }
  }

  return tmp_return;
}

void FUNC_Delta(int input_t, int input_tstep, int input_i, int input_j, int input_num[], int input_denom,
                int input_Delta[][N_AGE], STATE_SEIR *input_SEIR, STATE_SEIR *Delta_SEIR)
{
  int n_trans;

  n_trans = (int)(input_Delta[input_t - (259 - 1) * input_tstep][input_i] * ((double)input_num[input_j] / (double)input_denom) + 0.5);
  if (n_trans > (*input_SEIR).S[input_i])
  {
    n_trans = (*input_SEIR).S[input_i];
    (*input_SEIR).S[input_i] -= n_trans;
    (*Delta_SEIR).S[input_i] += n_trans;
  }
  else
  {
    (*input_SEIR).S[input_i] -= n_trans;
    (*Delta_SEIR).S[input_i] += n_trans;
  }
}
