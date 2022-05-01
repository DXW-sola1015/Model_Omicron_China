#include "head.h"

int main(int argc, char *argv[])
{
  char help[] = "PARAMTERS NOT ADEQUATE!!!";
  if (argc < 23)
  {
    fprintf(stderr, "%s\n", help);
    exit(1);
  }

  // SET ALL INPUT
  int SIM_MAX = atoi(argv[1]);
  int T_EPIDURATION = atoi(argv[2]);
  char *CITY_SELECTED = argv[3];

  int T_EPISTART = atoi(argv[4]);
  int N_initial_I = atoi(argv[5]);
  double BETA = atof(argv[6]);
  char *STRATEGY4BOOST = argv[7];
  char path_STRATEGY4BOOST[100];
  sprintf(path_STRATEGY4BOOST, "%s/VAC/VB2%s", CITY_SELECTED, STRATEGY4BOOST);
  char *TYPE4BOOST = argv[8];
  int BOOST_SUPPLY = atoi(argv[9]);
  double BOOST_COVERAGE = atof(argv[10]);
  double p_DRUG = atof(argv[11]);
  double rate4DRUG = atof(argv[12]);
  char *TYPE4CM = argv[13];
  char path_CM[100], path_CM4beta[100];
  char *TYPE4SUS = argv[14];
  char path_TYPE4SUS[100];
  sprintf(path_TYPE4SUS, "%s/POP/%s_sus", CITY_SELECTED, TYPE4SUS);
  double Kappa = atof(argv[15]);
  char *TYPE4GAMMA = argv[16];
  char path_TYPE4GAMMA[100];
  sprintf(path_TYPE4GAMMA, "%s/Gamma_%s", CITY_SELECTED, TYPE4GAMMA);
  int RNG_SEED = atoi(argv[17]);
  char *TYPE4VE = argv[18];
  char path_VE[100];
  sprintf(path_VE, "%s/VE/VE_%s_%s", CITY_SELECTED, TYPE4VE, TYPE4BOOST);
  char *WHE_60 = argv[19];
  char path_CONTRA[100];
  sprintf(path_CONTRA, "%s/POP/contraindication_%s", CITY_SELECTED, WHE_60);
  char *TYPE4OUTCOME = argv[20];
  char path_TYPE4OUTCOME[100];
  sprintf(path_TYPE4OUTCOME, "%s/PROB_%s/", CITY_SELECTED, TYPE4OUTCOME);
  char *path_SAVE = argv[21];
  double p_scale4symp = atof(argv[22]);
  int whe_delta = atoi(argv[23]);

  // TRANSMISSION RELATED PARAMETERS
  int T_STEP = 1, N[N_AGE], AGE_selected, whe_CS;
  int num_vac_trans = 0;
  double Beta, I_COUNT;

  // DECLARE ALL COMPARTMENTS
  STATE_SEIR STATE_S, STATE_V1, STATE_V1E, STATE_V1W, STATE_V2, STATE_V2E, STATE_V2W;
  STATE_SEIR STATE_V3, STATE_V3E, STATE_V3W, STATE_R, STATE_RDelta;

  // GENERATE ALL SAVE
  SAVE save_STATE_S, save_STATE_V1, save_STATE_V1E, save_STATE_V1W, save_STATE_V2, save_STATE_V2E, save_STATE_V2W;
  SAVE save_STATE_V3, save_STATE_V3E, save_STATE_V3W, save_STATE_R, save_STATE_RDelta;
  SAVE save_NEW_S, save_NEW_V1, save_NEW_V1E, save_NEW_V1W, save_NEW_V2, save_NEW_V2E, save_NEW_V2W;
  SAVE save_NEW_V3, save_NEW_V3E, save_NEW_V3W, save_NEW_R, save_NEW_RDelta;

  // VACCINE RELATED INPUT
  int VP_strategy[N_AGE][4], VB_strategy[N_AGE][1];
  int VP_index, VB_index, VP_dose, VB_dose, DOSE_REMAINED;
  int INPUT_VP_dose[T_EPISTART + T_EPIDURATION], INPUT_VB_dose[T_EPISTART + T_EPIDURATION];
  int V1_admin[N_AGE], V2_admin[N_AGE], VB_admin[N_AGE];

  char path_VP[100], path_VB[100], path_STRATEGY4VP[100];
  sprintf(path_VP, "%s/VAC/VP_ts", CITY_SELECTED);
  sprintf(path_VB, "%s/VAC/VB_ts", CITY_SELECTED);
  sprintf(path_STRATEGY4VP, "%s/VAC/VP_strategy", CITY_SELECTED);
  FUNC_read_col_int(path_VP, T_EPISTART + T_EPIDURATION, INPUT_VP_dose);
  FUNC_read_col_int(path_VB, T_EPISTART + T_EPIDURATION, INPUT_VB_dose);
  FUNC_read_matrix_int(path_STRATEGY4VP, N_AGE, 4, VP_strategy);
  FUNC_read_matrix_int(path_STRATEGY4BOOST, N_AGE, 1, VB_strategy);

  // POPULATION RELATED INPUT
  double RATE_contra[N_AGE], RATE_sus[SIM_MAX][N_AGE], CM[N_AGE][N_AGE], CM4Beta[N_AGE][N_AGE], Gamma[7];
  double Lambda[N_AGE], p_foi[N_AGE], p_V1[N_AGE], p_V1E[N_AGE], p_V1W[N_AGE];
  double p_V2[N_AGE], p_V2E[N_AGE], p_V2W[N_AGE];
  double p_V3[N_AGE], p_V3E[N_AGE], p_V3W[N_AGE];
  double p_R[N_AGE], p_RDelta[N_AGE];

  char path_N[100];
  sprintf(path_N, "%s/POP/population", CITY_SELECTED);
  FUNC_read_col_int(path_N, N_AGE, N);
  FUNC_read_col_double(path_CONTRA, N_AGE, RATE_contra);
  FUNC_read_matrix_double(path_TYPE4SUS, SIM_MAX, N_AGE, RATE_sus);

  // DISEASE PARAMETERS INPUT
  FUNC_read_col_double(path_TYPE4GAMMA, 7, Gamma);
  GAMMA day_Gamma = {Gamma[0] / (double)T_STEP, Gamma[1] / (double)T_STEP, Gamma[2] / (double)T_STEP,
                     Gamma[3] / (double)T_STEP, Gamma[4] / (double)T_STEP, Gamma[5] / (double)T_STEP,
                     Gamma[6] / (double)T_STEP};
  GAMMA p_Gamma = {1.0 - exp(-Gamma[0] / (double)T_STEP), 1.0 - exp(-Gamma[1] / (double)T_STEP),
                   1.0 - exp(-Gamma[2] / (double)T_STEP), 1.0 - exp(-Gamma[3] / (double)T_STEP),
                   1.0 - exp(-Gamma[4] / (double)T_STEP), 1.0 - exp(-Gamma[5] / (double)T_STEP),
                   1.0 - exp(-Gamma[6] / (double)T_STEP)};

  PROB_PROGRESS p_OUTCOME;
  double p_symp[N_AGE], p_hosp[N_AGE], p_icu[N_AGE], p_death_icu[N_AGE], p_death_hos[N_AGE];
  char path_SYMP[100], path_HOSP[100], path_ICU[100], path_DEATH_ICU[100], path_DEATH_HOSP[100];
  strcpy(path_SYMP, path_TYPE4OUTCOME);
  strcat(path_SYMP, "p_symp");
  strcpy(path_HOSP, path_TYPE4OUTCOME);
  strcat(path_HOSP, "p_hosp");
  strcpy(path_ICU, path_TYPE4OUTCOME);
  strcat(path_ICU, "p_icu");
  strcpy(path_DEATH_ICU, path_TYPE4OUTCOME);
  strcat(path_DEATH_ICU, "p_death_icu");
  strcpy(path_DEATH_HOSP, path_TYPE4OUTCOME);
  strcat(path_DEATH_HOSP, "p_death_h");
  FUNC_read_col_double(path_SYMP, N_AGE, p_symp);
  FUNC_read_col_double(path_HOSP, N_AGE, p_hosp);
  FUNC_read_col_double(path_ICU, N_AGE, p_icu);
  FUNC_read_col_double(path_DEATH_ICU, N_AGE, p_death_icu);
  FUNC_read_col_double(path_DEATH_HOSP, N_AGE, p_death_hos);

  for (int i = 0; i < N_AGE; i++)
  {
    p_OUTCOME.Symp[i] = p_scale4symp * p_symp[i];
    p_OUTCOME.Hosp[i] = p_hosp[i];
    p_OUTCOME.ICU[i] = p_icu[i];
    p_OUTCOME.Death_ICU[i] = p_death_icu[i];
    p_OUTCOME.Death_HOS[i] = p_death_hos[i];
    if (i <= 1)
    {
      p_OUTCOME.Drug[i] = 0.0;
    }
    else
    {
      p_OUTCOME.Drug[i] = p_DRUG;
    }
  }

  // VE INPUT
  VE VE_4_S, VE_4_V1, VE_4_V1E, VE_4_V1W, VE_4_V2, VE_4_V2E, VE_4_V2W;
  VE VE_4_V3, VE_4_V3E, VE_4_V3W, VE_4_R, VE_4_RDelta;
  double TMP_VE[7][12];
  FUNC_read_matrix_double(path_VE, 7, 12, TMP_VE);
  FUNC_read_VE(12, 0, TMP_VE, &VE_4_S);
  FUNC_read_VE(12, 1, TMP_VE, &VE_4_V1);
  FUNC_read_VE(12, 2, TMP_VE, &VE_4_V1E);
  FUNC_read_VE(12, 3, TMP_VE, &VE_4_V1W);
  FUNC_read_VE(12, 4, TMP_VE, &VE_4_V2);
  FUNC_read_VE(12, 5, TMP_VE, &VE_4_V2E);
  FUNC_read_VE(12, 6, TMP_VE, &VE_4_V2W);
  FUNC_read_VE(12, 7, TMP_VE, &VE_4_V3);
  FUNC_read_VE(12, 8, TMP_VE, &VE_4_V3E);
  FUNC_read_VE(12, 9, TMP_VE, &VE_4_V3W);
  FUNC_read_VE(12, 10, TMP_VE, &VE_4_R);
  FUNC_read_VE(12, 11, TMP_VE, &VE_4_RDelta);

  // DELTA COMPARMENT
  int N_Delta[123][N_AGE], AGE_SUM, AGE_EACH[8];
  char path_Delta[500];
  sprintf(path_Delta, "%s/Delta_inf", CITY_SELECTED);
  FUNC_read_matrix_int(path_Delta, 123, N_AGE, N_Delta);

  // INITIALIZE GSL RANDOM GENERATOR
  gsl_rng *RNG;
  RNG = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(RNG, RNG_SEED);

  // SAVE PATH AND DECLARATION
  char path_save1[500], path_save2[500], path_save3[500];
  sprintf(path_save1, "%s+new_I_Is_H_ICU_D_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta.txt", path_SAVE);
  sprintf(path_save2, "%s+state_S_Is_H_ICU_D_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta.txt", path_SAVE);
  sprintf(path_save3, "%s+vac_V1_V2_VB.txt", path_SAVE);

  FILE *f_save1;
  f_save1 = fopen(path_save1, "w");
  FILE *f_save2;
  f_save2 = fopen(path_save2, "w");
  FILE *f_save3;
  f_save3 = fopen(path_save3, "w");

  // SIM START
  for (int sim = 0; sim < SIM_MAX; sim++)
  {
    printf("SIM-%+3d: ", sim);

    sprintf(path_CM, "%s/CM/%s/cm_%d", CITY_SELECTED, TYPE4CM, sim);
    sprintf(path_CM4beta, "%s/CM/baseline/cm_%d", CITY_SELECTED, sim);
    FUNC_read_matrix_double(path_CM, N_AGE, N_AGE, CM);
    FUNC_read_matrix_double(path_CM4beta, N_AGE, N_AGE, CM4Beta);

    // INITIALIZE COMPARTMENT WITH 0
    for (int i = 0; i < N_AGE; i++)
    {
      V1_admin[i] = 0;
      V2_admin[i] = 0;
      VB_admin[i] = 0;
      FUNC_zero_SEIR(i, &STATE_S);
      FUNC_zero_SEIR(i, &STATE_V1);
      FUNC_zero_SEIR(i, &STATE_V1E);
      FUNC_zero_SEIR(i, &STATE_V1W);
      FUNC_zero_SEIR(i, &STATE_V2);
      FUNC_zero_SEIR(i, &STATE_V2E);
      FUNC_zero_SEIR(i, &STATE_V2W);
      FUNC_zero_SEIR(i, &STATE_V3);
      FUNC_zero_SEIR(i, &STATE_V3E);
      FUNC_zero_SEIR(i, &STATE_V3W);
      FUNC_zero_SEIR(i, &STATE_R);
      FUNC_zero_SEIR(i, &STATE_RDelta);

      FUNC_zero_SAVE(i, &save_STATE_S);
      FUNC_zero_SAVE(i, &save_STATE_V1);
      FUNC_zero_SAVE(i, &save_STATE_V1E);
      FUNC_zero_SAVE(i, &save_STATE_V1W);
      FUNC_zero_SAVE(i, &save_STATE_V2);
      FUNC_zero_SAVE(i, &save_STATE_V2E);
      FUNC_zero_SAVE(i, &save_STATE_V2W);
      FUNC_zero_SAVE(i, &save_STATE_V3);
      FUNC_zero_SAVE(i, &save_STATE_V3E);
      FUNC_zero_SAVE(i, &save_STATE_V3W);
      FUNC_zero_SAVE(i, &save_STATE_R);
      FUNC_zero_SAVE(i, &save_STATE_RDelta);

      FUNC_zero_SAVE(i, &save_NEW_S);
      FUNC_zero_SAVE(i, &save_NEW_V1);
      FUNC_zero_SAVE(i, &save_NEW_V1E);
      FUNC_zero_SAVE(i, &save_NEW_V1W);
      FUNC_zero_SAVE(i, &save_NEW_V2);
      FUNC_zero_SAVE(i, &save_NEW_V2E);
      FUNC_zero_SAVE(i, &save_NEW_V2W);
      FUNC_zero_SAVE(i, &save_NEW_V3);
      FUNC_zero_SAVE(i, &save_NEW_V3E);
      FUNC_zero_SAVE(i, &save_NEW_V3W);
      FUNC_zero_SAVE(i, &save_NEW_R);
      FUNC_zero_SAVE(i, &save_NEW_RDelta);

      STATE_S.C[i] = (int)(RATE_contra[i] * (double)N[i] + 0.5);
      STATE_S.S[i] = N[i] - STATE_S.C[i];

      for (int j = 0; j < N_AGE; j++)
      {
        CM[i][j] = CM[i][j] * RATE_sus[sim][i];
        CM4Beta[i][j] = CM4Beta[i][j] * RATE_sus[sim][i];
        CM4Beta[i][j] = ((Kappa * (1.0 - p_OUTCOME.Symp[j])) / day_Gamma.I2R + p_OUTCOME.Symp[j] / (day_Gamma.I2R * (1.0 - p_OUTCOME.Hosp[j]) + day_Gamma.Is2H * p_OUTCOME.Hosp[j])) * CM4Beta[i][j] * N[i] / N[j];
      }
    }

    Beta = BETA / FUNC_NGM(N_AGE, CM4Beta);

    VP_index = VB_index = 0;
    for (int t = 0; t < (T_EPISTART + T_EPIDURATION) * T_STEP; t++)
    {
      fprintf(f_save1, "%d %d", sim, t);
      fprintf(f_save2, "%d %d", sim, t);
      fprintf(f_save3, "%d %d", sim, t);

      // ALLOCATE INITIAL INFECTORS
      if (t == (T_EPISTART - 1) * T_STEP)
      {
        for (int i = 0; i < N_initial_I; i++)
        {
          AGE_selected = gsl_rng_uniform_int(RNG, N_AGE);
          whe_CS = gsl_rng_uniform_int(RNG, 2);
          if (whe_CS == 0)
          {
            while (STATE_S.C[AGE_selected] == 0)
            {
              AGE_selected = gsl_rng_uniform_int(RNG, N_AGE);
            }
            STATE_S.C[AGE_selected]--;
            STATE_S.Is_Drug_preR[AGE_selected]++;
            STATE_S.I[AGE_selected]++;
          }
          else
          {
            while (STATE_S.S[AGE_selected] == 0)
            {
              AGE_selected = gsl_rng_uniform_int(RNG, N_AGE);
            }
            STATE_S.S[AGE_selected]--;
            STATE_S.Is_Drug_preR[AGE_selected]++;
            STATE_S.I[AGE_selected]++;
          }
        }
      }

      // CALCULATE FOI
      for (int i = 0; i < N_AGE; i++)
      {
        Lambda[i] = 0;
        for (int j = 0; j < N_AGE; j++)
        {
          I_COUNT = STATE_S.I[j] + STATE_V1.I[j] + STATE_V1E.I[j] + STATE_V1W.I[j] +
                    STATE_V2.I[j] + STATE_V2E.I[j] + STATE_V2W.I[j] +
                    STATE_V3.I[j] + STATE_V3E.I[j] + STATE_V3W.I[j] +
                    STATE_R.I[j] + STATE_RDelta.I[j];
          Lambda[i] += CM[i][j] * Beta * I_COUNT / (double)N[j];
        }
        Lambda[i] = Lambda[i] / (double)T_STEP;
        p_foi[i] = 1.0 - exp(-VE_4_S.Inf * Lambda[i]);
        p_V1[i] = 1.0 - exp(-VE_4_V1.Inf * Lambda[i]);
        p_V1E[i] = 1.0 - exp(-VE_4_V1E.Inf * Lambda[i]);
        p_V1W[i] = 1.0 - exp(-VE_4_V1W.Inf * Lambda[i]);
        p_V2[i] = 1.0 - exp(-VE_4_V2.Inf * Lambda[i]);
        p_V2E[i] = 1.0 - exp(-VE_4_V2E.Inf * Lambda[i]);
        p_V2W[i] = 1.0 - exp(-VE_4_V2W.Inf * Lambda[i]);
        p_V3[i] = 1.0 - exp(-VE_4_V3.Inf * Lambda[i]);
        p_V3E[i] = 1.0 - exp(-VE_4_V3E.Inf * Lambda[i]);
        p_V3W[i] = 1.0 - exp(-VE_4_V3W.Inf * Lambda[i]);
        p_R[i] = 1.0 - exp(-VE_4_R.Inf * Lambda[i]);
        p_RDelta[i] = 1.0 - exp(-VE_4_RDelta.Inf * Lambda[i]);
      }

      // TRANSMISSION
      for (int i = 0; i < N_AGE; i++)
      {
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_S, VE_4_S, 1, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_S, &save_NEW_S, Kappa, p_foi, rate4DRUG, 1, 0);

        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V1, VE_4_V1, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V1, &save_NEW_V1, Kappa, p_V1, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V1E, VE_4_V1E, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V1E, &save_NEW_V1E, Kappa, p_V1E, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V1W, VE_4_V1W, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V1W, &save_NEW_V1W, Kappa, p_V1W, rate4DRUG, 1, 0);

        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V2, VE_4_V2, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V2, &save_NEW_V2, Kappa, p_V2, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V2E, VE_4_V2E, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V2E, &save_NEW_V2E, Kappa, p_V2E, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V2W, VE_4_V2W, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V2W, &save_NEW_V2W, Kappa, p_V2W, rate4DRUG, BOOST_COVERAGE, 1);

        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V3, VE_4_V3, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V3, &save_NEW_V3, Kappa, p_V3, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V3E, VE_4_V3E, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V3E, &save_NEW_V3E, Kappa, p_V3E, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_V3W, VE_4_V3W, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_V3W, &save_NEW_V3W, Kappa, p_V3W, rate4DRUG, 1, 0);

        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_R, VE_4_R, 0, &STATE_R, p_Gamma,
                                          p_OUTCOME, &save_STATE_R, &save_NEW_R, Kappa, p_R, rate4DRUG, 1, 0);
        num_vac_trans = FUNC_transmission(RNG, i, T_STEP, num_vac_trans, &STATE_RDelta, VE_4_RDelta, 0, &STATE_RDelta, p_Gamma,
                                          p_OUTCOME, &save_STATE_RDelta, &save_NEW_RDelta, Kappa, p_RDelta, rate4DRUG, 1, 0);
      }

      // VACCINATION
      // CHANGE VACCINE STRATEGY BY TIME
      if (t == (119 - 1) * T_STEP)
      {
        VP_index = 1;
      }
      else if (t == (244 - 1) * T_STEP)
      {
        VP_index = 2;
      }
      else if (t == (332 - 1) * T_STEP)
      {
        VP_index = 3;
      }

      VP_dose = (int)((double)INPUT_VP_dose[(int)((double)t / (double)T_STEP)] / (double)T_STEP);
      if (t >= (T_EPISTART - 1) * T_STEP)
      {
        VB_dose = (int)((double)BOOST_SUPPLY / (double)T_STEP);
      }
      else
      {
        VB_dose = (int)((double)INPUT_VB_dose[(int)((double)t / (double)T_STEP)] / (double)T_STEP);
      }

      // BOOSTER DOSE
      DOSE_REMAINED = FUNC_vac_allocate(VB_dose, 1, VB_index, VB_strategy, &STATE_V2W, &STATE_V3, VB_admin);
      // 2nd DOSE
      DOSE_REMAINED = FUNC_vac_allocate(VP_dose, 4, VP_index, VP_strategy, &STATE_V1W, &STATE_V2, V2_admin);
      // 1st DOSE
      VP_dose -= DOSE_REMAINED;
      if (VP_dose < 0)
      {
        VP_dose = 0;
      }
      DOSE_REMAINED = FUNC_vac_allocate(VP_dose, 4, VP_index, VP_strategy, &STATE_S, &STATE_V1, V1_admin);

      // ALLOCATE DELTA INFECTIONS
      if (whe_delta == 1)
      {
        if (t >= (259 - 1) * T_STEP && t <= (381 - 1) * T_STEP)
        {
          for (int i = 0; i < N_AGE; i++)
          {
            AGE_EACH[0] = STATE_S.S[i];
            AGE_EACH[1] = STATE_V1.S[i];
            AGE_EACH[2] = STATE_V1E.S[i];
            AGE_EACH[3] = STATE_V1W.S[i];
            AGE_EACH[4] = STATE_V2.S[i];
            AGE_EACH[5] = STATE_V2E.S[i];
            AGE_EACH[6] = STATE_V2W.S[i];
            AGE_EACH[7] = STATE_R.S[i];

            AGE_SUM = 0;
            for (int j = 0; j < 8; j++)
            {
              AGE_SUM += AGE_EACH[j];
            }

            if (AGE_SUM != 0)
            {
              FUNC_Delta(t, T_STEP, i, 0, AGE_EACH, AGE_SUM, N_Delta, &STATE_S, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 1, AGE_EACH, AGE_SUM, N_Delta, &STATE_V1, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 2, AGE_EACH, AGE_SUM, N_Delta, &STATE_V1E, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 3, AGE_EACH, AGE_SUM, N_Delta, &STATE_V1W, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 4, AGE_EACH, AGE_SUM, N_Delta, &STATE_V2, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 5, AGE_EACH, AGE_SUM, N_Delta, &STATE_V2E, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 6, AGE_EACH, AGE_SUM, N_Delta, &STATE_V2W, &STATE_RDelta);
              FUNC_Delta(t, T_STEP, i, 7, AGE_EACH, AGE_SUM, N_Delta, &STATE_R, &STATE_RDelta);
            }
          }
        }
      }

      // PRINT TO SAVE
      for (int i = 0; i < N_AGE; i++)
      {
        // new_I_Is_H_ICU_D_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta
        fprintf(f_save1, " %d %d %d %d", save_NEW_S.I[i], save_NEW_V1.I[i], save_NEW_V2.I[i], save_NEW_V3.I[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1E.I[i], save_NEW_V2E.I[i], save_NEW_V3E.I[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1W.I[i], save_NEW_V2W.I[i], save_NEW_V3W.I[i]);
        fprintf(f_save1, " %d %d", save_NEW_R.I[i], save_NEW_RDelta.I[i]);

        fprintf(f_save1, " %d %d %d %d", save_NEW_S.Is[i], save_NEW_V1.Is[i], save_NEW_V2.Is[i], save_NEW_V3.Is[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1E.Is[i], save_NEW_V2E.Is[i], save_NEW_V3E.Is[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1W.Is[i], save_NEW_V2W.Is[i], save_NEW_V3W.Is[i]);
        fprintf(f_save1, " %d %d", save_NEW_R.Is[i], save_NEW_RDelta.Is[i]);

        fprintf(f_save1, " %d %d %d %d", save_NEW_S.H[i], save_NEW_V1.H[i], save_NEW_V2.H[i], save_NEW_V3.H[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1E.H[i], save_NEW_V2E.H[i], save_NEW_V3E.H[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1W.H[i], save_NEW_V2W.H[i], save_NEW_V3W.H[i]);
        fprintf(f_save1, " %d %d", save_NEW_R.H[i], save_NEW_RDelta.H[i]);

        fprintf(f_save1, " %d %d %d %d", save_NEW_S.ICU[i], save_NEW_V1.ICU[i], save_NEW_V2.ICU[i], save_NEW_V3.ICU[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1E.ICU[i], save_NEW_V2E.ICU[i], save_NEW_V3E.ICU[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1W.ICU[i], save_NEW_V2W.ICU[i], save_NEW_V3W.ICU[i]);
        fprintf(f_save1, " %d %d", save_NEW_R.ICU[i], save_NEW_RDelta.ICU[i]);

        fprintf(f_save1, " %d %d %d %d", save_NEW_S.D[i], save_NEW_V1.D[i], save_NEW_V2.D[i], save_NEW_V3.D[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1E.D[i], save_NEW_V2E.D[i], save_NEW_V3E.D[i]);
        fprintf(f_save1, " %d %d %d", save_NEW_V1W.D[i], save_NEW_V2W.D[i], save_NEW_V3W.D[i]);
        fprintf(f_save1, " %d %d", save_NEW_R.D[i], save_NEW_RDelta.D[i]);

        save_STATE_S.S[i] = STATE_S.S[i] + STATE_S.C[i];
        save_STATE_V1.S[i] = STATE_V1.S[i] + STATE_V1.C[i];
        save_STATE_V1E.S[i] = STATE_V1E.S[i] + STATE_V1E.C[i];
        save_STATE_V1W.S[i] = STATE_V1W.S[i] + STATE_V1W.C[i];
        save_STATE_V2.S[i] = STATE_V2.S[i] + STATE_V2.C[i];
        save_STATE_V2E.S[i] = STATE_V2E.S[i] + STATE_V2E.C[i];
        save_STATE_V2W.S[i] = STATE_V2W.S[i] + STATE_V2W.C[i];
        save_STATE_V3.S[i] = STATE_V3.S[i] + STATE_V3.C[i];
        save_STATE_V3E.S[i] = STATE_V3E.S[i] + STATE_V3E.C[i];
        save_STATE_V3W.S[i] = STATE_V3W.S[i] + STATE_V3W.C[i];
        save_STATE_R.S[i] = STATE_R.S[i] + STATE_R.C[i];
        save_STATE_RDelta.S[i] = STATE_RDelta.S[i] + STATE_RDelta.C[i];

        // state_S_Is_H_ICU_D_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta
        fprintf(f_save2, " %d %d %d %d", save_STATE_S.S[i], save_STATE_V1.S[i], save_STATE_V2.S[i], save_STATE_V3.S[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1E.S[i], save_STATE_V2E.S[i], save_STATE_V3E.S[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1W.S[i], save_STATE_V2W.S[i], save_STATE_V3W.S[i]);
        fprintf(f_save2, " %d %d", save_STATE_R.S[i], save_STATE_RDelta.S[i]);

        fprintf(f_save2, " %d %d %d %d", save_STATE_S.Is[i], save_STATE_V1.Is[i], save_STATE_V2.Is[i], save_STATE_V3.Is[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1E.Is[i], save_STATE_V2E.Is[i], save_STATE_V3E.Is[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1W.Is[i], save_STATE_V2W.Is[i], save_STATE_V3W.Is[i]);
        fprintf(f_save2, " %d %d", save_STATE_R.Is[i], save_STATE_RDelta.Is[i]);

        fprintf(f_save2, " %d %d %d %d", save_STATE_S.H[i], save_STATE_V1.H[i], save_STATE_V2.H[i], save_STATE_V3.H[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1E.H[i], save_STATE_V2E.H[i], save_STATE_V3E.H[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1W.H[i], save_STATE_V2W.H[i], save_STATE_V3W.H[i]);
        fprintf(f_save2, " %d %d", save_STATE_R.H[i], save_STATE_RDelta.H[i]);

        fprintf(f_save2, " %d %d %d %d", save_STATE_S.ICU[i], save_STATE_V1.ICU[i], save_STATE_V2.ICU[i], save_STATE_V3.ICU[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1E.ICU[i], save_STATE_V2E.ICU[i], save_STATE_V3E.ICU[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1W.ICU[i], save_STATE_V2W.ICU[i], save_STATE_V3W.ICU[i]);
        fprintf(f_save2, " %d %d", save_STATE_R.ICU[i], save_STATE_RDelta.ICU[i]);

        fprintf(f_save2, " %d %d %d %d", save_STATE_S.D[i], save_STATE_V1.D[i], save_STATE_V2.D[i], save_STATE_V3.D[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1E.D[i], save_STATE_V2E.D[i], save_STATE_V3E.D[i]);
        fprintf(f_save2, " %d %d %d", save_STATE_V1W.D[i], save_STATE_V2W.D[i], save_STATE_V3W.D[i]);
        fprintf(f_save2, " %d %d", save_STATE_R.D[i], save_STATE_RDelta.D[i]);

        // vac_V1_V2_VB
        fprintf(f_save3, " %d %d %d", V1_admin[i], V2_admin[i], VB_admin[i]);
      }
      fprintf(f_save1, "\n");
      fprintf(f_save2, "\n");
      fprintf(f_save3, "\n");
    }
  }

  fclose(f_save1);
  fclose(f_save2);
  fclose(f_save3);

  return 0;
}