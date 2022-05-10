Modeling transmission of SARS-CoV-2 Omicron in China: instructions for
code running
================
Deng Xiaowei[1]
2022-05-10

This is the instruction for readers to repeat and understand the
construction for the model simulating the transmission of SARS-CoV-2
Omicron in China and its according health impact. The details for the
methodology and results refer to our paper published
(https://www.nature.com/articles/s41591-022-01855-7).

## Running platform and software version

If you want to repeat our demo, you must deploy **C** compiler and
[**GSL**](https://www.gnu.org/software/gsl/) (a library for C and C++).
The detialed version information refers to below.

    ## OS and hardware
    ## OS: Microsoft Windows 10 家庭中文版
    ## OS Version: 10.0.19044
    ## CPU: 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
    ## Memory size: 15.80 GB

    ## Soft version
    ## R version 4.1.0 (2021-05-18)
    ## cygwin (3.2.0-1)
    ## gcc-core (10.2.0-1)
    ## gsl (2.3-2)

## Linelist for **C** code, input and output

We separate our code, input, and output into three folders named in “C”,
“INPUT\_CN”, and “OUTPUT”.

Codes in folder “C” consist of three files, func.c, head.h, main.c.
“main.c” is the main body for the model, “func.c” includes all functions
used to read files, to process transmission, to allocated vaccine, and
“head.h” declares all relevant STRUCT and functions.

Most parameters used in our model are deposited in folder “INPUT\_CN”.
Detailed description refers to table below.

| Folder      | File                  | Dimentions(rows x columns) | Type   | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Limits                                                                                                                                                               |
|:------------|:----------------------|:---------------------------|:-------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CM/baseline | cm\_{i}               | 14x14                      | double | Contact matrix in baseline scenario setting                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | {i} in 0-199                                                                                                                                                         |
| POP         | contraindication\_{j} | 14x1                       | double | Proportion for people cannot be vacccinated                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | {j}=YES indicates limits for limits in vaccinating elderly according to real-world data                                                                              |
| POP         | {k}\_sus              | 200x14                     | double | Susceptibility for 200 max times simulation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | {k}=hetero indicates unequle susceptibility for 14 age-group                                                                                                         |
| POP         | population            | 14x1                       | int    | Population number                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |                                                                                                                                                                      |
| PROB\_SA    | p\_symp               | 14x1                       | double | Probability for developing symptom (adjusted)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |                                                                                                                                                                      |
| PROB\_SA    | p\_hosp               | 14x1                       | double | Probability for hospital addmission                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |                                                                                                                                                                      |
| PROB\_SA    | p\_icu                | 14x1                       | double | Probability for ICU addmission                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |                                                                                                                                                                      |
| PROB\_SA    | p\_death\_h           | 14x1                       | double | Probability for deaths occurring in hospitalized in common ward                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |                                                                                                                                                                      |
| PROB\_SA    | p\_death\_icu         | 14x1                       | double | Probability for deaths occurring in ICU                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |                                                                                                                                                                      |
| VAC         | VP\_strategy          | 14x4                       | int    | Age-group shift for primary vaccination                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |                                                                                                                                                                      |
| VAC         | VP\_ts                | 999x1                      | int    | Time series for supply of primary vaccination                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |                                                                                                                                                                      |
| VAC         | VB2{l}                | 14x1                       | int    | Target population for booster vaccination                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | {l}=18p indicates targetting at people above 18 years                                                                                                                |
| VAC         | VB\_ts                | 999x1                      | int    | Time series for supply of booster vaccination                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |                                                                                                                                                                      |
| VE          | VE\_{m}\_{n}\_{o}     | 7x12                       | double | VE (and transferring rate) for different outcome and vaccination status. Row from 1 to 7 refer to VE for infection, VE for onward transmission, VE for symptom, VE for hospitalization, VE for death, rate for vaccination transition, and rate for recovery. Column from 1 to 12 refer to susceptible population, people got the first shot, the first shot activate, people got the first shot and ready for the second shot, people got the second shot, the second shot activate, people got the second shot and ready for booster (waning status), people got booster vaccination, booster shot activate, the waning status after booster vaccination, people recovered from the Omicron infection, people recovered from the Delta infection | {m}=OPT indicates optimistic assumptions for VE; {n}=ONE indicates the relative waing rate for VE; {o}=Inactivated indicates the VE for inactivated booster schedule |
|             | Gamma\_{p}            | 7x1                        | double | Rate for status transition. Row 1-7 refer to rate for E2S, rate for I2R, rate for Is2H, rate for H2D, rate for H2R, rate for ICU2D, rate for ICU2R.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | {p}=L indicates longer generation time setting                                                                                                                       |

And, the folder “OUTPUT” consists all saved results. If you run the
code, three files will be created
(OPT=A\_baseline+new\_I\_Is\_H\_ICU\_D\_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta.txt,
OPT=A\_baseline+state\_S\_Is\_H\_ICU\_D\_byS-V1-V2-V3-V1E-V2E-V3E-V1W-V2W-V3W-R-RDelta.txt,
OPT=A\_baseline+vac\_V1\_V2\_VB.txt). The first two columns are index
for simulation and time. The file “new” indicates all variables are
daily new incidence, the order follows the compartment (I, Is, H, etc.),
then source (S, V1, V2, etc.), and then 14 age-group. The file “state”
indicates all variables are the temporary state, the following columns
are in similar order. The file “vac” refers to vaccination data, the
columns are ordered in vaccination status (V1, V2, etc.), then in 14
age-group.

## Compile **C** code

We implement the compiling process of **C** code in R using code below.

``` r
system("gcc -Wextra -Wall C/*.c -o SEIRS -lgsl")
```

    ## [1] 0

If you see 0 is returned, the compilation is successfully done. And you
may see some warning message as the ongoing chunk is showing. That is
the checkout we added in the function “FUNC\_transmission”, while we do
not need to check, the parameter is omitted and it will not hurt the
model running.

    ## [1] "C/func.c: In function 'FUNC_transmission':"                                  
    ## [2] "C/func.c:140:27: warning: unused parameter 'whe_output' [-Wunused-parameter]"
    ## [3] "  140 |                       int whe_output)"                               
    ## [4] "      |                       ~~~~^~~~~~~~~~"

## Execute the model code

As we deploy our environment for running **C** using
[**Cygwin**](http://www.cygwin.com/), we can execute the compiled **C**
program using syntax below in R.

``` r
n_sim = "10" # simulation times, max 200
t_simulate = "180" # time for simulation in days
input_folder = "INPUT_CN" # input folder for parameters
t_epistart = "457" # start date of epidemic, "457" indicates "2022-03-01"
n_initialinfectors = "20" # initial infectors introduced
R0 = "3.9" # reproduction number
booster_target = "18p" # target booster age-group
booster_regimen = "Inactivated" # booster vaccination type
booster_supply = "5000000" # booster supply since epidemic
booster_acceptance = "0.9" # rate for people accept to have booster shot
drug_uptake = "0" # proportion of drug uptake
drug_effect = "0.8" # effectiveness for drug
cm = "baseline" # contact matrix setting
sus_type = "hetero" # susceptibility setting
kappa = "1" # infectivity for asymptomatic cases
GT = "L" # rate for transmission setting, include generation time
RNG_seed = "1" # seed for random number generator, must be an int
VE_setting = "OPT_ONE" # VE setting, see table above for "INPUT_CN"
limits4elderly = "YES" # limits for vaccination for elderly
prob_folder = "SA" # folder name for probability parameter
save_prefix = "OUTPUT/OPT=A_baseline" # save path prefix
asymp_scale = "0.15" # adjusted ratio for asymptomatic cases
delta_infection = "0" # proportion for prior delta_infection

paste("./SEIRS", n_sim, t_simulate, input_folder, t_epistart,
      n_initialinfectors, R0, booster_target, booster_regimen,
      booster_supply, booster_acceptance, drug_uptake,
      drug_effect, cm, sus_type, kappa, GT, RNG_seed, 
      VE_setting, limits4elderly, prob_folder, save_prefix, 
      asymp_scale, delta_infection, sep = " ") -> tmp_command
print(tmp_command)
```

    ## [1] "./SEIRS 10 180 INPUT_CN 457 20 3.9 18p Inactivated 5000000 0.9 0 0.8 baseline hetero 1 L 1 OPT_ONE YES SA OUTPUT/OPT=A_baseline 0.15 0"

Once the code finishes, the following messages will be print on screen.

``` r
system(tmp_command, intern = T)
```

    ## [1] "SIM- +0: SIM- +1: SIM- +2: SIM- +3: SIM- +4: SIM- +5: SIM- +6: SIM- +7: SIM- +8: SIM- +9: "

[1] PhD candidate, School of Public Healh, Fudan University, Shanghai.
<sola1015@126.com>
