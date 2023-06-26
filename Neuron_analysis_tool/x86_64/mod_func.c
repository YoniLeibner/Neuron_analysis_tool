#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _beforestep_py_reg(void);
extern void _CaDynamics_E2_reg(void);
extern void _Ca_HVA_reg(void);
extern void _Ca_LVAst_reg(void);
extern void _epsp_reg(void);
extern void _g_total_reg(void);
extern void _Ih_human_linear_reg(void);
extern void _Ih_human_shifts_mul_add_reg(void);
extern void _Ih_reg(void);
extern void _Im_reg(void);
extern void _K_Pst_reg(void);
extern void _K_Tst_reg(void);
extern void _kv_reg(void);
extern void _na_reg(void);
extern void _Nap_Et2_reg(void);
extern void _NaTa_t_reg(void);
extern void _NaTs2_t_reg(void);
extern void _NMDA_reg(void);
extern void _ProbAMPA_reg(void);
extern void _ProbAMPANMDA2_ratio_reg(void);
extern void _ProbAMPANMDA_EMS_reg(void);
extern void _ProbAMPANMDA_reg(void);
extern void _ProbGABAAB_EMS_reg(void);
extern void _ProbGABAA_EMS_reg(void);
extern void _ProbGABAA_reg(void);
extern void _ProbNMDA_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mods/beforestep_py.mod\"");
    fprintf(stderr," \"mods/CaDynamics_E2.mod\"");
    fprintf(stderr," \"mods/Ca_HVA.mod\"");
    fprintf(stderr," \"mods/Ca_LVAst.mod\"");
    fprintf(stderr," \"mods/epsp.mod\"");
    fprintf(stderr," \"mods/g_total.mod\"");
    fprintf(stderr," \"mods/Ih_human_linear.mod\"");
    fprintf(stderr," \"mods/Ih_human_shifts_mul_add.mod\"");
    fprintf(stderr," \"mods/Ih.mod\"");
    fprintf(stderr," \"mods/Im.mod\"");
    fprintf(stderr," \"mods/K_Pst.mod\"");
    fprintf(stderr," \"mods/K_Tst.mod\"");
    fprintf(stderr," \"mods/kv.mod\"");
    fprintf(stderr," \"mods/na.mod\"");
    fprintf(stderr," \"mods/Nap_Et2.mod\"");
    fprintf(stderr," \"mods/NaTa_t.mod\"");
    fprintf(stderr," \"mods/NaTs2_t.mod\"");
    fprintf(stderr," \"mods/NMDA.mod\"");
    fprintf(stderr," \"mods/ProbAMPA.mod\"");
    fprintf(stderr," \"mods/ProbAMPANMDA2_ratio.mod\"");
    fprintf(stderr," \"mods/ProbAMPANMDA_EMS.mod\"");
    fprintf(stderr," \"mods/ProbAMPANMDA.mod\"");
    fprintf(stderr," \"mods/ProbGABAAB_EMS.mod\"");
    fprintf(stderr," \"mods/ProbGABAA_EMS.mod\"");
    fprintf(stderr," \"mods/ProbGABAA.mod\"");
    fprintf(stderr," \"mods/ProbNMDA.mod\"");
    fprintf(stderr," \"mods/SK_E2.mod\"");
    fprintf(stderr," \"mods/SKv3_1.mod\"");
    fprintf(stderr, "\n");
  }
  _beforestep_py_reg();
  _CaDynamics_E2_reg();
  _Ca_HVA_reg();
  _Ca_LVAst_reg();
  _epsp_reg();
  _g_total_reg();
  _Ih_human_linear_reg();
  _Ih_human_shifts_mul_add_reg();
  _Ih_reg();
  _Im_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _kv_reg();
  _na_reg();
  _Nap_Et2_reg();
  _NaTa_t_reg();
  _NaTs2_t_reg();
  _NMDA_reg();
  _ProbAMPA_reg();
  _ProbAMPANMDA2_ratio_reg();
  _ProbAMPANMDA_EMS_reg();
  _ProbAMPANMDA_reg();
  _ProbGABAAB_EMS_reg();
  _ProbGABAA_EMS_reg();
  _ProbGABAA_reg();
  _ProbNMDA_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
}
