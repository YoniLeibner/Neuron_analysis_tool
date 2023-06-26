/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__g_total
#define _nrn_initial _nrn_initial__g_total
#define nrn_cur _nrn_cur__g_total
#define _nrn_current _nrn_current__g_total
#define nrn_jacob _nrn_jacob__g_total
#define nrn_state _nrn_state__g_total
#define _net_receive _net_receive__g_total 
#define rates rates__g_total 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define zero_val _p[0]
#define g_total _p[1]
#define g_syn _p[2]
#define _g _p[3]
#define g_ref0	*_ppvar[0]._pval
#define _p_g_ref0	_ppvar[0]._pval
#define g_ref1	*_ppvar[1]._pval
#define _p_g_ref1	_ppvar[1]._pval
#define g_ref2	*_ppvar[2]._pval
#define _p_g_ref2	_ppvar[2]._pval
#define g_ref3	*_ppvar[3]._pval
#define _p_g_ref3	_ppvar[3]._pval
#define g_ref4	*_ppvar[4]._pval
#define _p_g_ref4	_ppvar[4]._pval
#define g_ref5	*_ppvar[5]._pval
#define _p_g_ref5	_ppvar[5]._pval
#define g_ref6	*_ppvar[6]._pval
#define _p_g_ref6	_ppvar[6]._pval
#define g_ref7	*_ppvar[7]._pval
#define _p_g_ref7	_ppvar[7]._pval
#define g_ref8	*_ppvar[8]._pval
#define _p_g_ref8	_ppvar[8]._pval
#define g_ref9	*_ppvar[9]._pval
#define _p_g_ref9	_ppvar[9]._pval
#define g_ref10	*_ppvar[10]._pval
#define _p_g_ref10	_ppvar[10]._pval
#define g_ref11	*_ppvar[11]._pval
#define _p_g_ref11	_ppvar[11]._pval
#define g_ref12	*_ppvar[12]._pval
#define _p_g_ref12	_ppvar[12]._pval
#define g_ref13	*_ppvar[13]._pval
#define _p_g_ref13	_ppvar[13]._pval
#define g_ref14	*_ppvar[14]._pval
#define _p_g_ref14	_ppvar[14]._pval
#define g_ref15	*_ppvar[15]._pval
#define _p_g_ref15	_ppvar[15]._pval
#define g_ref16	*_ppvar[16]._pval
#define _p_g_ref16	_ppvar[16]._pval
#define g_ref17	*_ppvar[17]._pval
#define _p_g_ref17	_ppvar[17]._pval
#define g_ref18	*_ppvar[18]._pval
#define _p_g_ref18	_ppvar[18]._pval
#define g_ref19	*_ppvar[19]._pval
#define _p_g_ref19	_ppvar[19]._pval
#define g_syn0	*_ppvar[20]._pval
#define _p_g_syn0	_ppvar[20]._pval
#define g_syn1	*_ppvar[21]._pval
#define _p_g_syn1	_ppvar[21]._pval
#define g_syn2	*_ppvar[22]._pval
#define _p_g_syn2	_ppvar[22]._pval
#define g_syn3	*_ppvar[23]._pval
#define _p_g_syn3	_ppvar[23]._pval
#define g_syn4	*_ppvar[24]._pval
#define _p_g_syn4	_ppvar[24]._pval
#define g_syn5	*_ppvar[25]._pval
#define _p_g_syn5	_ppvar[25]._pval
#define g_syn6	*_ppvar[26]._pval
#define _p_g_syn6	_ppvar[26]._pval
#define g_syn7	*_ppvar[27]._pval
#define _p_g_syn7	_ppvar[27]._pval
#define g_syn8	*_ppvar[28]._pval
#define _p_g_syn8	_ppvar[28]._pval
#define g_syn9	*_ppvar[29]._pval
#define _p_g_syn9	_ppvar[29]._pval
#define g_syn10	*_ppvar[30]._pval
#define _p_g_syn10	_ppvar[30]._pval
#define g_syn11	*_ppvar[31]._pval
#define _p_g_syn11	_ppvar[31]._pval
#define g_syn12	*_ppvar[32]._pval
#define _p_g_syn12	_ppvar[32]._pval
#define g_syn13	*_ppvar[33]._pval
#define _p_g_syn13	_ppvar[33]._pval
#define g_syn14	*_ppvar[34]._pval
#define _p_g_syn14	_ppvar[34]._pval
#define g_syn15	*_ppvar[35]._pval
#define _p_g_syn15	_ppvar[35]._pval
#define g_syn16	*_ppvar[36]._pval
#define _p_g_syn16	_ppvar[36]._pval
#define g_syn17	*_ppvar[37]._pval
#define _p_g_syn17	_ppvar[37]._pval
#define g_syn18	*_ppvar[38]._pval
#define _p_g_syn18	_ppvar[38]._pval
#define g_syn19	*_ppvar[39]._pval
#define _p_g_syn19	_ppvar[39]._pval
#define g_syn20	*_ppvar[40]._pval
#define _p_g_syn20	_ppvar[40]._pval
#define g_syn21	*_ppvar[41]._pval
#define _p_g_syn21	_ppvar[41]._pval
#define g_syn22	*_ppvar[42]._pval
#define _p_g_syn22	_ppvar[42]._pval
#define g_syn23	*_ppvar[43]._pval
#define _p_g_syn23	_ppvar[43]._pval
#define g_syn24	*_ppvar[44]._pval
#define _p_g_syn24	_ppvar[44]._pval
#define g_syn25	*_ppvar[45]._pval
#define _p_g_syn25	_ppvar[45]._pval
#define g_syn26	*_ppvar[46]._pval
#define _p_g_syn26	_ppvar[46]._pval
#define g_syn27	*_ppvar[47]._pval
#define _p_g_syn27	_ppvar[47]._pval
#define g_syn28	*_ppvar[48]._pval
#define _p_g_syn28	_ppvar[48]._pval
#define g_syn29	*_ppvar[49]._pval
#define _p_g_syn29	_ppvar[49]._pval
#define g_syn30	*_ppvar[50]._pval
#define _p_g_syn30	_ppvar[50]._pval
#define g_syn31	*_ppvar[51]._pval
#define _p_g_syn31	_ppvar[51]._pval
#define g_syn32	*_ppvar[52]._pval
#define _p_g_syn32	_ppvar[52]._pval
#define g_syn33	*_ppvar[53]._pval
#define _p_g_syn33	_ppvar[53]._pval
#define g_syn34	*_ppvar[54]._pval
#define _p_g_syn34	_ppvar[54]._pval
#define g_syn35	*_ppvar[55]._pval
#define _p_g_syn35	_ppvar[55]._pval
#define g_syn36	*_ppvar[56]._pval
#define _p_g_syn36	_ppvar[56]._pval
#define g_syn37	*_ppvar[57]._pval
#define _p_g_syn37	_ppvar[57]._pval
#define g_syn38	*_ppvar[58]._pval
#define _p_g_syn38	_ppvar[58]._pval
#define g_syn39	*_ppvar[59]._pval
#define _p_g_syn39	_ppvar[59]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  0;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_g_total", _hoc_setdata,
 "rates_g_total", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "g_total_g_total", "S/cm2",
 "g_syn_g_total", "S",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"g_total",
 "zero_val_g_total",
 0,
 "g_total_g_total",
 "g_syn_g_total",
 0,
 0,
 "g_ref0_g_total",
 "g_ref1_g_total",
 "g_ref2_g_total",
 "g_ref3_g_total",
 "g_ref4_g_total",
 "g_ref5_g_total",
 "g_ref6_g_total",
 "g_ref7_g_total",
 "g_ref8_g_total",
 "g_ref9_g_total",
 "g_ref10_g_total",
 "g_ref11_g_total",
 "g_ref12_g_total",
 "g_ref13_g_total",
 "g_ref14_g_total",
 "g_ref15_g_total",
 "g_ref16_g_total",
 "g_ref17_g_total",
 "g_ref18_g_total",
 "g_ref19_g_total",
 "g_syn0_g_total",
 "g_syn1_g_total",
 "g_syn2_g_total",
 "g_syn3_g_total",
 "g_syn4_g_total",
 "g_syn5_g_total",
 "g_syn6_g_total",
 "g_syn7_g_total",
 "g_syn8_g_total",
 "g_syn9_g_total",
 "g_syn10_g_total",
 "g_syn11_g_total",
 "g_syn12_g_total",
 "g_syn13_g_total",
 "g_syn14_g_total",
 "g_syn15_g_total",
 "g_syn16_g_total",
 "g_syn17_g_total",
 "g_syn18_g_total",
 "g_syn19_g_total",
 "g_syn20_g_total",
 "g_syn21_g_total",
 "g_syn22_g_total",
 "g_syn23_g_total",
 "g_syn24_g_total",
 "g_syn25_g_total",
 "g_syn26_g_total",
 "g_syn27_g_total",
 "g_syn28_g_total",
 "g_syn29_g_total",
 "g_syn30_g_total",
 "g_syn31_g_total",
 "g_syn32_g_total",
 "g_syn33_g_total",
 "g_syn34_g_total",
 "g_syn35_g_total",
 "g_syn36_g_total",
 "g_syn37_g_total",
 "g_syn38_g_total",
 "g_syn39_g_total",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 4, _prop);
 	/*initialize range parameters*/
 	zero_val = 0;
 	_prop->param = _p;
 	_prop->param_size = 4;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 60, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _g_total_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 4, 60);
  hoc_register_dparam_semantics(_mechtype, 0, "pointer");
  hoc_register_dparam_semantics(_mechtype, 1, "pointer");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
  hoc_register_dparam_semantics(_mechtype, 3, "pointer");
  hoc_register_dparam_semantics(_mechtype, 4, "pointer");
  hoc_register_dparam_semantics(_mechtype, 5, "pointer");
  hoc_register_dparam_semantics(_mechtype, 6, "pointer");
  hoc_register_dparam_semantics(_mechtype, 7, "pointer");
  hoc_register_dparam_semantics(_mechtype, 8, "pointer");
  hoc_register_dparam_semantics(_mechtype, 9, "pointer");
  hoc_register_dparam_semantics(_mechtype, 10, "pointer");
  hoc_register_dparam_semantics(_mechtype, 11, "pointer");
  hoc_register_dparam_semantics(_mechtype, 12, "pointer");
  hoc_register_dparam_semantics(_mechtype, 13, "pointer");
  hoc_register_dparam_semantics(_mechtype, 14, "pointer");
  hoc_register_dparam_semantics(_mechtype, 15, "pointer");
  hoc_register_dparam_semantics(_mechtype, 16, "pointer");
  hoc_register_dparam_semantics(_mechtype, 17, "pointer");
  hoc_register_dparam_semantics(_mechtype, 18, "pointer");
  hoc_register_dparam_semantics(_mechtype, 19, "pointer");
  hoc_register_dparam_semantics(_mechtype, 20, "pointer");
  hoc_register_dparam_semantics(_mechtype, 21, "pointer");
  hoc_register_dparam_semantics(_mechtype, 22, "pointer");
  hoc_register_dparam_semantics(_mechtype, 23, "pointer");
  hoc_register_dparam_semantics(_mechtype, 24, "pointer");
  hoc_register_dparam_semantics(_mechtype, 25, "pointer");
  hoc_register_dparam_semantics(_mechtype, 26, "pointer");
  hoc_register_dparam_semantics(_mechtype, 27, "pointer");
  hoc_register_dparam_semantics(_mechtype, 28, "pointer");
  hoc_register_dparam_semantics(_mechtype, 29, "pointer");
  hoc_register_dparam_semantics(_mechtype, 30, "pointer");
  hoc_register_dparam_semantics(_mechtype, 31, "pointer");
  hoc_register_dparam_semantics(_mechtype, 32, "pointer");
  hoc_register_dparam_semantics(_mechtype, 33, "pointer");
  hoc_register_dparam_semantics(_mechtype, 34, "pointer");
  hoc_register_dparam_semantics(_mechtype, 35, "pointer");
  hoc_register_dparam_semantics(_mechtype, 36, "pointer");
  hoc_register_dparam_semantics(_mechtype, 37, "pointer");
  hoc_register_dparam_semantics(_mechtype, 38, "pointer");
  hoc_register_dparam_semantics(_mechtype, 39, "pointer");
  hoc_register_dparam_semantics(_mechtype, 40, "pointer");
  hoc_register_dparam_semantics(_mechtype, 41, "pointer");
  hoc_register_dparam_semantics(_mechtype, 42, "pointer");
  hoc_register_dparam_semantics(_mechtype, 43, "pointer");
  hoc_register_dparam_semantics(_mechtype, 44, "pointer");
  hoc_register_dparam_semantics(_mechtype, 45, "pointer");
  hoc_register_dparam_semantics(_mechtype, 46, "pointer");
  hoc_register_dparam_semantics(_mechtype, 47, "pointer");
  hoc_register_dparam_semantics(_mechtype, 48, "pointer");
  hoc_register_dparam_semantics(_mechtype, 49, "pointer");
  hoc_register_dparam_semantics(_mechtype, 50, "pointer");
  hoc_register_dparam_semantics(_mechtype, 51, "pointer");
  hoc_register_dparam_semantics(_mechtype, 52, "pointer");
  hoc_register_dparam_semantics(_mechtype, 53, "pointer");
  hoc_register_dparam_semantics(_mechtype, 54, "pointer");
  hoc_register_dparam_semantics(_mechtype, 55, "pointer");
  hoc_register_dparam_semantics(_mechtype, 56, "pointer");
  hoc_register_dparam_semantics(_mechtype, 57, "pointer");
  hoc_register_dparam_semantics(_mechtype, 58, "pointer");
  hoc_register_dparam_semantics(_mechtype, 59, "pointer");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 g_total /ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Neuron_analysis_tool/Neuron_analysis_tool/mods/g_total.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates();
 
static int  rates (  ) {
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  );
 hoc_retpushx(_r);
}

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   }

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {
   g_total = g_ref0 + g_ref1 + g_ref2 + g_ref3 + g_ref4 + g_ref5 + g_ref6 + g_ref7 + g_ref8 + g_ref9 + g_ref10 + g_ref11 + g_ref12 + g_ref13 + g_ref14 + g_ref15 + g_ref16 + g_ref17 + g_ref18 + g_ref19 ;
   g_syn = g_syn0 + g_syn1 + g_syn2 + g_syn3 + g_syn4 + g_syn5 + g_syn6 + g_syn7 + g_syn8 + g_syn9 + g_syn10 + g_syn11 + g_syn12 + g_syn13 + g_syn14 + g_syn15 + g_syn16 + g_syn17 + g_syn18 + g_syn19 + g_syn20 + g_syn21 + g_syn22 + g_syn23 + g_syn24 + g_syn25 + g_syn26 + g_syn27 + g_syn28 + g_syn29 + g_syn30 + g_syn31 + g_syn32 + g_syn33 + g_syn34 + g_syn35 + g_syn36 + g_syn37 + g_syn38 + g_syn39 ;
   }
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Neuron_analysis_tool/Neuron_analysis_tool/mods/g_total.mod";
static const char* nmodl_file_text = 
  ":Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones\n"
  ":Comment: add comment here\n"
  "\n"
  "NEURON	{\n"
  "	SUFFIX g_total\n"
  "	RANGE g_total, g_syn, zero_val\n"
  "	POINTER g_ref0, g_ref1, g_ref2, g_ref3, g_ref4, g_ref5, g_ref6, g_ref7, g_ref8, g_ref9, g_ref10, g_ref11, g_ref12, g_ref13, g_ref14, g_ref15, g_ref16, g_ref17, g_ref18, g_ref19\n"
  "	POINTER g_syn0, g_syn1, g_syn2, g_syn3, g_syn4, g_syn5, g_syn6, g_syn7, g_syn8, g_syn9, g_syn10, g_syn11, g_syn12, g_syn13, g_syn14, g_syn15, g_syn16, g_syn17, g_syn18, g_syn19\n"
  "	POINTER g_syn20, g_syn21, g_syn22, g_syn23, g_syn24, g_syn25, g_syn26, g_syn27, g_syn28, g_syn29, g_syn30, g_syn31, g_syn32, g_syn33, g_syn34, g_syn35, g_syn36, g_syn37, g_syn38, g_syn39\n"
  "}\n"
  "\n"
  "UNITS	{\n"
  "	(S) = (siemens)\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "}\n"
  "\n"
  "PARAMETER	{\n"
  "     zero_val = 0.0\n"
  "}\n"
  "\n"
  "ASSIGNED	{\n"
  "	g_total	(S/cm2)\n"
  "	g_syn	(S)\n"
  "	g_ref0\n"
  "	g_ref1\n"
  "	g_ref2\n"
  "	g_ref3\n"
  "	g_ref4\n"
  "	g_ref5\n"
  "	g_ref6\n"
  "	g_ref7\n"
  "	g_ref8\n"
  "	g_ref9\n"
  "	g_ref10\n"
  "	g_ref11\n"
  "	g_ref12\n"
  "	g_ref13\n"
  "	g_ref14\n"
  "	g_ref15\n"
  "	g_ref16\n"
  "	g_ref17\n"
  "	g_ref18\n"
  "	g_ref19\n"
  "	g_syn0\n"
  "	g_syn1\n"
  "	g_syn2\n"
  "	g_syn3\n"
  "	g_syn4\n"
  "	g_syn5\n"
  "	g_syn6\n"
  "	g_syn7\n"
  "	g_syn8\n"
  "	g_syn9\n"
  "	g_syn10\n"
  "	g_syn11\n"
  "	g_syn12\n"
  "	g_syn13\n"
  "	g_syn14\n"
  "	g_syn15\n"
  "	g_syn16\n"
  "	g_syn17\n"
  "	g_syn18\n"
  "	g_syn19\n"
  "	g_syn20\n"
  "	g_syn21\n"
  "	g_syn22\n"
  "	g_syn23\n"
  "	g_syn24\n"
  "	g_syn25\n"
  "	g_syn26\n"
  "	g_syn27\n"
  "	g_syn28\n"
  "	g_syn29\n"
  "	g_syn30\n"
  "	g_syn31\n"
  "	g_syn32\n"
  "	g_syn33\n"
  "	g_syn34\n"
  "	g_syn35\n"
  "	g_syn36\n"
  "	g_syn37\n"
  "	g_syn38\n"
  "	g_syn39\n"
  "}\n"
  "\n"
  "STATE	{ \n"
  "}\n"
  "\n"
  "BREAKPOINT	{\n"
  "	g_total = g_ref0+g_ref1+g_ref2+g_ref3+g_ref4+g_ref5+g_ref6+g_ref7+g_ref8+g_ref9+g_ref10+g_ref11+g_ref12+g_ref13+g_ref14+g_ref15+g_ref16+g_ref17+g_ref18+g_ref19\n"
  "    g_syn = g_syn0+g_syn1+g_syn2+g_syn3+g_syn4+g_syn5+g_syn6+g_syn7+g_syn8+g_syn9+g_syn10+g_syn11+g_syn12+g_syn13+g_syn14+g_syn15+g_syn16+g_syn17+g_syn18+g_syn19+g_syn20+g_syn21+g_syn22+g_syn23+g_syn24+g_syn25+g_syn26+g_syn27+g_syn28+g_syn29+g_syn30+g_syn31+g_syn32+g_syn33+g_syn34+g_syn35+g_syn36+g_syn37+g_syn38+g_syn39\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL{\n"
  "}\n"
  "\n"
  "PROCEDURE rates(){\n"
  "}\n"
  ;
#endif
