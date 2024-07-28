/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
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
 
#define nrn_init _nrn_init__Ih_human_shifts_mul_add
#define _nrn_initial _nrn_initial__Ih_human_shifts_mul_add
#define nrn_cur _nrn_cur__Ih_human_shifts_mul_add
#define _nrn_current _nrn_current__Ih_human_shifts_mul_add
#define nrn_jacob _nrn_jacob__Ih_human_shifts_mul_add
#define nrn_state _nrn_state__Ih_human_shifts_mul_add
#define _net_receive _net_receive__Ih_human_shifts_mul_add 
#define rates rates__Ih_human_shifts_mul_add 
#define states states__Ih_human_shifts_mul_add 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gIhbar _p[0]
#define gIhbar_columnindex 0
#define m_alpha_shift_v _p[1]
#define m_alpha_shift_v_columnindex 1
#define m_tau_add _p[2]
#define m_tau_add_columnindex 2
#define m_tau_mul _p[3]
#define m_tau_mul_columnindex 3
#define dist_thresh _p[4]
#define dist_thresh_columnindex 4
#define dist _p[5]
#define dist_columnindex 5
#define max_gbar_val _p[6]
#define max_gbar_val_columnindex 6
#define gIhbar_exp _p[7]
#define gIhbar_exp_columnindex 7
#define ihcn _p[8]
#define ihcn_columnindex 8
#define gIh _p[9]
#define gIh_columnindex 9
#define m _p[10]
#define m_columnindex 10
#define mInf _p[11]
#define mInf_columnindex 11
#define mTau _p[12]
#define mTau_columnindex 12
#define mAlpha _p[13]
#define mAlpha_columnindex 13
#define mBeta _p[14]
#define mBeta_columnindex 14
#define Dm _p[15]
#define Dm_columnindex 15
#define v _p[16]
#define v_columnindex 16
#define _g _p[17]
#define _g_columnindex 17
 
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
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Ih_human_shifts_mul_add", _hoc_setdata,
 "rates_Ih_human_shifts_mul_add", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define ehcn ehcn_Ih_human_shifts_mul_add
 double ehcn = -45;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ehcn_Ih_human_shifts_mul_add", "mV",
 "gIhbar_Ih_human_shifts_mul_add", "S/cm2",
 "m_alpha_shift_v_Ih_human_shifts_mul_add", "mV",
 "max_gbar_val_Ih_human_shifts_mul_add", "S/cm2",
 "gIhbar_exp_Ih_human_shifts_mul_add", "S/cm2",
 "ihcn_Ih_human_shifts_mul_add", "mA/cm2",
 "gIh_Ih_human_shifts_mul_add", "S/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ehcn_Ih_human_shifts_mul_add", &ehcn_Ih_human_shifts_mul_add,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Ih_human_shifts_mul_add",
 "gIhbar_Ih_human_shifts_mul_add",
 "m_alpha_shift_v_Ih_human_shifts_mul_add",
 "m_tau_add_Ih_human_shifts_mul_add",
 "m_tau_mul_Ih_human_shifts_mul_add",
 "dist_thresh_Ih_human_shifts_mul_add",
 "dist_Ih_human_shifts_mul_add",
 "max_gbar_val_Ih_human_shifts_mul_add",
 "gIhbar_exp_Ih_human_shifts_mul_add",
 0,
 "ihcn_Ih_human_shifts_mul_add",
 "gIh_Ih_human_shifts_mul_add",
 0,
 "m_Ih_human_shifts_mul_add",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	gIhbar = 1e-05;
 	m_alpha_shift_v = -20;
 	m_tau_add = 20;
 	m_tau_mul = 5;
 	dist_thresh = 0.0031;
 	dist = 0;
 	max_gbar_val = 0.001;
 	gIhbar_exp = 1e-05;
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Ih_human_shifts_mul_add_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 18, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Ih_human_shifts_mul_add /ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Neuron_analysis_tool/Neuron_analysis_tool/mods/Ih_human_shifts_mul_add.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dm = ( mInf - m ) / mTau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mTau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mTau)))*(- ( ( ( mInf ) ) / mTau ) / ( ( ( ( - 1.0 ) ) ) / mTau ) - m) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
    if ( - 145.0 + m_alpha_shift_v > v ) {
     v = - 145.0 ;
     }
   mAlpha = 0.001 * 6.43 * ( v + 154.9 + m_alpha_shift_v ) / ( exp ( ( v + 154.9 + m_alpha_shift_v ) / 11.9 ) - 1.0 ) ;
   mBeta = 0.001 * 193.0 * exp ( ( v + m_alpha_shift_v ) / 33.1 ) ;
   mInf = mAlpha / ( mAlpha + mBeta ) ;
   mTau = 1.0 / ( mAlpha + mBeta ) * m_tau_mul + m_tau_add ;
     return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  m = m0;
 {
   rates ( _threadargs_ ) ;
   m = mInf ;
   max_gbar_val = 0.015 ;
   gIhbar_exp = ( - 0.8696 + 2.087 * exp ( dist * dist_thresh ) ) * gIhbar ;
   if ( gIhbar_exp > max_gbar_val ) {
     gIhbar_exp = max_gbar_val ;
     }
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gIh = gIhbar_exp * m ;
   ihcn = gIh * ( v - ehcn ) ;
   }
 _current += ihcn;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/ems/elsc-labs/segev-i/yoni.leibner/PycharmProjects/Neuron_analysis_tool/Neuron_analysis_tool/mods/Ih_human_shifts_mul_add.mod";
static const char* nmodl_file_text = 
  ":Comment :\n"
  ":Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006\n"
  "\n"
  "NEURON	{\n"
  "	SUFFIX Ih_human_shifts_mul_add\n"
  "	NONSPECIFIC_CURRENT ihcn\n"
  "	RANGE gIhbar, gIh, ihcn, m_alpha_shift_v, m_tau_add, m_tau_mul, dist, dist_thresh, gIhbar_exp, max_gbar_val\n"
  "}\n"
  "\n"
  "UNITS	{\n"
  "	(S) = (siemens)\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "}\n"
  "\n"
  "PARAMETER	{\n"
  "	gIhbar = 0.00001 (S/cm2) \n"
  "	ehcn =  -45.0 (mV)\n"
  "	m_alpha_shift_v = -20 (mV)\n"
  "	m_tau_add = 20\n"
  "	m_tau_mul = 5\n"
  "	dist_thresh = 0.0031\n"
  "	dist = 0\n"
  "	max_gbar_val = 0.001 (S/cm2)\n"
  "	gIhbar_exp = 0.00001 (S/cm2)\n"
  "}\n"
  "\n"
  "ASSIGNED	{\n"
  "	v	(mV)\n"
  "	ihcn	(mA/cm2)\n"
  "	gIh	(S/cm2)\n"
  "	mInf\n"
  "	mTau\n"
  "	mAlpha\n"
  "	mBeta\n"
  "}\n"
  "\n"
  "STATE	{ \n"
  "	m\n"
  "}\n"
  "\n"
  "BREAKPOINT	{\n"
  "	SOLVE states METHOD cnexp\n"
  "	gIh = gIhbar_exp*m\n"
  "	ihcn = gIh*(v-ehcn)\n"
  "}\n"
  "\n"
  "DERIVATIVE states	{\n"
  "	rates()\n"
  "	m' = (mInf-m)/mTau\n"
  "}\n"
  "\n"
  "INITIAL{\n"
  "	rates()\n"
  "	m = mInf\n"
  "	max_gbar_val = 0.015\n"
  "	gIhbar_exp = (-0.8696 + 2.087 * exp(dist*dist_thresh)) * gIhbar\n"
  "	if(gIhbar_exp > max_gbar_val){\n"
  "            gIhbar_exp = max_gbar_val\n"
  "    }\n"
  "}\n"
  "\n"
  "PROCEDURE rates(){\n"
  "	UNITSOFF\n"
  "        if(-145 + m_alpha_shift_v > v){\n"
  "            v = -145\n"
  "        }\n"
  "		mAlpha =  0.001*6.43*(v+154.9 + m_alpha_shift_v)/(exp((v+154.9 + m_alpha_shift_v)/11.9)-1)\n"
  "		mBeta  =  0.001*193*exp((v+m_alpha_shift_v)/33.1)\n"
  "		mInf = mAlpha/(mAlpha + mBeta)\n"
  "		mTau = 1/(mAlpha + mBeta) * m_tau_mul + m_tau_add\n"
  "	UNITSON\n"
  "}\n"
  ;
#endif
