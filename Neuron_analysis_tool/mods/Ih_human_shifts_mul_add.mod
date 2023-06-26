:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih_human_shifts_mul_add
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn, m_alpha_shift_v, m_tau_add, m_tau_mul, dist, dist_thresh, gIhbar_exp, max_gbar_val
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gIhbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)
	m_alpha_shift_v = -20 (mV)
	m_tau_add = 20
	m_tau_mul = 5
	dist_thresh = 0.0031
	dist = 0
	max_gbar_val = 0.001 (S/cm2)
	gIhbar_exp = 0.00001 (S/cm2)
}

ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	gIh	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIh = gIhbar_exp*m
	ihcn = gIh*(v-ehcn)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
	max_gbar_val = 0.015
	gIhbar_exp = (-0.8696 + 2.087 * exp(dist*dist_thresh)) * gIhbar
	if(gIhbar_exp > max_gbar_val){
            gIhbar_exp = max_gbar_val
    }
}

PROCEDURE rates(){
	UNITSOFF
        if(-145 + m_alpha_shift_v > v){
            v = -145
        }
		mAlpha =  0.001*6.43*(v+154.9 + m_alpha_shift_v)/(exp((v+154.9 + m_alpha_shift_v)/11.9)-1)
		mBeta  =  0.001*193*exp((v+m_alpha_shift_v)/33.1)
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta) * m_tau_mul + m_tau_add
	UNITSON
}
