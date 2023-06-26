:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
:Comment: add comment here

NEURON	{
	SUFFIX g_total
	RANGE g_total, g_syn, zero_val
	POINTER g_ref0, g_ref1, g_ref2, g_ref3, g_ref4, g_ref5, g_ref6, g_ref7, g_ref8, g_ref9, g_ref10, g_ref11, g_ref12, g_ref13, g_ref14, g_ref15, g_ref16, g_ref17, g_ref18, g_ref19
	POINTER g_syn0, g_syn1, g_syn2, g_syn3, g_syn4, g_syn5, g_syn6, g_syn7, g_syn8, g_syn9, g_syn10, g_syn11, g_syn12, g_syn13, g_syn14, g_syn15, g_syn16, g_syn17, g_syn18, g_syn19
	POINTER g_syn20, g_syn21, g_syn22, g_syn23, g_syn24, g_syn25, g_syn26, g_syn27, g_syn28, g_syn29, g_syn30, g_syn31, g_syn32, g_syn33, g_syn34, g_syn35, g_syn36, g_syn37, g_syn38, g_syn39
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
     zero_val = 0.0
}

ASSIGNED	{
	g_total	(S/cm2)
	g_syn	(S)
	g_ref0
	g_ref1
	g_ref2
	g_ref3
	g_ref4
	g_ref5
	g_ref6
	g_ref7
	g_ref8
	g_ref9
	g_ref10
	g_ref11
	g_ref12
	g_ref13
	g_ref14
	g_ref15
	g_ref16
	g_ref17
	g_ref18
	g_ref19
	g_syn0
	g_syn1
	g_syn2
	g_syn3
	g_syn4
	g_syn5
	g_syn6
	g_syn7
	g_syn8
	g_syn9
	g_syn10
	g_syn11
	g_syn12
	g_syn13
	g_syn14
	g_syn15
	g_syn16
	g_syn17
	g_syn18
	g_syn19
	g_syn20
	g_syn21
	g_syn22
	g_syn23
	g_syn24
	g_syn25
	g_syn26
	g_syn27
	g_syn28
	g_syn29
	g_syn30
	g_syn31
	g_syn32
	g_syn33
	g_syn34
	g_syn35
	g_syn36
	g_syn37
	g_syn38
	g_syn39
}

STATE	{ 
}

BREAKPOINT	{
	g_total = g_ref0+g_ref1+g_ref2+g_ref3+g_ref4+g_ref5+g_ref6+g_ref7+g_ref8+g_ref9+g_ref10+g_ref11+g_ref12+g_ref13+g_ref14+g_ref15+g_ref16+g_ref17+g_ref18+g_ref19
    g_syn = g_syn0+g_syn1+g_syn2+g_syn3+g_syn4+g_syn5+g_syn6+g_syn7+g_syn8+g_syn9+g_syn10+g_syn11+g_syn12+g_syn13+g_syn14+g_syn15+g_syn16+g_syn17+g_syn18+g_syn19+g_syn20+g_syn21+g_syn22+g_syn23+g_syn24+g_syn25+g_syn26+g_syn27+g_syn28+g_syn29+g_syn30+g_syn31+g_syn32+g_syn33+g_syn34+g_syn35+g_syn36+g_syn37+g_syn38+g_syn39
}


INITIAL{
}

PROCEDURE rates(){
}
