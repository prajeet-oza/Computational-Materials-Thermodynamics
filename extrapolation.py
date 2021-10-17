
def G_AB(xA, xB):
	value = xA*xB*(-5500+6000*(xB-xA))
	return value

def G_BC(xB, xC):
	value = xB*xC*(-8200-3300*(xC-xB))
	return value

def G_AC(xA, xC):
	value = xA*xC*(5330-8350*xC)
	return value

def Muggianu(mfA, mfB, mfC):
	GE_Mug = 0

	xA_AB = (1+mfA-mfB)/2
	xB_AB = (1+mfB-mfA)/2
	xA_AC = (1+mfA-mfC)/2
	xC_AC = (1+mfC-mfA)/2
	xB_BC = (1+mfB-mfC)/2
	xC_BC = (1+mfC-mfB)/2

	f_AB = mfA*mfB/(xA_AB*xB_AB)
	f_BC = mfB*mfC/(xB_BC*xC_BC)
	f_AC = mfA*mfC/(xA_AC*xC_AC)

	GE_Mug += f_AB*G_AB(xA_AB, xB_AB) + f_AC*G_AC(xA_AC, xC_AC) + \
				f_BC*G_BC(xB_BC, xC_BC)
	return GE_Mug

def Kohler(mfA, mfB, mfC):
	GE_Koh = 0

	xA_AB = (1+mfA-mfB)/2 + (1-mfA-mfB)*(mfA-mfB)/2/(mfA+mfB)
	xB_AB = (1+mfB-mfA)/2 + (1-mfB-mfA)*(mfB-mfA)/2/(mfA+mfB)
	xA_AC = (1+mfA-mfC)/2 + (1-mfA-mfC)*(mfA-mfC)/2/(mfA+mfC)
	xC_AC = (1+mfC-mfA)/2 + (1-mfC-mfA)*(mfC-mfA)/2/(mfA+mfC)
	xB_BC = (1+mfB-mfC)/2 + (1-mfB-mfC)*(mfB-mfC)/2/(mfB+mfC)
	xC_BC = (1+mfC-mfB)/2 + (1-mfC-mfB)*(mfC-mfB)/2/(mfB+mfC)

	f_AB = mfA*mfB/(xA_AB*xB_AB)
	f_BC = mfB*mfC/(xB_BC*xC_BC)
	f_AC = mfA*mfC/(xA_AC*xC_AC)

	GE_Koh += f_AB*G_AB(xA_AB, xB_AB) + f_AC*G_AC(xA_AC, xC_AC) + \
				f_BC*G_BC(xB_BC, xC_BC)
	return GE_Koh

def Colinet(mfA, mfB, mfC):
	GE_Col = 0

	x1A_AB = (1+mfA-mfB)/2 + (1-mfA-mfB)/2
	x2A_AB = (1+mfA-mfB)/2 - (1-mfA-mfB)/2
	x1B_AB = (1+mfB-mfA)/2 + (1-mfB-mfA)/2
	x2B_AB = (1+mfB-mfA)/2 - (1-mfB-mfA)/2
	x1A_AC = (1+mfA-mfC)/2 + (1-mfA-mfC)/2
	x2A_AC = (1+mfA-mfC)/2 - (1-mfA-mfC)/2
	x1C_AC = (1+mfC-mfA)/2 + (1-mfC-mfA)/2
	x2C_AC = (1+mfC-mfA)/2 - (1-mfC-mfA)/2
	x1B_BC = (1+mfB-mfC)/2 + (1-mfB-mfC)/2
	x2B_BC = (1+mfB-mfC)/2 - (1-mfB-mfC)/2
	x1C_BC = (1+mfC-mfB)/2 + (1-mfC-mfB)/2
	x2C_BC = (1+mfC-mfB)/2 - (1-mfC-mfB)/2

	f1_AB = mfA*mfB/(x1A_AB*x1B_AB)
	f1_BC = mfB*mfC/(x1B_BC*x1C_BC)
	f1_AC = mfA*mfC/(x1A_AC*x1C_AC)

	f2_AB = mfA*mfB/(x2A_AB*x2B_AB)
	f2_BC = mfB*mfC/(x2B_BC*x2C_BC)
	f2_AC = mfA*mfC/(x2A_AC*x2C_AC)

	GE_Col += 0.5*(f1_AB*G_AB(x1A_AB, x1B_AB) + f1_AC*G_AC(x1A_AC, x1C_AC) + \
				f1_BC*G_BC(x1B_BC, x1C_BC))
	GE_Col += 0.5*(f2_AB*G_AB(x2A_AB, x2B_AB) + f2_AC*G_AC(x2A_AC, x2C_AC) + \
				f2_BC*G_BC(x2B_BC, x2C_BC))
	return GE_Col

mfB = 0.38
mfC = 0.42
mfA = 1-mfB-mfC

print('\n GE_ABC, Muggianu scheme: ', Muggianu(mfA, mfB, mfC))
print('\n GE_ABC, Kohler scheme: ', Kohler(mfA, mfB, mfC))
print('\n GE_ABC, Colinet scheme: ', Colinet(mfA, mfB, mfC))
print('\n with units: J/mol\n')