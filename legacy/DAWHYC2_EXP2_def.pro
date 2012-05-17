pro DAWHYC2_EXP2_def

close,/all

;*************************************************************************************************************************

dir_salida = '/home/salazar/JFSDoc/CalculosDAWHYC_JFSDoc/ResDAWHYC/DAWHYC2_Astronomia'

;PARÁMETROS FIJOS
cp = 3.e13	;erg*cm-2*K-1
S = 2.89e13	;erg*cm-2*año-1
sigma = 1789.	;erg*cm-2*año-1*K-4
gamma = 0.3	;año-1
pp = 1.	;adimensional
Topt = 295.5	;K
q = 20.	;K
Ab = 0.25	;adimensinal
Aw = 0.75	;adimensinal
Agf = 0.50	;adimensinal
grad_T = (-0.0065)	;K*m-1, Trenberth 95, p.10

fi = 0.1
mm = 0.35
nn = 0.1
Pmax = (1./mm)^(1/nn) ;=36251 [mm/año], cuando ac=1.0
EVPmax = Pmax
;EVPmax = 1511. ;[mm/año] corresponde a Ts=26 grados centígrados
;Pmax = EVPmax	;[mm/año]
;ac_crit = mm*Pmax^nn

;*********************
;experimento otro: cambiando estos parámetros
Ec = 1.	;emisividad de las nubes
Es = 1.	;emisividad de la superficie
;*********************

deltat = 0.01 ;[años] tamaño de paso temporal para dTsdt y dadt
t_integracion = 1 ;[año] cada cuanto se guardan valores de las variables
t_modelacion = 10000	;[años] período de modelación por cada L
niter = t_integracion/deltat	;# de iteraciones en los RK4

Lini = 1.992
Lfin = 4.004	;hasta 1.604 para que la tabla quede hasta 1.800 exactamente
deltaL = 0.004
niterL = (Lfin-Lini)/deltaL

;;zc_vector = [1000.,2000.,3000.,4000.,5000.,6000.,7000.,8000.]
zc_vector = [6000.]

nzc = n_elements(zc_vector)

FOR zczc = 0, nzc-1 DO BEGIN 	;zc

zc = zc_vector(zczc)
;;Ac = 1.- fi*(zc/1000.)
Ac = 0.7

resultados2 = dblarr(19,niterL+1)

for ii = 0, niterL do begin

	L = deltaL*ii + Lini
print, ii, ' de ', niterL
	;*********************************************************************************************************************
	;PROCESAMIENTO

	;valores iniciales
	aclouds = 0.01	;adimensional, 0 para reproducir modelo original, 0.01 para iniciar con area de nubes
	awhite = 0.01	;adimensional
	ablack = 0.01	;adimensional
	Ts=295.5	;temperatura en la superficie, valor inicial para rk4

	d = t_modelacion ;numero de años en el eje de las abscisas - iteraciones de t - dimension de los vectores de resultados
	time=fltarr(d+1)
	temperature=dblarr(d+1)
	white_temperature=dblarr(d+1)
	black_temperature=dblarr(d+1)
	white_area=dblarr(d+1)
	black_area=dblarr(d+1)
	clouds_area=dblarr(d+1)
	evap=dblarr(d+1)
	prec=dblarr(d+1)

	FOR t=0,t_modelacion,t_integracion DO BEGIN

		;runge-kutta4 para resolver Ts
		xx = pp - awhite - ablack
		As = xx*Agf + ablack*Ab + awhite*Aw

		for j=1,niter do begin

			;actualización temperatura
			k1=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*(Ts+grad_T*zc)^4-sigma*Es*Ts^4)
			k2=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*((Ts+k1/2)+grad_T*zc)^4-sigma*Es*(Ts+k1/2)^4)
			k3=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*((Ts+k2/2)+grad_T*zc)^4-sigma*Es*(Ts+k2/2)^4)
			k4=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*((Ts+k3)+grad_T*zc)^4-sigma*Es*(Ts+k3)^4)
			Ts = Ts+deltat*(k1/6+k2/3+k3/3+k4/6)

			;actualización areas
			Tc=Ts+zc*grad_T
			Tlb=q*(As-Ab)+Ts
			Tlw=q*(As-Aw)+Ts

			if Ts gt 277 then begin
				;if Ts lt 299 then begin
					Tanual = Ts - 273.	;(°C)
					I = 12.*(Tanual/5.)^1.5
					a = (6.7e-7)*I^3. - (7.7e-5)*I^2. + (1.8e-2)*I + 0.49
					EVP = 12.*16*(10.*(Ts - 273.)/I)^a
					E = min([1.,EVP/EVPmax])
				;endif else begin
					;E = 1.
				;endelse
			endif else begin
				E = 0.
			endelse

			;if aclouds lt ac_crit then P = (1./Pmax)*(aclouds/mm)^(1./nn) else P = 1.
			P = (1./Pmax)*(aclouds/mm)^(1./nn)

			;nubes
		   	k1c=(1-aclouds)*E-aclouds*P
			k2c=(1-(aclouds+k1c/2))*E-(aclouds+k1c/2)*P
			k3c=(1-(aclouds+k2c/2))*E-(aclouds+k2c/2)*P
			k4c=(1-(aclouds+k3c))*E-(aclouds+k3c)*P
			aclouds=aclouds+deltat*(k1c/6+k2c/3+k3c/3+k4c/6)

			;margaritas
			if (Tlw gt 278) and (Tlw lt 313) then Bw=1-0.003265*(Topt-Tlw)^2 else Bw=0
			if (Tlb gt 278) and (Tlb lt 313) then Bb=1-0.003265*(Topt-Tlb)^2 else Bb=0

			k1w=awhite*(xx*Bw-gamma)
			k2w=(awhite+k1w/2)*(xx*Bw-gamma)
			k3w=(awhite+k2w/2)*(xx*Bw-gamma)
			k4w=(awhite+k3w)*(xx*Bw-gamma)
			awhite=awhite+deltat*(k1w/6+k2w/3+k3w/3+k4w/6)

			k1b=ablack*(xx*Bb-gamma)
			k2b=(ablack+k1b/2)*(xx*Bb-gamma)
			k3b=(ablack+k2b/2)*(xx*Bb-gamma)
			k4b=(ablack+k3b)*(xx*Bb-gamma)
			ablack=ablack+deltat*(k1b/6+k2b/3+k3b/3+k4b/6)

			xx = pp - awhite - ablack
			As = xx*Agf + ablack*Ab + awhite*Aw

		endfor

		;print, awhite, ablack, aclouds, Ts, I, a, EVP, E
		time(t)=t
		temperature(t)=Ts-273
		white_temperature(t)=Tlw-273
		black_temperature(t)=Tlb-273
		white_area(t)=awhite
		black_area(t)=ablack
		clouds_area(t)=aclouds
		evap(t)=E
		prec(t)=P

	ENDFOR

	;*********************************************************************************************************************
	;TABLAS

	resultados = fltarr(9,d+1)
	resultados(0,*) = time
	resultados(1,*) = temperature
	resultados(2,*) = white_temperature
	resultados(3,*) = black_temperature
	resultados(4,*) = white_area
	resultados(5,*) = black_area
	resultados(6,*) = clouds_area
	resultados(7,*) = evap*(1-clouds_area)
	resultados(8,*) = prec*clouds_area

	resultados = resultados(*,5001:10000)				;to get rid numerical noise

	resultados2(0,ii) = L
	resultados2(1,ii) = min(resultados(1,*))	;Ts
	resultados2(2,ii) = mean(resultados(1,*))	;Ts
	resultados2(3,ii) = max(resultados(1,*))	;Ts
	resultados2(4,ii) = min(resultados(4,*))	;aw
	resultados2(5,ii) = mean(resultados(4,*))	;aw
	resultados2(6,ii) = max(resultados(4,*))	;aw
	resultados2(7,ii) = min(resultados(5,*))	;ab
	resultados2(8,ii) = mean(resultados(5,*))	;ab
	resultados2(9,ii) = max(resultados(5,*))	;ab
	resultados2(10,ii) = min(resultados(6,*))	;ac
	resultados2(11,ii) = mean(resultados(6,*))	;ac
	resultados2(12,ii) = max(resultados(6,*))	;ac
	resultados2(13,ii) = min(resultados(7,*))	;E
	resultados2(14,ii) = mean(resultados(7,*))	;E
	resultados2(15,ii) = max(resultados(7,*))	;E
	resultados2(16,ii) = min(resultados(8,*))	;P
	resultados2(17,ii) = mean(resultados(8,*))	;P
	resultados2(18,ii) = max(resultados(8,*))	;P

	;renovando condiciones iniciales = CAGÁNDOLA JUAN!!!! ya lo arreglé, por eso esto queda comentado... JA!!!
	;valores iniciales
	;aclouds = mean(resultados(6,*))	;adimensional, 0 para reproducir modelo original, 0.01 para iniciar con area de nubes
	;awhite = mean(resultados(4,*))	;adimensional
	;ablack = mean(resultados(5,*))	;adimensional
	;Ts=mean(resultados(1,*))	;temperatura en la superficie, valor inicial para rk4

endfor

;*************************************************************************************************************************
;RESULTADOS

;TABLAS
openw, 1, dir_salida+'/DAWHYC2_EXP2_L1940_'+strcompress(string(Ac))+'_'+strcompress(string(fix(zc/1000.)))+'_activa.txt'
printf, 1, 'zc='+strcompress(string(fix(zc)))
printf, 1, 'Ac='+strcompress(string(Ac))
printf, 1, format = '(9x,"L", 2x, "Tsmin(C)", 2x, "Tsmed(C)",2x, "Tsmax(C)",5x, "awmin",5x, "awmed",5x, "awmax",5x, "abmin",5x, "abmed",5x, "abmax",5x, "acmin",5x, "acmed",5x, "acmax",6x, "Emin",6x, "Emed",6x, "Emax",6x, "Pmin",6x, "Pmed",6x, "Pmax")'
printf, 1, resultados2, format = '(19f10.4)'
close, 1

ENDFOR	;zc

PRINT, 'Ya terminé!!!'

stop
end




























