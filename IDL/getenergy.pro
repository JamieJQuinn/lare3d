FUNCTION getenergy, wkdir=wkdir


IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

file = wkdir + '/en.dat'
file_info = fstat(100)
IF (file_info.open EQ 1) THEN close,100
openr,100,file
file_info = fstat(100)
fileint = assoc(100,lonarr(1),0,/packed)

; set the size of the variables in bytes
prec = reform(fileint[0])

; set the number of variables (including time) in the energy data
; if you change this you also need to alter the structure below
nv = reform(fileint[1])

; calculate the number of outputs (the -8 is due to the prec/columns output)
outs = (file_info.size - 8) / nv / prec

energy_mask = assoc(100,(prec EQ 4) ? fltarr(nv,outs) : dblarr(nv,outs) ,8,/packed) ; again 8 offset due to prec/columns
data = energy_mask[0]

energy = {points: outs, time: reform(data(0,*)), j_max: reform(data(1,*)), $
         en_b: reform(data(2,*)), en_ke: reform(data(3,*)), en_int: reform(data(4,*)), $
         heating_visc: reform(data(5,*)), heating_ohmic: reform(data(6,*)), eta_perp: reform(data(7,*))}

close,100

RETURN, energy

END
