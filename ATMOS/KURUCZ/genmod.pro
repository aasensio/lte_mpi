pro genmod, in, out
	openr,2,in
	c = ''
	for i = 0, 3 do readf,2,c

; Read abundances
	abund = fltarr(92)
	
	readf,2,c
	res = strsplit(c,' ',/extract)
	abund[0:1] = [res[6],res[8]]

	abund_scale = res[2]

	from = 2
	to = from + 5
	for i = 0, 14 do begin
		readf,2,c
		res = strsplit(c,' ',/extract)		
		abund[from:to] = [res[3],res[5],res[7],res[9],res[11],res[13]]
		from = from + 6
		to = to + 6
	endfor

	for i = 0, 2 do readf,2,c
	temp = dblarr(10,72)
	readf,2,temp
	close,2

	h = reform(temp[0,*])
	T = reform(temp[1,*])
	n_e = reform(temp[3,*])
	vturb = reform(temp[6,*])
	vmacro = fltarr(72)

	openw,2,out
	printf,2,'# Model  cmass   T[K]     ne [cm^-3]    vturb [cm/s]'
	printf,2,"'CMASS'"
	printf,2,'MODIFIED_ABUNDANCES'  ; /NORMAL_ABUNDANCES
	printf,2,'ABUNDANCE_SCALE'
	printf,2,abund_scale
	printf,2,'ABUNDANCE_CHANGE'
	printf,2,abund
	printf,2,72
	for i = 0, 71 do begin
		printf,2,h[i],T[i],n_e[i],vturb[i],vmacro[i]
	endfor

	close,2

	stop
end