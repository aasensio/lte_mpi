; In the original database, J is in format F5.2, which is not enough
; for molecular lines. This program reads the Kurucz format, expands
; the size of the J values and rewrites it
pro kurucz_expand_Jsize

; Read linelists
	n_atomic = file_lines('kurucz.list')
	print, 'Reading atomic linelist...'
	str_atomic = strarr(n_atomic)
	openr,2,'kurucz.list'
	readf,2,str_atomic
	close,2

	str_atomic_new = str_atomic

	for i = 0L, n_atomic-1 do begin
		s = str_atomic[i]

		s1 = strmid(s, 0, 11+7+6+12)
		s2 = strmid(s, 11+7+6+12, 5)
		s3 = strmid(s, 11+7+6+12+5, 1+10+12)
		s4 = strmid(s, 11+7+6+12+5+1+10+12, 100)

		str_atomic_new[i] = s1+' '+s2+s3+' '+s4

		if (i / 1000 eq i / 1000.d0) then begin
			print, i,n_atomic
		endif
		
	endfor

	openw,2,'kurucz_expanded.list',width=200
	for i = 0L, n_atomic-1 do begin
		printf,2,str_atomic_new[i]
	endfor
	close,2

	stop
	
end