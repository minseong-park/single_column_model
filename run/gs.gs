'reinit'
'set display white'
'clear'
'set grads off'

taxis=off
var = q


if ( taxis = on )
'sdfopen output.nc'
'set z 1 40'
'set t 1 last'
'd 'var
endif

if ( taxis = off )
'sdfopen output.nc'

'set t last'
'q dims'
st=sublin(result,5);last=subwrd(st,9)
t=1
while(t<=last)
'c'
'set t 't
say t
'set z 1 50'
'set vrange -10 400'
'd 'var
'gxprint x.png'
'!convert x.png -trim -quality 100 't'.png'
'!rm -f x.png'

t=t+1
endwhile


endif
