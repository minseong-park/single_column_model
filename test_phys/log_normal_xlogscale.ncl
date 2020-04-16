rb_0    = 1.       ; fisrt radius     (boundary value), unit : micrometer

rb_last = 10000.   ; end radius       (boundary value), unit : micrometer

rb_num  = 1000     ; number of radius (boundary value)

rb_ratio    = (rb_last / rb_0)^( 1. / (rb_num-1) )  ; common ratio (boundary value)

rb_array    = new(rb_num, "float")                  ; array of radius (boundary value)

   r_num    = rb_num-1                              ; r_num = rb_num - 1
 r_array    = new(r_num , "float")                  ; array of radius (center value)

mu    = 2000      ; parameter ( mean )

sigma = 0.5       ; parameter ( standard deviation )

num_droplet = 1000 ; total number of droplets


print("log(mu)  : "+log(mu) )
print("  sigma  : "+sigma   )


;=========== Geometric sequence =========


do i = 0, rb_num -1

    rb_array(i) = rb_0 * ( ( rb_ratio )^(i) )  

end do


;======== Calculate  center value =======


do i = 0, r_num -1

     r_array(i) = rb_array(i) + rb_array(i+1) / 2.

end do


;========= Log-normal distribution ======


PD  = new(rb_num, "float")                  ; PDF (Probability density function)

do i = 0, rb_num -1

 PD(i) = 1. / ( sqrt( 2. * 4.0 * atan(1.0) ) * sigma * rb_array(i) ) *( exp( -1 * (  log( rb_array(i) ) - log( mu ) )^2  / ( 2 * ( sigma )^2 )  ) )

end do


;== Calculate probabilities =============
;== and number of droplets with a PDF ===


PB      = new(rb_num-1, "float")             ; Probability
N       = new(rb_num-1, "float")             ; number of droplets
dN_dlnr = new(rb_num-1, "float")             ; dN/dlnr


do i = 0, rb_num-2 ; 0~998

     PB(i)  = ( ( rb_array(i+1) - rb_array(i) ) * ( PD(i+1) + PD(i) ) ) / 2.
      N(i)  = PB(i)*num_droplet
dN_dlnr(i)  = N(i) / ( log(rb_array(i+1) ) - log(rb_array(i)) )

end do

; print("PB sum : "+sum(PB))
; print(" N sum : "+sum(N))
; print(dN_dlnr)

;=============== Drawing ================


w1 = gsn_open_wks("x11", "PDF"     )
w2 = gsn_open_wks("x11", "N"       )
w3 = gsn_open_wks("x11", "dN_dlnr" )


res1 = True
res1@tiMainString   = "Log normal distribution (PDF)"
res1@tiYAxisString  = "PDF"
res1@tiXAxisString  = "r (micrometer)"
res1@gsnLeftString  = "Mu = "+mu+""
res1@gsnRightString = "Sigma = "+sigma+""

p = gsn_csm_xy(w1, rb_array, PD, res1)

res2 = True
res2@tiMainString   = "Log normal distribution (N)"
res2@tiYAxisString  = "N (number of droplets)"
res2@tiXAxisString  = "r (micrometer)"
res2@gsnLeftString  = "Total number of droplets : "+num_droplet+""

p = gsn_csm_xy(w2, rb_array(0:rb_num-2), PB, res2)

res3 = True
res3@tiMainString   = "Log normal distribution (dN/dlnr)"
res3@tiYAxisString  = "dN/dlnr"
res3@tiXAxisString  = "r (micrometer, Log scale)"


ln_r = log(rb_array(0:rb_num-2))
p = gsn_csm_xy(w3, ln_r, dN_dlnr, res3 )



