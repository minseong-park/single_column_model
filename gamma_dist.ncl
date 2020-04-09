rb_0 = 1.           ; first radius of boundary

rb_r = 1.05         ; raito of radius of boundary

rb_num = 100      ; number of radius of boundary

n0 = 1000           ; total number of droplet

shape   = 8.

scale   = 5.


rb_array = new(rb_num  , "float") ; array of boundary radius
; r_array = new(rb_num-1, "float") ; array of center radius
 PDF     = new(rb_num   , "float") ; PDF
  PB     = new(rb_num-1 , "float") ; PB
   N     = new(rb_num-1 , "float") ; Number of droplets
dN_dlnr  = new(rb_num-1 , "float") ; dN/dlnr

do i = 0, rb_num-1
 rb_array(i) = rb_0*(rb_r)^(i) ; array of boundary radius
end do

;do i = 0, rb_num-2
; r_array(i) = ( rb_array(i) + rb_array(i+1) ) / 2. ; array of center radius
;end do


;== PDF ==


do i = 0, rb_num-1

 PDF(i) = 1. * (  (rb_array(i))^(shape-1) ) * exp(-1 * rb_array(i) / scale ) * ( 1. / ( gamma(shape) * ( scale )^(shape) ) )

end do


;== PB & N ==


do i = 0, rb_num-2

 PB(i) = ( PDF(i+1) + PDF(i) ) * ( ( rb_array(i+1) ) - ( rb_array(i) ) ) / 2.
  N(i) = PB(i)*n0

end do

print("PB sum : "+sum(PB))

print(" N sum : "+sum(N))

;== dN/dlnr ==

do i = 0, rb_num-2

 dN_dlnr(i) = N(i) / ( log(rb_array(i+1)) - log(rb_array(i)) )

end do


w = gsn_open_wks("x11", "Gamma distribution")

res1 = True
res1@tiMainString  = "Gamma distribution ( PDF )"
res1@tiYAxisString = "PDF"
res1@tiXAxisString = "r ( micrometer )"
res1@gsnLeftString = "shape : "+shape+""
res1@gsnRightString= "scale : "+scale+""

p = gsn_csm_xy(w, rb_array, PDF, res1)



res2 = True
res2@tiMainString  = "Gamma distribution ( PB )"
res2@tiYAxisString = "PB"
res2@tiXAxisString = "r ( micrometer )"

p = gsn_csm_xy(w, rb_array(0:rb_num-2), PB, res2)



res3 = True
res3@tiMainString  = "Gamma distribution ( N )"
res3@tiYAxisString = "N"
res3@tiXAxisString = "r"
res3@gsnLeftString = "Total number of droplets : "+n0+""

p = gsn_csm_xy(w, rb_array(0:rb_num-2), N, res3)


;res4 = True
;res4@tiMainString  = "Gamma distribution ( dN/dlnr )"
;res4@tiYAxisString = "dN/dlnr"
;res4@tiXAxisString = "r ( micrometer, log scale )"

;log_x = log(rb_array(0:rb_num-2))

;p = gsn_csm_xy(w, log_x, dN_dlnr ,res4)



