; does up the ramp non-linear flux correction
;   nam is string root filename
;   nr is number of up the ramp frames
;   md is mean dark image with high s/n reference pixels
;   limt <0 |limt| is data min to consider (slow mode), limt>0 |limt| is data max value to consider
;   a,b,q are intercept (dark), slope (flux), and non-linear (quadratic) flux terms
;   in lsq fit first frame uses independent variable = 0
pro rampfix,nam,nr,md,limt,a,b,q
nx = 2048
d=fltarr(nx,nx)
a=fltarr(nx,nx)
b=fltarr(nx,nx)
q=fltarr(nx,nx)
; compute inverse arrays needed for quadratic fits
c = dblarr(5,nr+1)
arr = dblarr(3,3)
arrt = dblarr(3,3,nr+1)
da = dblarr(nx,nx,3)
; first elemnt not invertible
arrt(*,*,0) = 1.
for ii = 1,nr do begin
  for l = 0,4 do begin
    c(l,ii) = c(l,ii-1) + float(ii-1)^l
  endfor
  arr(0,0)=c(0,ii)
  arr(0,1)=c(1,ii)
  arr(0,2)=c(2,ii)
  arr(1,1)=c(2,ii)
  arr(1,2)=c(3,ii)
  arr(2,2)=c(4,ii)
  arr(1,0)=arr(0,1)
  arr(2,0)=arr(0,2)
  arr(2,1)=arr(1,2)
  if ii gt 2 then begin
    ari = invert(arr)
    arrt(*,*,ii) = ari(*,*)
  endif
endfor
;stop
islow = 0
if limt lt 0 then islow = 1
fl = abs(limt)
m = intarr(nx,nx)
mt = m      ;loop count for good pixels
f = m
f(*,*) = 0  ;count of good pixels up the ramp
; get just ref. pixels
m(*,*) = 1
m(4:nx-5,4:nx-5) = 0
refpix = where(m eq 1)
m(*,*) = 1
for ii = 0,nr-1 do begin
  fn = string(ii,format='(i5.5)')
  frnm = nam + fn +'.fits'
  ; dummy data inserted here
  print,frnm
  ;d(*,*) = 50000.-float(ii)*10 + float(ii)^2*0.1
  fxread,frnm,d
  ; correct reference pixels for median removed mean dark reference pixels (pixel-pixel variations)
  d(refpix) = d(refpix) - md(refpix)  ; xxx remove this line to turn off ref pix-pix correction
  ; find 'saturated' pixels
  if islow eq 1 then bdpix = where(d(*,*) lt limt, nb)
  if islow eq 0 then bdpix = where(d(*,*) gt limt, nb)
  ; correct all data by median of left/rigth reference pixels
  for iy = 4,nx-5 do begin
    corr = (median(d(0:3,iy))+ median(d(nx-4:nx-1,iy)))/2.
    d(4:nx-5,iy) = d(4:nx-5,iy) - corr
  endfor
  mt(*,*) = 1
  if nb gt 0 then mt(bdpix) = 0
  m = m*mt
  ; where f=0 implies zero good frame pixels
  ; f=1 is one good pixel...(not an index from 0)
  f = f + m
  ll = 1.
  ;stop
  for i = 0,2 do begin
    da(*,*,i) = da(*,*,i) + d(*,*)*m*ll
    ll = ll*float(ii)
  endfor
endfor
; now evaluate for all pixels quadratic fits
for ix = 0,nx-1 do begin
  for iy = 0,nx-1 do begin
    pi = f(ix,iy)                   ;pixel count index
    if f(ix,iy) gt 2 then begin       ; only pixels with 3 or more ramp values
      a(ix,iy)=arrt(0,0,pi)*da(ix,iy,0)+arrt(0,1,pi)*da(ix,iy,1)+arrt(0,2,pi)*da(ix,iy,2)
      b(ix,iy)=arrt(1,0,pi)*da(ix,iy,0)+arrt(1,1,pi)*da(ix,iy,1)+arrt(1,2,pi)*da(ix,iy,2)
      q(ix,iy)=arrt(2,0,pi)*da(ix,iy,0)+arrt(2,1,pi)*da(ix,iy,1)+arrt(2,2,pi)*da(ix,iy,2)
    endif
    if f(ix,iy) eq 2 then begin       ; only pixels with 3 or more ramp values
      a(ix,iy)=da(ix,iy,0)-da(ix,iy,1)
      b(ix,iy)=-1.*da(ix,iy,0)+2*da(ix,iy,1)
      q(ix,iy)= 0
    endif
  endfor
endfor
; adjust slope (flux) for negative slow-mode values
if islow eq 1 then b = -b
;stop
end
; get mean dark adjusted reference pixels
;   md is output
;   nam is file name root
;   nr is number of frames in ramp
; 
pro dkrefpix,nam,nr,md
nx = 2048
md = fltarr(nx,nx)
  for ii = 0,nr-1 do begin
    fn = string(ii,format='(i5.5)')
    frnm = nam + fn +'.fits'
    print,frnm
    fxread,frnm,d
    md = md + d
  endfor
  md = md/nr
  for iy = 0, nx-1 do begin
    md(0:3,iy) = md(0:3,iy) - median(md(0:3,iy))
    md(nx-4:nx-1,iy) = md(nx-4:nx-1,iy) - median(md(nx-4:nx-1,iy))
    md(4:nx-5,iy) = 0
  endfor
  end
  ; Find sigma statistics
  ; nam is ramp filename root
  ; nr is number of frames
  ; bpix is bad pixel mask (1 is bad pix)
  ; a,b,q are returned from rampfix
  ; limt defined as in rampfix
  ; sigs is output standard deviation by pixel
  ;
  pro frmsigs,nam,nr,md,bpix,limt,a,b,q,sigs
  islow = 0
  if limt lt 0 then islow = 1
  fl = abs(limt)
  nx = 2048
  ; cumulative sig and pixel count and ref pix mask
  sigs = fltarr(nx,nx)
  f = fltarr(nx,nx)
  mask = fltarr(nx,nx)
  mask(4:nx-5,4:nx-5)=1
  refpix = where(mask lt 0.9)
  for ii = 0,nr-1 do begin
    fn = string(ii,format='(i5.5)')
    frnm = nam + fn +'.fits'
    print,frnm
    fxread,frnm,d
  ; trim reference pixels
  ; correct reference pixels for median removed mean dark reference pixels (pixel-pixel variations)
    d(refpix) = d(refpix) - md(refpix)    ; xxx remove this line to turn off ref pix-pix corrections
    if islow eq 1 then good = where(d(*,*) ge limt, ng)
    if islow eq 0 then good = where(d(*,*) le limt, ng)
  ; correct all data by median of left/rigth reference pixels
    for iy = 4,nx-5 do begin
      corr = (median(d(0:3,iy))+ median(d(nx-4:nx-1,iy)))/2.
      d(4:nx-5,iy) = d(4:nx-5,iy) - corr
    endfor
  ; find 'saturated' pixels

  ; recall that for islow=1 slope must be inverted
    if islow eq 0 then v = a + b*ii + q*float(ii)^2
    if islow eq 1 then v = a - b*ii + q*float(ii)^2
    if ng gt 0 then begin
      sigs(good) = (v(good)-d(good))^2 + sigs(good)
      f(good) = f(good) + 1
    endif 
  endfor
  good = where(f gt 0,ng)
  if ng gt 0 then sigs(good) = sqrt(sigs(good)/f(good))
  ;stop
  end
    
    
