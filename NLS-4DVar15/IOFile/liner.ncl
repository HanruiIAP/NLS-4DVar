

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************

f = new((/101,2/), float)

f         = asciiread ("./rmse.txt",(/101,2/),"float")
;f(0,:)    = asciiread ("./01.txt",40,"float")
;f(1,:)    = asciiread ("./lorenz96.txt",40,"float")

x = ispan(1, 101, 1)

;exit
;************************************************
; plotting parameters
;************************************************
wks   = gsn_open_wks ("png","liner")                  ; open workstation

resources                   = True                      ; plot mods desired


resources@xyLineColors      = (/"red"/)          ; change line color

resources@xyDashPatterns    =(/0/)

resources@trXMinF                =1
resources@trXMaxF                =40


plot  =gsn_csm_xy (wks,x,transpose(f),resources) ; create plot

end













