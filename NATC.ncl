begin
da   = asciiread("./NA_IBtracks_TCG.csv", -1, "string")
tcinf = da(1:)

yr   = tointeger(str_get_field(tcinf, 2, ","))
mo   = tointeger(str_get_field(tcinf, 3, ","))
lon  = tofloat  (str_get_field(tcinf, 4, ","))
lat  = tofloat  (str_get_field(tcinf, 5, ","))

year            = ispan(1979, 2014, 1)
year@units      = "yr"
TCG           = year 
TCG(:)        = 0

; print(mo)
mo_s = 7
mo_e = 10 
Lons = 260
Lone = 360
Lats = 0
Late = 30

do i = 0, dimsizes(year)-1
    as = ind((yr.eq.year(i)).and.(mo.ge.mo_s).and.(mo.le.mo_e).and.((lon+360).ge.Lons).and.((lon+360).le.Lone).and.(lat.ge.Lats).and.(lat.le.Late) )
    if (.not.ismissing(as(0))) then
        TCG(i) = dimsizes(as)
    end if  
    delete(as)
end do

if (fileexists("./ALTC.txt")) then
    system("rm -rf "+"./ALTC.txt")   
end if   
asciiwrite("./ALTC.txt", TCG)
end 
