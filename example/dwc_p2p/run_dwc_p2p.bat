set field_path=..\data
set field_in=mag_plane_8
set field_out=%field_in%_dwc_p2p.grd
set height_dwc=8
..\..\bin\conti2d %field_path%\%field_in%.grd -G%field_out% -E5 -H0 -T%height_dwc% -D+L10  

pause