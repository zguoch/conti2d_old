set field_path=..\data
set field_in=mag_topo
set field_out=%field_in%_dwc_s2p.vtk
set topofile=%field_path%\topo.grd
..\..\bin\conti2d %field_path%\%field_in%.grd -G%field_out% -E0 -H0 -T%topofile% -D+L1500 

pause