set field_path=../data
set field_in=mag_plane_0
set field_out=%field_in%_uwc_p2p.vtk
set height_uwc=8
D:\SoftwareProject\conti2d\windows\conti2d\x64\Release\conti2d %field_path%/%field_in%.grd -G%field_out% -H%height_uwc%


pause