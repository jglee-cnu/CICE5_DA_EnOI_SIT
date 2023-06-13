module spherical_dists
contains

subroutine haversine_formula(lon1,lat1,lon2,lat2,dist)
implicit none
real,intent(in)::lon1,lon2,lat1,lat2
real,intent(out)::dist
real,parameter::pi=3.141592,mean_earth_radius=6371.0088
real::lonr1,lonr2,latr1,latr2
real::delangl,dellon,dellat,a

lonr1=lon1*(pi/180.);lonr2=lon2*(pi/180.)
latr1=lat1*(pi/180.);latr2=lat2*(pi/180.)
dellon=lonr2-lonr1
dellat=latr2-latr1
a=(sin(dellat/2))**2+cos(latr1)*cos(latr2)*(sin(dellon/2))**2
delangl=2*asin(sqrt(a)) !2*asin(sqrt(a))
dist=delangl*mean_earth_radius

end subroutine

subroutine great_circle_distance(lon1,lat1,lon2,lat2,dist)
implicit none
real,intent(in)::lon1,lon2,lat1,lat2
real,intent(out)::dist
real,parameter::pi=3.141592,mean_earth_radius=6371.0088
real::lonr1,lonr2,latr1,latr2
real::delangl,dellon

lonr1=lon1*(pi/180.);lonr2=lon2*(pi/180.)
latr1=lat1*(pi/180.);latr2=lat2*(pi/180.)
dellon=lonr2-lonr1
delangl=acos(sin(latr1)*sin(latr2)+cos(latr1)*cos(latr2)*cos(dellon))
dist=delangl*mean_earth_radius

end subroutine

subroutine vincenty_formula(lon1,lat1,lon2,lat2,dist)
implicit none
real,intent(in)::lon1,lon2,lat1,lat2
real,intent(out)::dist
real,parameter::pi=3.141592,mean_earth_radius=6371.0088
real::lonr1,lonr2,latr1,latr2
real::delangl,dellon,nom,denom

lonr1=lon1*(pi/180.);lonr2=lon2*(pi/180.)
latr1=lat1*(pi/180.);latr2=lat2*(pi/180.)
dellon=lonr2-lonr1
nom=sqrt((cos(latr2)*sin(dellon))**2. + (cos(latr1)*sin(latr2)-sin(latr1)*cos(latr2)*cos(dellon))**2.)
denom=sin(latr1)*sin(latr2)+cos(latr1)*cos(latr2)*cos(dellon)
delangl=atan2(nom,denom)
dist=delangl*mean_earth_radius

end subroutine
end module
