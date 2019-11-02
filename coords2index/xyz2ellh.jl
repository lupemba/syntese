
"""
    xyz2ellh(X,Y,Z,semi_major_axis=6378137.,flattening=1/298.257223563)

    Go from xyz elipsoidal coordinates to latitude, longitude, height
    (see B.R. Bowring, "The accuracy of geodetic latitude and height equations",
    Survey Review, v28 #218, October 1985 pp.202-206).
"""
function xyz2ellh(X,Y,Z,semi_major_axis=6378137.,flattening=1/298.257223563)
    e2 = flattening*(2-flattening)


    elat=1.e-12
    eht=1.e-5
    p=sqrt(X^2+Y^2)
    lat=atan(Z,p./(1-e2))
    h=0
    dh=1
    dlat=1


    while (dlat>elat) | (dh>eht)
      lat0=lat
      h0=h
      v=semi_major_axis/sqrt(1-e2*sin(lat)*sin(lat))
      h=p*cos(lat)+Z*sin(lat)-(semi_major_axis^2)/v  # Bowring formula
      lat=atan(Z, p*(1-e2*v/(v+h)))
      dlat=abs(lat-lat0)
      dh=abs(h-h0)
    end
    lon=atan(Y,X)
    return lat, lon, h
end
