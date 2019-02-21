cd('/Users/matthis/Dropbox/UTexas/KateBerkeley2018/2018-04-05-S03')
gpx = gpxread('activity_2603878236.gpx');

figure(1)
geoshow(gpx.Latitude, gpx.Longitude)


figure(2)
subplot(311)
plot(gpx.Latitude,'-o')

subplot(312)
plot(gpx.Longitude,'-o')

subplot(313)
plot(gpx.Elevation,'-o')

webmap('World Topographic Map')
wmline(gpx)


gpx