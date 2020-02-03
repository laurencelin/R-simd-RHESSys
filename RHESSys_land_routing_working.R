library(methods)

source("~/Dropbox/LIB_Rscript/LIB_misc.r")
source("~/Dropbox/LIB_Rscript/LIB_dailytimeseries3.r")
#############################################################
#	assume all patches are grid, 1m, 10m, 30m
#	assume all zones are patches ~ so patch should not be too small
#	assume one stratam per patch
#	patch variable is stored as a vector ~ consistent to GIS readout but without NA cells
#	climate uses only tmax, tmin, precip
#############################################################
	
	#------------------------------------------------------------------------------------------
	SimpsonsYsum = function(x){
		len = length(x); # must be odd number and >= 3
		scale = c(1,rep(c(4,2),len)[1:(len-2)],1); 
		return <- sum(scale*x)
	}#function
	
	TrapezoidalYsum = function(x){
		return <- sum(x) - 0.5*x[1] - 0.5*x[length(x)]
	}
	
	pow = function(a,b){return <- a^b}
	
	# (this part can be read in from stationList.csv)	
	# https://stackoverflow.com/questions/17288197/reading-a-csv-file-organized-horizontally
	# modified by Lin June 23 2018
	read.tcsv = function(file, header=T, sep=",", vskip=0, hskip=0, ...) {
		n = max(count.fields(file, sep=sep), na.rm=TRUE)
		x = readLines(file)
		.splitvar = function(x, sep, n) {
			var = unlist(strsplit(x, split=sep))
		    length(var) = n
		    return(var)
		}#function
		x = do.call(cbind, lapply(x[(vskip+1):length(x)], .splitvar, sep=sep, n=n))
		x = apply(x[(hskip+1):dim(x)[1],], 1, paste, collapse=sep) 
		out = read.csv(text=x, sep=sep, header=header, skip=0, ...)
		return(out)
	}#function	
	
	#------------------------------------------------------------------------------------------
	G_STD = 9.80665 # (m2/s) standard gravitational accel.
	P_STD = 101325.0 # (Pa)standard pressure at 0.0 m elevation
	T_STD = 288.15 # (K) standard temp at 0.0 m elevation
	MA = 28.9644e-3 # (kg/mol) molecular weight of air
	MW = 18.0148e-3 # (kg/mol) molecular weight of water
	CP = 1010.0 # (J/kg*K) specific heat of air
	LR_STD = 0.0065 # (-K/m) standard temperature lapse rate
	RR = 8.3143 # (m3 PA / mol K) gas law constant
	SBC = 5.67e-8 # (W/m2*K4) Stefan Boltzmann COnstant
	EPS = 0.6219 # (MW/MA) unitless ratio of molecular weights
	HVAP = 42.7        # Heat of vaporization, kJ/mol 
	KELVIN = 273.16
	PI = 3.14159265359
	seconds_per_day = 86400
	ICE_DENSITY = 917.0  # (kg/m3) density of ice 
	RAD2PAR = 0.5 # (DIM) ratio PAR/ SWtotal
	EPAR = 4.55 # (umol/J) PAR photon energy ratio
	SECPERRAD = 13750.9871 # seconds per radian of hour angle 
	LITTER_ALBEDO = 0.1  # changed from 0.02 to 0.15 based on Oke 1987 
	WATER_ALBEDO = 0.05 # average liquid water albedo for solar angle of 60 (Dingman) 
	PARTICLE_DENSITY = 2.65 # soil particle density g/cm3 (Dingman) 
	DtoR = pi/180;
	RtoD = 1/DtoR
	Kdown2PAR = 1000 * RAD2PAR * EPAR; ## why 1000? is it # KJ/m2 -> # J/m2 ?
	Io_array = c(1445.0,1431.0,1410.0,1389.0,1368.0,1354.0,1354.0,1375.0,1403.0,1424.0,1438.0,1445.0)
	declination_array = c(-23.0,-22.0,-21.0,-19.0,-17.0,-15.0,-12.0,-9.0,-6.0,-3.0,0.0,3.0,6.0,9.0,12.0,14.0,17.0,19.0,21.0,22.0,
		23.0,23.5,23.5,23.0,21.5,20.0,18.0,16.0,14.0,12.0,9.0,6.0,3.0,0.0,-3.0,-6.0,-9.0,-12.0,-15.0,-17.0,-19.0,-21.0,-22.0,-23.0,-23.5,-23.5)
	declination_array_unique = unique(declination_array)
	declinationUnique2declination = match(declination_array ,declination_array_unique)
	declination_array_unique = declination_array_unique*DtoR
	air_mass_array = c(2.90,3.05,3.21,3.39,3.69,3.82,4.07,4.37,4.72,5.12,5.60,6.18,6.88,7.77,8.90,10.39,12.44,15.36,19.79,26.96,30.00);
	
	

		
	#----------------------- static information ------- should be develop to read from worldfile.csv, landuseList.csv, vegeationList.csv, soilList.csv, stationList.csv
	period = seq.Date(from=as.Date('2000-10-1'), to=as.Date('2002-12-31'),by='day')
	
	# leapYearQ = ifelse(Y%%100==0, Y%%400==0, Y%%4==0)
	station = read.tcsv(file='defs/station_list.csv',stringsAsFactors=F)
	station$max_min_frac = 1/(station$max_snow_temp - station$min_rain_temp)
		#station$stationID
		#station$seriesName	
		#range(worldfile$zoneBaseID)
		
	## ... this is easy for read in and stored; but not cool with BTS and worldfile.zone2station
	climateDTS = list()
	for(i in 1:dim(station)[1]){
		climateDTS[[i]]=read.tcsv(paste('defs/',station$seriesName[i],sep=''), stringsAsFactors=F );
		#climateDTS[[i]]$date = as.Date(climateDTS[[i]]$date,'%d/%m/%Y')	
	}#i
	
	
		
	# just idead here		
	what = read.csv(paste('defs/',station$seriesName[1],sep=''), stringsAsFactors=F)
	what$date = as.Date(paste(what$year,what$month,what$day,sep='-')) 
	match(period,what$date)
		
	
	veg = read.tcsv(file='defs/veg_list.csv', vskip=4, hskip=1, stringsAsFactors=F)
	soil = read.tcsv(file='defs/soil_list.csv', vskip=4, hskip=1, stringsAsFactors=F) 
		#idealy, remove snow parameter in soil; DOM production rate; 

	# single catchment per worldfile
	worldfile = read.csv('worldfiles_workflows/world1000RZsoil3b_average4.csv')
	worldfile.key = cbind(1:dim(worldfile)[2],colnames(worldfile))
	worldfile.patch_index = tapply(seq_len(dim(worldfile)[1]),INDEX=match(worldfile$patchID,unique(worldfile$patchID)),FUN=function(x){x[1]})
	worldfile.zone_index = tapply(seq_len(dim(worldfile)[1]),INDEX=match(worldfile$zoneID,unique(worldfile$zoneID)),FUN=function(x){x[1]})
	worldfile.hill_index = tapply(seq_len(dim(worldfile)[1]),INDEX=match(worldfile$hillID,unique(worldfile$hillID)),FUN=function(x){x[1]})
	worldfile.patch_num = length(worldfile.patch_index)
	worldfile.zone_num = length(worldfile.zone_index)
	worldfile.hill_num = length(worldfile.hill_index)
	worldfile.patch2zone = match(worldfile$zoneID[worldfile.patch_index],worldfile$zoneID[worldfile.zone_index])
	worldfile.patch2hill = match(worldfile$hillID[worldfile.patch_index],worldfile$hillID[worldfile.hill_index])
	#worldfile.patch2soil
	worldfile.patch2soil = 
	worldfile.zone2station = match(worldfile$zoneBaseID[worldfile.zone_index],station$stationID)
	
		# location -> angles 
		location = data.frame(latitude = worldfile$basinLatitude[1]*DtoR) # = catchment_location_Matrix
		location$cos_latitude = cos(location$latitude)
		location$sin_latitude = sin(location$latitude)
		location$tan_latitude = tan(location$latitude)
		
		zone = data.frame( cos_aspect=cos(worldfile$zoneAspect[worldfile.zone_index]*DtoR) ) #--remove?
			zone$sin_aspect = sin(worldfile$zoneAspect[worldfile.zone_index]*DtoR) #--remove?
			zone$cos_slope = cos(worldfile$zoneSlope[worldfile.zone_index]*DtoR) #--remove?
			zone$cos_halfslope = cos(worldfile$zoneSlope[worldfile.zone_index]*DtoR*0.5) #--remove?
			zone$sin_slope = sin(worldfile$zoneSlope[worldfile.zone_index]*DtoR) #--remove?
## Dec 2018, how are the e_horizon and w_horizon included in the calculation?

		zone$pa = P_STD * (1-LR_STD*worldfile$zoneZ[worldfile.zone_index]/T_STD)^(G_STD*MA/LR_STD/RR) 			
		zone$z_delta = worldfile$zoneZ[worldfile.zone_index] - station$z[worldfile.zone2station]  #station elevation
		zone$precipAdjust = (1+ zone$z_delta* station$lapse_rate_precip_default[worldfile.zone2station]) 
		zone$tAdjustDry = zone$z_delta* station$lapse_rate[worldfile.zone2station]
		zone$tAdjustWet = zone$z_delta* station$wet_lapse_rate[worldfile.zone2station] ## <<-- simply this
		zone$kdownA = zone$cos_slope*location$cos_latitude - zone$cos_slope*zone$cos_aspect*location$sin_latitude
		zone$kdownB = zone$cos_aspect*zone$sin_slope*location$cos_latitude + zone$cos_slope*location$sin_latitude
		zone$kdownC = zone$sin_slope*zone$sin_aspect		
		#zone$tsoil = (station_timeseries_tmax[1,hillCSID] -zone$tAdjustDry + station_timeseries_tmin[1,hillCSID] -zone$tAdjustDry)*0.5
			
	
	
		
	#-----------------------
		# simulation daily
		catchmentATM = dailyTimeSeries(seq.Date(from=as.Date('2000-1-1'), to=as.Date('2002-12-31') ,by="day"))
		catchmentATM$declinUniqueID = (declinationUnique2declination[1+catchmentATM$doy %/% 8]-1)*24+1
		## this is daylength for U.S.A. by NOAA (vailded for certain long period of time)
		# NOAA calculation https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
	             	# century starts from 2000 Jan 1st
	    AJday = catchmentATM$doy + (catchmentATM$year-1)*365.25 + 1721423.25
		## ---- -from NOAA
		JulianCentury =(AJday-2451545)/36525  # (YYYY-1)*365.25+MMdd +0.5+1721423.25  -timezone/24  (time zone 0 is Greenwich)
		GeomMean_long_sun_deg = (280.46646+ JulianCentury*(36000.76983 + JulianCentury*0.0003032))%%360
		GeomMean_Anom_sun_deg = 357.52911+ JulianCentury*(35999.05029 - 0.0001537* JulianCentury)
		Sun_eq_Ctr = sin(GeomMean_Anom_sun_deg*DtoR)*(1.914602-JulianCentury*(0.004817+0.000014* JulianCentury)) + 
			sin((2* GeomMean_Anom_sun_deg)*DtoR)*(0.019993-0.000101* JulianCentury)+sin((3* GeomMean_Anom_sun_deg)*DtoR)*0.000289
		Sun_True_long_deg = GeomMean_long_sun_deg + Sun_eq_Ctr
		Sun_app_long_deg = Sun_True_long_deg-0.00569-0.00478*sin( (125.04-1934.136*JulianCentury)*DtoR )
		mean_Obliq_Ecliptic_deg = 23+(26+((21.448-JulianCentury*(46.815+ JulianCentury*(0.00059-JulianCentury*0.001813))))/60)/60
		Obliq_Corr_deg = mean_Obliq_Ecliptic_deg + 0.00256*cos((125.04-1934.136* JulianCentury)*DtoR)
		Sun_delin_deg = asin( sin(Obliq_Corr_deg*DtoR)*sin(Sun_app_long_deg*DtoR) )/DtoR
		HA_sunrise_deg= acos( cos(90.833*DtoR)/(location$cos_latitude*cos(Sun_delin_deg*DtoR)) - location$tan_latitude*tan(Sun_delin_deg*DtoR) )/DtoR 
		catchmentATM$dayhours = (HA_sunrise_deg - 2.076*sqrt(220)/60 )*8/60  #220 m is the basin average DEM 
		catchmentATM$daylength = catchmentATM$dayhours * 3600
		
		## 27 case
		specialNum = 1/2.9 - 1.0e-7
		D24H = rep(seq_len(24),27)
		hour_angle = (D24H*3600-43200)*0.0041667*DtoR
		sqrtatm_trans_estimate = sqrt((station$sea_level_clear_sky_trans[worldfile.zone2station]+
					worldfile$zoneZ[worldfile.zone_index]*station$atm_trans_lapse_rate[worldfile.zone2station])*
					(1-exp(-station$trans_coeff1[worldfile.zone2station])))
		
		
					
		catchment24HATM = data.frame( cos_sza = rep(0,27*24) )
		catchment24HATM$optical_air_mass = rep(0,27*24)
		catchment24HATM$kdown_diffuse_flatsqrtA = rep(0,27*24)
		catchment24HATM$Kdown_directA = rep(0,27*24)
		catchment24HATM$Kdown_directB = rep(0,27*24)
		catchment24HATM$Kdown_directC = rep(0,27*24)
		
		catchment24HATM$cos_sza = cos(declination_array_unique[rep(1:27,each=24)])*location$cos_latitude*cos(hour_angle) + 
					  		      sin(declination_array_unique[rep(1:27,each=24)])*location$sin_latitude; # 24-vector
		
			cond1 = catchment24HATM$cos_sza > 0
			cond2 = specialNum>catchment24HATM$cos_sza & cond1
			cond3 = specialNum<=catchment24HATM$cos_sza & cond1
		catchment24HATM$optical_air_mass[cond3] = 1.0/(catchment24HATM$cos_sza[cond3] + 1.0e-7);
			ML = acos(catchment24HATM$cos_sza[cond2]) %/% 0.0174533 - 69+1; 
		catchment24HATM$optical_air_mass[cond2] = air_mass_array[ML]
		catchment24HATM$kdown_diffuse_flatsqrtA[cond1] = sqrt(catchment24HATM$optical_air_mass[cond1]) 
		catchment24HATM$Kdown_directA[cond1] = cos(declination_array_unique[rep(1:27,each=24)[cond1]])*cos(hour_angle[cond1])*catchment24HATM$optical_air_mass[cond1]
		catchment24HATM$Kdown_directB[cond1] = sin(declination_array_unique[rep(1:27,each=24)[cond1]])*catchment24HATM$optical_air_mass[cond1]
		catchment24HATM$Kdown_directC[cond1] = cos(declination_array_unique[rep(1:27,each=24)[cond1]])*cos(hour_angle[cond1])*catchment24HATM$optical_air_mass[cond1]
		
		i=1; zone$KDDbase1 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=2; zone$KDDbase2 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=3; zone$KDDbase3 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=4; zone$KDDbase4 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=5; zone$KDDbase5 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=6; zone$KDDbase6 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=7; zone$KDDbase7 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=8; zone$KDDbase8 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=9; zone$KDDbase9 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=10; zone$KDDbase10 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		
		i=11; zone$KDDbase11 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=12; zone$KDDbase12 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=13; zone$KDDbase13 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=14; zone$KDDbase14 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=15; zone$KDDbase15 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=16; zone$KDDbase16 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=17; zone$KDDbase17 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=18; zone$KDDbase18 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=19; zone$KDDbase19 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=20; zone$KDDbase20 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		
		i=21; zone$KDDbase21 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=22; zone$KDDbase22 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=23; zone$KDDbase23 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=24; zone$KDDbase24 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=25; zone$KDDbase25 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=26; zone$KDDbase26 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		i=27; zone$KDDbase27 = tapply(zone$kdownA%x%catchment24HATM$Kdown_directA[(24*(i-1)+1):(24*i)] + zone$kdownB%x%catchment24HATM$Kdown_directB[(24*(i-1)+1):(24*i)] - zone$kdownC%x%catchment24HATM$Kdown_directC[(24*(i-1)+1):(24*i)],rep(1:worldfile.zone_num,each=24),function(x){sum(x[x>0])})
		# what = rowSums(zone[,14:40])
		# plot(what,type='l')
		# zone$test <- NULL
		# zone[1:10,]
		
		
	timemax_BTS = dim(period)[1]; #daily
	timemax_STS = 24; #hourly
	checkKdown = matrix(NA,365,4)						 
	system.time({		
		hillCSID = 1
		for(BTS in seq_len(365)){ 
			#daily
		
			## zone / patch [certain day but spatial by patch]
			zone$skyfallwater = climateDTS[[c(1,1)]]$precip_m[BTS]*zone$precipAdjust #worldfile.zone2station
				zone_MatrixCOND = zone$skyfallwater>0
			zone$tmax = station_timeseries_tmax[BTS,hillCSID]
			zone$tmin = station_timeseries_tmin[BTS,hillCSID]
				zone$tmax[zone_MatrixCOND] = zone$tmax[zone_MatrixCOND] - zone$tAdjustWet[zone_MatrixCOND]
				zone$tmin[zone_MatrixCOND] = zone$tmin[zone_MatrixCOND] - zone$tAdjustWet[zone_MatrixCOND]
				zone_MatrixCOND = !zone_MatrixCOND
				zone$tmax[zone_MatrixCOND] = zone$tmax[zone_MatrixCOND] - zone$tAdjustDry[zone_MatrixCOND]
				zone$tmin[zone_MatrixCOND] = zone$tmin[zone_MatrixCOND] - zone$tAdjustDry[zone_MatrixCOND]

			zone$tavg = (zone$tmax + zone$tmin)*0.5;
			zone$tday = station$temcf[worldfile.zone2station]*(zone$tmax - zone$tavg)+zone$tavg;
			zone$tnight = (zone$tday+zone$tmin)*0.5;
			zone$tdewpoint = zone$tmin
				zone_MatrixCOND = station$pptmin[worldfile.zone2station] < zone$skyfallwater
			zone$Delta_T = zone$tmax - zone$tmin; zone$Delta_T[zone_MatrixCOND] = zone$Delta_T[zone_MatrixCOND]*0.75;
			
			
				zone$tmp2 = (station$max_snow_temp[worldfile.zone2station]-zone$tavg)*station$max_min_frac[worldfile.zone2station]
				zone$tmp2[zone$tmp2<0] = 0
			zone$tmp1 = zone$tavg-station$[hillCSID,'min_rain_temp']
			zone_Matrix[zone$tmp1]<=0,'tmp2'] = 1
			zone$snow = zone$tmp2 * zone$skyfallwater
			zone$precip = (1-zone$tmp2) * zone$skyfallwater
			
			zone$es = 613.75*exp( (17.502*zone$tavg)/( 240.97+zone$tavg) ); #---->>> log(a) = log(613)+17*tavg/(240+tvag)
			zone$edewpoint = 613.75*exp( (17.502*zone$tdewpoint)/( 240.97+zone$tdewpoint) );#---->>> log(b)
			zone$vpd = zone$es - zone$edewpoint; zone_Matrix[zone$vpd<0,'vpd']=0; ## log(a-b) = log(a) + log(1-b/a) = log(a) + log(1-rhu)
			zone$rhu = zone$edewpoint/zone$es; ## log(b/a) = log(b) - log(a) --> b/a = exp[log(b) - log(a)]
			
			zone$cloud_fraction = exp(-station_char[hillCSID,'trans_coeff1']*(zone$Delta_T']^station_char[hillCSID,'trans_coeff2']))
			zone$atm_trans = (station_char[hillCSID,'sea_level_clear_sky_trans']+
				zone$z*station_char[hillCSID,'atm_trans_lapse_rate'])*(1-zone$cloud_fraction']);
			zone$sqrtatm_trans = sqrt(zone$atm_trans)
			
			zone$tsoil = 0.9*zone$tsoil + 0.1*zone$tavg
			
			zone$kdown_direct = 0
			zone$kdown_diffure_flat =0
			
		
			#-------------------------------------- hourly loop -------------------------------------------------#
			declinID = catchmentATM[BTS,'declinUniqueID'];
			declinIDname = paste('KDDbase', declinID,sep='')
			
			# basin [all 24 hours of a certain day but aspatial by patch ]
			# catchment24HATM[,'cos_sza'] = catchmentATM[ BTS,'cos_declin']*catchment_location_Matrix['cos_latitude']* catchment24HATM[,'cos_hour_angle'] + 
					  					  # catchmentATM[ BTS,'sin_declin']*catchment_location_Matrix['sin_latitude']; # 24-vector
			# catchment24HATM[,'cos_declin_cos_hourangle']= catchmentATM[ BTS,'cos_declin']*catchment24HATM[,'cos_hour_angle']; 
			# catchment24HATM[,'cos_declin_sin_hourangle']= catchmentATM[ BTS,'cos_declin']*catchment24HATM[,'sin_hour_angle'];
				# optical_air_mass = 1.0/(catchment24HATM[,'cos_sza'] + 1.0e-7); 
				# ML = acos(catchment24HATM[,'cos_sza'])%/%0.0174533-69; ML[ML<1]=1; ML[ML>21]=21;
				# optical_air_mass[optical_air_mass>2.9] = air_mass_array[ML[optical_air_mass>2.9]] # 24-vector. -- cannot simplified more 
			# catchment24HATM[,'optical_air_mass'] = optical_air_mass
			# catchment24HATM[catchment24HATM[,'cos_sza']<=0,c('cos_declin_cos_hourangle','cos_declin_sin_hourangle','optical_air_mass')]=0;
			
					# # these are different by day because of 'cos_declin' and we only need the positive part. we can take out lot of OP
			# catchment24HATM[,'kdown_diffuse_flatA'] = catchment24HATM[,'cos_sza']*sqrt(catchment24HATM[,'optical_air_mass'])
			# catchment24HATM[,'kdown_diffuse_flatB'] = catchment24HATM[,'cos_sza']*catchment24HATM[,'optical_air_mass']
			# catchment24HATM[,'Kdown_directA'] = catchment24HATM[,'cos_declin_cos_hourangle']*catchment24HATM[,'optical_air_mass']
			# catchment24HATM[,'Kdown_directB'] = catchmentATM[BTS ,'sin_declin']*catchment24HATM[,'optical_air_mass']
			# catchment24HATM[,'Kdown_directC'] = -catchment24HATM[,'cos_declin_sin_hourangle']*catchment24HATM[,'optical_air_mass']
			#### from this line above, we should have pre-calculate the 46 cases
			
			for(STS in D24H[catchment24HATM[declinID,,'cos_sza']>0] ){
				
				# basin -- 
				
				# zone / patch [certain hour of a certain day but spatial by patch ]
				# zone$tmp1'] = zone$kdownA']*catchment24HATM[STS,'Kdown_directA']+
										# zone$kdownB']*catchment24HATM[STS,'Kdown_directB']+
										# zone$kdownC']*catchment24HATM[STS,'Kdown_directC']
				# zone$kdown_direct'] = zone$kdown_direct']+ifelse(zone$tmp1']>0,zone$tmp1'],0)
				
				
				# zone$tmp1'] = catchment24HATM[declinID, STS,'kdown_diffuse_flatA']*zone$sqrtatm_trans'] -
										# catchment24HATM[declinID, STS,'kdown_diffuse_flatB']*zone$atm_trans']	
										
				zone$tmp1 = catchment24HATM[declinID, STS,'kdown_diffuse_flatsqrtA']*zone$sqrtatm_trans
				zone_MatrixCOND = zone$tmp1< 1
				zone_Matrix[zone_MatrixCOND,'kdown_diffure_flat'] = zone_Matrix[zone_MatrixCOND,'kdown_diffure_flat'] + 
					catchment24HATM[declinID,STS,'cos_sza']*(1-zone_Matrix[zone_MatrixCOND,'tmp1'])*zone_Matrix[zone_MatrixCOND,'tmp1']	
				
			}#STS (very slow here)
			
						
			
			#basin
				#snowpack calculations
			
			
			#zone / patch [daily; spatial by patch ]
			
			zone$kdown_direct_flat'] = sum(catchment24HATM[declinID,,'cos_sza']*catchment24HATM[declinID,,'optical_air_mass'])*
				Io_array[catchmentATM[BTS,'month']] * zone$atm_trans']*3.6
			zone$kdown_direct'] = zone_Matrix[,declinIDname]*Io_array[catchmentATM[BTS,'month']]*zone$atm_trans'] *3.6
			
			zone$kdown_diffure_flat'] = zone$kdown_diffure_flat']*Io_array[catchmentATM[BTS,'month']]*3.6 #<<-----
			# zone$kdown_diffure_flat'] = (sum(catchment24HATM[declinID,,'kdown_diffuse_flatA'])*zone$sqrtatm_trans'] -
												  # sum(catchment24HATM[declinID,,'kdown_diffuse_flatB'])*zone$atm_trans'])*
												  # Io_array[catchmentATM[BTS,'month']]*3.6 #<<----- problem
												  
												  # zone$sqrtatm_trans']*catchment24HATM[declinID,,'kdown_diffuse_flatsqrtA']
												  
			zone$kdown_diffure'] = zone$kdown_diffure_flat'] * zone$cos_halfslope'] *  zone$cos_halfslope']
			
			# radrat = (Kdown_direct+ Kdown_diffuse)/(kdown_direct_flat+ Kdown_diffuse_flat); radrat[radrat<=0]=1
			zone$PAR_direct'] = zone$kdown_direct'] * Kdown2PAR;
			zone$PAR_diffuse'] = zone$kdown_diffure'] * Kdown2PAR;
			
			 
			
			
			
			# zone$Ldown_night'] = zone$cloud_fraction'] * (SBC * (86400-catchmentATM[BTS,'daylength'])*0.001) * pow(zone$tnight']+273,4) + (1.0 - zone$cloud_fraction']) * ( 59.38 + 113.7*pow((zone$tnight']+273)/273.16,6) + 96.96 * pow(4650*(zone$edewpoint']/1000)/(zone$tnight']+273)/25,0.5)) * (86400-catchmentATM[BTS,'daylength'])*0.001;

			# zone$Ldown_day'] = zone$cloud_fraction'] * (SBC * catchmentATM[BTS,'daylength']/1000) * pow(zone$tday']+273,4)
					# + (1.0 - zone$cloud_fraction']) * ( 59.38 + 113.7*pow((zone$tday']+273)/273.16,6)
					# + 96.96 * pow(4650*(zone$edewpoint']/1000)/(zone$tday']+273)/25,0.5)) * catchmentATM[BTS,'daylength']/1000;
	
			# zone$Ldown'] = zone$Ldown_day'] + zone$Ldown_night'];
		
			
			# zone[0].metv.swavgfd = (zone[0].Kdown_direct + zone[0].Kdown_diffuse)
				# / zone[0].metv.dayl * 1000.0  ;
			# zone[0].metv.ppfd = (zone[0].PAR_direct + zone[0].PAR_diffuse)
				# / zone[0].metv.dayl ;
				
			
			##-------------------------------- vegetation
			# gpp
			# mr, gr
			# some patches do not have veg
			# coverf -- for light, rain intercept, dep (i.e. above ground canopy cover % of the patch area)
			# gapf -- light direct through from canopy
			# ksat0(%impervious) -- for infiltration
				
				
			checkKdown[BTS,] = c(
				mean(zone$kdown_direct']),	
				mean(zone$kdown_direct_flat']),
				mean(zone$kdown_diffure']),
				mean(zone$kdown_diffure_flat'])
				)#c
				
				
				
		}#BTS
	
		## note Ldown = long wave; kdown = shortwave; 
		
		
	})#system.time
		




	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
