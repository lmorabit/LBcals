##########################################
## script to read in relevant catalogs  ##
## to find appropriate LB calibrators   ##
##########################################

require('astro')
require('celestial')
require('stringr')
require('argparser')

p <- arg_parser("Find LB calibrators around a target")
p <- add_argument(p,"--mJIVE_off",flag=TRUE,help="If flag is set, do not search mJIVE catalog")
p <- add_argument(p,"--VLSS_off",flag=TRUE,help="If flag is set, do not search VLSS catalog")
p <- add_argument(p,"--VLSSr_off",flag=TRUE,help="If flag is set, do not search VLSSr catalog")
p <- add_argument(p,"--WENSS_off",flag=TRUE,help="If flag is set, do not search WENSS catalog")
p <- add_argument(p,"--TGSS_off",flag=TRUE,help="If flag is set, do not search TGSS catalog")
p <- add_argument(p,"--RA_deg",type="double",help="Central RA in degrees")
p <- add_argument(p,"--Dec_deg",type="double",help="Central Dec in degrees")
p <- add_argument(p,"--radius",type="double",help="Radius to search in degrees")
p <- add_argument(p,"--filestem",type="character",help="Output filestem (default LBcals_)",default="LBcals_")

args <- parse_args(p)

## parse the flags
if (args$mJIVE_off) fmJIVE <- FALSE else fmJIVE <- TRUE
if (args$VLSS_off) fVLSS <- FALSE else fVLSS <- TRUE
if (args$VLSSr_off) fVLSSr <- FALSE else fVLSSr <- TRUE
if (args$WENSS_off) fWENSS <- FALSE else fWENSS <- TRUE
if (args$TGSS_off) fTGSS <- FALSE else fTGSS <- TRUE

catalog_location <<- './'

search_catalogs <- function(tgtRA, tgtDE, search_rad, filestem, mJIVE=TRUE, vlss=TRUE, vlssr=TRUE, wenss=TRUE, gmrt=TRUE) {

	if(mJIVE) {
		print('Searching the mJIVE catalog ...')
	        mJIVE_cat <- paste( catalog_location, 'mJIVE/catalog.dat', sep="" )
		widths <- c(8,3,3,7,4,3,6,8,8,6,6,6,2,7,6,3,3,10,9,4,3,9,8,6,8,8,6,6,6,6,6,6,6,2)
		colnames <- c('MJV','RAh','RAm','RAs','DEd','DEm','DEs','Fpk','Fint','BMaj','BMin','BPA','l_Vpk','Vpk','e_Vpk','RAVh','RAVm','RAVs','e_RAVs','DEVd','DEVm','DEVs','e_DEVs','V/Fpk','Vint','e_Vint','V/Fint','DMaj','e_DMaj','DMin','e_DMin','DPA','e_DPA','C')
	        Mjive <- read.fwf(mJIVE_cat, widths=widths, col.names=colnames, na.strings="NA", header=FALSE)
		colindex <- c(1,2,3,4,5,6,7,8,9)
		Mjive_slim <- Mjive[,colindex]
		## convert sexigesimal to degrees
		MjiveRA <<- ( Mjive_slim$RAh + (Mjive_slim$RAm + Mjive_slim$RAs/60.)/60. ) * 15.
	        MjiveDE <<- ( Mjive_slim$DEd + (Mjive_slim$DEm + Mjive_slim$DEs/60.)/60. )
	
		# find the distances
        	distances <- cenang( tgtRA, tgtDE, MjiveRA, MjiveDE, units = "deg", method = "vincenty")
		# index the data
                distance_index <- which(distances < search_rad)
		Mjive_filtered <- Mjive_slim[distance_index,]
                distance_filtered <- distances[ distance_index ]
		# now sort it by distance
		Mjive_sorted <- Mjive_filtered[order(distance_filtered),]
		distance_sorted <- distance_filtered[order(distance_filtered)]
		distance_index_sorted <- distance_index[order(distance_filtered)]
        	print(paste('there are ',length(distance_sorted),' sources',sep=""))
                if(length(distance_sorted) > 0) {
			colstoprint <- c('MJV','RAh','RAm','RAs','DEd','DEm','DEs','Fpk','Fint','distance','\n',sep=" ")
                        outfile <- paste(filestem,"_mJIVEmatches.dat",sep="")
                        cat(colstoprint, file=outfile)
                        for(i in seq_along(distance_sorted)) {
                                ss <- paste(c(Mjive_sorted[i,],distance_sorted[i],"\n"),sep=" ")
                                cat(ss, file=outfile,append=TRUE)
                        }
                }

	}

	if(vlss) {
		print('Reading in VLSS catalog ...')
	        VLSS_cat <- paste( catalog_location, 'VLSS/catalog.dat',sep="" )
        	widths <- c(2,3,6,4,3,5,10,2,5,2,5,9,8,7,5,7,5,9,7,5,7)
	        colnames <- c('RAh','RAm','RAs','DEd','DEm','DEs','Si','l_MajAx','MajAx','l_MinAx','MinAx','PA','Field','Xpos','Ypos','e_RAs','e_DEs','e_Si','e_MajAx','e_MinAx','e_PA')
        	vlssDF <- read.fwf(VLSS_cat, widths=widths, col.names=colnames, na.strings="NA", header=FALSE)
		colindex <- c(1,2,3,4,5,6,7)
		vlssDF_slim <- vlssDF[,colindex]
	        ## convert sexigesimal to degrees
        	vlssRA <- ( vlssDF_slim$RAh + (vlssDF_slim$RAm + vlssDF_slim$RAs/60.)/60. ) * 15.
	        vlssDE <- ( vlssDF_slim$DEd + (vlssDF_slim$DEm + vlssDF_slim$DEs/60.)/60. )

                # find the distances
                distances <- cenang( tgtRA, tgtDE, vlssRA, vlssDE, units = "deg", method = "vincenty")
                # index the data
                distance_index <- which(distances < search_rad)
                vlssDF_filtered <- vlssDF_slim[distance_index,]
                distance_filtered <- distances[ distance_index ]
                # now sort it by distance
                vlssDF_sorted <- vlssDF_filtered[order(distance_filtered),]
                distance_sorted <- distance_filtered[order(distance_filtered)]
                distance_index_sorted <- distance_index[order(distance_filtered)]
                print(paste('there are ',length(distance_sorted),' sources',sep=""))
                if(length(distance_sorted) > 0) {
			colstoprint <- c('RAh','RAm','RAs','DEd','DEm','DEs','Si','distance','\n',sep=" ")
                        outfile <- paste(filestem,"_VLSSmatches.dat",sep="")
                        cat(colstoprint, file=outfile)
                        for(i in seq_along(distance_sorted)) {
                                ss <- paste(c(vlssDF_sorted[i,],distance_sorted[i],"\n"),sep=" ")
                                cat(ss, file=outfile,append=TRUE)
                        }
                }

	}

	if(vlssr) {
		print('Reading in VLSSr catalog ...')
		VLSSr_cat <- paste( catalog_location, 'VLSSr/catalog.dat', sep="")
	        colnames <- c('RAdeg','DEdeg','Sp','MajAx','MinAx','PA','e_Sp','Rrms','Rpk','RFlx','Xpos','Ypos','Field','JD')
        	vlssrDF <- read.table(VLSSr_cat, col.names=colnames, header=FALSE, sep="|")
		colindex <- c(1,2,3)
		vlssrDF_slim <- vlssrDF[,colindex]
		# get RA and deg
		vlssrRA <<- vlssrDF_slim$RAdeg
		vlssrDE <<- vlssrDF_slim$DEdeg

                # find the distances
                distances <- cenang( tgtRA, tgtDE, vlssrRA, vlssrDE, units = "deg", method = "vincenty")
                # index the data
                distance_index <- which(distances < search_rad)
                vlssrDF_filtered <- vlssrDF_slim[distance_index,]
                distance_filtered <- distances[ distance_index ]
                # now sort it by distance
                vlssrDF_sorted <- vlssrDF_filtered[order(distance_filtered),]
                distance_sorted <- distance_filtered[order(distance_filtered)]
                distance_index_sorted <- distance_index[order(distance_filtered)]
                print(paste('there are ',length(distance_sorted),' sources',sep=""))
                if(length(distance_sorted) > 0) {
			colstoprint <- c('RAdeg','DEdeg','Sp','distance','\n')
			colindex <- c(1,2,3)
                        outfile <- paste(filestem,"_VLSSrmatches.dat",sep="")
                        cat(colstoprint, file=outfile)
                        for(i in seq_along(distance_sorted)) {
                                ss <- paste(c(vlssrDF_sorted[i,colindex],distance_sorted[i],"\n"),sep=" ")
                                cat(ss, file=outfile,append=TRUE)
                        }
                }

	}

        if(wenss) {
                print('Reading in WENSS catalog ...')
                wenssMain_cat <- paste( catalog_location, 'WENSS/main.dat', sep="")
		widths <- c(15,1,3,3,6,4,3,5,4,3,6,4,3,5,3,3,7,8,5,4,4,5,10)
		colnames <- c('Name','f_Name','RA1950h','RA1950m','RA1950s','DE1950d','DE1950m','DE1950s','RAh','RAm','RAs','DEd','DEm','DEs','flg1','flg2','Speak','Sint','MajAxis','MinAxis','PA','Nse','Frame')
                wenssDF <- read.fwf(wenssMain_cat, widths=widths, col.names=colnames, na.strings="NA", header=FALSE)
		colindex <- c(9,10,11,12,13,14,17,18)
		wenssDF_slim <- wenssDF[,colindex]
                # get RA and deg
                wenssRA <<- ( wenssDF$RAh + (wenssDF$RAm + wenssDF$RAs/60.)/60. ) * 15.
                wenssDE <<- ( wenssDF$DEd + (wenssDF$DEm + wenssDF$DEs/60.)/60. )

                # find the distances
                distances <- cenang( tgtRA, tgtDE, wenssRA, wenssDE, units = "deg", method = "vincenty")
                # index the data
                distance_index <- which(distances < search_rad)
                wenssDF_filtered <- wenssDF_slim[distance_index,]
                distance_filtered <- distances[ distance_index ]
                # now sort it by distance
                wenssDF_sorted <- wenssDF_filtered[order(distance_filtered),]
                distance_sorted <- distance_filtered[order(distance_filtered)]
                distance_index_sorted <- distance_index[order(distance_filtered)]
                print(paste('there are ',length(distance_sorted),' sources',sep=""))
                if(length(distance_sorted) > 0) {
			colstoprint <- c('RAh','RAm','RAs','DEd','DEm','DEs','Speak','Sint','distance','\n')
                        outfile <- paste(filestem,"_WENSSmatches.dat",sep="")
                        cat(colstoprint, file=outfile)
                        for(i in seq_along(distance_sorted)) {
                                ss <- paste(c(wenssDF_sorted[i,],distance_sorted[i],"\n"),sep=" ")
                                cat(ss, file=outfile,append=TRUE)
                        }
                }

        }

        if(gmrt) {
                print('Reading in TGSS catalog ...')
                gmrt_cat <- paste( catalog_location, 'TGSS/TGSSADR1_7sigma_catalog.tsv',sep="")
		gmrt_dat <- read.csv(gmrt_cat, header=TRUE, sep="\t" )
                # get RA and deg
                gmrtRA <<- gmrt_dat$RA
                gmrtDE <<- gmrt_dat$DEC

                colindex <- c(2, 4, 6, 8)
		gmrt_colnames <- names(gmrt_dat)
		gmrt_colnames_slim <- gmrt_colnames[colindex]
                gmrt_dat_slim <- gmrt_dat[,colindex]
		## reorder for output
		gmrt_dat_slim <- gmrt_dat_slim[,c(1,2,4,3)]

                # find the distances
                distances <- cenang( tgtRA, tgtDE, gmrtRA, gmrtDE, units = "deg", method = "vincenty")
                # index the data
                distance_index <- which(distances < search_rad)
                gmrt_filtered <- gmrt_dat_slim[distance_index,]
                distance_filtered <- distances[ distance_index ]
                # now sort it by distance
                gmrt_sorted <- gmrt_filtered[order(distance_filtered),]
                distance_sorted <- distance_filtered[order(distance_filtered)]
                distance_index_sorted <- distance_index[order(distance_filtered)]
                print(paste('there are ',length(distance_sorted),' sources',sep=""))
                if(length(distance_sorted) > 0) {
                        colstoprint <- c('RA','DE','Speak','Sint','distance','\n')
                        outfile <- paste(filestem,"_TGSSmatches.dat",sep="")
                        cat(colstoprint, file=outfile)
                        for(i in seq_along(distance_sorted)) {
                                ss <- paste(c(gmrt_sorted[i,],distance_sorted[i],"\n"),sep=" ")
                                cat(ss, file=outfile,append=TRUE)
                        }
                }

        }



	return()
}

combine_output <- function(filestem) {

	print("Combining catalogues and formatting for Northstar ...")

	outfile <- paste(filestem,'_Northstar.cat',sep="")
	calnum <- 1
	mJIVEfile <- paste(filestem,"_mJIVEmatches.dat",sep="")
	if (file.exists(mJIVEfile)) {
		mJIVE_data <- read.table(mJIVEfile, header=TRUE)
		resolved_ratio <- mJIVE_data$Fpk / mJIVE_data$Fint
		for (i in seq_along(mJIVE_data$MJV)) {
			if ( resolved_ratio > 0.5 && resolved_ratio < 2.0 ) {
				calname <- paste('Cal',calnum,sep="")
				calnum <- calnum + 1
				ra_form <- paste(mJIVE_data$RAh[i], str_pad(mJIVE_data$RAm[i],2,pad="0"), str_pad(format(mJIVE_data$RAs[i],nsmall=2),5,pad="0"), sep=":")
				de_form <- paste(mJIVE_data$DEd[i], str_pad(mJIVE_data$DEm[i],2,pad="0"), str_pad(format(mJIVE_data$DEs[i],nsmall=2),5,pad="0"), sep=":")
				ss <- paste(calname, ra_form, de_form, 'J2000\n',sep="\t")
				cat(ss, file=outfile, append=TRUE)
			}
		}		
	}

	vlssfile <- paste(filestem,"_VLSSmatches.dat",sep="")
	if (file.exists(vlssfile)) {
		vlss_data <- read.table(vlssfile, header=TRUE)
		for (i in seq_along(vlss_data$RAdeg)) {
			calname <- paste('Cal',calnum,sep="")
			calnum <- calnum + 1
			ra_form <- deg2hms(vlss_data$RAdeg[i], type='cat', sep=":")
			de_form <- substring(deg2dms(vlss_data$DEdeg[i], type='cat', sep=":"),2)
			ss <- paste(calname, ra_form, de_form, 'J2000\n',sep="\t")
                        cat(ss, file=outfile, append=TRUE)
		}
	}

        vlssrfile <- paste(filestem,"_VLSSrmatches.dat",sep="")
        if (file.exists(vlssrfile)) {
	        vlssr_data <- read.table(vlssrfile, header=TRUE)
                for (i in seq_along(vlssr_data$RAdeg)) {
        	        calname <- paste('Cal',calnum,sep="")
                        calnum <- calnum + 1
                        ra_form <- deg2hms(vlssr_data$RAdeg[i], type='cat', sep=":")
                        de_form <- substring(deg2dms(vlssr_data$DEdeg[i], type='cat', sep=":"),2)
                        ss <- paste(calname, ra_form, de_form, 'J2000\n',sep="\t")
                        cat(ss, file=outfile, append=TRUE)
		}
        }

        wenssfile <- paste(filestem,"_WENSSmatches.dat",sep="")
        if (file.exists(wenssfile)) {
        	wenss_data <- read.table(wenssfile, header=TRUE)
		for (i in seq_along(wenss_data$RAh)) {
			calname <- paste('Cal',calnum,sep="")
                        calnum <- calnum + 1
			ra_form <- paste(wenss_data$RAh[i], str_pad(wenss_data$RAm[i],2,pad="0"), str_pad(format(wenss_data$RAs[i],nsmall=2),5,pad="0"), sep=":")
			de_form <- paste(wenss_data$DEd[i], str_pad(wenss_data$DEm[i],2,pad="0"), str_pad(format(wenss_data$DEs[i],nsmall=2),5,pad="0"), sep=":")
			ss <- paste(calname, ra_form, de_form, 'J2000\n',sep="\t")
                        cat(ss, file=outfile, append=TRUE)
		}			
        }

	gmrtfile <- paste(filestem,"_TGSSmatches.dat",sep="")
        if (file.exists(gmrtfile)) {
	        gmrt_data <- read.table(gmrtfile, header=TRUE)
                resolved_ratio <- gmrt_data$Speak / gmrt_data$Sint
                for (i in seq_along(gmrt_data$RA)) {
                        if ( resolved_ratio[i] > 0.3 && resolved_ratio[i] < 2.6 ) {
				if ( gmrt_data$Speak[i] > 0.3 ) {
	        	 	        calname <- paste('Cal',calnum,sep="")
        	                	calnum <- calnum + 1
	        	                ra_form <- deg2hms(gmrt_data$RA[i], type='cat', sep=":")
        	        	        de_form <- substring(deg2dms(gmrt_data$DE[i], type='cat', sep=":"),2)
                	        	ss <- paste(calname, ra_form, de_form, 'J2000\n',sep="\t")
	                        	cat(ss, file=outfile, append=TRUE)
				}
			}	
                }
        }

	print(paste('Number of sources written to catalog:',calnum-1,sep=" "))

	return()
}


search_catalogs( args$RA_deg, args$Dec_deg, args$radius, args$filestem, mJIVE=fmJIVE, vlss=fVLSS, vlssr=fVLSSr, wenss=fWENSS, gmrt=fTGSS )
combine_output( args$filestem )
