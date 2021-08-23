library(usmap)
library(ggplot2)
library(av)

TRI_WIDTH <- 40000
TRI_MAX_HEIGHT <- 2000000
TRI_COLOR <- "red"
TRI_ALPHA <- 0.3
AVG_WINDOW <- 14

BOUND_TRI_X <- c(2640000, 2600000, 2620000, -2140000, -2100000, -2120000)
BOUND_TRI_Y <- c(800000, 800000, 1600000, -1900000, -1900000, -2700000)

FIPS_COLNAME <- "FIPS"
LAT_COLNAME <- "Lat"
LONG_COLNAME <- "Long_"
POP_COLNAME <- "Population"

CASES_PER <- 1000
PNG_DIR <- "./c19_maps_png/"
PNG_EXT <- ".png"

FRAMERATE = 20
END_FRAMES <- 10
OUT_MOVIE_NAME <- "output.mp4"


# https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv
c19_csv <- read.csv("./time_series_covid19_confirmed_US.csv")
# https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv
fips_pop <- read.csv("./UID_ISO_FIPS_LookUp_Table.csv")

state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
                 "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia",
                 "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
                 "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
                 "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
                 "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
                 "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina",
                 "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia",
                 "Washington", "West Virginia", "Wisconsin", "Wyoming")

# Trim out entries with no FIPS code
c19_csv_trim <- c19_csv[which(!is.na(c19_csv[,FIPS_COLNAME])),]
fips_pop_trim <- fips_pop[which(!is.na(fips_pop[,FIPS_COLNAME])),]

# Also trim by state name, because our map doesn't have PR, Guam, USVI, etc.
c19_csv_trim <- c19_csv_trim[which(c19_csv_trim[,"Province_State"] %in% state_names),]
fips_pop_trim <- fips_pop_trim[which(fips_pop_trim[,"Province_State"] %in% state_names),]

# Also trim away entries that have no Latitude/Longitude coordinates
c19_csv_trim <- c19_csv_trim[which(!(c19_csv_trim[,LAT_COLNAME]==0 & c19_csv_trim[,LONG_COLNAME]==0)),]
fips_pop_trim <- fips_pop_trim[which(!(is.na(fips_pop_trim[,LAT_COLNAME]) | is.na(fips_pop_trim[,LAT_COLNAME]))),]

# Add rownames for indexing
rownames(fips_pop_trim) <- fips_pop_trim[,FIPS_COLNAME]
rownames(c19_csv_trim) <- c19_csv_trim[,FIPS_COLNAME]

# Make fips_pop_trim match c19_csv_trim, sorted by FIPS index
fips_pop_trim <- fips_pop_trim[as.character(sort(c19_csv_trim[,FIPS_COLNAME])),] # Note that fips_pop_trim is being indexed by c19_csv_trim values here
c19_csv_trim <- c19_csv_trim[as.character(sort(c19_csv_trim[,FIPS_COLNAME])),]

# Add Population numbers to time_series dataset
c19_csv_trim <- cbind(Population=fips_pop_trim[as.character(c19_csv_trim[,FIPS_COLNAME]),POP_COLNAME], c19_csv_trim)

# Need to transform Lat/Long coordinates for mapping with R package usmap. . .
fips_coords <- data.frame(
  Long_ = c19_csv_trim[,LONG_COLNAME],
  Lat = c19_csv_trim[,LAT_COLNAME],
  FIPS = c19_csv_trim[,FIPS_COLNAME]
)
trans_coords <- usmap_transform(fips_coords)
# Can lose ordering?  Did in previous version.  Re-order by FIPS code to be safe
trans_coords <- trans_coords[order(trans_coords$FIPS),]

# Separate time series dataset from the rest of array, to calculate weekly new cases
total_cases <- c19_csv_trim[,(which(colnames(c19_csv_trim) == "Combined_Key")+1):ncol(c19_csv_trim)]
all_dates <- colnames(total_cases) #grab this for later
# Add n weeks of zeros to left side of matrix (corresponding to dates before COVID-19 data started being collected)
total_cases <- cbind(matrix(data=0, ncol=AVG_WINDOW, nrow=nrow(total_cases)), total_cases)

# Subtract from 2 weeks ago, to get new cases from the past 2 weeks
new_cases <- total_cases[,15:ncol(total_cases)] - total_cases[,1:(ncol(total_cases)-AVG_WINDOW)]
# Add n weeks of zeros again, this time for moving average
new_cases <- cbind(matrix(data=0, ncol=AVG_WINDOW, nrow=nrow(new_cases)), new_cases)


# Define moving average function
mov_avg <- function(x, n){filter(as.numeric(x), rep(1/n, n), sides=1)}
# Apply moving average function.  Note that apply() results in transposed matrix, so we transpose it back
avg_new_cases <- t(apply(new_cases, 1, mov_avg, n=AVG_WINDOW))[,(AVG_WINDOW+1):ncol(new_cases)]

# Fix dates
all_dates <- as.Date(gsub("X", "", gsub("\\.", "/", all_dates)), format = "%m/%d/%y")
# Add back in dates
colnames(avg_new_cases) <- all_dates

# Get number of new cases as a fraction of total population.  Cases per n (typically n = 1000).
avg_new_cases_per <- avg_new_cases/c19_csv_trim[,POP_COLNAME]*CASES_PER

# Normalize for plotting
norm_avg_new_cases_per <- avg_new_cases_per/max(avg_new_cases_per)

# Set up triangle coordinates to plot

trans_long_colname <- paste0(LONG_COLNAME,".1")
trans_lat_colname <- paste0(LAT_COLNAME,".1")

triang_coords <- c()
triang_coords <- cbind(triang_coords, x1=trans_coords[,trans_long_colname]-TRI_WIDTH)
triang_coords <- cbind(triang_coords, x2=trans_coords[,trans_long_colname]+TRI_WIDTH)
triang_coords <- cbind(triang_coords, x3=trans_coords[,trans_long_colname])
triang_coords <- cbind(triang_coords, y1=trans_coords[,trans_lat_colname])
triang_coords <- cbind(triang_coords, y2=trans_coords[,trans_lat_colname])
triang_coords <- cbind(triang_coords, y3=NA) # This will vary by time

# Create png output directory
if (!dir.exists(PNG_DIR)){
  dir.create(PNG_DIR)
}


filter_start <- FALSE # Use this to only start writing images when the cases number is higher than filter

d=1
while (d <= ncol(avg_new_cases_per)) {
  
  cat(paste0(d, "\n"))
  
  triang_coords[,"y3"] = trans_coords[,trans_lat_colname] + norm_avg_new_cases_per[,d]*TRI_MAX_HEIGHT
  filtered_triangs <- triang_coords[which(avg_new_cases_per[,d] > 1),]
  
  if (filter_start) {
    p = plot_usmap() +
      # Two invisible triangles for bottom right and top left, to make sure that
      # map doesn't move on plot with different triangle sizes
      geom_polygon(data=data.frame(x=BOUND_TRI_X, y=BOUND_TRI_Y, group=c(-1,-1,-1,-2,-2,-2)),
                   aes(x=x, y=y, group=group), fill="white", alpha=0) +
      annotate("text", x = -2200000, y = 1600000, hjust = 0,
               label = "Weekly New Coronavirus Cases", size=6) +
      annotate("text", x = -2200000, y = 1200000, hjust = 0,
               label = paste0("Source: https://github.com/CSSEGISandData/COVID-19/\n", all_dates[d]), size=4)
    
      
    if(length(filtered_triangs) > 0) {
      if(length(filtered_triangs) == 6) {
        # this is to deal with R array-to-vector aggravation
        tri_df = data.frame(x=filtered_triangs[c("x1", "x2", "x3")],
                            y=filtered_triangs[c("y1", "y2", "y3")],
                            group=rep(1,3))
      } else {
        tri_df = data.frame(x=c(as.vector(t(filtered_triangs[,c("x1", "x2", "x3")]))),
                            y=c(as.vector(t(filtered_triangs[,c("y1", "y2", "y3")]))),
                            group=rep(1:nrow(filtered_triangs),each=3))
      }
      p = p + geom_polygon(data=tri_df, aes(x=x, y=y, group=group),
                           fill=TRI_COLOR, alpha=TRI_ALPHA)
    }

    legend_x <- c(900000, 1050000, 1200000, 1350000)
    legend_x <- as.vector(rbind(legend_x, legend_x+TRI_WIDTH*2, legend_x+TRI_WIDTH))
    legend_bottom <- -2000000
    legend_scale = max(avg_new_cases_per)
    legend_y <- c(legend_bottom+TRI_MAX_HEIGHT*1/legend_scale, legend_bottom+TRI_MAX_HEIGHT*2/legend_scale,
                  legend_bottom+TRI_MAX_HEIGHT*5/legend_scale, legend_bottom+TRI_MAX_HEIGHT*10/legend_scale)
    legend_y <- as.vector(rbind(rep(legend_bottom, 4), rep(legend_bottom, 4), legend_y))
    p <- p + geom_polygon(data=data.frame(x=legend_x, y=legend_y, group=c(1,1,1,2,2,2,3,3,3,4,4,4)),
                          aes(x=x, y=y, group=group), fill=TRI_COLOR, alpha=0.5)
    p <- p + annotate("text", x=legend_x[1], y=legend_y[1], hjust=-0.2, vjust=1.5, label = "1", size=4)
    p <- p + annotate("text", x=legend_x[4], y=legend_y[4], hjust=-0.2, vjust=1.5, label = "2", size=4)
    p <- p + annotate("text", x=legend_x[7], y=legend_y[7], hjust=-0.2, vjust=1.5, label = "5", size=4)
    p <- p + annotate("text", x=legend_x[10], y=legend_y[10], hjust=0.2, vjust=1.5, label = "10", size=4)
    p <- p + annotate("text", x=legend_x[4], y=legend_y[4], hjust=0.27, vjust=3, label = "Cases per", size=4)
    p <- p + annotate("text", x=legend_x[1], y=legend_y[1], hjust=0.11, vjust=4.5, label = "1,000 people", size=4)
    
    ggsave(paste0(PNG_DIR, all_dates[d], PNG_EXT), p)
  }
  
  if (length(filtered_triangs) > 0 & !(filter_start)) {
    filter_start <- TRUE
    if (d > 1) {
      d <- d - 2
    } else {
      d <- d - 1 # don't want to go back to negative index value if first date has cases over filter level
    }
  }
  d <- d + 1
}

png_files <- list.files(PNG_DIR, pattern = paste0("\\", PNG_EXT, "$"), full.names=TRUE)

av_encode_video(png_files[c(1:length(png_files), rep(length(png_files), END_FRAMES))], OUT_MOVIE_NAME, framerate = FRAMERATE)
