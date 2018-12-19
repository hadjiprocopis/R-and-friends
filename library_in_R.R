# in case of crash. last change was legend.direction

# if you get:
# Error: Continuous variable () supplied to discrete scale_hue.
# then use scale_.._continuous
#
# this file is good for R2.11, 2.14 requires some adjustments to ggplot2 things - but the most important
# is that that ggplot can not detach legends so... stick with this version for now.

# sets the range between the specified parameter tomin:tomax inclusive
normalise_vector <- function(
	datain=NULL,
	tomin=0.0,
	tomax=1.0,
	fromin=min(datain), # these are either real min/max values if known so not to recalc, or min/max vals of superset data this data is part of
	fromax=max(datain)
){
	l = (tomin - tomax) / (fromin - fromax)
	return( l * (datain-fromin) + tomin )
}
zscore_vector <- function(
	datain=NULL,
	dmean=mean(datain), # specify meand and sd if you already know them so we don't recalculate them in here
	dsd=sd(datain)
){
	return( (datain - dmean) / dsd )
}

process_column_specs_and_read_input_files <- function(
	inputFilesStr="", # a string of input files, separated by param_ifs2
	inputFilesLabelsStr="", # a string of file labels, separated by param_ifs2, can be used for labeling input files (instead of using the filename), if empty, then it will be same with the filename
	doColumnsStr="", # a string of column numbers to do. Each file has its own string separated by param_ifs2. Each file string has its own column numbers separated by param_ifs.
	dontDoColumnsStr="", # a string of 'dont' columns
	param_ifs1, # the IFS (char or string) for separating individual column numbers
	param_ifs2, # the IFS for separating the groups of col numbers per filename
	inputFilesHaveHeaderStr, # True or false if have header in input file - this is again a string but each header spec (T or F) is separated by param_ifs2
	defaultInputFilesHaveHeader=FALSE,
	inputFilesIFSStr, # the column separator in the input files - this is again a string but each header spec (T or F) is separated by param_ifs2
	defaultInputFilesIFS='\t',
	return_back_data_from_files=TRUE # input files must be read but whether their data is returned back is controlled by this - use FALSE if a lot of huge data is opened at the same time, but it is not necessary (i.e. sequential processing) - then all data will be nulled and you need to re-read the files
){
	inputFiles=unlist(strsplit(inputFilesStr, param_ifs2)); num_input_files=length(inputFiles)
	inputFilesLabels=unlist(strsplit(inputFilesLabelsStr, param_ifs2)); for(i in 1:num_input_files){ if( inputFilesLabels[i] == "" ){ inputFilesLabels[i] = inputFiles[i] } }
	doColumns_per_file=unlist(strsplit(doColumnsStr, param_ifs2)); if( length(doColumns_per_file) == 0 ){ doColumns_per_file=c("") }
	dontDoColumns_per_file=unlist(strsplit(dontDoColumnsStr, param_ifs2)); if( length(dontDoColumns_per_file) == 0 ){ dontDoColumns_per_file = c("") }
	oridata=list();
	actual_num_cols=vector(mode='integer', length=num_input_files);
	columns_to_do=list(); column_names=list()
	num_rows=vector(mode='integer', length=num_input_files);
	num_columns_to_do=vector(mode='integer', length=num_input_files);
	num_paddings=vector(mode='integer', length=num_input_files)
	inputFilesIFS=unlist(strsplit(inputFilesIFSStr, param_ifs2)); if( length(inputFilesIFS) == 0 ){ inputFilesIFS=c("") }; inputFilesIFS[inputFilesIFS==""] = defaultInputFilesIFS
	inputFilesHaveHeader=rep(FALSE, times=num_input_files)
	crap=unlist(strsplit(inputFilesHaveHeaderStr, param_ifs2)); i=1; for(cc in crap){if( cc=="TRUE" ){ inputFilesHaveHeader[i] = T } else { if( cc=="" ){inputFilesHaveHeader[i] = defaultInputFilesHaveHeader}}; i = i + 1}
	if( return_back_data_from_files == FALSE ){ nrows_from_file = 1 } else { nrows_from_file = -1 } # read only the header if dont need to return data back
	i=1; for(inputFile in inputFiles){
		write(paste("library_in_R.R::process_input_files_and_column_specs : reading data from '", inputFile, "', header=", inputFilesHaveHeader[i], ", ifs='", inputFilesIFS[i], "' ...", sep=''), stderr())
		tryCatch({oridata[[i]] = read.table(file=inputFile, header=inputFilesHaveHeader[i], sep=inputFilesIFS[i], nrows=nrows_from_file, check.names=FALSE, fill=T)}, error=function(ex){write(paste("library_in_R.R::process_input_files_and_column_specs : error reading file '", inputFile, "'\n",ex,sep=""),stderr()); return(NULL)})
		write("library_in_R.R::process_input_files_and_column_specs : done.", stderr())
		actual_num_cols[i]=ncol(oridata[[i]])
		if( is.null(doColumns_per_file) || (length(doColumns_per_file)==0) || (doColumns_per_file[i] == "") ){
			# no columns specified, do all of them
			columns_to_do[[i]] = seq(from=1, to=actual_num_cols[i])
		} else { columns_to_do[[i]] = as.integer(unlist(strsplit(doColumns_per_file[i], param_ifs1))) }
		if( !(is.na(dontDoColumns_per_file[i]) || (dontDoColumns_per_file[i]=="")) ){ columns_to_do[[i]] = setdiff(columns_to_do[[i]], as.integer(unlist(strsplit(dontDoColumns_per_file[i], param_ifs1)))) }

		num_columns_to_do[i] = length(columns_to_do[[i]]); num_rows[i] = nrow(oridata[[i]])
		if( num_rows[i]==0 ){ write(paste("library_in_R.R::process_input_files_and_column_specs : file '", inputFile, "' : no rows left to do", sep=''), stderr()); return(NULL) }
		if( num_columns_to_do[i]==0 ){ write(paste("library_in_R.R::process_input_files_and_column_specs : file '", inputFile, "' : no columns left to do", sep=''), stderr()); return(NULL) }
#		if( inputFilesHaveHeader[i] == TRUE ){ column_names[[i]] = names(oridata[[i]])[columns_to_do[[i]]] } else { column_names[[i]] = paste("col:", columns_to_do[[i]], sep='') }
		if( inputFilesHaveHeader[i] == TRUE ){ column_names[[i]] = names(oridata[[i]]) } else { column_names[[i]] = paste("col:", seq(from=1, to=actual_num_cols[i]), sep='') }
		col_nam_str = ""; for(j in 1:num_columns_to_do[i]){ if( length(column_names[[i]]) > 0 ){ col_nam_str=paste(col_nam_str, "'", column_names[[i]][columns_to_do[[i]][j]], "'::", columns_to_do[[i]][j], ",", sep='') } else { col_nam_str=paste(col_nam_str, "'NA'::", columns_to_do[[i]][j], ",", sep='') } }
		write(paste("library_in_R.R::process_input_files_and_column_specs : input file: '", inputFile, "', columns to do: ", col_nam_str, sep=''), stderr())
		num_paddings[i] = as.integer(log(actual_num_cols[i])/log(10)+1)
		if( return_back_data_from_files == FALSE ){ oridata[[i]] = NULL } # don't return the data back...
		i = i + 1
	}

	return(list(
		inputFiles=inputFiles,
		num_input_files=num_input_files,
		oridata=oridata,
		actual_num_cols=actual_num_cols,
		columns_to_do=columns_to_do,
		num_rows=num_rows,
		num_columns_to_do=num_columns_to_do,
		column_names=column_names,
		file_labels=inputFilesLabels,
		num_paddings=num_paddings,
		inputFilesHaveHeader=inputFilesHaveHeader,
		inputFilesIFS=inputFilesIFS
	));
}

# breaks is a vector of length numbins+1, the first element is min(dat) and last element is max(dat)
# this is what hist and hist.su output ($breaks)
discretize_data <- function(breaks, dat){
	num_items = length(dat);
	ret=vector(mode='integer', length=num_items)
	num_breaks=length(breaks)
	avrg = vector(mode='double', length=num_breaks-1)
	for(i in 1:(num_breaks-1)){ avrg[i] = (breaks[i]+breaks[i+1])/2.0 }
	for(j in 1:num_items){
		for(i in 2:num_breaks){
			if( dat[j] < breaks[i] ){ ret[j] = avrg[i-1]; break }
		}
	}
	return(ret)
} 
# as above but instead of output the mean value of the bin which is the hist$mids,
# it ASSUMES mids are sorted!!! in ascending order
# it outputs the id of the bin, e.g. 1, 2, etc.
# labels start from 1 <<< 
# breaks start from the values min value or less and end with max value or more (like output hist$breaks)
# returns NULL if error
label_data <- function(
	breaksin,
	dat
){
	num_items = length(dat);
	num_breaksin=length(breaksin)
	print(paste("label_data : for ", num_items, " items and ", num_breaksin, " breaks ...", sep=''), stderr())
	#print(paste(breaksin))
	ret = NULL
	if( num_breaksin == 1 ){
		# if breaks==1, we got a small problem...
		ret = rep(x=0,times=num_items)
	} else {
		ret=vector(mode='integer', length=num_items)
		for(j in 1:num_items){
			adat = dat[j]
			for(i in 2:num_breaksin){
				#print(paste("checking ", adat, " < ", breaksin[i]))
				if( adat <= breaksin[i] ){
					ret[j] = i-1;
					break
				}
			}
			if( ret[j] == 0 ){
				write(paste("label_data : can not find appropriate bin for value '", adat, "', breaks are: (", paste(breaksin, ",", sep=''), ").", sep=''), stderr());
				return(NULL)
			}
		}
	}
print(ret)
	return(ret)
}
# calculate bins(breaks) for vector optionally specify: numbins or breaks or nothing or equal height
# returns NULL if error
calculate_histogram_breaks <- function(values, numbins=NA, equal_height_bins=FALSE, density_num_points=512, density_kernel='biweight', dens=NULL, use.ks=TRUE, calculate.error=FALSE){
	breaks = NULL; den = dens; error = NA
	spn = diff(range(values, na.rm=T)); mini=min(values, na.rm=T); maxi=max(values, na.rm=T)
	if( is.na(numbins) ){
#		if( is.null(den) ){ den = density(values,n=density_num_points,na.rm=T,kernel=density_kernel) }; numbins = den$bw
		ah=hist(values,breaks="FD", plot=F); numbins=length(ah$mids)
	}
	if( (numbins == 1) || (mini==maxi) ){
# this is the case where no variation in data and therefore we create 2 bins otherwise ...
		mini = mini - 0.02*abs(mini)
		maxi = maxi + 0.02*abs(maxi)
		spn = maxi - mini
	}
	if( equal_height_bins == FALSE ){
		breaks = seq(from=mini, to=maxi, by=(spn/numbins));
		# this extends the breaks slightly left and right to include ...
		breaks[1] = mini - 0.02*abs(mini); breaks[length(breaks)] = maxi + 0.02*abs(maxi)
		return(list(breaks=breaks,density=NULL,error=NA))
	} else if( equal_height_bins == "bins" ){
		breaks=vector(mode='double',length=(numbins+1))
		numvalues = length(values)
		sortvalues = sort(values)
		avecounts = as.integer(numvalues/numbins+0.5) # counts per bin average
		j = 1; for(i in seq(1,numvalues,avecounts)){ breaks[j] = sortvalues[i]; j = j + 1 }
		breaks[1] = mini - 0.02*abs(mini); breaks[length(breaks)] = maxi + 0.02*abs(maxi)
		return(list(breaks=breaks,density=NULL,error=error))
	}
	write(paste("calculate_histogram_breaks : don't know what to do with this equal_height_bins value : '",equal_height_bins,"', can be FALSE, 'bin', 'density'. But 'density' is obsolete.", sep=""), stderr())
	return(NULL)
}
# values: a matrix of N-dimensional vectors 
# breaks : a list of N vectors of possibly different length, length=number of breaks
# each vector represents the breaks in a given dimension 1..N
# for one dimension, this works as:
# values = matrix(...)
# breaks = l(c(...))
# colprefix = is what to prefix the columns (1..N) holding the bin number. e.g. with default is label1, label2, ...labelN
# rangeprefix = is what to prefix the columns (1..N) holding the range for the corresponding bin number (label).
calculate_histogram <-function(values, breaks, colprefix="binlabel", rangeprefix="binrange"){
	ndims <- ncol(values); nrows <- nrow(values)
	nbreaks <- lapply(breaks, length)
	if( length(breaks) != ndims ){ write(paste(sep='',Sys.getenv("HOME"), "/bin/library_in_R.R::calculate_histogram : number of columns of matrix 'values' (",ndims,") must be the same as the number of vectors in the list 'breaks' (",length(breaks),"). Can't continue", sep=""), stderr()); return }
# this is the structure of the returned matrix, add new cols here but also change the dims of the matrix (ncol) below
# the matrix is created with one dummy row (removeme)
	counts <- matrix(ncol=1+1+ndims+ndims+ndims,dimnames=list(
		c('removeme'), # the row name
	# these are the col names:
		c('counts', 'probabilities',
			paste(colprefix,seq(1,ndims,1),sep=''), # the histogram labels (i.e. bin number) for each input col
			paste(paste(rangeprefix,"_from",sep=''), seq(1,ndims,1),sep=''), # the FROM range for the bin (from:to) for each input col
			paste(paste(rangeprefix,"_to",sep=''), seq(1,ndims,1),sep='') # the FROM range for the bin (from:to) for each input col
		)
	)) # first is a string of indices - one to each dimension - to the bin, second is the count, e.g. "1,2,3" 5 means box i=1,j=2,k=3 (i,j,k=dimensions) count=5
	# we will label each dimension to its own breaks, then use the labels as index to array of counts
	labels <- list()
	for(i in 1:ndims){ labels[[i]] <- label_data(breaksin=breaks[[i]], dat=unlist(values[i])) }
	ind = vector(mode='integer', length=ndims)
	ranges_from = vector(mode='numeric', length=ndims)
	ranges_to = vector(mode='numeric', length=ndims)
	strind = "" ; num_adds=1
#browser()
	for(r in 1:nrows){
		for(i in 1:ndims){
			ind[i] = labels[[i]][r]
			ranges_from[i] = breaks[[i]][ind[i]]
			ranges_to[i] = breaks[[i]][ind[i]+1]
		}
		strind = paste(ind, collapse=",", sep="")
		if( strind %in% rownames(counts) ){
			counts[strind,'counts'] = counts[strind,'counts'] + 1
		} else {
			num_adds = num_adds + 1
			counts=rbind(
				counts, c(
				# these are the default values (when creating a new empty row) for each col of output, length is important, if you add more columns (see above) you need to add more things here
					1, # counts
					NA, # probabilities
					ind, # labels
					ranges_from, # ranges from
					ranges_to # ranges to
				)
			)
			rownames(counts)[num_adds] <- strind
		}
	}
	counts[,'probabilities'] = counts[,'counts'] / sum(counts[,'counts'], na.rm=T)
	return(counts[-1,,drop=FALSE]) # shitty R drops dimensions automatically,... so 1 row matrix becomes shit
	# the result is a sparse matrix, i.e. those bins with zero count are not included!!!
}
# the returned histogram from calculate_histogram does not include zero-count bins - this functions corrects this if you need it.
# a_histogram : a matrix representing a histogram previously calculated using calculate_histogram function
# numbins : a vector of bin-numbers per input dimension of the histogram.
unsparse_histogram <- function(a_histogram, numbins){
	ret <- a_histogram
	labels = list()
	for(i in 1:length(numbins)){ labels[[i]] = seq(1,numbins[i],1) }
	all_combinations = expand.grid(labels)
	extra_columns = ncol(a_histogram)-2-length(numbins)
	rnames = rownames(a_histogram); num_adds=length(rnames)+1
	for(arow in 1:nrow(all_combinations)){
		comb_label = paste(all_combinations[arow,], collapse=',')
		if( !(comb_label %in% rnames) ){
			if( extra_columns > 0 ){
				ret = rbind(ret, c(0,0.0,unlist(all_combinations[arow,]),rep(x=NA,times=extra_columns)))
			} else { 
				ret = rbind(ret, c(0,0.0,unlist(all_combinations[arow,])))
			}
			rownames(ret)[num_adds] <- comb_label
			num_adds = num_adds + 1
		}
	}
	return(ret)
}

# Entropy of a set of values in a vector representing probabilities (i.e. sum=1, >=0). will be ignored if prob=0
# this stupid R returns a list() whenever all elements are 0 and sum() fails.
#simple_entropy <- function(a_vec){ sum(sapply(a_vec[a_vec>0], function(p){-p*log10(p)}), na.rm=T) }
# so use the for-loop, it looks same time...
#simple_entropy <- function(a_vec){ sum=0; for(i in a_vec[a_vec>0]){ sum = sum - i*log10(i) }; return(sum) }
shannon_entropy <- function(p){
	if (min(p) < 0 || sum(p) <= 0){ return(NA) }
	p.norm <- p[p>0]/sum(p)
	return( -sum(log2(p.norm)*p.norm) )
}
simple_entropy <- function(p){
	if (min(p) < 0 || sum(p) <= 0){ return(NA) }
	p.norm <- p[p>0]/sum(p)
	return( -sum(log10(p.norm)*p.norm) )
}

plot_heatmap <- function(
# the dataframe with 3/6 entries: x, y and values. optional cols: frame_color, alpha, extra_values
# for your own safety, these columns contain numerical data, for frame_color and alpha, ggplot will determine
# the actual colors and alphas given the variation in these columns
# the last entry 'extra_values' is extra data to be printed on each cell using its own sprintf format (see below)
	datain=NULL,

	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	background=element_blank(),
	cell.frame_color=alpha('black',0.25), # the border around each cell
	cell.fill=NA, # the fill color of each cell this is usually controlled by the input data col 'values'
	fontsize=20,
	axis.label.x=NULL,
	axis.label.y=NULL,
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.x.angle=0, # or 90 for rotating 90 degrees
	axis.labels.x.fontsize=0.8*fontsize,
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y.angle=0,
	axis.labels.y.fontsize=0.8*fontsize,

	box.labels.bgcol="green", # how to draw the boxes
	box.labels.bgfill="skyblue",
	box.labels.bgalpha=0.8,
	box.labels.fontsize=3.5, # fontsize of the box , set it to zero in order not to print any labels
# you can also put anything in this format as it goes through sprintf e.g. "\nxxx=%.2f\nahahah\txxx"
	box.labels.sprintfFormat="%.2f", # how to print the 'values' column of the data to the cells (this assumes floats) 
	box.labels.extra_values.sprintfFormat="%s", # format for the extra_values of the data (good luck with the string!!! maybe data.frame(newdata, stringsAsFactors=F) can help.
#	box.labels.palette.start=alpha('cornflowerblue',0.2), # the 3 colors for the gradient of the boxes
#	box.labels.palette.middle=alpha('white',0.2),
#	box.labels.palette.end=alpha('red',0.2),
	box.labels.palette.start='blue', # the 3 colors for the gradient of the boxes
#	box.labels.palette.middle='white',
	box.labels.palette.middle=NULL,
	box.labels.palette.end='red',
#	box.labels.middle=mean(datain$values), # the midpoint of the range where the mid color will be mapped
	box.labels.middle=NULL,
	box.labels.breaks=NULL, # a c(1,2,3) representing the legend values (i.e. here 1, 2 and 3), default is auto
	plot.title=NULL,
	legend.fontsize=12,
	legend.direction="horizontal", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
	legend.title=NULL,
	legend.show=TRUE
){
# the tickmarks are in the centre of each box, if don't like it: change +- 0.5 below
	if( is.null(datain) ){ write("library_in_R.R::plot_heatmap : datain is required and is a dataframe with 'x', 'y', 'values' and optional 'extra_values' (for extra labels) and 'color' for frame-color if required (NA for those not required)", stderr()); return(NULL) }
	datain_colnames = colnames(datain)
	if( ! "x" %in% datain_colnames ){ write("library_in_R.R::plot_heatmap : datain requires 'x' as column.", stderr()); return(NULL); }
	if( ! "y" %in% datain_colnames ){ write("library_in_R.R::plot_heatmap : datain requires 'y' as column.", stderr()); return(NULL); }
	if( ! "values" %in% datain_colnames ){ write("library_in_R.R::plot_heatmap : datain requires 'values' as column.", stderr()); return(NULL); }

	suppressMessages(require(ggplot2))
	p <-ggplot()
	aes_string = "aes(xmin=x-0.5,xmax=x+0.5,ymin=y-0.5,ymax=y+0.5,fill=values"
	other_string = ""
	if( "frame_color" %in% colnames(datain) ){ aes_string = paste(aes_string, ",colour=frame_color", sep='') }
	if( "alpha" %in% colnames(datain) ){ aes_string = paste(aes_string, ",alpha=alpha", sep='') }
	aes_string = paste(aes_string, ")", sep='')
#write(paste("aesstring: ", aes_string), stderr())
	if( !("frame_color" %in% colnames(datain)) && !is.null(cell.frame_color) ){
		p <- p + geom_rect(
			data=datain,
			eval(parse(text=aes_string)),
			colour=cell.frame_color # for this we have to make this if-stat
		)
	} else { 
		p <- p + geom_rect(
			data=datain,
			eval(parse(text=aes_string))
		)
	}
	p <- p + scale_alpha_continuous(guide=F)+scale_color_continuous(guide=F)
#	p <- p + geom_rect(data=datain,eval(parse(text=aes_string)),color=cell.frame_color)+scale_alpha_continuous(legend=F)+scale_color_continuous(legend=F)
	if( is.null(axis.labels.x) ){
		# don't show x-axis or its labels or tickmarks (show x-axis title)
		# unfortunately hiding tickmarks is for both X,Y axes
		p <- p + theme(axis.text.x=element_blank(),axis.ticks=element_blank())+scale_x_continuous(breaks=NULL)
	} else {
		p <- p + theme(axis.text.x=element_text(size=axis.labels.x.fontsize,angle=axis.labels.x.angle,hjust=1,vjust=1))+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)
	}
	if( is.null(axis.labels.y) ){
		# don't show y-axis or its labels or tickmarks (show y-axis title)
		p <- p + theme(axis.text.y=element_blank(),axis.ticks=element_blank())+scale_y_continuous(breaks=NULL)
	} else {
		p <- p + theme(axis.text.y=element_text(size=axis.labels.y.fontsize,angle=axis.labels.y.angle,hjust=1))+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)
	}

	p <- p + theme(panel.background=background,panel.grid.major=gridline.major,panel.grid.minor=gridline.minor)+
	xlab(axis.label.x)+ylab(axis.label.y)+
	theme(axis.title.x=element_text(size=fontsize),axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5)) # magic numbers

	# no variation in values (which is the counts basically for each x,y coords) causes some problem if legend.direction="horizontal"
	if( sd(datain$values) < 0.00001 ){ legend.direction="vertical" }

	if( !is.null(axis.limits.x) ){ p <- p + xlim(axis.limits.x) }
	if( !is.null(axis.limits.y) ){ p <- p + ylim(axis.limits.y) }

	p <- p + theme(title=plot.title,plot.title=element_text(size=fontsize),legend.direction=legend.direction)+
	scale_alpha_continuous(guide=FALSE)
	dummy_str = "scale_fill_gradient2(name=legend.title,low=box.labels.palette.start,high=box.labels.palette.end,breaks=box.labels.breaks,guide=legend.show"
	if( !is.null(box.labels.middle) ) dummy_str = paste(dummy_str, ",midpoint=box.labels.middle",sep='')
	if( !is.null(box.labels.palette.middle) ) dummy_str = paste(dummy_str, ",box.labels.palette.middle",sep='')
	dummy_str = paste(dummy_str, ")", sep='')
	p <- p + eval(parse(text=dummy_str))

#	scale_fill_continuous(name=legend.title)+scale_alpha_continuous(name=legend.title)+
#	scale_fill_gradient2(name=legend.title,limits=c(box.labels.start,box.labels.end),low=box.labels.palette.start,mid=box.labels.palette.middle,high=box.labels.palette.end,breaks=box.labels.breaks)+#	scale_fill_continuous(name=legend.title)+scale_alpha_continuous(name=legend.title)+
#c(rgb(colorRamp(c(box.labels.palette.start,box.labels.palette.end))(seq(0,1,length=box.labels.palette.num_shades)), maxColorValue=255)))+

# superhack: aes can not read variables from environment so do that to add mxc
#p <- p + geom_text(eval(parse(text=paste('aes(x=(breaksFrom+breaksTo)/2,y=counts+',(0.1*max(data.boxes$counts)),',label=counts)', sep=''))),size=box.height_labels.fontsize_multiplier*fontsize)

	if( box.labels.fontsize > 0.0 ){
		if( "extra_values" %in% colnames(datain) ){
			p <- p + geom_text(data=datain,eval(parse(text=paste('aes(x=x,y=y,label=sprintf("', box.labels.sprintfFormat, box.labels.extra_values.sprintfFormat, '",values, extra_values))', sep=''))),size=box.labels.fontsize,bgcol=box.labels.bgcol,bgalpha=box.labels.bgalpha)
#aes(x=x+0.5,y=y+0.5,label=sprintf(paste(box.labels.sprintfFormat,"\n%s",sep=''),values, extra_values)),size=box.labels.fontsize,bgcol=box.labels.bgcol,bgalpha=box.labels.bgalpha)
		} else {
			p <- p + geom_text(data=datain,eval(parse(text=paste('aes(x=x,y=y,label=sprintf("', box.labels.sprintfFormat, '",values))', sep=''))),size=box.labels.fontsize,bgcol=box.labels.bgcol,bgalpha=box.labels.bgalpha)
		}
	}
	return(p);
}
# will plot one more density functions of the given data in the dataframe - the density is calculated here using ggplot
# if you have calculated density already then use plot_points (see below).
plot_density <- function(
	datain, # the dataframe with 2 columns: datax and class, datax is the data and name is the class the data belongs to.
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	density.kernel='biweight',

	density.alpha=0.25,
	background=element_blank(),

	histogram.show=FALSE, # add optional histogram with...
	histogram.binwidth=0, # binwith==0: automatic calculation or set it >0
	histogram.alpha=0.2,
	histogram.fill='blue',

	fontsize=20,
	axis.label.x="data",
	axis.label.y="density",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	plot.title=NULL,

	legend.title=NULL,
	legend.fontsize=12,
	legend.direction="horizontal", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
	legend.position="right",
	legend.show=TRUE
){
	suppressMessages(require(ggplot2))
	data_class_names <- unique(datain[,'class']); num_classes <- length(data_class_names)

	p <- ggplot(datain,aes(x=datax,colour=class,fill=class))+
	geom_density(alpha=density.alpha,kernel=density.kernel)+
	theme(panel.grid.major=gridline.major,panel.grid.minor=gridline.minor,panel.background=background)+
	theme(axis.text.x=element_text(size=fontsize*0.8)) + # these are magic numbers for just and pointsize
	theme(axis.text.y=element_text(size=fontsize*0.8,hjust=1))+
	theme(axis.title.x=element_text(size=fontsize))+
	theme(axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5))+
	xlab(axis.label.x)+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)+
	ylab(axis.label.y)+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)+
# resort to a combination of these limit shits. coord on x-axis does not work...
	coord_cartesian(ylim=axis.limits.y)+
#	xlim(axis.limits.x)+ylim(axis.limits.y)+ # put the limits after the scale_x...!!!!, alternative is coord_cartesian(ylim=axis.limits.y,xlim=axis.limits.x)
	theme(title=plot.title,plot.title=element_text(size=fontsize),legend.position=legend.position,legend.text=element_text(size=legend.fontsize,hjust=0))
	if( !is.null(axis.limits.x) ){ p <- p + xlim(axis.limits.x) }
	if( histogram.show == TRUE ){
		if( histogram.binwidth > 0 ){ p <- p + geom_histogram(aes(y = ..density..),binwidth=histogram.binwidth,alpha=histogram.alpha,fill=histogram.fill) }
		else { p <- p + geom_histogram(aes(y = ..density..),alpha=histogram.alpha,fill=histogram.fill) }
	}
	if( num_classes > 1 ){
		p <- p + scale_color_brewer(guide=FALSE)+
#		scale_fill_brewer(name=legend.title,legend=legend.show,legend.direction=legend.direction)+
		scale_fill_brewer(name=legend.title,guide=legend.show)
	}
	return(p);
}
# plot multiple y-values with a common x-value each as pair as a point, e.g. x=1, y1=1,y2=2, or just x=1,y=1, etc.
plot_points <- function(
	datain, # a list of dataframes each with a x and y
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	class.labels=NULL, # data ontains 'class' col which may contain string labels or integers
			   # if it's integers, then supply this vector of labels to be placed in the legend
			   # if class.labels is not null, it overwrites the class cols
	class.different_colors=TRUE, # FALSE (for no different colors), TRUE (for different colors, automatic setup) or a vector of colors for manual setup
	class.different_shapes=FALSE, # specify if each class is drawn with different shape
	points.size=0.5,
	points.alpha=1.0,
	background=element_blank(),
	point.color=c(alpha('black',0.25)), # the color of each point
	fontsize=20,
	axis.label.x="x",
	axis.label.y="y",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	plot.title=NULL,

	legend.title=NULL,
	legend.show=TRUE,
	legend.direction="horizontal" # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
){
	suppressMessages(require(ggplot2))
	num_classes <- length(datain)
	if( is.null(class.labels) ){ all_labels <- 1:num_classes } else { all_labels <- class.labels }
	str = "a_list<-list("; for(i in 1:num_classes){ str = paste(str,"df",all_labels[i],"=datain[[",i,"]],",sep='') }; str=substr(str, 1, nchar(str)-1); str=paste(str,")",sep='')
	eval(parse(text=str))
	newData = melt(a_list, id.vars='x')
	p <- ggplot()+
	geom_line(data=newData,aes(x=x,y=value,color=L1),size=points.size,alpha=points.alpha)+
	theme(panel.grid.major=gridline.major,panel.grid.minor=gridline.minor,panel.background=background)+
	theme(axis.text.x=element_text(size=fontsize*0.8)) + # these are magic numbers for just and pointsize
	theme(axis.text.y=element_text(size=fontsize*0.8,hjust=1))+
	theme(axis.title.x=element_text(size=fontsize))+
	theme(axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5))+
#	coord_cartesian(ylim=axis.limits.y,xlim=axis.limits.x)+
	xlab(axis.label.x)+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)+
	ylab(axis.label.y)+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)+
	theme(title=plot.title,plot.title=element_text(size=fontsize))

	if( !is.null(axis.limits.x) ){ p <- p + xlim(axis.limits.x) }
	if( !is.null(axis.limits.y) ){ p <- p + ylim(axis.limits.y) }

	if( !is.null(class.different_colors) ){
		p <- p + scale_color_manual(breaks=paste("df",all_labels,sep=''),labels=all_labels,values=class.different_colors,name=legend.title,guide=legend.show)
	} else { scale_color_manual(guide=legend.show) }
	return(p);
}
# like histogram bars
plot_boxes <- function(
	data.boxes, # the dataframe with 3 entries: breaksFrom,breaksTo, counts (e.g. (0,1,2),(1,2,3) and 5,7,8: between 0,1->5, 1,2->7, 2,3->8)
	data.lines=NULL, # optional other data to be plotted, this is a list of lists, each list contains: a dataframe, ($data) of (x,y) pairs
			# optional colour name ($color), the legend name is taken from the list's $, e.g. data.lines$fuck, data.lines$fick
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2),
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	background=element_blank(),
	box.color.outline=alpha('black',0.75),
	box.color.fill=alpha('#D55E00',0.75),
	fontsize=20,
	axis.label.x="x",
	axis.label.y="y",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	box.height_labels.show=TRUE,
	box.height_labels.fontsize_multiplier=0.25, # multiply fontsize by this parameter to get fontsize of height labels, may be you want them smaller
	plot.title=NULL,
	legend.title=NULL,
	legend.direction="horizontal" # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
){
	suppressMessages(require(ggplot2))
	p <- ggplot(data=data.boxes)
	if( box.height_labels.show == TRUE ){
		# superhack: aes can not read variables from environment so do that to add mxc
		p <- p + geom_text(eval(parse(text=paste('aes(x=(breaksFrom+breaksTo)/2,y=counts+',(0.1*max(data.boxes$counts)),',label=counts)', sep=''))),size=box.height_labels.fontsize_multiplier*fontsize)
	}
	# this is totally stupid but...
	if( is.null(box.color.outline) ){
		if( is.null(box.color.fill) ){
			p <- p +geom_rect(aes(xmin=breaksFrom,xmax=breaksTo,ymin=0,ymax=counts))
		} else { 
			p <- p +geom_rect(aes(xmin=breaksFrom,xmax=breaksTo,ymin=0,ymax=counts),fill=box.color.fill)
		}
	} else { 
		if( is.null(box.color.fill) ){
			p <- p +geom_rect(aes(xmin=breaksFrom,xmax=breaksTo,ymin=0,ymax=counts),colour=box.color.outline)
		} else { 
			p <- p +geom_rect(aes(xmin=breaksFrom,xmax=breaksTo,ymin=0,ymax=counts),colour=box.color.outline,fill=box.color.fill)
		}
	}
	p <- p + theme(panel.background=background,panel.grid.major=gridline.major,panel.grid.minor=gridline.minor)+
	xlab(axis.label.x)+ylab(axis.label.y)+
	theme(axis.title.x=element_text(size=fontsize),axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5)) # magic numbers

	if( is.null(axis.labels.x) ){
		# don't show x-axis or its labels or tickmarks (show x-axis title)
		# unfortunately hiding tickmarks is for both X,Y axes
		p <- p + theme(axis.text.x=element_blank(),axis.ticks=element_blank())+scale_x_continuous(breaks=NULL)
	} else {
		p <- p + theme(axis.text.x=element_text(size=fontsize*0.8))+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)
	}
	if( is.null(axis.labels.y) ){
		# don't show y-axis or its labels or tickmarks (show y-axis title)
		p <- p + theme(axis.text.y=element_blank(),axis.ticks=element_blank())+scale_y_continuous(breaks=NULL)
	} else {
		p <- p + theme(axis.text.y=element_text(size=fontsize*0.8))+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)
	}

#	p <- p + coord_cartesian(ylim=axis.limits.y,xlim=axis.limits.x)+
	if( !is.null(axis.limits.x) ){ p <- p + xlim(axis.limits.x) }
	if( !is.null(axis.limits.y) ){ p <- p + ylim(axis.limits.y) }
	p <- p + scale_fill_continuous(name=legend.title)+
	scale_alpha_continuous(guide=FALSE)+
	theme(title=plot.title,plot.title=element_text(size=fontsize))

	if( !is.null(data.lines) ){
		rnames = names(data.lines)
		for(a_line_name in names(data.lines)){
			if( 'color' %in% names(a_line) ){
				p <- p +
				geom_line(data=a_line$data,aes(x=x,y=y),color=a_line$color)
			} else {
				p <- p +
				geom_line(data=a_line$data,aes(x=x,y=y))+scale_colour_brewer("clarity")
			}
		}
	}
	return(p)
}
# like histogram bars
plot_barplot <- function(
# the input data is just a number for each barplot, its height (and optionally errrobar height in another dataframe below)
# the valuesin df contains 3 cols. First col is 'cat1', second col is 'cat2', third col is 'value' and optional 4th col is 'errorbar'
# corresponding to the 'value'.
# cat1: is the category by which barplots are either stacked together or side-by-side, all those values with same 'cat1' value will be stacked together.
# cat2: is the category by which stacked barplots(by 'cat1') will be drawn along the x-axis.
# for example, cat1 is dataset name and cat2 is columnname in the dataset. This will stack data from each dataset and then
# stacked bars will be ploted along the x-axis for each columnname.
	# there are 3 columns in the dataframe of values.
	# each row represents a dataset and each column the same column in each dataset, the value is what we want to plot.
	# so, each bar in the plot will be of one column stacked or side-by-side over the many datasets 
	valuesin,
	stacked=FALSE, # stacked or side-by-side bars?
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2),
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	background=element_blank(),
	bar.color.outline=alpha('black',0.75),
	bar.color.alpha=0.75,
	fontsize=20,
	axis.label.x=NULL,
	axis.label.y=NULL,
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.x.angle=0, # or 90 for rotating 90 degrees
	axis.labels.x.fontsize=0.8*fontsize,
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y.angle=0,
	axis.labels.y.fontsize=0.8*fontsize,

	errorbars.show=TRUE,

	bar.height_labels.show=TRUE,
	bar.height_labels.fontsize_multiplier=0.25, # multiply fontsize by this parameter to get fontsize of height labels, may be you want them smaller
	plot.title=NULL,
	legend.title=NULL,
	
	legend.direction="vertical", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
	legend.position="right",
	legend.show=TRUE
){
	suppressMessages(require(ggplot2))
	p <- ggplot(data=valuesin, aes(x=cat2,y=value,fill=cat1))
	gbar_str="geom_bar(";
	
	if( stacked == T ){ gbar_str=paste(gbar_str, "position='fill'",sep='') } else { gbar_str=paste(gbar_str, "position=position_dodge()", sep='') }
	if( is.null(bar.color.outline) == F ){ gbar_str=paste(gbar_str, ",colour=bar.color.outline", sep='') }
	if( is.null(bar.color.alpha) == F ){  gbar_str=paste(gbar_str, ",alpha=bar.color.alpha", sep='') }
	gbar_str=paste(gbar_str, ")", sep='')
#write(paste("bin/library_in_R.R::plot_barplot : gbar_str='", gbar_str, "'", sep=''), stderr())
	p <- p + eval(parse(text=gbar_str))

	if( ('errorbar' %in% colnames(valuesin)) && (errorbars.show==TRUE) ){
		if( stacked == T ){
#			p <- p + geom_segment(aes(x=cat2,y=value,xend=cat2,yend=value+errorbar)) 
#		        geom_point(aes(x=cat2,y=value+errorbar),shape = "|",show_guide = FALSE)
		} else {
			p <- p + geom_errorbar(aes(ymin=value-errorbar,ymax=value+errorbar), width=0.2, position=position_dodge(.9))
		}
	}
                
	p <- p + theme(panel.background=background,panel.grid.major=gridline.major,panel.grid.minor=gridline.minor,legend.direction=legend.direction)+
	xlab(axis.label.x)+ylab(axis.label.y)+
	theme(axis.title.x=element_text(size=fontsize),axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5)) # magic numbers
	if( is.null(axis.labels.x) ){
		# don't show x-axis or its labels or tickmarks (show x-axis title)
		# unfortunately hiding tickmarks is for both X,Y axes
		p <- p + theme(axis.text.x=element_blank(),axis.ticks=element_blank())+scale_x_discrete(breaks=NULL)
	} else {
		p <- p + theme(axis.text.x=element_text(size=axis.labels.x.fontsize,angle=axis.labels.x.angle,hjust=1,vjust=1))+scale_x_discrete(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)
	}
	if( is.null(axis.labels.y) ){
		# don't show y-axis or its labels or tickmarks (show y-axis title)
		p <- p + theme(axis.text.y=element_blank(),axis.ticks=element_blank())+scale_y_continuous(breaks=NULL)
	} else {
		p <- p + theme(axis.text.y=element_text(size=axis.labels.y.fontsize,angle=axis.labels.y.angle,hjust=1))+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)
	}

	if( !is.null(axis.limits.x) ){ p <- p + xlim(axis.limits.x) }
	if( !is.null(axis.limits.y) ){ p <- p + ylim(axis.limits.y) }
	p <- p + scale_fill_discrete(name=legend.title,guide=legend.show)+
	theme(title=plot.title,plot.title=element_text(size=fontsize),legend.position=legend.position)
	return(p)
}

plot_scatterplot <- function(
	data, # the dataframe with 3 entries: x, y and class, class says which class the (x,y) pair belongs to for annotation
	      # OPTIONALLY there is a 4th entry: label for labelling each point with a text label - if this is missing no labels
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	class.labels=NULL, # data ontains 'class' col which may contain string labels or integers
			   # if it's integers, then supply this vector of labels to be placed in the legend
			   # if class.labels is not null, it overwrites the class cols
	class.colors.manual=NULL, # a manual color scale for each class in the data, NULL will calculate this scale itself using the min and max colors below
	class.colors.scale.begin='red',
	class.colors.scale.end='green',
	class.different_colors=TRUE, # FALSE (for no different colors), TRUE (for different colors, automatic setup) or a vector of colors for manual setup
	class.different_shapes=FALSE, # specify if each class is drawn with different shape
	points.size=2,
	points.alpha=0.75,
	background=element_blank(),
	fontsize=20,
	axis.label.x="x",
	axis.label.y="y",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	axis.log10scale.x=FALSE, # set axes to logscale (log_10)
	axis.log10scale.y=FALSE,
	axis.logscale.x=FALSE, # set axes to logscale (log_e)
	axis.logscale.y=FALSE,
	plot.title=NULL,

	legend.title=NULL,
	legend.fontsize=12,
	legend.show=TRUE,
	legend.direction="horizontal", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically

	correlation.method="", # calculate and print a correlation coefficient between x and y using this method (pearson or spearman)
	correlation.test.method=correlation.method,
	add.bestfit.line=FALSE, # add best-fit lines and curves to show the trend of the correlation
	add.bestfit.curve=FALSE, # add best-fit lines and curves to show the trend of the correlation
	render=TRUE,
	lines.horizontal=list(), # this takes 3 params to draw a horizontal line, e.g. (c(12,"red","dashed"),c(13,"blue","dashed"))
	lines.vertical=list(), # as above
	# if box.labels.fontsize > 0 then it draws a label/class next to each point (random pos around point) and a line joining them
	# this will remove the legends of class color etc.
	box.labels.fontsize=2.5, # fontsize of the labels if there is a column 'label' in data
	box.labels.linewidth=0.25, # the lines joining labels to points, if 0, then no lines drawn
	box.labels.leaway=0.5, # the labels/class are drawn around the point but at random pos, this is the stdev of this randomness, the higher the further away from point
	box.labels.alpha=0.8, # alpha of labels
	box.labels.linealpha=0.8 # alpha of lines joining the labels
#	seed=-1 # set the random generator seed so that labels are drawn at the same place all the time
){
	suppressMessages(require(ggplot2))
	all_classes = unique(as.vector(data$class)); num_classes = length(all_classes)
	if( is.null(class.labels) ){ all_labels = all_classes } else { all_labels = class.labels }
	if( correlation.method != "" ){
		corres = list(length=num_classes)
		for(i in 1:num_classes){
			aclass=all_classes[i]
			corres[i] = list(calculate_correlation(data$x[data$class==aclass], data$y[data$class==aclass], correlation.method=correlation.method, correlation.test.method=correlation.test.method))
			all_labels[i] = paste(all_labels[i], ", R: ", corres[[i]]$correlation, ", p-value:", corres[[i]]$p.value, sep="")
		}
	}
	has_text_labels = 'label' %in% colnames(data)

	if( class.different_colors == TRUE ){
		if( class.different_shapes == FALSE ){ p <- ggplot(data=data,aes(x=x,y=y,color=class)) }
		else { p <- ggplot(data=data,aes(x=x,y=y,color=class,shape=class)) }
	} else {
		if( class.different_shapes == FALSE ){ p <- ggplot(data=data,aes(x=x,y=y)) }
		else { p <- ggplot(data=data,aes(x=x,y=y,shape=class)) }
	}			
	if( class.different_colors == TRUE ){
		if( class.different_shapes == FALSE ){
			p <- p +
				geom_point(size=points.size,alpha=points.alpha)+
				scale_shape_manual(guide=F)
		} else {
			p <- p +
				geom_point(size=points.size,alpha=points.alpha)
		}
		if( has_text_labels == F ){
			p <- p +
				scale_shape_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title)
			if( num_classes > 1 ){
				if( is.null(class.colors.manual) && (num_classes>1) ){
#					p <- p + scale_color_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title,low=class.colors.scale.begin,high=class.colors.scale.end)
#					p <- p + scale_color_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title,legend=legend.show)
# this stupid shit seems to draw no colors when breaks=c(1:num_classes),labels=all_labels in any even labels are ok below so no problem
					p <- p + scale_color_discrete(name=legend.title)
				} else {
#					p <- p + scale_color_manual(breaks=c(1:num_classes),labels=all_labels,values=class.colors.manual,legend=legend.show)
					p <- p + scale_color_manual(values=class.colors.manual,guide=legend.show)
				}
			}
		}
	}
	p <- p +
	theme(panel.grid.major=gridline.major,panel.grid.minor=gridline.minor,panel.background=background, legend.direction=legend.direction)+
	theme(axis.text.x=element_text(size=fontsize*0.8)) + # these are magic numbers for just and pointsize
	theme(axis.text.y=element_text(size=fontsize*0.8,hjust=1))+
	theme(axis.title.x=element_text(size=fontsize))+
	theme(axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5))+
	coord_cartesian(ylim=axis.limits.y,xlim=axis.limits.x)+
	xlab(axis.label.x)+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)+
	ylab(axis.label.y)+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)+
	scale_alpha_continuous(guide=FALSE)+theme(legend.position="bottom",legend.text=element_text(size=legend.fontsize,hjust=0),title=plot.title,plot.title=element_text(size=fontsize))
	if( (add.bestfit.line == TRUE) && (add.bestfit.curve == FALSE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),method=lm,alpha=0.35,fill=alpha('#D55E00',0.105))
	} else if( (add.bestfit.line == FALSE) && (add.bestfit.curve == TRUE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),alpha=0.35,fill=alpha('#D55E00',0.105))
	} else if( (add.bestfit.line == TRUE) && (add.bestfit.curve == TRUE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),method=lm,alpha=0.35,fill=alpha('#D55E00',0.105),se=T)+
			 geom_smooth(data=data,aes(x=x,y=y,color=class),alpha=0.35,fill=alpha('#D55E00',0.105),se=F)
	}
	if( axis.logscale.x ){ p <- p + scale_x_log() }
	if( axis.logscale.y ){ p <- p + scale_y_log() }
	if( axis.log10scale.x ){ p <- p + scale_x_log10() }
	if( axis.log10scale.y ){ p <- p + scale_y_log10() }

	if( has_text_labels == T ){
		oldseed <- NULL; if( exists(".Random.seed") ){ oldseed <- get(".Random.seed",.GlobalEnv) }; .Random.seed = 1234 # set this seed fixed so that labels are drawn at the same place
		# if labels are to be drawn on each point, then find the position randomly around the point and also connect it with a line to the point
		N = nrow(data)
		data$randx = rnorm(n=N, mean=data$x, sd=box.labels.leaway)
		data$randy = rnorm(n=N, mean=data$y, sd=box.labels.leaway)
#p <- ggplot(data=data) + geom_text(aes(x=randx,y=randy,label=label,color=class),size=box.labels.fontsize,alpha=box.labels.alpha)
		p <- p + geom_text(data=data,aes(x=randx,y=randy,label=label,color=class),size=box.labels.fontsize,alpha=box.labels.alpha)
		if( box.labels.linewidth > 0.0 ){
			p <- p + geom_segment(data=data,aes(x=x,y=y,xend=randx,yend=randy,color=class),size=box.labels.linewidth,alpha=box.labels.linealpha)
		}
		if( !is.null(oldseed) ){ assign(".Random.seed", oldseed, .GlobalEnv) }
	}
	if( legend.show == F ){ p <- p + theme(legend.position = "none") } # what can i say? nothing else works

	# aes superhack for intercept lines
	if( !is.null(lines.horizontal) ){ for(al in lines.horizontal){ p <- p + geom_hline(eval(parse(text=paste('aes(yintercept=',al[1],')',sep=''))),colour=al[2],linetype=al[3]) } }
	if( !is.null(lines.vertical) ){ for(al in lines.vertical){ p <- p + geom_vline(eval(parse(text=paste('aes(xintercept=',al[1],')',sep=''))),colour=al[2],linetype=al[3]) } }
#	legend <- p + theme(keep="legend_box",title=plot.title,legend.text=element_text(size=legend.fontsize,hjust=0),legend.justification="center")
#	p <- p + theme(legend.position = "none") # remove the legend and put it underneath
	return(p)
}

old_plot_scatterplot <- function(
	data, # the dataframe with 3 entries: x, y and class, class says which class the (x,y) pair belongs to for annotation
	gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	class.labels=NULL, # data ontains 'class' col which may contain string labels or integers
			   # if it's integers, then supply this vector of labels to be placed in the legend
			   # if class.labels is not null, it overwrites the class cols
	class.colors.manual=NULL, # a manual color scale for each class in the data, NULL will calculate this scale itself using the min and max colors below
	class.colors.scale.begin='red',
	class.colors.scale.end='green',
	class.different_colors=TRUE, # FALSE (for no different colors), TRUE (for different colors, automatic setup) or a vector of colors for manual setup
	class.different_shapes=FALSE, # specify if each class is drawn with different shape
	points.size=2,
	points.alpha=0.75,
	background=element_blank(),
	fontsize=20,
	axis.label.x="x",
	axis.label.y="y",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	plot.title=NULL,
	legend.title=NULL,
	legend.fontsize=12,
	legend.show=TRUE,
	legend.direction="horizontal", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
	correlation.method="", # calculate and print a correlation coefficient
	correlation.test.method=correlation.method,
	add.bestfit.line=FALSE, # add best-fit lines and curves to show the trend
	add.bestfit.curve=FALSE, # add best-fit lines and curves to show the trend
	render=TRUE
){
	suppressMessages(require(ggplot2))
	all_classes = unique(as.vector(data$class)); num_classes = length(all_classes)
	if( is.null(class.labels) ){ all_labels = all_classes } else { all_labels = class.labels }
	if( correlation.method != "" ){
		corres = list(length=num_classes)
		for(i in 1:num_classes){
			aclass=all_classes[i]
			corres[i] = list(calculate_correlation(data$x[data$class==aclass], data$y[data$class==aclass], correlation.method=correlation.method, correlation.test.method=correlation.test.method))
			all_labels[i] = paste(all_labels[i], ", R: ", corres[[i]]$correlation, ", p-value:", corres[[i]]$p.value, sep="")
		}
	}

	if( class.different_colors == TRUE ){
		if( class.different_shapes == FALSE ){ p <- ggplot(data=data,aes(x=x,y=y,color=class)) }
		else { p <- ggplot(data=data,aes(x=x,y=y,color=class,shape=class)) }
	} else {
		if( class.different_shapes == FALSE ){ p <- ggplot(data=data,aes(x=x,y=y)) }
		else { p <- ggplot(data=data,aes(x=x,y=y,shape=class)) }
	}			
	if( class.different_colors == TRUE ){
		if( class.different_shapes == FALSE ){ p <- p + geom_point(size=points.size,alpha=points.alpha)+scale_shape_manual(guide=F) }
		else { p <- p + geom_point(size=points.size,alpha=points.alpha)+scale_shape_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title) }
		if( is.null(class.colors.manual) ){
#			p <- p + scale_color_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title,low=class.colors.scale.begin,high=class.colors.scale.end)
			if( num_classes > 1 ){ p <- p + scale_color_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title) }
		} else {
			p <- p + scale_color_manual(breaks=c(1:num_classes),labels=all_labels,values=class.colors.manual,guide=legend.show)
		}
	} else {
		if( class.different_shapes == FALSE ){ p <- p + geom_point(size=points.size,alpha=points.alpha)+scale_shape_manual(guide=F) }
		else {
			p <- p + geom_point(size=points.size,alpha=points.alpha)+scale_shape_discrete(breaks=c(1:num_classes),labels=all_labels,name=legend.title)
		}
	}
	p <- p +
	theme(panel.grid.major=gridline.major,panel.grid.minor=gridline.minor,panel.background=background)+
	theme(axis.text.x=element_text(size=fontsize*0.8)) + # these are magic numbers for just and pointsize
	theme(axis.text.y=element_text(size=fontsize*0.8,hjust=1))+
	theme(axis.title.x=element_text(size=fontsize))+
	theme(axis.title.y=element_text(size=fontsize,angle=90,hjust=0.5))+
	coord_cartesian(ylim=axis.limits.y,xlim=axis.limits.x)+
	xlab(axis.label.x)+scale_x_continuous(breaks=axis.labels.x$breaks,labels=axis.labels.x$labels)+
	ylab(axis.label.y)+scale_y_continuous(breaks=axis.labels.y$breaks,labels=axis.labels.y$labels)+
	theme(legend.position="bottom",legend.text=element_text(size=legend.fontsize,hjust=0),title=plot.title)
	if( (add.bestfit.line == TRUE) && (add.bestfit.curve == FALSE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),method=lm,alpha=0.35,fill=alpha('#D55E00',0.105))
	} else if( (add.bestfit.line == FALSE) && (add.bestfit.curve == TRUE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),alpha=0.35,fill=alpha('#D55E00',0.105))
	} else if( (add.bestfit.line == TRUE) && (add.bestfit.curve == TRUE) ){
		p <- p + geom_smooth(data=data,aes(x=x,y=y,color=class),method=lm,alpha=0.35,fill=alpha('#D55E00',0.105),se=T)+
			 geom_smooth(data=data,aes(x=x,y=y,color=class),alpha=0.35,fill=alpha('#D55E00',0.105),se=F)
	}
#	legend <- p + theme(keep="legend_box",title=plot.title,legend.text=element_text(size=legend.fontsize,hjust=0),legend.justification="center")
#	p <- p + theme(legend.position = "none") # remove the legend and put it underneath
	return(p)
}

# plot multiple y-values with a common x-value each as pair as a point, e.g. x=1, y1=1,y2=2, or just x=1,y=1, etc.
plot_scatterplot_and_density_of_many_groups <- function(
	data, # dataframe contains x, y and class - class could be filename and x:col1, y:col2
	class.labels=NULL, # data contains 'class' col which may contain string labels or integers
			   # if it's integers, then supply this vector of labels to be placed in the legend
			   # if class.labels is not null, it overwrites the class cols
	scatterplot.background=element_blank(),
	scatterplot.gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	scatterplot.gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	density_top.background=element_blank(),
	density_top.gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	density_top.gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),
	density_right.background=element_blank(),
	density_right.gridline.major=element_line(colour=alpha('blue',0.15),size=0.2), # grid color and alpha
	density_right.gridline.minor=element_line(colour=alpha('purple',0.05),size=0.2),

	class.different_colors=TRUE, # FALSE (for no different colors), TRUE (for different colors, automatic setup) or a vector of colors for manual setup
	class.different_shapes=FALSE, # specify if each class is drawn with different shape
	class.colors.scale.begin='red',
	class.colors.scale.end='green',
	scatterplot.points.size=1.5,
	scatterplot.points.alpha=0.35,
	density_right.points.size=0.5,
	density_right.points.alpha=1.0,
	density_top.points.size=0.5,
	density_top.points.alpha=1.0,
	fontsize=20,
	axis.label.x="x",
	axis.label.y="y",
	axis.limits.x=NULL,
	axis.limits.y=NULL,
	axis.labels.x=data.frame(breaks=NULL,labels=NULL),
	axis.labels.y=data.frame(breaks=NULL,labels=NULL),
	plot.title=NULL,
	legend.title=NULL,
	legend.fontsize=12,
	legend.show=TRUE,
	legend.direction="horizontal", # horizontal or vertical legend (the gradient) to the colors, this is for space-management basically
	correlation.method="", # calculate and print a correlation coefficient
	correlation.test.method=correlation.method,
	add.bestfit.line=FALSE, # add best-fit lines and curves to show the trend
	add.bestfit.curve=FALSE, # add best-fit lines and curves to show the trend
	render=TRUE
){
	suppressMessages(require(ggplot2))
	p <- plot_scatterplot(
		data=data,
		axis.label.x=axis.label.x,
		axis.label.y=axis.label.y,
		axis.labels.x=axis.labels.x,
		axis.labels.y=axis.labels.y,
		axis.limits.x=axis.limits.x,
		axis.limits.y=axis.limits.y,
		class.labels=class.labels,

		class.colors.scale.begin=class.colors.scale.begin,
		class.colors.scale.end=class.colors.scale.end,
		class.different_colors=class.different_colors,

		class.different_shapes=class.different_shapes,
		background=scatterplot.background,gridline.major=scatterplot.gridline.major,gridline.minor=scatterplot.gridline.minor,
		points.size=scatterplot.points.size,points.alpha=scatterplot.points.alpha,
		fontsize=fontsize,

		correlation.method=correlation.method,
		correlation.test.method=correlation.method,
		add.bestfit.line=add.bestfit.line, # add best-fit lines and curves to show the trend
		add.bestfit.curve=add.bestfit.curve, # add best-fit lines and curves to show the trend
	)
#	a_legend <- p + theme(keep="legend_box",legend.text=element_text(size=legend.fontsize,hjust=0),legend.direction=legend.direction)
# opts(keep=".." is deprecated so use this hack:
	a_legend <- extract_legend_from_plot(p)

	p1 <- p + theme(legend.position = "none") # remove legend
	suppressMessages(library(KernSmooth))
	# if the number of bins is too small - OR MOST LIKELY data has no variation,
	# bkde will fail (Error in if (L == 0) warning("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'") :   missing value where TRUE/FALSE needed
	dat <- data[data[,'class']==1,'x']; rang <- range(dat)
	if( rang[1] == rang[2] ){ write(paste("plot_scatterplot_and_density_of_many_groups : no variation for x-col :  ", rang[1], " to ", rang[2], sep=''), stderr()); return(NULL); }
	denF1x <- tryCatch(bkde(dat, range.x=rang), error=function(e){write("plot_scatterplot_and_density_of_many_groups : too little bins or data has no variation (x-col)", stderr()); return(NULL); })

	dat <- data[data[,'class']==1,'y']; rang <- range(dat);
	if( rang[1] == rang[2] ){ write(paste("plot_scatterplot_and_density_of_many_groups : no variation for y-col :  ", rang[1], " to ", rang[2], sep=''), stderr()); return(NULL); }
	denF1y <- tryCatch(bkde(dat, range.x=rang), error=function(e){write("plot_scatterplot_and_density_of_many_groups : too little bins or data has no variation (y-col)", stderr()); return(NULL); })

	if( 2 %in% data[,'class'] ){
		dat <- data[data[,'class']==2,'x']; rang <- range(dat); denF2x <- bkde(dat, range.x=rang)
		dat <- data[data[,'class']==2,'y']; rang <- range(dat); denF2y <- bkde(dat, range.x=rang);
		m=max(denF1x$y,denF2x$y); s=seq(0,m,length=4); density_top.ylabels = data.frame(breaks=s,labels=c(0,format(s[-1],scientific=T,digits=2)))
		m=max(denF1y$y,denF2y$y); s=seq(0,m,length=3); density_right.ylabels = data.frame(breaks=s,labels=c(0,format(s[-1],scientific=T,digits=2)))
		l1 <- list(data.frame(x=denF1x$x, y=denF1x$y), data.frame(x=denF2x$x,y=denF2x$y))
		l2 <- list(data.frame(x=denF1y$x, y=denF1y$y), data.frame(x=denF2y$x,y=denF2y$y))
	} else {
		denF2x <- NULL; denF2y <- NULL
		m=max(denF1x$y); s=seq(0,m,length=4); density_top.ylabels = data.frame(breaks=s,labels=c(0,format(s[-1],scientific=T,digits=2)))
		m=max(denF1y$y); s=seq(0,m,length=3); density_right.ylabels = data.frame(breaks=s,labels=c(0,format(s[-1],scientific=T,digits=2)))
		l1 <- list(data.frame(x=denF1x$x, y=denF1x$y))
		l2 <- list(data.frame(x=denF1y$x, y=denF1y$y))
	}
	p2 <- plot_points(
		data=l1,
		legend.show=FALSE,
		axis.label.x='',
		axis.label.y='',
		axis.labels.x=axis.labels.x,
		axis.labels.y=density_top.ylabels,
		axis.limits.x=axis.limits.x,
		class.different_colors=c('red','green'),
		background=density_top.background,gridline.major=density_top.gridline.major,gridline.minor=density_top.gridline.minor,
		points.size=density_top.points.size,points.alpha=density_top.points.alpha,
		fontsize=fontsize
	)
	p3 <- plot_points(
		data=l2,
		legend.show=FALSE,
		axis.label.x='',
		axis.label.y='',
		axis.labels.x=axis.labels.x,
		axis.labels.y=density_right.ylabels,
		axis.limits.x=axis.limits.x,
		class.different_colors=c('red','green'),
		background=density_right.background,gridline.major=density_right.gridline.major,gridline.minor=density_right.gridline.minor,
		points.size=density_right.points.size,points.alpha=density_right.points.alpha,
		fontsize=fontsize
	); p3 <- p3 + coord_flip()
	if( render == TRUE ){
		if( is.null(plot.title) ){ 
			multiplot(ncol=2,nrow=3,widths=c(7,3),heights=c(3,7,2.5),list(p1,2,1), list(p2,1,1), list(p3,2,2), legend=list(a_legend,3,1:2))
		} else {
			multiplot(ncol=2,nrow=4,widths=c(7,3),heights=c(1.5,3,7,2.5),legend=list(a_legend,4,1),list(p1,3,1), list(p2,2,1), list(p3,3,2))
			grid.text(plot.title, vp=viewport(layout.pos.row=1,layout.pos.col=1:2),just=c('centre','top'))
		}
	}
	return(list(scatterplot=p1,density_top=p2,density_right=p3,legend=a_legend))
}
multiplot <- function(
	ncol,	# number of columns
	nrow,	# number of rows
	widths=NULL,	# widths as a vector of widths for each column, e.g. c(1,2)
	heights=NULL,# hights for each row e.g. c(1,2)
	legend=NULL, # a list (like the ones below) which contains a legend element (can extract from a ggplot like: legend <- p + theme(keep = "legend_box",title='title')
	... # one or more lists: e.g. list(p1,1,2),list(p2,1,3) **** first is ROW, second is COL *** each with the ggplot object and x,y coords to be placed on the (nrow,ncol) matrix
) {
	# plot a scatter plot or heatmap with marginal histograms/densities on bottom, right etc. leave null for not having one
	source(paste(Sys.getenv("HOME"),"/PROJECTS/ICR/bin/library_align.R",sep=''))
	suppressMessages(library(grid))
	if( is.null(widths) ){ widths = rep(1,times=ncol) }
	if( is.null(heights) ){ heights = rep(1,times=nrow) }
	grid_layout <- grid.layout(nrow=nrow, ncol=ncol, widths=widths, heights=heights)
	grid.newpage()
	pushViewport( viewport( layout=grid_layout, width=1, height=1 ) )
	align.plots(grid_layout, ...)
	if( !is.null(legend) ){ print(legend[[1]], vp=viewport(layout.pos.row=legend[[2]],layout.pos.col=legend[[3]])) }
}

# a ggplot plot title formatter use it as:
# theme(title=title, plot.title=format_title(width=unit(20, "cm")))
# it produces the formatting basically
ggplot_string_formatter <- function (width=unit(1, "npc"), family = "", face = "plain", colour ="black", size = 10,hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 1.1)
{
	vj <- vjust
	hj <- hjust
	structure(function(label, x = hj, y = vj, ..., vjust = vj,
	hjust = hj, default.units = "npc") {
		textGrob(splitString(label, width), x, y, hjust = hjust, vjust= vjust, ...,default.units = default.units, gp = gpar(fontsize = size,col = colour, fontfamily = family, fontface = face,lineheight = lineheight), rot = angle)
	}, class = "theme", type = "title", call = match.call())
}

# will calculate correlation (pearson or spearman), do a test to get a p-value for this correlation and
# best-fit a line to the data (will return intercept and slope)
calculate_correlation <- function(
	x,
	y,
	correlation.method="pearson",
	correlation.test.method="pearson"
){
	# correlation value
	corellis = cor(x, y, method=correlation.method)
	# test of significance
	signi = cor.test(x, y, method=correlation.test.method)

	# fit a line to the two data
	corelifit = lm(x~y)
	return(list(
		correlation=corellis,
		p.value=signi$p.value,
		lm.intercept=corelifit$coefficients[1],
		lm.slope=corelifit$coefficients[2]
	))
	# the line can be plot using abline(coef=c(intercept,slope)) etc.
}
# opts(keep="legend_box" ... is deprecated
# use this to get the legend of any plot
extract_legend_from_plot <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

