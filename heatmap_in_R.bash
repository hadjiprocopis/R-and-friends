#!/bin/bash

set -o pipefail
set -e

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

appname=`basename "$0"`
bindir=$( cd "$( dirname "$0" )" && pwd )
now=`date`
rand=$$
png_width=1024
png_height=768
png_fontsize=18
axis_labels_x_angle=45
axis_labels_y_angle=45
extra_options=""
print_values=0
value_to_replace_NA_with=0

function display_help_exit {
	echo "Usage : $1 ... "
	echo "This script will read a 2D data from input_file and make a heat map (no dendrograms or anything). The format of the input data: first row and first column contain element names, e.g.:"
	echo "XXX A B C"
	echo "X 1 2 3"
	echo "Y 3 4 5"
	echo "Z 10 11 12"
	echo "it will then plot the value of each row->col (e.g.X->A) as a color shade"
	echo "NOTE : if there are NA values in the input they will be replaced by default value which is '${value_to_replace_NA_with}'. This value can be changed via the -N option."
	echo " -i input_file"
	echo "        defines the input file which each column represents a difference and the first column is the gene name"
	echo " -o output_basename"
	echo "		where the results go as '.png' file and '.dentrogram.png' (with dentrograms) file"
	echo "-3 png_width:png_height:png_fontsize : output png options, including the fontsize of all labels etc."
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file, default is white space"
	echo "        which is '' (i.e. nothing in between the quotes)"
	echo " -N value_to_replace_NA_with : specify a value to replace all instances of NA in the input file."
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs)"
	echo ""
	echo " -x xlabel : x-axis label (name), default is empty."
	echo " -y ylabel : y-axis label (name), default is empty"
	echo " -t title : default is empty"
	echo " -r : do not order rows and columns using hierarchical clustering (like R's heatmap function does). By default sorting rows and columns is ON so that similar values are clustered together."
	echo " -s : data is symmetric, only if -r is not used, and i don't think it matters ..."
	echo " -m values_cutoff_min : all data values below this will have a smaller alpha"
	echo " -M values_cutoff_max : all data values greater than this will have a larger alpha"
	echo " -X angle : the labels on the x-axis will be rotated by this amount, default is 45 (degrees) but usually it is either 0 (horizontal) or 90 (vertical)"
	echo " -Y angle : same, default is 45 (degrees) etc."
	echo " -F printf_format : the values will be printed in the middle of each cell using this format, e.g. '%.2f' or 'v=%.2f'"
	echo " -b color : draw a border around each cell of this color (can be 'black' or 'alpha("black", 0.2)' in this case check quotes)"
	echo " -c color : the cells will have a gradient of colors STARTING from this color."
	echo " -C color : the cells will have a gradient of colors ENDING at this color."
	echo " -g color : if your data is around a middle value, e.g. -1 to 0 to +1, then this color will be mapped to the middle value (which you set using -G) - this is not always necessary, unless you really want that the middle value has this exact color."
	echo " -G value : if you specify a middle color, then you may want to specify a middle value to map this color to."
	echo " -S min:max : scale the input to be between min and max."
	echo " -L fontsize : box labels fontsize, set it to zero if you don't want any labels to be printed in the boxes."
	echo "	      Example:"
	echo "${appname} -i file -o out -3 2048:1436:18 -t fufufu -x aa -y bb -X 90 -Y 0 -c white -C red -g blue"
	echo "program by Andreas Hadjiprocopis"
	echo "andreashad2 at gmail.com"
	echo "Institute of Cancer Research"
	echo "London, 2012"
	echo "Free for everyone to copy, use, modify / quote author"
	echo ""
	exit 1
}

inputFile=""
outputBasename=""
separator=""
separatorOut="\t"
appname="$0"
rand=$$
Rscript_name="${rand}.R"
delete_Rscript_when_done=1
order_using_R_heatmap=1
heatmap_symmetric="FALSE"
scale_from="NULL"
scale_to="NULL"
heatmapplus_title="" # for the heatmap() only

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

while getopts "i:o:R:Hm:M:x:y:l:t:X:Y:F:b:c:C:g:G:rsS:L:e:3:PN:" OPTION; do
	case $OPTION in
		i)
			inputFile="$OPTARG"
			;;
		3)
			png_width=`echo "${OPTARG}" | cut -d':' -f1`
			png_height=`echo "${OPTARG}" | cut -d':' -f2`
			png_fontsize=`echo "${OPTARG}" | cut -d':' -f3`
			;;
		o)
			outputBasename="$OPTARG"
			;;
		e)
			separator=",sep='$OPTARG'"
			;;
		N)
			value_to_replace_NA_with="$OPTARG"
			;;
		R)
			delete_Rscript_when_done=0
			Rscript_name="$OPTARG"
			;;
		H)
			display_help_exit "$appname"
			;;
# other options
		r)
			order_using_R_heatmap=0
			;;
		s)
			heatmap_symmetric="TRUE"
			;;
		S)
			scale_from=`echo "${OPTARG}" | cut -d':' -f1`
			scale_to=`echo "${OPTARG}" | cut -d':' -f2`
			;;
		P)
			print_values=1
			;;
		m)
			values_cutoff_min=`echo "${OPTARG}" | cut -d':' -f1`
			values_cutoff_min_replacement=`echo "${OPTARG}" | cut -d':' -f2`
			;;
		M)
			values_cutoff_max=`echo "${OPTARG}" | cut -d':' -f1`
			values_cutoff_max_replacement=`echo "${OPTARG}" | cut -d':' -f2`
			;;
		x)
			extra_options+=",axis.label.x='$OPTARG'"
			;;
		y)
			extra_options+=",axis.label.y='$OPTARG'"
			;;
		l)
			extra_options+=",box.labels.midpoint='$OPTARG'"
			;;
		L)
			extra_options+=",box.labels.fontsize=$OPTARG"
			;;
		t)
			extra_options+=",plot.title='$OPTARG'"
			heatmapplus_title=",main='$OPTARG'"
			;;
		X)
			axis_labels_x_angle="$OPTARG"
			;;
		Y)
			axis_labels_y_angle="$OPTARG"
			;;
		F)
			extra_options+=",box.labels.sprintfFormat='$OPTARG'"
			;;
		b)
			extra_options+=",cell.color='$OPTARG'"
			;;
		c)
			extra_options+=",box.labels.palette.start='$OPTARG'"
			;;
		C)
			extra_options+=",box.labels.palette.end='$OPTARG'"
			;;
		g)
# this is useful if yur valeus are centred around a number (e.g. -1 0 1), then you specify the color for middle as well
# as extremes and also specify middle value with -G
			extra_options+=",box.labels.palette.middle='$OPTARG'"
			;;
		G)
			extra_options+=",box.labels.middle=$OPTARG"
			;;
		*)
			echo "${appname} : unrecognised option $OPTION"
			display_help_and_exit
			;;
	esac
done

if [ ! -f "$inputFile" ]; then
	echo "$0 : input file $inputFile does not exist."
	exit 1
fi
if [ "$Rscript_name" = "" ]; then
	echo "$0 : script file name is not defined (use the -R option to do it)."
	exit 1
fi
if [ "$outputBasename" = "" ]; then
	echo "$0 : you need to specify an output file using the '-o' option."
	exit 1
fi

cat<<EOC>"$Rscript_name"
source('${bindir}/library_in_R.R')

# R script by $appname on $now
write("$appname : input file='$inputFile', output file='$outputBasename', R file='$Rscript_name'", append=FALSE, stderr());

write("$appname : reading input file $inputFile", stderr())
# keep column 'Probe Set' as char and do not replace spaces or hyphens in the header column names
datatmp <- tryCatch({
	# returns the last thing got
	as.matrix(read.table("$inputFile" $separator, row.names=1, header=T, as.is=1, check.names=FALSE))
	}, error=function(ex){
		write(paste("$appname : error reading file $inputFile\n",ex,sep=""),stderr())
		return(NULL) # or null if error
	}
); if( is.null(datatmp) ){ quit(status=1) }

rnames = rownames(datatmp); cnames = colnames(datatmp); nR = nrow(datatmp); nC = ncol(datatmp)
write(paste("$appname : done, read ", nC, " columns and ", nR, " rows from input file '${inputFile}'.", sep=''), stderr())
NAspresent = datatmp[is.na(datatmp)]
if( length(NAspresent) > 0 ){
	value_to_replace_NA_with = ${value_to_replace_NA_with}
	write(paste("$appname : there are ", length(NAspresent), " NA values in the input. They will be replaced by value '", value_to_replace_NA_with, "'.", sep=""), stderr())
	datatmp[is.na(datatmp)] = value_to_replace_NA_with
}
minv = min(datatmp)
maxv = max(datatmp)
if( minv == maxv ){
	write("$appname : there is no variation in the input file '${inputFile}', will disable ordering of heatmap as it will complain. I will not do any scaling either...", stderr())
} else {
	# dont do scaling or heatmap ordering if no variation in data (above)
	if( !is.null($scale_from) && !is.null($scale_to) ){
		sf=$scale_from; st=$scale_to;
		write(paste("$appname : scaling input data between $scale_from and $scale_to (currently min: ", minv, ", max: ", maxv, ") ...", sep=''), stderr())
		l = (st-sf) / (maxv-minv)
		datatmp = l * (datatmp - sf) -minv
		write("$appname : done scaling.", stderr())
	}
	# make a heatmap and then use the ordering of the heatmap to the data
	if( $order_using_R_heatmap == 1 ){
		suppressMessages(library(heatmap.plus))
		png(filename="$outputBasename.dentrogram.png", width=${png_width}, height=${png_height}, pointsize=${png_fontsize})
		aheat <- heatmap.plus(
			datatmp # no comma here
			${heatmapplus_title},
			symm=${heatmap_symmetric}
		) # this plots something, dont know how to turn it off all we need is the ordering.
		dev.off()
		newdata = matrix(nrow=nR, ncol=nC)
		for(i in 1:nR){
			for(j in 1:nC){
				newdata[i,j] = datatmp[aheat\$rowInd[i],aheat\$colInd[j]]
			}
		}
		colnames(newdata) <- cnames[aheat\$colInd]
		rownames(newdata) <- rnames[aheat\$rowInd]
		cnames = colnames(newdata); rnames = rownames(newdata);
		datatmp = newdata; newdata = NULL
	}
} # no variation

have_cutoff=0
cutoff_min=NA; cutoff_min_replacement=NA;
cutoff_max=NA; cutoff_max_replacement=NA;
if( "$values_cutoff_min" != "" ){ cutoff_min = $values_cutoff_min+0; cutoff_min_replacement = $values_cutoff_min_replacement+0; have_cutoff=1 }
if( "$values_cutoff_max" != "" ){ cutoff_max = $values_cutoff_max+0; cutoff_max_replacement = $values_cutoff_max_replacement+0; have_cutoff=1 }

# now put it into heatmap format dataframe of 'x', 'y', 'values' and optional 'frame_color', 'alpha', 'extra_values'
num_outRows = nR * nR # the input is a matrix of rows, the total items is nR * nC
datain_m = matrix(ncol=3, nrow=nR*nC)
arow = 1; for(r in 1:nR){
	for(c in 1:nC){
		datain_m[arow, 1] = c # x
		datain_m[arow, 2] = r # y
		datain_m[arow, 3] = datatmp[r,c]
		arow = arow + 1
	}
}; datatmp = NULL
datain = data.frame(datain_m, stringsAsFactors=F); colnames(datain) = c('x', 'y', 'values')

if( ${print_values} == 1 ){
	for(arow in 1:num_outRows){
		r = datain[arow, 'x']; c = datain[arow, 'y']; v = datain[arow, 'values']
		write(paste(rnames[r], ' :: ', cnames[c], ' = ', v, sep=''), stderr())
	}
}

if( have_cutoff == 1 ){
	if( is.na(cutoff_min) ) cutoff_min = min(datain_m[,3])
	if( is.na(cutoff_max) ) cutoff_max = max(datain_m[,3])
	datain[,'alpha'] = rep(2, num_outRows)
	for(i in 1:num_outRows){
		if( datain_m[i, 3] < cutoff_min ){ datain[i, 'alpha'] = 0 }
		else if( datain_m[i, 3] > cutoff_max ){ datain[i, 'alpha'] = 1 }
	}
}
datain_m = NULL;
p <- plot_heatmap(
	datain = datain,
	axis.labels.x = data.frame(breaks=rep(1:nC), labels=cnames),
	axis.labels.y = data.frame(breaks=rep(1:nR), labels=rnames),
	axis.labels.x.angle = ${axis_labels_x_angle},
	axis.labels.y.angle = ${axis_labels_y_angle},
	legend.show=F, legend.direction='vertical',
	fontsize=${png_fontsize} # NO COMMA HERE and leave it last before the \${} options
	${extra_options}
)
if( is.null(p) ){ write("${appname} : call to plot_heatmap has failed.", stderr()); quit(status=1) }
png(filename="$outputBasename.png", width=${png_width}, height=${png_height}, pointsize=${png_fontsize})
print(p)
dev.off()

write(paste("$appname : heatmap image written in '$outputBasename.png'", sep=""), stderr())
EOC

echo "$appname : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
	echo "$appname : error executing R script ($Rscript_name)";
#	cat /tmp/rlog$rand.txt; echo "$appname : error executing R script ($Rscript_name)";
	echo "Command line:"
	echo "$0 ${*}"
	exit 1
fi
echo "Command line:"
echo "$0 ${*}"
rm -f /tmp/rlog$rand.txt /tmp/$rand.*
if [ "${delete_Rscript_when_done}" -eq 1 ]; then rm -f "$Rscript_name"; fi
exit 0

