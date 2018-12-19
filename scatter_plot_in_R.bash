#!/bin/bash

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

appname=`basename "$0"`
bindir=$( cd "$( dirname "$0" )" && pwd )
now=`date`
rand=$$
_internal_arrays_ifs_="_::_" # stupid bash arrays!!!
param_ifs=':'
defaultInputFilesIFS='\t'
defaultInputFilesHaveHeader='FALSE'
seed=$$
png_width=1024
png_height=768
png_fontsize=18
num_density_points=1024
num_threads=1

alpha="1.0"

function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file [-c col_num ... -c ... | -C col_num ... -C ...] [-z col_name ... -z ... | -Z col_name ... -Z ...] [-1] [-e input_data_separator] [-E output_data_separator] [-R R_output_script_name] [-2 seed] [-3 png_width:png_height:png_fontsize] [-p num_threads]"
	echo "other params: "
	echo "This script will do a scatter plot of the first column versus the second column given labelling the points with optional third column for each input file."
	echo " -i inputFile"
	echo "	defines the input file."
	echo " -o outputBasename"
	echo "	defines the output basename. The column number and '.png' '.info.txt' '.labels.txt' '.discrete.txt' may be output files (the first for the images, for each column), the second for histogram info (bins, mids, etc) per column, and the other two may be produced if -I or -L are specified, see below."
	echo " -c col_num"
	echo "	choose the column number you wish to operate as a comma-separated list of column numbers (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
	echo " -C col_num"
	echo "	choose the column number you DO NOT wish to operate as a comma-separated list of column numbers (starting from 1) - repeate -C options if necessary, all columns will be processed unless this option is specified"
	echo " -z col_name"
	echo "	choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (define with option -1), default is '${param_ifs}' or repeat -z options"
	echo " -Z col_name"
	echo "	choose the column name you DO NOT wish to operate as a list of column names or repeate -Z options. Use unique separator (define with option -1), default is '${param_ifs}' or repeat -Z options"
	echo " -1"
	echo "	include columns not selected for processing to the output, default is not to print un-processed columns"
	echo " -2 random_number_generator_seed"
	echo "	useful to repeat experiments, otherwise seed is selected at random using the process id"
	echo " -3 png_width:png_height:png_fontsize : the size of the output image if any in pixels and font size, default: ${png_width}:${png_height}:${png_fontsize}"
	echo " -5 : do not show data legends"
# standard features
	echo " -p num_threads"
	echo "  if there are mutliple data to be processed, and num_threads > 1, then this processing will be done in parallel. Default num_threads is ${num_threads}."
	echo " -h"
	echo "	optional, that the first row is not numerical but it holds column names, these names will be copied to the output. If the script fails, this could be the reason (i.e. calculations with column names)."
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file. If the script fails, this could be the reason (i.e. spaces are used as col separator but in actual terms they are tab separated)."
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs. If this option is not used, then the output script name will be deleted."
	echo "	Otherwise, the R script will not be deleted and you can  re-run it again using R < script_name"
	echo ""
	echo "Other options:"
	echo "other params: [-T title] [-X xlabel] [-Y ylabel] [-x min:max] [-y min:max] [-A alpha] [-w function] [-W function]"
	echo "  -A alpha : alpha value for drawing the boxes, see also -j, -J, default alpha value is ${alpha} (from 0 to 1)"
	echo "  -T title : title (%cnumber% replaces column number, %cname% replaces column name)"
	echo "  -x min:max : the range of values of x-axis, useful for drawing histograms on the same scale for comparisons or in a movie"
	echo "  -y min:max : same for y-axis"
	echo "  -X label : the label for the x-axis, default is the column header name if a header exists or just 'col:N' for the Nth column"
	echo "  -Y label : the label for the y-axis, default depends on the histogram mode, i.e. density, counts, probability etc."
	echo "  -w function : specify an R function to do the pre-processing of the x-column, e.g. '{x=log(3);return(log10(v/x))}', where v is the variable with the data."
	echo "  -W function : specify an R function for the y-column."
	echo "  -d 'intercept,"color","linetype" : draw a vertical line at specified intercept with specified color and linetype. The last 2 must be quoted and understood by ggplot2. The intercept is on scale AFTER any trasnforms of the input wiht -w/-W. Example -d '12,\"red\",\"dashed\"'."
	echo "  -D 'intercept,"color","linetype" : same for the Y-axis - i.e. horizontal line"
	echo "For a VOLCANO plot, do : -W '{return(-log10(v))}' which plots x against -log10(y) assuming y is a p-value and x is a fold change, yeah whatever."
	echo "  -l fontsize : put a label around each point with its label/class - you must then supply 3 columns for each file (x,y) and the optional label. The fontsize is the size for the font to draw the label."
	echo ""
	echo "EXAMPLE:"
	echo "program by Andreas Hadjiprocopis"
	echo "andreashad2 at gmail.com"
	echo "Institute of Cancer Research"
	echo "London, 2010"
	echo "Free for everyone to copy, use, modify / quote author"
	echo ""
	exit 1
}
function column_names_to_column_numbers {
	local   _inputfile="$1"; shift # must have header
	local   _ifs="$1"; shift
	local   _column_names=${@}
	local   _rand=$$
	local   _i _ret="" _ifs2=""
	if [ "${_ifs}" != "" ]; then _ifs2="--ifs '${_ifs}'"; fi
	# dont echo anything apart from the last one before return (unless > /dev/stderr)
#	echo extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit
cat << EOC | sh
	extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit >& /dev/null
EOC
	if [ $? -ne 0 ]; then
		echo "$appname : ${FUNCNAME}, line ${LINENO} : call to the following command has failed"
		echo extract_columns.pl --input \"${_inputfile}\" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit
		return 1
	fi
	oldIFS="${IFS}"; IFS=${param_ifs}
	for _a_col_name in ${_column_names[*]}; do
#	      echo "checking for ${_a_col_name} (/tmp/${_rand}.cn2cn)" > /dev/stderr
		_a_col=`awk -F '\t' -v n="${_a_col_name}" '{if($4==n){printf ($2+1)"'"${param_ifs}"'"}}' /tmp/${_rand}.cn2cn`
		if [ "${_a_col}" == "" ]; then
#			echo "$appname : ${FUNCNAME}, line ${LINENO} : could not find column name '${_a_col_name}' in the header of the input file '${_inputfile}'."
#			return 1
			continue # baby
		fi
#	       echo "found ${_a_col}" > /dev/stderr
		_ret+="${param_ifs}${_a_col%${param_ifs}}"
	done
	IFS="${oldIFS}"
	echo "${_ret#${param_ifs}}"
	rm -f /tmp/${_rand}.cn2cn
	return 0
}
inputFilesStr=""
inputFilesLabelsStr=""
inputFilesHaveHeaderStr=""
inputFilesIFSStr=""; currentInputFileIFS=${defaultInputFilesIFS}
doColumnsStr=""
dontDoColumnsStr=""
_doColumnsStr=""
_dontDoColumnsStr=""
outputFile=""
# for standard arguments
separatorOut='\t'
header="FALSE"
Rscript_name="${rand}.R"
outputPrecision="5"
include_discarded_columns_in_output=0
delete_Rscript_when_done=1
x_log10_scale='FALSE'; y_log10_scale='FALSE'
x_log_scale='FALSE'; y_log_scale='FALSE'
limitsX=""; limitsY=""; xlabel=""; ylabel=""; title=""; legend_title=""
x_column_preprocess_function="{return(v)}"
y_column_preprocess_function="{return(v)}"
label_each_point="4.0" # no labels if 0
legend_show="TRUE" # control if want to remove legends
horizontal_lines="list("
vertical_lines="list("

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

lastInputFile=""
while getopts "c:C:z:Z:i:o:R:Hhe:E:p:12:3:4:5b:B:T:X:Y:x:y:A:F:w:W:l:jJkKd:D:" OPTION
do
	case $OPTION in
		i)
			if [ "${inputFilesStr}" == "" ]; then
				inputFilesStr="${OPTARG}"
			else
				inputFilesStr+="${_internal_arrays_ifs_}$OPTARG"
				doColumnsStr+=${_doColumnsStr#${param_ifs}}${_internal_arrays_ifs_}; _doColumnsStr=""
				dontDoColumnsStr+=${_dontDoColumnsStr#${param_ifs}}${_internal_arrays_ifs_}; _dontDoColumnsStr=""
				inputFilesHaveHeaderStr+=${_internal_arrays_ifs_}
				inputFilesIFSStr+=${_internal_arrays_ifs_}
				inputFilesLabelsStr+=${_internal_arrays_ifs_}
			fi
			lastInputFile="$OPTARG"
			;;
                o)
                        if [ "$outputFile" != "" ]; then
                                echo "$appname : you can specify only one output folder"
                                exit 1
                        fi
                        outputFile=$OPTARG
                        ;;
		F)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column labels specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			inputFilesLabelsStr+="${OPTARG}"
			;;
		c)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column names specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			_doColumnsStr="${_doColumnsStr}${param_ifs}$OPTARG"
			;;
		C)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column names specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			_dontDoColumnsStr="${_dontDoColumnsStr}${param_ifs}$OPTARG"
			;;
		z)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column names specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			crap=$(column_names_to_column_numbers "${lastInputFile}" "${currentInputFileIFS}" "${OPTARG}")
			if [ $? -ne 0 ]; then echo "${appname} : could not find column name '${OPTARG}' in input file '${lastInputFile}', got: ${crap}"; exit 1; fi
			if [ "${crap}" != "" ]; then if [ "${_doColumnsStr}" != "" ] ; then _doColumnsStr+=${param_ifs}; fi; _doColumnsStr+=${crap}; fi
			;;
		Z)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column names specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			crap=$(column_names_to_column_numbers "${lastInputFile}" "${currentInputFileIFS}" "${OPTARG}")
			if [ $? -ne 0 ]; then echo "${appname} : could not find column name '${OPTARG}' in input file '${lastInputFile}', got: ${crap}"; exit 1; fi
			if [ "${crap}" != "" ]; then if [ "${_dontDoColumnsStr}" != "" ] ; then _dontDoColumnsStr+=${param_ifs}; fi; _dontDoColumnsStr+=${crap}; fi
			;;
		h)
			# in case there is a header (column names at the first row of file)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : setting default 'have header' to TRUE (for all files not having corresponding -h)."; defaultInputFilesHaveHeader="TRUE"
			else inputFilesHaveHeaderStr+="TRUE"; fi
			;;
		e)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : setting default input files IFS to '${OPTARG}' (for all files not having corresponding -e)."; defaultInputFilesIFS="${OPTARG}"
			else inputFilesIFSStr+="${OPTARG}"; fi
			currentInputFileIFS="${OPTARG}"
			;;
		# standard arguments
		2)
			seed=$OPTARG
			;;
		3)
			png_width=`echo "${OPTARG}" | cut -d':' -f1`
			png_height=`echo "${OPTARG}" | cut -d':' -f2`
			png_fontsize=`echo "${OPTARG}" | cut -d':' -f3`
			;;
                4)
                        outputPrecision=$OPTARG
                        ;;
		5)
			legend_show="FALSE"
			;;
		R)
			delete_Rscript_when_done=0
			Rscript_name="$OPTARG"
			;;
		H)
			display_help_exit "$appname"
			;;
		E)
			separatorOut=$OPTARG
			;;
		p)
			num_threads=$OPTARG
			;;
		# other options
		T)
			main_title="$OPTARG"
			;;
		X)
			xlabel="${OPTARG}"
			;;
		Y)
			ylabel="${OPTARG}"
			;;
		x)
			minx=`echo "$OPTARG"|cut -d':' -f1`
			maxx=`echo "$OPTARG"|cut -d':' -f2`
			limitsX=",limits=c($minx,$maxx)"
			;;
		y)
			miny=`echo "$OPTARG"|cut -d':' -f1`
			maxy=`echo "$OPTARG"|cut -d':' -f2`
			limitsY=",limits=c($miny,$maxy)"
			;;
		A)
			alpha=$OPTARG
			;;
		w)
			x_column_preprocess_function="$OPTARG"
			;;
		W)
			y_column_preprocess_function="$OPTARG"
			;;
		l)
			label_each_point="$OPTARG"
			;;
		j)
			x_log10_scale='TRUE'
			;;
		J)
			y_log10_scale='TRUE'
			;;
		k)
			x_log_scale='TRUE'
			;;
		K)
			y_log_scale='TRUE'
			;;
		d)
			vertical_lines+="c($OPTARG)," # requires 2,"red","dashed" - the intercept is on the scale after -w -W transforms
			;;
		D)
			horizontal_lines+="c($OPTARG),"
			;;
	esac
done
doColumnsStr+=${_doColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
dontDoColumnsStr+=${_dontDoColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
inputFilesHaveHeaderStr+=${_internal_arrays_ifs_}
inputFilesIFSStr+=${_internal_arrays_ifs_}
inputFilesLabelsStr+=${_internal_arrays_ifs_}
horizontal_lines=${horizontal_lines%,}")"; vertical_lines=${vertical_lines%,}")"

if [ "${inputFilesStr}" == "" ]; then
	echo "$0 : an input file must be specified using -i."
	exit 1
fi
if [ "$Rscript_name" = "" ]; then
	echo "$0 : script file name is not defined (use the -R option to do it)."
	exit 1
fi
if [ "$outputFile" = "" ]; then
	echo "$0 : you need to specify an output file using the '-o' option."
	exit 1
fi
if [ "${binspec_hist}" == "" ]; then
	binspec_hist="breaks=${numbins},"
	binspec_ggplot="binwidth=binwidth,"
fi

density_cmd="den=density(adfr\$x, na.rm=T, n=${num_density_points}, kernel='biweight')"

cat << EOC > "$Rscript_name"
set.seed($seed); write("$appname : starting R, seed is $seed", stderr());

suppressMessages(library(ggplot2))

# optional 
source("${bindir}/library_in_R.R")

# precision and when this precision is not enough
options(digits=$outputPrecision, scipen=100)

do_processing <- function(col_ids, oridata, num_padding, column_names, file_labels, include_column_numbers_in_output_file){
	output_basename=paste("${outputFile}",sep="")
	f1 = col_ids[1] # the file index
	x = col_ids[2]; y = col_ids[3]; cl = col_ids[4]; la = col_ids[5]
	p <- NULL
	a_title=NULL; if( "${main_title}" != "" ){ a_title=gsub("%cnumber1%", x, gsub("%cname1%", column_names[x], gsub("%cnumber2%", y, gsub("%cname2%", column_names[y], "${main_title}")))) }
	else { a_title=paste(file_labels[x], sep='') }
	if( "${xlabel}" == "" ){ xlabel=paste(column_names[[f1]][x], sep=''); } else { xlabel = "${xlabel}" }
	if( "${ylabel}" == "" ){ ylabel=paste(column_names[[f1]][y], sep=''); } else { ylabel = "${ylabel}" }
	crap = oridata[[f1]][complete.cases(oridata[[f1]]),]
	da = data.frame(
		x=unlist(lapply(crap[,x], x_column_preprocess_function)),
		y=unlist(lapply(crap[,y], y_column_preprocess_function)),
		stringsAsFactors=F
	)
	if( cl > 0 ){ da\$class = crap[,cl] } else { da\$class = rep(1, times=nrow(crap)) }
	# optional label
	if( la > 0 ){ da\$label = crap[,la]; colnames(da) <- c('x','y', 'class', 'label') } else { colnames(da) <- c('x','y', 'class') } 
	p <- plot_scatterplot(
		data=da,
		axis.label.x=xlabel,
		axis.label.y=ylabel,
		box.labels.fontsize=${label_each_point},
		axis.log10scale.x=${x_log10_scale},
		axis.log10scale.y=${y_log10_scale},
		axis.logscale.x=${x_log_scale},
		axis.logscale.y=${y_log_scale},
		legend.show=${legend_show}, # don't need legend if have labels and 1 class
		lines.horizontal=${horizontal_lines},
		lines.vertical=${vertical_lines},
		plot.title=a_title,
		render=F
	)
	a_filename=paste(sep='',output_basename, ".png");
	png(filename=a_filename, width=${png_width}, height=${png_height}, pointsize=${png_fontsize})
	print(p)
	dev.off()

	write(paste("$appname : col ", x, " vs ", y, " done.", sep=''), stderr())

	return(TRUE)
} # end do_processing function

# R script by $appname on $now
inp_datas <- process_column_specs_and_read_input_files(
	inputFilesStr="${inputFilesStr}",
	inputFilesLabelsStr="${inputFilesLabelsStr}",
	doColumnsStr="${doColumnsStr}",
	dontDoColumnsStr="${dontDoColumnsStr}",
	param_ifs1="${param_ifs}",
	param_ifs2="${_internal_arrays_ifs_}",
	inputFilesHaveHeaderStr="${inputFilesHaveHeaderStr}",
	inputFilesIFSStr="${inputFilesIFSStr}",
	defaultInputFilesHaveHeader=${defaultInputFilesHaveHeader}, # this is boolean!!! dont quote
	defaultInputFilesIFS="${defaultInputFilesIFS}",
	return_back_data_from_files=TRUE
); if( is.null(inp_datas) ){ write("${appname} : call to process_column_specs_and_read_input_files has failed.", stderr()); quit(status=1) }

x_column_preprocess_function <- function(v)${x_column_preprocess_function}
y_column_preprocess_function <- function(v)${y_column_preprocess_function}

column_pairs = list() # this is a mixture of file numbers and column in that file number, e.g. 1 2 3 4 (file 1, col2 with fil3, col4)
for(f1 in 1:inp_datas\$num_input_files){
	# wr are expecting at least 2 columns: x and y. Optional third is class/labels
	if( inp_datas\$num_columns_to_do[f1] == 2 ){
		# just x, y
		column_pairs[[f1]] = c(f1, inp_datas\$columns_to_do[[f1]][1], inp_datas\$columns_to_do[[f1]][2], -1, -1)
	} else if( inp_datas\$num_columns_to_do[f1] == 3 ){
		# we have class
		column_pairs[[f1]] = c(f1, inp_datas\$columns_to_do[[f1]][1], inp_datas\$columns_to_do[[f1]][2], inp_datas\$columns_to_do[[f1]][3], -1)
	} else if( inp_datas\$num_columns_to_do[f1] == 4 ){
		# we have label column
		column_pairs[[f1]] = c(f1, inp_datas\$columns_to_do[[f1]][1], inp_datas\$columns_to_do[[f1]][2], inp_datas\$columns_to_do[[f1]][3], inp_datas\$columns_to_do[[f1]][4])
	} else { write(paste("${appname} : at least two column names/numbers are expected for each input file, with an optional third, the spec for input file '", inp_datas\$inputFiles[f1], "' had ", inp_datas\$num_columns_to_do[f1], " items.", sep=''), stderr()); quit(status=1) }
}
num_processes = length(column_pairs)

stime = NULL; res = NULL
if( num_processes > 1 ){
	include_column_numbers_in_output_file=TRUE
	if( ${num_threads} > 1 ){
		suppressMessages(library(multicore))
		write(paste(sep='',"${appname} : parallelising ", num_processes, " processes over ${num_threads} threads..."), stderr())
		stime <- system.time({ res<-mclapply(column_pairs, do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file, mc.cores=${num_threads}) })
	} else {
		write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
		stime <- system.time({ res<-lapply(column_pairs,  do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file) })
	}
} else {
	include_column_numbers_in_output_file=FALSE
	write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
	stime <- system.time({ res<-lapply(column_pairs,  do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file) })
}
failed=F; for(i in 1:length(res)){ if( (is.logical(res[[i]])==FALSE) || (res[[i]]==FALSE) ){ write(paste("${appname} : processing failed for process ", i, " (col spec:", paste(column_pairs[[i]],sep='',collapse=','),"). Message was:\n", res[[i]], sep=''), stderr()); failed=T } }
if( failed == T ){
	write(paste(sep='',"\n${appname} : failed."), stderr())
	quit(status=1)
}
write(paste(sep='',"\n${appname} : done all, with times user:", stime[1], ", system:", stime[2], ", elapsed:", stime[3]), stderr())
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

