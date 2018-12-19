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

# the first mode is the default
declare -a known_histogram_modes=('counts' 'probability' 'density'); histogram_mode=${known_histogram_modes[0]}
declare -a known_density_kernels=('biweight' 'gaussian' 'rectangular' 'triangular' 'epanechnikov' 'cosine' 'optcosine'); density_kernel=${known_density_kernels[0]}
declare -a known_density_estimators=('ks.kde' 'base.density'); density_estimator=${known_density_estimators[0]}

function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file [-c col_num ... -c ... | -C col_num ... -C ...] [-z col_name ... -z ... | -Z col_name ... -Z ...] [-1] [-e input_data_separator] [-E output_data_separator] [-R R_output_script_name] [-2 seed] [-3 png_width:png_height:png_fontsize] [-p num_threads]"
	echo "other params: [-b numbins | -B actual_bins] [-m histogram_mode] [-N num_density_points] [-T title] [-X xlabel] [-Y ylabel] [-x min:max] [-y min:max] [-A alpha] [-l] [-Q] [-P] [-k density_estimation_kernel] [-K select_density_estimator]"
	echo "This script will calculate the histogram of each column specified (separately)"
	echo " -i inputFile"
	echo "	defines the input file."
	echo " -o outputBasename"
	echo "	defines the output basename. The column number and '.png' '.info.txt' '.labels.txt' '.discrete.txt' may be output files (the first for the images, for each column), the second for histogram info (bins, mids, etc) per column, and the other two may be produced if -I or -L are specified, see below."
	echo " -c col_num"
	echo "	choose the column number you wish to operate as a comma-separated list of column numbers (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
	echo " -C col_num"
	echo "	choose the column number you DO NOT wish to operate as a comma-separated list of column numbers (starting from 1) - repeate -C options if necessary, all columns will be processed unless this option is specified"
	echo " -z col_name"
	echo "	choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (-1), default is '${param_ifs}' or repeat -z options"
	echo " -Z col_name"
	echo "	choose the column name you DO NOT wish to operate as a list of column names or repeate -Z options. Use unique separator (-1), default is '${param_ifs}' or repeat -Z options"
	echo " -1"
	echo "	include columns not selected for processing to the output, default is not to print un-processed columns"
	echo " -2 random_number_generator_seed"
	echo "	useful to repeat experiments, otherwise seed is selected at random using the process id"
	echo " -3 png_width:png_height:png_fontsize : the size of the output image if any in pixels and font size, default: ${png_width}:${png_height}:${png_fontsize}"
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
	echo "other params: [-b numbins | -B actual_bins] [-m histogram_mode] [-N num_density_points] [-T title] [-X xlabel] [-Y ylabel] [-x min:max] [-y min:max] [-A alpha] [-l] [-Q] [-I] [-L]"
	echo "  -b numbins : the number of bins > 1"
	echo "  -B actual_bins : the actual bins (breaks) as a comma-separated string, e.g. A,B,C means two bins one between A-B and the other between B-C. Specifying the actual bins works only for histogram mode 'counts' and 'probability' - for example if your values are all 0 or 1, then use '-B -0.5,0.5,1.5'"
	echo "  -m histogram mode : how to calculate the histogram. Known modes are: ${known_histogram_modes[@]}. Default mode is ${histogram_mode}"
	echo "  -N num_density_points : the number of points during the calculation of density if density is needed (e.g. if histogram mode is NOT counts OR equal-height bins is requested), default is ${num_density_points} points. In general, the larger the number of points the more accurate the density BUT data can not be invented..."
	echo "  -Q : Calculate variable bins (breaks) so that each bin contains the same number of points - more or less. The method is based on estimating the density and then using the density for this calculation, therefore the result may be correct with density but when it comes to actual counts, it might have slight errors - i.e. not all bins have the same height"
	echo "  -l : pass each input column via a logarithmic function. The transform is log(column+1.001*min(column)) which does not require that values are positive."
	echo "  -P : do not save histogram as images, just save the bins etc into the '.info.txt' file per column. Default is to save histograms to images."
	echo "  -k : kernel to be used for the density estimation if needed. Known kernels (R, ?density): ${known_density_kernels[@]}. Default kernel is '${density_kernel}'."
	echo "  -K : select the density estimator. Known estimators are ${known_density_estimators[@]}. Default is '${density_estimator}'."
	echo "  DISCRETISE DATA OPTIONS: discretise each column by assigning each value to its corresponding bin then either:"
	echo "  -I : save the bin mids for each input column value. Output file postfix: '.N.discrete.txt', where N is the column number"
	echo "  -L : save the bin ID (first bin has ID=1, etc.) for each input column value. Output file postfix: '.N.labels.txt', where N is the column number"
	echo "  AESTHETICS OPTIONS:"
	echo "  -A alpha : alpha value for drawing the boxes, see also -j, -J, default alpha value is ${alpha} (from 0 to 1)"
	echo "  -T title : title (%cnumber% replaces column number, %cname% replaces column name)"
	echo "  -x min:max : the range of values of x-axis, useful for drawing histograms on the same scale for comparisons or in a movie"
	echo "  -y min:max : same for y-axis"
	echo "  -X label : the label for the x-axis, default is the column header name if a header exists or just 'col:N' for the Nth column"
	echo "  -Y label : the label for the y-axis, default depends on the histogram mode, i.e. density, counts, probability etc."
	echo "  -j histogram_fill_color : color to fill the histogram, default is ${histogram_fill_color}, alpha (-A) is also relevant."
	echo "  -J histogram_border_color : border color of the histogram, default is ${histogram_border_color}"
	echo ""
	echo "EXAMPLE:"
	echo "create_histogram_in_R.bash -i all1 -o gyy -R rr -m counts -I -L -b 30 -Q -k 'gaussian' -c 1"
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

limitsX="NULL"; limitsY="NULL"; xlabel=""; ylabel=""; title=""; legend_title=""
do_density=0
assign_data_to_bins_and_save_it_to_file=0
label_input_data_and_dump_it=0
numbins=""
actualBins=""
binwidth_cmd="data_range=diff(range(adfr\$x,na.rm=T)); numbins=as.integer(0.5+(data_range/den\$bw)); binwidth=data_range/numbins;" # this just sets the binwidth automatically by calling the R's density function
binspec_hist=""
binspec_ggplot=""
equal_height_bins=0
produce_images=1
histogram_fill_color='white'
histogram_border_color='darkgreen'

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

lastInputFile=""
while getopts "c:C:z:Z:i:o:R:Hhe:E:p:12:3:4:b:B:G:T:X:Y:x:y:A:gm:ILN:QPk:K:F:j:J:" OPTION
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
			_doColumnsStr="${_doColumnsStr}${param_ifs}$OPTARG"
			;;
		C)
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
		b)
			numbins=${OPTARG}
			binwidth_cmd="data_range=diff(range(adfr\$x,na.rm=T)); numbins=${numbins}; binwidth=data_range/numbins;"
			;;
		B)
			# the actual bins (may be un-equal bins) as a comma-separated string(e.g. '1,2,3') specified:
			binwidth_cmd="data_range=diff(range(adfr\$x,na.rm=T)); binwidth=NA; numbins=NA;"
			binspec_hist="breaks=c(${OPTARG}),"
			binspec_ggplot="breaks=c(${OPTARG}),position='dodge',"
			;;
		Q)
			equal_height_bins=1
			;;
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
			limitsX="c($minx,$maxx)"
			;;
		y)
			miny=`echo "$OPTARG"|cut -d':' -f1`
			maxy=`echo "$OPTARG"|cut -d':' -f2`
			limitsY="c($miny,$maxy)"
			;;
		A)
			alpha=$OPTARG
			;;
		g)
			log_preprocess=1
			;;
		m)
			histogram_mode=$OPTARG
			found=0; for hm in ${known_histogram_modes[@]}; do if [ "${histogram_mode}" == "${hm}" ]; then found=1; break; fi; done
			if [ ${found} -eq 0 ]; then echo "$appname : specified histogram mode '$histogram_mode' is now known, known modes are: ${known_histogram_modes[@]}"; exit 1; fi
			;;
		k)
			density_kernel=$OPTARG
			found=0; for dk in ${known_density_kernels[@]}; do if [ "${density_kernel}" == "${dk}" ]; then found=1; break; fi; done
			if [ ${found} -eq 0 ]; then echo "$appname : specified density kernel '$density_kernel' is not known, known kernels are: ${known_density_kernels[@]}"; exit 1; fi
			;;
		K)
			density_estimator=$OPTARG
			found=0; for dk in ${known_density_estimators[@]}; do if [ "${density_estimator}" == "${dk}" ]; then found=1; break; fi; done
			if [ ${found} -eq 0 ]; then echo "$appname : specified density kernel '$density_estimator' is not known, known kernels are: ${known_density_estimators[@]}"; exit 1; fi
			;;
		I)
			assign_data_to_bins_and_save_it_to_file=1
			;;
		L)
			label_input_data_and_dump_it=1
			;;
		N)
			num_density_points=$OPTARG
			;;
		P)
			produce_images=0
			;;
		j)
			histogram_fill_color=$OPTARG
			;;
		J)
			histogram_border_color=$OPTARG
			;;
	esac
done
doColumnsStr+=${_doColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
dontDoColumnsStr+=${_dontDoColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
inputFilesHaveHeaderStr+=${_internal_arrays_ifs_}
inputFilesIFSStr+=${_internal_arrays_ifs_}
inputFilesLabelsStr+=${_internal_arrays_ifs_}

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

do_processing <- function(in_col_id, oridata, num_padding, column_names, file_labels, input_filenames, include_column_numbers_in_output_file){
	if( include_column_numbers_in_output_file == TRUE ){
		output_basename = paste(sep='.', "${outputFile}", sprintf(paste("%0", num_padding, "d", sep=''), in_col_id))
	} else {
		output_basename=paste("${outputFile}",sep="")
	}

	p <- NULL
	# logarithmic scale by taking the log of the input col
	a_title=NULL; if( "${main_title}" != "" ){ a_title=gsub("%cnumber1%", in_col_id, gsub("%cname1%", column_names[in_col_id], "${main_title}")) }
	else { a_title=paste(file_labels[in_col_id], sep='') }
	if( "${xlabel}" == "" ){ xlabel=paste(column_names[in_col_id], "::", in_col_id, paste(" (", file_labels[1], ")", sep=''), sep=''); } else { xlabel = "${xlabel}" }
	crap=oridata[!is.na(oridata[,in_col_id]),in_col_id] # do't count NA's AT ALL (not even averaging
	if( "${log_preprocess}" == "1" ){ adfr <- data.frame(x=log(crap-min(crap)+1)); xlabel=paste("log, ", xlabel, sep='') } else { adfr <- data.frame(x=crap)  }
	crap=NULL
	${density_cmd}
	${binwidth_cmd} # sets a var 'binwidth' to the binwidth either auto or specified, except when bins are variable

	write(paste("$appname : col ", in_col_id, " doing histogram with number of bins ", numbins, ", binwidth is ", binwidth, " data range: ", data_range, sep=""), stderr())

	breaks=NULL; if( (${equal_height_bins} == 1) && !is.na(numbins) ){
		if( "${histogram_mode}" != "counts" ){ write("$appname : warning, equal-height bins are applicable only for histogram mode 'counts' (and not for '${histogram_mode}'), will skip their calculation now...", stderr()); next }

		# use density from ks package: use.ks=T, it's better (?)
		chb=calculate_histogram_breaks(values=adfr\$x, numbins=numbins, equal_height_bins="bins", density_kernel='${density_kernel}', use.ks=T, calculate.error=T)
		breaks = chb\$breaks; den = chb\$density;
		write(paste("$appname : equal-height-bins, error: ", chb\$error, ", breaks:(",length(breaks),") : ", paste(breaks, collapse=",", sep=""), sep=""), stderr())
	}
	if( "${histogram_mode}" == "density" ){
		if( "${ylabel}" == "" ){ylabel="density"}else{ylabel="${ylabel}"}
		p <- ggplot(data=na.omit(adfr), aes(x=x))+geom_histogram(aes(y=..density..),colour="${histogram_border_color}",${binspec_ggplot} fill="${histogram_fill_color}",alpha=${alpha})+geom_density(colour='red',na.rm=T,n=${num_density_points})
	} else if( "${histogram_mode}" == "counts" ){
		if( "${ylabel}" == "" ){ylabel="count"}else{ylabel="${ylabel}"}
		if( !is.null(breaks) ){
			p <- ggplot(data=na.omit(adfr), aes(x=x))+geom_histogram(colour="${histogram_border_color}",breaks=breaks,position='dodge', fill="${histogram_fill_color}",alpha=${alpha})
		} else { 
			p <- ggplot(data=na.omit(adfr), aes(x=x))+geom_histogram(colour="${histogram_border_color}",${binspec_ggplot} fill="${histogram_fill_color}",alpha=${alpha})
		}
	} else if( "${histogram_mode}" == "probability" ){
		if( "${ylabel}" == "" ){ylabel="probability"}else{ylabel="${ylabel}"}
# from http://stackoverflow.com/questions/5033240/plot-probability-with-ggplot2-not-density
#		p <- ggplot(adfr, aes(x=x))+geom_histogram(aes(y = ..density..),colour="${histogram_border_color}",binwidth=density(adfr\$x,na.rm=T)\$bw,fill="${histogram_fill_color}",alpha=${alpha})+geom_density()
		p <- ggplot(data=na.omit(adfr), aes(x=x))+
			geom_histogram(aes(y=..density..),colour="${histogram_border_color}",${binspec_ggplot} fill="${histogram_fill_color}",alpha=${alpha})
#			geom_density(colour='red',na.rm=T,n=${num_density_points})
	}
	if( ${produce_images} == 1 ){
		p <- p +opts(axis.text.x=theme_text(size=${png_fontsize}*0.8)) + # these are magic numbers for just and pointsize
			opts(axis.text.y=theme_text(size=${png_fontsize}*0.8,hjust=1)) +
			opts(axis.title.x=theme_text(size=${png_fontsize})) +
			opts(axis.title.y=theme_text(size=${png_fontsize},angle=90,hjust=0.5)) +
			opts(panel.grid.major=theme_line(colour=alpha('blue',0.15),size=0.2),panel.grid.minor=theme_line(colour=alpha('purple',0.05),size=0.2),panel.background=theme_blank())+
			xlab(xlabel)+ylab(ylabel)+
# resort to a combination of these limit shits. coord on x-axis does not work...
			opts(legend.position="none")+
			xlim($limitsX)+
			ylim($limitsY)

# IT'S A SHIT!!!!! setting x and y limits AAAAAAFFAFAFAFAF

#			coord_cartesian(ylim=$limitsY)
#ORI			coord_cartesian(ylim=$limitsY)+
# scale_x_continuous(name=xlabel $limitsX) + scale_y_continuous(name=ylabel $limitsY)

#		if( "$limitsX" != "NULL" ){ p <- p + xlim($limitsX) }

		if( "${main_title}" != "" ){ p <- p + opts(title=a_title) } else { p <- p + opts(title=NULL) }
		a_filename = paste(output_basename, ".png", sep='')
		png(filename=a_filename, width=${png_width}, height=${png_height}, pointsize=${png_fontsize})
		print(p)
		dev.off()
	}

	# ok, now if we need to save, extract the breaks
	a_filename = paste(output_basename, ".info.txt", sep='')
	info = ggplot_build(p)\$data[[1]][[1]]
	FILEH <- file(a_filename, open="w")
		br=c(info\$xmin, tail(info\$xmax, 1))
		mi=info\$x # mids (of each bin) - number of bins is the length of these (not the breaks!)
		co=info\$count # counts (per bin)
		de=info\$density # density (per bin)
		write(paste("filename  : (", input_filenames[1], ")", sep=""), FILEH)
		write(paste("column    : (", in_col_id, ")", sep=""), FILEH)
		write(paste("num bins  : (", length(mi), ")", sep=""), FILEH)
		write(paste("bins      : (", paste(br, collapse=", "), ")", sep=""), FILEH)
		write(paste("mids      : (", paste(mi, collapse=", "), ")", sep=""), FILEH)
		write(paste("counts    : (", paste(co, collapse=", "), ")", sep=""), FILEH)
		write(paste("density   : (", paste(de, collapse=", "), ")", sep=""), FILEH)
	close(FILEH)
	
	# write out if required
	if( "${assign_data_to_bins_and_save_it_to_file}" == "1" ){
		# set all cols to the same (max) num rows, otherwise we can't save
		data_out_discrete = discretize_data(br, adfr\$x)
		if( "${header}" == "TRUE" ){ data_out_discrete = c(nams[in_col_id],data_out_discrete) }
		a_filename = paste(output_basename, ".discrete.txt", sep='')
		write.table(data_out_discrete, file=a_filename, sep="$separatorOut", col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE, na='')
		data_out_discrete=NULL
	}
	if( "${label_input_data_and_dump_it}" == "1" ){
		# set all cols to the same (max) num rows, otherwise we can't save
		data_out_labels = label_data(br, adfr\$x)
		if( "${header}" == "TRUE" ){ data_out_labels = c(nams[in_col_id],data_out_labels) }
		a_filename = paste(output_basename, ".labels.txt", sep='')
		write.table(data_out_labels, file=a_filename, sep="$separatorOut", col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE, na='')
		data_out_labels=NULL
	}
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

num_processes = inp_datas\$num_columns_to_do[1]
stime = NULL; res = NULL
if( num_processes > 1 ){
	include_column_numbers_in_output_file=TRUE
	if( ${num_threads} > 1 ){
		suppressMessages(library(multicore))
		write(paste(sep='',"${appname} : parallelising ", num_processes, " processes over ${num_threads} threads..."), stderr())
		stime <- system.time({ res <- mclapply(inp_datas\$columns_to_do[[1]], do_processing, inp_datas\$oridata[[1]], inp_datas\$num_paddings[1], inp_datas\$column_names[[1]], inp_datas\$file_labels[[1]], inp_datas\$inputFiles, include_column_numbers_in_output_file, mc.cores=${num_threads}) })
	} else {
		write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
		stime <- system.time({ res <- lapply(inp_datas\$columns_to_do[[1]],  do_processing, inp_datas\$oridata[[1]], inp_datas\$num_paddings[1], inp_datas\$column_names[[1]], inp_datas\$file_labels[[1]], inp_datas\$inputFiles, include_column_numbers_in_output_file) })
	}
} else {
	include_column_numbers_in_output_file=FALSE
	write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
	stime <- system.time({ res <- lapply(inp_datas\$columns_to_do[[1]],  do_processing, inp_datas\$oridata[[1]], inp_datas\$num_paddings[1], inp_datas\$column_names[[1]], inp_datas\$file_labels[[1]], inp_datas\$inputFiles, include_column_numbers_in_output_file) })
}
failed=F; for(i in 1:length(res)){ if( (is.logical(res[[i]])==FALSE) || (res[[i]]==FALSE) ){ write(paste("${appname} : processing failed for process ", i, "(col spec:", paste(inp_datas\$columns_to_do[[1]][i],sep='',collapse=','),"). Message was:\n", res[[i]], sep=''), stderr()); failed=T } }
if( failed == T ){
	write(paste(sep='',"\n${appname} : failed.", sep=''), stderr())
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

