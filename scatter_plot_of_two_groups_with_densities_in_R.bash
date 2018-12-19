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
png_fontsize=20
points_size=1.5
points_alpha=0.35
num_threads=1
		 
function display_help_exit {
	echo "Usage : $1 -i input_file -c x:y | -z nx:ny [-i input_file -c x:y | -z nx:ny ...] -o output_file [-e input_data_separator] [-E output_data_separator] [-R R_output_script_name] [-2 seed] [-3 png_width:png_height:png_fontsize]"
	echo "specify one or more files of which two columns 'x' and 'y' (using -c x:y) will be used to create a scatter plot of column 'x' against 'y'. More files will be plotted on the same scatter plot with different colors - hope that the scales match."
	echo "other params: [-x min:max] [-y min:max] [-T title] [-t binwidth] [-X xlabel] [-Y ylabel] [-a points_alpha] [-s points_size] [-r pearson|spearman|...] [-b] [-B]"
	echo "This script will take N input files (N>=1) and plot pairwise columns, each file a different color, each column will have its own histogram."
	echo "-x min:max : specify range of the x-axis from min to max"
	echo "-y min:max : same for y-axis"
	echo "-T title : title of the plot"
	echo "-X xlabel : label of the x-axis"
	echo "-Y ylabel : label of the y-axis"
	echo "-F label1 [-F label2 ...] : specify as many labels as input files, each label will be used to label data from that file, if not present, the name of the file will be used instead"
	echo "-g legend_title : legend title"	
	echo "-t binwidth : plot a histogram as well for each column and use this binwidth"
	echo "-s points_size : the size of the points in the scatterplot, default is ${points_size}"
	echo "-a points_alpha: the alpha of the points in the scatterplot, default is ${points_alpha}"
	echo "-r method: find correlation, do a correlation test and best-fit a line on the scatter plot using the specified method, e.g. pearson or spearman. By default no correlation will be calculated."
	echo "-b : add a best-fit line on each scatterplot with confidence region"
	echo "-B : add a best-fit curve on each scatterplot with confidence region"
	echo "This script will plot densities of all the columns specified"
	echo " -o outputFile"
	echo "	defines the output file."
	echo " -2 random_number_generator_seed"
	echo "	useful to repeat experiments, otherwise seed is selected at random using the process id"
	echo " -3 png_width:png_height:png_fontsize : the size of the output image if any in pixels, default: ${png_width}:${png_height}:${png_fontsize}"
# standard features
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
	echo "example:"
	echo "scatter_plot_of_two_groups_with_densities_in_R.bash -i d1 -c 1,2 -i d2 -c 1,2 -o ii -R rr -h -t fdu -g dgdgdgdg -F cel1 -F cel2"
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
inputFilesLabelsStr="" # in legends and titles, use these labels instead of the filename
inputFilesHaveHeaderStr=""
inputFilesIFSStr=""; currentInputFileIFS=${defaultInputFilesIFS}
doColumnsStr=""
dontDoColumnsStr=""
_doColumnsStr=""
_dontDoColumnsStr=""

outputFile=""
# for standard arguments
separatorOut='\t'
Rscript_name="${rand}.R"
outputPrecision="5"
include_discarded_columns_in_output=0
delete_Rscript_when_done=1

limitsX="NULL"; limitsY="NULL"
xlabel=""; ylabel=""
plot_histogram_too=0
main_title=""
legend_title=""
last_input_file=""
use_different_colors='TRUE'
use_different_shapes='FALSE'
correlation_method=""
add_best_fit_line='FALSE'
add_best_fit_curve='FALSE'

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

lastInputFile=""
while getopts "c:z:i:o:R:Hhe:E:p:2:3:4:x:y:X:Y:t:a:T:g:F:C:Z:Ss:r:bB" OPTION; do
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
		o)
			if [ "$outputFile" != "" ]; then
				echo "$appname : you can specify only one output folder"
				exit 1
			fi
			outputFile=$OPTARG
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
		X)
			xlabel="${OPTARG}"
			;;
		Y)
			ylabel="${OPTARG}"
			;;
		F)
			if [ "${input_file_labels}" == "" ]; then input_file_labels="${OPTARG}"; else input_file_labels="${input_file_labels}${param_ifs}${OPTARG}"; fi
			;;
		T)
			main_title="${OPTARG}"
			;;
		g)
			legend_title="+labs(colour='${OPTARG}')"
			;;
		t)
			plot_histogram_too="$OPTARG"
			;;
		a)
			points_alpha=$OPTARG
			;;
		s)
			points_size=$OPTARG
			;;
		r)
			correlation_method=$OPTARG
			;;
		C)
			use_different_colors='FALSE'
			;;
		S)
			use_different_shapes='TRUE'
			;;
		b)
			add_best_fit_line='TRUE'
			;;
		B)
			add_best_fit_curve='TRUE'
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

cat << EOC > "$Rscript_name"
# R script by $appname on $now
set.seed($seed); write("$appname : starting R, seed is $seed", stderr());
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
source("${bindir}/library_in_R.R")
# precision and when this precision is not enough
options(digits=$outputPrecision, scipen=100)

do_processing <- function(columns_to_do_in, oridata_in, num_paddings, column_names, file_labels, include_column_numbers_in_output_file=TRUE){
	the_files = 1:columns_to_do_in[1]
	c1 = columns_to_do_in[2]; c2 = columns_to_do_in[3]
	num_input_files = length(the_files)
	if( include_column_numbers_in_output_file == TRUE ){
		output_basename = sprintf(paste("${outputFile}.%0", num_paddings[1], "d", ".%0", num_paddings[1], "d", sep=''), c1, c2)
	} else { output_basename=paste("${outputFile}",sep="") }
	write(paste("do_processing : i am called for columns: ", c1, ", ", c2, ", output basename: '", output_basename, "'.", sep=""), stderr())

	alldata = NULL; details=NULL; a_title=NULL; colnames_per_file_x=c(); colnames_per_file_y=c(); class.labels=c()
	if( "${details}" != "" ){ details="${details}" }; if( "${main_title}" != "" ){ a_title="${main_title}" }
	for(i in 1:num_input_files){
		datafr=data.frame(x=oridata_in[[the_files[i]]][!is.na(oridata_in[[the_files[i]]][,c1]),c1],
				  y=oridata_in[[the_files[i]]][!is.na(oridata_in[[the_files[i]]][,c2]),c2],
				  class=toString(i)
		)
		if( is.null(alldata) ){ alldata=datafr } else { alldata=rbind(alldata, datafr) }
		colnames_per_file_x[i] = column_names[[the_files[i]]][c1]
		colnames_per_file_y[i] = column_names[[the_files[i]]][c2]
		class.labels[i]=paste(file_labels[the_files[i]]," (",colnames_per_file_x[i],":",colnames_per_file_y[i],")",sep='')
		if( "${details}" != "" ){ details=gsub("%cnumber1%", c1, gsub("%cnumber2%", c2, gsub(sprintf("%%fname%d%%", i), file_labels[the_files[i]], gsub("%cname1%", column_names[[the_files[i]]][c1], gsub("%cname2%", column_names[[the_files[i]]][c2], details))))) }
		if( "${main_title}" != "" ){ a_title=gsub("%cnumber1%", c1, gsub("%cnumber2%", c2,gsub(sprintf("%%fname%d%%", i), file_labels[the_files[i]],gsub("%cname1%", column_names[[the_files[i]]][c1], gsub("%cname2%", column_names[[the_files[i]]][c2], a_title))))) }
	}
	if( !is.null(details) ){ details=paste(details, "\n", sep="") }
	if( is.null(a_title) ){ a_title=paste(unique(file_labels), collapse=' vs ', sep='') }
	if( "${xlabel}" == "" ){ xlabel = paste(colnames_per_file_x,collapse='/') } else { xlabel = "${xlabel}" }
	if( "${ylabel}" == "" ){ ylabel = paste(colnames_per_file_y,collapse='/') } else { ylabel = "${ylabel}" }
	a_filename=paste(sep='',output_basename, ".png");
write(a_title, stderr())
	png(filename=a_filename, width=${png_width}, height=${png_height}, pointsize=${png_fontsize})
	# the scatter plot
	ret=plot_scatterplot_and_density_of_many_groups(
		data=alldata,
		axis.label.x=xlabel,
		axis.label.y=ylabel,
		axis.limits.x=${limitsX},
		axis.limits.y=${limitsY},
		class.labels=class.labels,
		scatterplot.points.size=${points_size},
		scatterplot.points.alpha=${points_alpha},
		class.different_colors=${use_different_colors},
		class.colors.scale.begin='red',
		class.colors.scale.end='green',
		class.different_shapes=${use_different_shapes},
		fontsize=${png_fontsize},
		correlation.method="${correlation_method}",
		correlation.test.method="${correlation_method}",
		add.bestfit.line=${add_best_fit_line},
		add.bestfit.curve=${add_best_fit_curve},
		plot.title=a_title,
		legend.fontsize=15,
		render=TRUE
	)
	if( is.null(ret) == T ){ write(paste("do_processing : error call to plot_scatterplot_and_density_of_many_groups has failed (usually it fails if no variation in data) , for columns: ", c1, ",", c2, ", output basename: '", output_basename, "'.", sep=""), stderr()); return(FALSE); }
	dev.off()
	write(paste("do_processing : finished, for columns: ", c1, ",", c2, ", output basename: '", output_basename, "'.", sep=""), stderr())
	return(TRUE)
} # end of do_processing function

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

# find all pairwise combinations of columns to do from each input file all mixed up, e.g. f1:1,2 f2:3,4
column_pairs = list() # this is a mixture of file numbers and column in that file number, e.g. 1 2 3 4 (file 1, col2 with fil3, col4)
i=1; l = -1
for(f1 in 1:inp_datas\$num_input_files){
	if( l == -1 ){ l = inp_datas\$num_columns_to_do[f1] } else { if( l != inp_datas\$num_columns_to_do[f1] ){ write(paste("${appname} : all files must have the same number of columns (or specified columns to be processed), the columns of file '", inp_datas\$inputFiles[f1], "' is not ", l, " like the previous files.", sep=''), stderr()); quit(status=1) } }
	for(p1 in 1:inp_datas\$num_columns_to_do[f1]){
		for(p2 in 1:inp_datas\$num_columns_to_do[f1]){
			if( p2 <= p1 ){ next }
			column_pairs[[i]] = c(inp_datas\$num_input_files, inp_datas\$columns_to_do[[f1]][p1], inp_datas\$columns_to_do[[f1]][p2])
			i = i + 1
		}
	}
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
#debug(do_processing)
#do_processing(column_pairs[[1]],inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, T)
EOC

echo "$appname : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
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
