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
num_threads=1

inputFilesStr=""
inputFilesHaveHeaderStr=""
inputFilesIFSStr=""; currentInputFileIFS=${defaultInputFilesIFS}
_doColumnsStr=""
_dontDoColumnsStr=""
doColumnsStr=""
dontDoColumnsStr=""
inputFilesLabelsStr="" # in legends and titles, use these labels instead of the filename
columnLabelsStr="" # use these labels instead of column headers
outputFile=""
# for standard arguments
separatorOut='\t'
Rscript_name="${rand}.R"
outputPrecision="5"
include_discarded_columns_in_output=0
delete_Rscript_when_done=1
limitsX="NULL"; limitsY="NULL"
xlabel=""; main_title=""; sub_title=""; alpha="0.15,"
plot_histogram_too="-1"

function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file [-c col_num -c | -z col_name -z col_name] [-C col_num_not_to_do | -Z col_name_not_to_do] [-e input_data_separator] [-E output_data_separator] [-R R_output_script_name] [-2 seed] [-3 png_width:png_height:png_fontsize] [-4 output_file_precision] [-p num_threads] [-1 params_ifs]"
	echo "This script will ..."
	echo "STANDARD OPTIONS:"
	echo " -i inputFile"
	echo "	defines the input file."
	echo " -o outputFile"
	echo "	defines the output file."
	echo " -c col_num"
	echo "	choose the column number you wish to operate (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
	echo " -z col_name"
	echo "	choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (-1), default is '${param_ifs}' or repeat -z options"
	echo " -C col_num"
	echo "	choose the column number NOT to operate (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
	echo " -Z col_name"
	echo "	choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (-1), default is '${param_ifs}' or repeat -Z options"
	echo " -1 params_ifs : dont use, default is '${params_ifs}'"
	echo " -2 random_number_generator_seed"
	echo "	useful to repeat experiments, otherwise seed is selected at random using the process id"
	echo " -3 png_width:png_height:png_fontsize : the size of the output image if any in pixels and font size, default: ${png_width}:${png_height}:${png_fontsize}"
	echo " -4 precision : the number of digits in writing out results, R may not honor this."
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
	echo " -t binwidth"
	echo "        plot histogram too with specified bin width if binwidth>0 OR AUTOMATIC binwidth calculation if binwidth==0"
	echo ""
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
if [ "$*" = "" ]; then
	display_help_exit $appname
fi

lastInputFile=""
while getopts "i:o:F:c:C:z:Z:L:1:2:3:4:R:Hhe:E:p:x:y:T:a:t" OPTION; do
	case $OPTION in
		# other options here:
		t)
			plot_histogram_too=0
			;;
		# standard arguments, nothing to change here
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
		L)
			if [ "${lastInputFile}" == "" ]; then echo "${appname} : an input file for which the column names specified here apply to must be given (-i) BEFORE THIS option."; exit 1; fi
			if [ "${columnLabelsStr}" != "" ]; then columnLabelsStr+=${param_ifs}; fi; columnLabelsStr+=${OPTARG}
			# this is a dirty hack
			xlabel=${OPTARG}
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
		1)
			param_ifs=$OPTARG
			;;
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
		h)
			# in case there is a header (column names at the first row of file)
			header="TRUE"
			;;
		E)
			separatorOut=$OPTARG
			;;
		p)
			num_threads=$OPTARG
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
		T)
			main_title="${OPTARG}"
			;;
		a)
			alpha="${OPTARG}"
			;;
		# leave this LAST
		*)
			echo "${appname} : unknown option '${OPTION} ${OPTARG}'."; exit 1
			;;
	esac
done
doColumnsStr+=${_doColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
dontDoColumnsStr+=${_dontDoColumnsStr#${param_ifs}}${_internal_arrays_ifs_}
inputFilesHaveHeaderStr+=${_internal_arrays_ifs_}
inputFilesIFSStr+=${_internal_arrays_ifs_}
inputFilesLabelsStr+=${_internal_arrays_ifs_}
columnLabelsStr+=${_internal_arrays_ifs_}

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
set.seed($seed); write("$appname : starting R, seed is $seed, loading libraries ...", stderr());

# R script by $appname on $now
suppressMessages(library(flexmix)) # put it before ggplot2 as it gives problesm about signature 'empty'
suppressMessages(library(mclust))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

source("${bindir}/library_in_R.R")

# depending on what you need to do here, you will get a spec of column numbers and file numbers in 'columns_to_do_in'
# e.g. 1 2 1 3, which means file 1, col 2 and file 2 col 3. This depends entirely on the column spec enumeration, after this function
# so most of this code should be adjusted.
do_processing <- function(columns_to_do_in, oridata, num_paddings, column_names, file_labels, include_column_numbers_in_output_file=TRUE){
	the_file = columns_to_do_in[1]; the_column = columns_to_do_in[2]; the_file_label = file_labels[the_file]
	datain_col = oridata[[the_file]][,the_column]
	if( (sd(datain_col) < 0.000001) ){ write(paste("${appname} : column: ", the_column, " file: '",the_file_label, "', has no variation, will not continue.", sep=""), stderr()); return(TRUE); }

	if( include_column_numbers_in_output_file == TRUE ){
		output_basename = paste(sep='.', "${outputFile}", sprintf(paste("%0", num_paddings[the_file], "d", sep=''), the_column))
	} else {
		output_basename=paste("${outputFile}",sep="")
	}
	pngoutfilename=paste(sep='', output_basename, ".png");

	write(paste("${appname} : do_processing() : i am called for column: ", the_column, " of file '", the_file_label, "', png output: '", pngoutfilename, "'.", sep=""), stderr())

	a_file_label = paste(" (", file_labels[the_file], ")", sep='');
	if( "${xlabel}" == "" ){ xlabel=paste(column_names[[the_file]][the_column], "::", the_column, a_file_label, sep=''); } else { xlabel = "${xlabel}" }
	a_title=NULL; if( "${main_title}" != "" ){ a_title=gsub("%cnumber1%", the_column, gsub("%cname1%", xlabel, "${main_title}")) }
	else {
		a_title=NULL
	}

	if( "$plot_histogram_too" != "-1" ){ showh = TRUE } else { showh = FALSE }

	dataframe_toplot=data.frame(datax=datain_col,class=1)
	p <- plot_density(
		data=dataframe_toplot,
		axis.label.x=xlabel,
		axis.limits.x=${limitsX},
		axis.limits.y=${limitsY},
		fontsize=${png_fontsize},
		density.alpha=${alpha},
		density.kernel='biweight',
		plot.title=a_title,
		histogram.show=showh,
		histogram.binwidth=$plot_histogram_too,
		histogram.alpha=0.2,
		histogram.fill='blue',
		legend.show=FALSE
	)
	png(filename=pngoutfilename, width=${png_width}, height=${png_height},pointsize=${png_fontsize})
	print(p)
	dev.off()
	# precision and when this precision is not enough
	write(paste("do_processing : column: ", the_column, " file: '",the_file_label, "' : finished, output basename: '", output_basename, "'.", sep=""), stderr())
	return(TRUE)
} # end of do_processing function

write("${appname} : processing input params and reading data ... ", stderr())
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

# The following loop creates all processes by enumerating the input data columns to be processed, you need to change this:
# find all pairwise combinations of columns to do from each input file all mixed up, e.g. f1:1,2 f2:3,4
column_pairs = list() # this is a mixture of file numbers and column in that file number, e.g. 1 2 3 4 (file 1, col2 with fil3, col4)
i=1;
for(f1 in 1:inp_datas\$num_input_files){
	# straightforward 1 file, all columns to do, 2nd file, all columns to do.
	for(p1 in 1:inp_datas\$num_columns_to_do[f1]){
		column_pairs[[i]] = c(f1, inp_datas\$columns_to_do[[f1]][p1])
		i = i + 1
	}
}
# end change

num_processes = length(column_pairs)
stime = NULL; res=NULL
if( num_processes > 1 ){
	include_column_numbers_in_output_file=TRUE
	if( ${num_threads} > 1 ){
		suppressMessages(library(multicore))
		write(paste(sep='',"${appname} : parallelising ", num_processes, " processes over ${num_threads} threads..."), stderr())
		stime <- system.time({ res <- mclapply(column_pairs, do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file, mc.cores=${num_threads}) })
	} else {
		write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
		stime <- system.time({ res <- lapply(column_pairs,  do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file) })
	}
} else {
	include_column_numbers_in_output_file=FALSE
	write(paste(sep='',"${appname} : running sequentially ", num_processes, " processes..."), stderr())
	stime <- system.time({ res <- lapply(column_pairs,  do_processing, inp_datas\$oridata, inp_datas\$num_paddings, inp_datas\$column_names, inp_datas\$file_labels, include_column_numbers_in_output_file) })
}
failed=F; for(i in 1:length(res)){ if( (is.logical(res[[i]])==FALSE) || (res[[i]]==FALSE) ){ write(paste("${appname} : processing failed for process ", i, "(", paste(column_pairs[[i]],sep='',collapse=','),"). Message was:\n", res[[i]], sep=''), stderr()); failed=T } }
if( failed == T ){
	write(paste(sep='',"\n${appname} : failed."), stderr())
	quit(status=1)
}
write(paste(sep='',"\n${appname} : done all, with times user:", stime[1], ", system:", stime[2], ", elapsed:", stime[3]), stderr())
EOC

echo "$appname : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
#	cat /tmp/rlog$rand.txt; echo "$appname : error executing R script ($Rscript_name)";
	echo "$appname : error executing R script ($Rscript_name)";
	echo "Command line:"
	echo "$0 ${*}"
	exit 1
fi
echo "Command line:"
echo "$0 ${*}"
rm -f /tmp/rlog$rand.txt /tmp/$rand.* 
if [ "${delete_Rscript_when_done}" -eq 1 ]; then rm -f "$Rscript_name"; fi
exit 0
