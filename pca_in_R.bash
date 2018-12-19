#!/bin/bash

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

appname=`basename "$0"`
now=`date`
rand=$$
param_ifs=':'
                 
knownPreprocessingTransformations=(	'none' 'doScale=FALSE;doCenter=FALSE',
					'scale' 'doScale=TRUE;doCenter=FALSE'
					'center' 'doScale=FALSE;doCenter=TRUE'
					'scale_and_center' 'doScale=TRUE;doCenter=TRUE'
				)

function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file"
	echo "This script will read an excel (xlsx) file and output each sheet onto"
	echo "a separate file in an output folder"
	echo " -i inputFile"
	echo "        defines the input file."
	echo " -o outputFile"
	echo "        Output folder name, all input's sheets will be saved in"
	echo "	      this folder. Hopefully, the sheet's name will not mess up"
	echo "	      with the file system (e.g question marks, stars etc"
        echo " -c col_num"
        echo "        choose the column number you wish to operate as a comma-separated list of column numbers (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
        echo " -C col_num"
        echo "        choose the column number you DO NOT wish to operate as a comma-separated list of column numbers (starting from 1) - repeate -C options if necessary, all columns will be processed unless this option is specified"
        echo " -z col_name"
        echo "        choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (-1), default is '${param_ifs}' or repeat -z options"
        echo " -Z col_nam"
        echo "        choose the column number you DO NOT wish to operate as a list of column names or repeate -Z options. Use unique separator (-1), default is '${param_ifs}' or repeat -Z options"
        echo " -1 param_separator"
        echo "        character to separate the column names when you give more than one names via -n or -N. Default separator for this is '${param_ifs}'."
	echo " -l save_loadings_file"
	echo "	      File to save the loadings (the PCA transformation) into."
	echo "	      Then data can be transformed without re-calculating PCA"
	echo " -L load_loadings_file"
	echo "	      File to load the loadings (the PCA transformation) from."
	echo "	      PCA will not be calculated but the matrix in this file"
	echo "	      will transform the data. The -v option will be ignored,"
	echo "	      you will have to remove the columns you don't want manually"
	echo " -s summary_file"
	echo "        Save a summary of the PCA onto the specified file"
	echo " -d N"
	echo "        Remove the *LAST (i.e. less important)* N principal components irrespective of variance"
	echo " -D N"
	echo "        Remove the *FIRST (i.e. most important)* N principal components irrespective of variance"
	echo " -v proportion_of_variance"
	echo "        optional, a number from 0 to 1 denoting"
	echo "	      the proportion of variance you want to capture at the output file - default is 1,"
	echo "	      after the PCA transform. Use 1 and you will get all the components."
	echo "	      Use 0.75 and you will get the first X components whose total variance"
	echo "        proportion is 75% of the total (i.e. explain the data by 75% or the"
	echo "	      effect of the discarded columns will be a loss of 25% in the information"
	echo "	      your data conveys). Usually you use a number like 0.99 or 0.9 and see"
	echo "	      how many components you can elininate."
	echo " -V variance_cumvariance_file"
	echo "        Save or Load (depending on context) the variance and cumulative variance for"
	echo "        each principal component of the PCA in specified file as a vector."
	echo " -F"
	echo "        If the number of rows in the data file is less than the number of columns, then"
	echo "        the number of principal components produced will be equal to the number of rows,"
	echo "        in which case the program will not proceed unless you specify the -F option here"
	echo "        THIS OPTIONS IS NOT IMPLEMENTED : does it make sense to have less observations than variables?"
	echo "        NO, so you need to have as many PCs as observations, which I can't be bother"
	echo " -u"
	echo "        Flag denoting to do reverse PCA on some already-pcaed data by applying"
	echo "        then inverse of the loadings matrix to the PCAed data."
	echo "        For this, it is required that a loadings matrix is supplied (i.e. use -L file)"
	echo " -t preprocessingTransformationName"
	echo "       one of the following strings:"
	for((i=0;i<${#knownPreprocessingTransformations[@]};i+=2)); do echo "          ${knownPreprocessingTransformations[i]} -> ${knownPreprocessingTransformations[i+1]}"; done
	echo "       The aim of these transformations is to scale the data so as to have zero mean and unit standard deviation"
	echo " NOTE:"
	echo "       if you want to use correlation matrix instead of covariance, then apply -t scale_and_center"

# standard features
	echo " -h"
	echo "        optional, that the first row is not numerical but it holds column names"
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file, default is white space"
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs"
	echo ""
	echo "	      Example:"
	echo "$appname -i a.txt -o abc -h -c 16 -d 'gamma' -R out.R -b 5"
	echo "read column 16 from 'a.txt' and write to output basename 'abc'"
	echo "The distribution to test is 'gamma' and the number of bins when"
	echo "constructing histograms is 5"
	echo "-h means that the input file's first row is not numerical but it consists of row names "
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
#       echo extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit > /dev/stderr
cat << EOC | sh
	extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit > /dev/null
EOC
	if [ $? -ne 0 ]; then
		echo "$appname : ${FUNCNAME}, line ${LINENO} : call to the following command has failed"
		echo extract_columns.pl --input \"${_inputfile}\" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit
		return 1
	fi
	oldIFS="${IFS}"; IFS=${param_ifs}
	for _a_col_name in ${_column_names[*]}; do
#	      echo "checking for ${_a_col_name} (/tmp/${_rand}.cn2cn)" > /dev/stderr
		_a_col=`awk -F '\t' -v n="${_a_col_name}" '{if($4==n){print ($2+1); exit 2}}' /tmp/${_rand}.cn2cn`
		if [ $? -ne 2 ]; then
#		       echo "$appname : ${FUNCNAME}, line ${LINENO} : could not find column name '${_a_col_name}' in the header of the input file '${_inputfile}'."
#		       return 1
			continue # baby
		fi
#	       echo "found ${_a_col}" > /dev/stderr
		_ret+="${param_ifs}${_a_col}"
	done
	IFS="${oldIFS}"
	echo "${_ret#${param_ifs}}"
	rm -f /tmp/${_rand}.cn2cn
	return 0
}
inputFile=""
outputFile=""
# how much variance to keep, 1.0 is max but this can be achieved by first N components (i.e. last add zero variance),
# so put a large number here if we want to retain all components unless specified otherwise)
proportion_of_variance="2.0"
save_loadings_file=""
load_loadings_file=""
doReversePCA=0
# for standard arguments
separator=''; Rseparator=''
separatorOut='\t'
header="FALSE"
Rscript_name="${rand}.R"
outputPrecision="5"
doColumns=""
dontDoColumns=""
doColumnsNames=""
dontDoColumnsNames=""
do_only_selected_columns=0
preprocessingTransformation="doScale=FALSE;doCenter=FALSE"
include_discarded_columns_in_output=0
removeNFirstComponents=0
removeNLastComponents=0
force_pca_even_if_numrows_less_numcols=0

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

while getopts "c:C:z:Z:1:2v:l:L:i:o:R:Hhe:E:p:t:ud:D:F" OPTION
do
	case $OPTION in
                1)
                        param_ifs="$OPTARG" # ifs for separating column names ONLY
                        ;;
                c)
                        if [ "${doColumns}" == "" ]; then doColumns="${OPTARG}"; else doColumns="${doColumns}${param_ifs}$OPTARG"; fi
                        do_only_selected_columns=1
                        ;;
                C)
                        if [ "${dontDoColumns}" == "" ]; then dontDoColumns="${OPTARG}"; else dontDoColumns="${dontDoColumns}${param_ifs}$OPTARG"; fi
                        do_only_selected_columns=1
                        ;;
                z)
                        if [ "${doColumnsNames}" == "" ]; then doColumnsNames="${OPTARG}"; else doColumnsNames="${doColumnsNames}${param_ifs}$OPTARG"; fi
                        do_only_selected_columns=1
                        ;;
                Z)
                        if [ "${dontDoColumnsNames}" == "" ]; then dontDoColumnsNames="${OPTARG}"; else dontDoColumnsNames="${dontDoColumnsNames}${param_ifs}$OPTARG"; fi
                        do_only_selected_columns=1
                        ;;
		2)
			include_discarded_columns_in_output=1
			;;
		v)
			proportion_of_variance="$OPTARG"
			;;
		l)
			save_loadings_file="$OPTARG"
			;;
		L)
			load_loadings_file="$OPTARG"
			;;
		i)
			if [ "$inputFile" != "" ]; then
				echo "$appname : you can specify only one input file"
				exit 1
			fi
			inputFile="$OPTARG"
			;;
		o)
			if [ "$outputFile" != "" ]; then
				echo "$appname : you can specify only one output folder"
				exit 1
			fi
			outputFile=$OPTARG
			;;
		# standard arguments
		R)
			if [ "$Rscript_name" != "" ]; then
				echo "$appname : you can specify output R script name once"
				exit 1
			fi
			Rscript_name="$OPTARG"
			;;
		H)
			display_help_exit "$appname"
			;;
		h)
			# in case there is a header (column names at the first row of file)
			header="TRUE"
			;;
		e)
			separator="$OPTARG"; Rseparator=",sep='${separator}'"
			;;
		E)
			separatorOut=$OPTARG
			;;
		d)
			removeNLastComponents=$OPTARG
			;;
		D)
			removeNFirstComponents=$OPTARG
			;;
		p)
			outputPrecision=$OPTARG
			;;
		u)
			doReversePCA=1
			;;
		F)
			echo "-F is not implemented, you need to fix it."
			exit 1
			force_pca_even_if_numrows_less_numcols=1
			;;
		t)
			found=0;
			for((preID=0;preID<${#knownPreprocessingTransformations[@]};preID+=2)); do if [ $OPTARG = ${knownPreprocessingTransformations[preID]} ]; then preprocessingTransformation=${knownPreprocessingTransformations[preID+1]}; found=1; break; fi; done
			if [ $found -eq 0 ]; then echo "$appname : preprocessing transformation $OPTARG is not known, I only know about these ${knownPreprocessingTransformations[@]}"; exit 1; fi
			echo "$appname : doing preprocessing: ${preprocessingTransformation}"
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
if [ "$outputFile" = "" ]; then
	echo "$0 : you need to specify an output file using the '-o' option."
	exit 1
fi
if [[ $doReversePCA -eq 1 && "${load_loadings_file}" == "" ]]; then
	echo "$0 : reverse PCA (aka un-PCA , the -u flag) requires that loadings file of the PCA be supplied with -L option."
	exit 1
fi
if [[ "${save_loadings_file}" == "" && "${load_loadings_file}" == "" ]]; then
	echo "$0 : you need to specify a loadings file"
	exit 1
fi
if [[ "${load_loadings_file}" != "" && ${do_only_selected_columns} -eq 1 ]]; then
	echo "$0 : you can not load a loadings file and select specific columns, pca will be applied to the columns with names such us pca(...)"
	exit 1
fi
asked_proportion_of_variance=`echo "$proportion_of_variance<=1.0"|bc -l`
if [[ ${asked_proportion_of_variance} -eq 1 && $doReversePCA -eq 1 ]]; then
	echo "$0 : can not eliminate principal components (-v) based on variance and do reverse PCA at the same time (-u)."
	exit 1
fi
if [[ ${asked_proportion_of_variance} -eq 1 && "$removeNLastComponents" != "0" ]]; then
	echo "$0 : can not eliminate principal components based on variance and eliminate N last components at the same time (-d)."
	exit 1
fi
if [[ ${asked_proportion_of_variance} -eq 1 && "$removeNFirstComponents" != "0" ]]; then
	echo "$0 : can not eliminate principal components based on variance and eliminate N first components at the same time (-D)."
	exit 1
fi
if [[ (${do_only_selected_columns} -eq 1) && (${include_discarded_columns_in_output} -eq 1) && ("$removeNFirstComponents" != "0" || "$removeNLastComponents" != "0" || ${asked_proportion_of_variance} -eq 1) ]]; then
	echo "$0 : including discarded columns to the output (-2) is not allowed when you asked to remove first (-d) or last (-D) components or remove components based on variance (-v). You will make a pie of an output."
	exit 1
fi
#if [[ ("${doColumns}" != "" || "${doColumnsNames}" != "") && ("${dontDoColumns}" != "" || "${dontDoColumnsNames}" != "") ]]; then
#        echo "$0 : either -c or -C options but not both can be used (either specify columns to do, or columns not to do - both makes no sense)"
#        exit 1
#fi

if [ "${doColumnsNames}" != "" ]; then
        crap=$(column_names_to_column_numbers "${inputFile}" "${separator}" ${doColumnsNames})
        if [ $? -ne 0 ]; then echo ${crap}; exit 1; fi
        if [ "${crap}" != "" ]; then
                doColumns+="${param_ifs}${crap}"
                echo "$appname : extra column numbers to do via their column names : ${crap}"
        fi
fi
if [ "${dontDoColumnsNames}" != "" ]; then
        crap=$(column_names_to_column_numbers "${inputFile}" "${separator}" ${dontDoColumnsNames})
        if [ $? -ne 0 ]; then echo ${crap}; exit 1; fi
	if [ "${crap}" != "" ]; then
                dontDoColumns+="${param_ifs}${crap}"
                echo "$appname : extra column numbers not to do via their column names : ${crap}"
        fi
fi
doColumns=${doColumns#${param_ifs}}; doColumns=${doColumns##${param_ifs}}; dontDoColumns=${dontDoColumns#${param_ifs}}; dontDoColumns=${dontDoColumns##${param_ifs}}

# IN PCA it is the same to scale data yourself and then apply PCA or ask prcomp to apply PCA on data
#   but to do scaling first.
#   In the first case you apply the loadings as t(pca$rotation %*% t(scaled_data)) where
#   scaled_data = scale(data, center=T, scale=T) and pca=prcomp(data)
#   In the second case you apply the loadings to the transformed data as:
#   t(pca$rotation %*% t(scale(data, scale=pca$scale, center=pca$center))) and pca=prcomp(data,scale=T,center=T)

cat << EOC > "$Rscript_name"
# R script by $appname on $now
write("$appname : reading data from '$inputFile'", stderr())
tryCatch({oridata = as.matrix(read.table("$inputFile", header=$header $Rseparator, check.names=FALSE))}, error=function(ex){write(paste("$appname : error reading file $inputFile\n",ex,sep=""),stderr())});
write("$appname : done.", stderr())
actual_num_cols=ncol(oridata)
if( "${doColumns}" == "" ){
# no columns specified, do all of them
	columns_to_do = seq(from=1, to=actual_num_cols)
} else { columns_to_do = as.integer(unlist(strsplit("${doColumns}", "${param_ifs}"))) }
if( "${dontDoColumns}" != "" ){ columns_to_do = setdiff(columns_to_do, as.integer(unlist(strsplit("${dontDoColumns}", "${param_ifs}")))) }

num_columns_to_do = length(columns_to_do); num_rows = nrow(oridata)
if( num_rows==0 ){ write("$appname : no rows left to do", stderr()) }
if( num_columns_to_do==0 ){ write("$appname : no columns left to do", stderr()) }
if( (num_columns_to_do==0) || (num_rows==0) ){ quit(status=0) }

proportion_of_variance = $proportion_of_variance

$preprocessingTransformation # this will set variables doScale doCenter to T or F depending on what user selected (-t)

# the loadings and all other files contain as many items as whatever columns were PCAed
# the file .cols contains the columns which were PCAed as a vector
# if an input with more columns than max(.cols), then it assumes that it has to select only those PCAed, pca them (or reverse pca)
# and then return the result merged with all those input columns filtered out, keeping the same order
# if the input contains as many cols as the size of .cols then all of it will be used without any filtering out.
if( "$load_loadings_file" != "" ){
	# read from file the loadings
	# just read a fucking vector incheck.names=F, quote='"', header=T
	tryCatch({columns_to_do=read.table("${load_loadings_file}.cols", row.names=c(1), check.names=F, quote='"', header=F)[,1]}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}.cols\n",ex,sep=""),stderr())});
	m = max(columns_to_do); l = length(columns_to_do)
	if( l == actual_num_cols ){ columns_to_do = 1:actual_num_cols } else {
		if( m != actual_num_cols ){
			write(paste("$app_name : the number of dimensions (columns) of the input data is ", actual_num_cols, ".",
				"The previous PCA was done using ", length(columns_to_do), " dimensions (PCs) on a data with at least ", m, " dimensions.",
                	       	" Your input data must have either ", m, " (implies mixed PCAed and non-PCAed columns, filtering will be used from file '${load_loadings_file}.cols') or ",
				l, " (implies purely PCAed columns, no filtering at all).", sep=""), stderr())
			quit(status=1)
		}
	}

	tryCatch({loadings=as.matrix(read.table(file="$load_loadings_file", row.names=c(1), check.names=F, quote='"'))}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}\n",ex,sep=""),stderr())});
	colnames_of_pca_columns = rownames(loadings)
	tryCatch({preprocessing_transformation_matrix_scale=read.table(file="${load_loadings_file}.scale", row.names=c(1), check.names=F, quote='"', header=F)[,1]}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}.scale\n",ex,sep=""),stderr())});
	tryCatch({preprocessing_transformation_matrix_center=read.table(file="${load_loadings_file}.center", row.names=c(1), check.names=F, quote='"', header=F)[,1]}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}.center\n",ex,sep=""),stderr())});
	tryCatch({variances=read.table("${load_loadings_file}.var", row.names=c(1), check.names=F, quote='"', header=F)[,1]}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}.var\n",ex,sep=""),stderr())});
	tryCatch({cum_variances=read.table("${load_loadings_file}.cumvar", row.names=c(1), check.names=F, quote='"', header=F)[,1]}, error=function(ex){write(paste("$appname : error reading file ${load_loadings_file}.cumvar\n",ex,sep=""),stderr())});
	write("$appname : variance and cum.variance info for each component read from file ${load_loadings_file}.var/.cumvar", stderr())
	data = oridata[, columns_to_do]
	if ( $doReversePCA == 1 ){
		# un-pca, reverse it, first invert the loading matrix and then apply it to the data just like with the normal loadings matrix
		loadings = solve(loadings)
		# we have to un-scale and un-center, don't use scale because the order is the same, do it one step at the
		# time to reverse the original scale/center
		new_data=scale(
			scale(data %*% loadings, center=F, scale=1/preprocessing_transformation_matrix_scale),
			scale=F, center=-preprocessing_transformation_matrix_center
		)
	} else {
		transformed_data = scale(data, scale=preprocessing_transformation_matrix_scale, center=preprocessing_transformation_matrix_center)
		new_data = transformed_data %*% loadings
	}
} else {
	# do pca from scratch, no loadings file
	data = oridata[, columns_to_do]; names_of_pca_columns = colnames(data)
	if( (num_rows <= 2) || (nrow(data) < ncol(data)) ){
		if( ${force_pca_even_if_numrows_less_numcols} == 0 ){
			write(paste("$appname : sorry, but the number of data rows (", nrow(data), ") must be at least equal to the number of columns (", ncol(data), ") unless you specify the -F option, will not proceed with pca, but have not failed either."), stderr());
			quit(status=0)
		}
		# quit anyway, no sense to have less observations than fucking variables.
	}

	PC_colnames = paste("PC(", seq(from=1, to=num_columns_to_do, by=1), ")", sep="")
	pca_colnames = paste("pca(", colnames(oridata[, columns_to_do]), ")", sep="")

	if( doCenter == TRUE ){ center_param = TRUE; } else { center_param = rep(0.0, num_columns_to_do) }
	if( doScale == TRUE ){
		scale_param = sd(data)
		zero_variance_columns = which(scale_param==0.0);
		if( length(zero_variance_columns) > 0 ){
			# we have to make sure that those cols with no variation do not give an error when divided by stdev, set stdev = 1, crap but hack
			cn = colnames(data);
			write(paste("$appname : warning! column ", zero_variance_columns, " (",cn[zero_variance_columns],") has no variation", sep=""), stderr());
			scale_param[zero_variance_columns] = 1.0; # if a stdev is zero then set it to 1
		}
	} else { scale_param = rep(1.0, num_columns_to_do) }
	pca_data = prcomp(data, center=center_param, scale=scale_param, retx=TRUE)
	summar=summary(pca_data); sink(stderr()); print(summar)
	preprocessing_transformation_matrix_scale = pca_data\$scale
	preprocessing_transformation_matrix_center = pca_data\$center
	# it means we have not applied centering, create a zero vector and save it : is(pca_data\$scale, "logical") == TRUE )
	loadings = pca_data\$rotation
	variances = summar\$importance[2,]
	cum_variances = summar\$importance[3,]
	new_data = pca_data\$x

	# write auxiliary data to files
	write.table(loadings, file="$save_loadings_file", append=FALSE, row.names=pca_colnames, col.names=PC_colnames)
	write.table(summar\$importance[2,], file="${save_loadings_file}.var", row.names=PC_colnames, col.names=F)
	write.table(summar\$importance[3,], file="${save_loadings_file}.cumvar", row.names=PC_colnames, col.names=F)
	write.table(columns_to_do, file="${save_loadings_file}.cols", row.names=pca_colnames, col.names=F)
	write.table(preprocessing_transformation_matrix_scale, file="${save_loadings_file}.scale", row.names=pca_colnames, col.names=F)
	write.table(preprocessing_transformation_matrix_center, file="${save_loadings_file}.center", row.names=pca_colnames, col.names=F)
	write("$appname : loadings written in file '$save_loadings_file'. Scale and center vectors in .scale/.center. Columns which were PCAed in .cols. Variance and cumulative variance of PCs in .var/.cumvar. PCA summary in .pca_summary.", stderr());
	sink(file="${save_loadings_file}.pca_summary"); print(summar); sink(stderr());
}

if( (proportion_of_variance > 0) && (proportion_of_variance < 1) && ($doReversePCA != 1) ){
	# keep only those PCs which represent the specified cum variance

	num_components_initial = num_columns_to_do
	num_components_to_keep_total = 0
	index=0; num_components_to_keep_from = $removeNFirstComponents + 1; proportion_of_variance_kept = 0.0
	for(val in variances){
		index = index + 1
		# remove the first (most important)/last N components first and then see if we can find the variance specified, if not, tough luck
		if( index <= $removeNFirstComponents ) next;
		if( index > (num_columns_to_do-$removeNLastComponents) ) break;

		proportion_of_variance_kept = proportion_of_variance_kept + val
		num_components_to_keep_total = num_components_to_keep_total + 1
		if( proportion_of_variance_kept > proportion_of_variance ){ break }
	}
	num_components_to_keep_to = num_components_to_keep_from + num_components_to_keep_total - 1
	if( num_components_to_keep_total == 0 ){ write("$appname : error zero principal components were selected (check variance constraints or number of components selected using options -d -D and -v)", stderr()); quit(status=1); }
	if( num_components_initial > num_components_to_keep_total ){ save_data = new_data[,num_components_to_keep_from:num_components_to_keep_to]; write(paste("$appname : selecting ", num_components_to_keep_total, " components: ",num_components_to_keep_from," to ", num_components_to_keep_to," (incl.) which encompass ", (100*proportion_of_variance_kept), " % of the variance of the original data - initially there were ", num_components_initial, " components.", sep=""), stderr()) } else { save_data = new_data }

	# merge with the un-pcaed columns if requested but PCs will go at the end
	PC_colnames = paste("PC(", num_components_to_keep_from:num_components_to_keep_total, ")", sep="")
	pca_colnames = paste("pca(", colnames(oridata)[columns_to_do], ")", sep="")
	save_data = cbind(oridata[,-columns_to_do], save_data)
	colnames(save_data) <- c(colnames(oridata)[-columns_to_do], PC_colnames)
} else {
	if( $include_discarded_columns_in_output ){ save_data = oridata; save_data[,columns_to_do] = new_data; colnames(save_data)[columns_to_do] <- colnames(new_data) } else { save_data = new_data }
}

# precision and when this precision is not enough
options(digits=$outputPrecision, scipen=100)
write.table(save_data, file="$outputFile", sep="$separatorOut", col.names=$header, row.names=FALSE, append=FALSE, quote=FALSE)

write("$appname : done with precision $outputPrecision", stderr())
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
rm -f /tmp/rlog$rand.txt /tmp/$rand.* "$Rscript_name"
exit 0
