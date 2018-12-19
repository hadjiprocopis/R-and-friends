#!/bin/bash

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

# name, description, default distance metric
availableClusteringMethods=(	'kmeans_hw' 'Hartigan-Wong (1979) (default)' 'n/a'
				'kmeans_lloyd'	'Lloyd (1957) kmeans' 'n/a'
				'kmeans_forgy' 'Forgy (1965) kmeans flavour' 'n/a'
				'kmeans_macqueen' 'MacQueen (1967) kmeans flavour' 'n/a'
				'pam' 'robust but really slow and intensive clustering (R package: cluster). Its output is mostly different to the other kmeans-based methods.' 'euclidean'
				'amap' 'kmeans from the amap R package, its output is comparable to kmeans_* methods' 'euclidean'
			)
# algorith it applies to, name, description
availableDistanceMetrics=(	'amap'	'euclidean'	'euclidean (default)'
				'amap'	'maximum'	'maximum'
				'amap'	'manhattan'	'manhattan'
				'amap'	'canberra'	'canberra'
				'amap'	'binary'	'binary'
				'amap'	'pearson'	'pearson'
				'amap'	'correlation'	'correlation'
				'amap'	'spearman'	'spearman'
				'amap'	'kendall'	'kendall'

				'pam'	'euclidean'	'euclidean (default)'
				'pam'	'manhattan'	'manhattan'
			)

function display_help_exit {
	echo "Usage : $1 -i input_file -o outputFile -c col [-r from_row:to_row] [-x minX:maxX] [-y minY:maxY] [-h] [-t max_iterations] [-f] [-e input_separator] [-E output_separator]"
	echo "This script will read a N-dimensional (=N columns) data from input_file and cluster it"
	echo " -i input_file"
	echo "        defines the input file."
	echo " -o output_file"
	echo "		where the results go"
	echo " -S seed"
	echo "          seed to the random number generator"
	echo " -c col"
	echo "	      choose the column number you wish to fit (use -c option repeatedly for more columns)"
	echo " -h"
	echo "        optional, that the first row is not numerical but it holds column names"
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file, default is white space"
	echo "        which is '' (i.e. nothing in between the quotes)"
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs)"
	echo " -s"
	echo "        scale data (unit variance and zero mean)"
	echo " -t max_iterations"
	echo "        maximum number of iterations to do (does not apply to 'pam' method which does not take any iterations input)"
	echo " -f"
	echo "        experiment clustering for number of clusters equal to 2 to N (specify N with -N option) and plot them"
	echo "        from http://www.statmethods.net/advstats/cluster.html:"
	echo "        A plot of the within groups sum of squares by number of clusters extracted can help determine the appropriate number of clusters. The analyst looks for a bend in the plot similar to a scree test in factor analysis."
	echo " -m method"
	echo "        method to use:"
	for((i=0;i<${#availableClusteringMethods[@]};i+=2)); do echo ${availableClusteringMethods[i]} :  ${availableClusteringMethods[i+1]}; done
	echo " -M distance_metric"
	echo "        metric to use (for some algorithms):"
	for((i=0;i<${#availableDistanceMetrics[@]};i+=2)); do echo ${availableDistanceMetrics[i]} :  ${availableDistanceMetrics[i+1]}; done
	echo ""
	echo "	      Example:"
	echo "$app_name -i a.txt -o abc -h"
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

inputFile=""
outputFile=""
col1=""
separator=""
separatorOut="\t"
header="FALSE"
app_name="$0"
rand="$$"
scaleData=""
Rscript_name="${rand}.R"
findNumberOfClustersAndExit="FALSE"
clusteringMethod=""
distanceMetric=""
numClusters=0
limits=""
Seed="${rand}"
maxIterations=100

if [ "$*" = "" ]; then
	display_help_exit $app_name
fi

app_name=`basename $0`
now=`date`

while getopts "i:o:c:hR:e:E:HsfN:S:m:M:t:" OPTION
do
	case $OPTION in
		i)
			inputFile="$OPTARG"
			;;
		o)
			outputFile=$OPTARG
			;;
		R)
			Rscript_name="$OPTARG"
			;;
		c)
			col1="$col1,$OPTARG"
			;;
		H)
			display_help_exit "$app_name"
			;;
		h)
			# in case there is a header (column names at the first row of file)
			header="TRUE"
			;;
		s)
			scaleData="TRUE"
			;;
		f)
			findNumberOfClustersAndExit="TRUE"
			;;
		e)
			separator="$OPTARG"
			;;
		m)
			clusteringMethod="$OPTARG"
			;;
		M)
			distanceMetric="$OPTARG"
			;;
		t)
			maxIterations="$OPTARG"
			;;
		N)
			numClusters="$OPTARG"
			;;
		E)
			separatorOut="$OPTARG"
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
if [ "$numClusters" -le 1 ]; then
	echo "$0 : you need to specify the number of clusters using the '-N' option - and this number must be an integer greater than 1."
	exit 1
fi
suc1=-1
for((cm=0;cm<${#availableClusteringMethods[@]};cm+=3)); do if [ "${availableClusteringMethods[cm]}" == "$clusteringMethod" ]; then suc1=${cm}; break; fi; done
if [ $suc1 -lt 0 ]; then echo "$0 : invalid clustering method ${clusteringMethod}, you need to specify a clustering method (-m) from the following"; for((i=0;i<${#availableClusteringMethods[@]};i+=3)); do echo ${availableClusteringMethods[i]} :  ${availableClusteringMethods[i+1]}; done; exit 1; fi
if [ "$distanceMetric" == "" ]; then distanceMetric=${availableClusteringMethods[suc1+2]}; fi

suc2=-1
for((cm=0;cm<${#availableDistanceMetrics[@]};cm+=3)); do if [ "${availableDistanceMetrics[cm]}" == "${clusteringMethod}" ] && [ "${availableDistanceMetrics[cm+1]}" == "${distanceMetric}" ]; then suc2=${cm}; break; fi; done
if [ $suc2 -lt 0 ]; then echo "$0 : invalid distance metric ${distanceMetric} specified for clustering method ${clusteringMethod}, will set it to default (${availableClusteringMethods[suc1+2]})."; distanceMetric=${availableClusteringMethods[suc1+2]}; fi

echo "$0 : using ${clusteringMethod} clustering method and ${distanceMetric} distance metric"

col1=${col1/,/}
if [ ${Seed} = "" ]; then Seed="${rand}"; fi

cat << EOC > "$Rscript_name"
# R script by $app_name on $now
library(cluster)
write(paste("$app_name : random number generator seed : ", $Seed, sep=""), file="$outputFile", append=FALSE);
set.seed($Seed)

tryCatch({ data = read.table("$inputFile", sep="$separator", header=$header)}, error=function(ex){write(paste("$app_name : error reading file $inputFile\n",ex,sep=""),stderr())});
dataAllDimensions = as.matrix(data)
dataSelectedColumn = as.matrix(data[,$col1])
num_rows_data = length(dataSelectedColumn)

dataIn = dataSelectedColumn
if( "$scaleData" != "" ){
	dataIn = scale(dataSelectedColumn, scale=T, center=T)
}

if( "$clusteringMethod" == "pam" ){
	if( "$findNumberOfClustersAndExit" == "TRUE" ){
		wss <- (nrow(dataIn)-1)*sum(apply(dataIn,2,var))
		for (i in 2:$numClusters){ wss[i] <- sum(pam(dataIn, $numClusters, metric="${distanceMetric}")\$withinss) }
		png(filename="$outputFile.png")
		plot(1:$numClusters, wss, type="b", xlab="Number of Clusters (method:$clusteringMethod, dist:$distanceMetric, iters:$maxIterations)", ylab="Within groups sum of squares")
		dev.off()
		write("$app_name : a plot of within groups sum of squares by number of clusters (max $numClusters) extracted is written in $outputFile.png", stderr())
		quit(status=0)
	}

	# K-Means Cluster Analysis
	write("$app_name : executing pam clustering with ${distanceMetric} distance metric, for $numClusters clusters", stderr())
	fit <- pam(dataIn, $numClusters, metric="${distanceMetric}")
	write("$app_name : done!", stderr())
	# get cluster means
	aggregate(dataIn,by=list(fit\$clustering),FUN=mean)
	# append cluster assignment
	dataOut <- data.frame(dataIn, fit\$clustering)
	centers = fit\$medoids
	sizes = fit\$size
	write.table(file="$outputFile.centers", fit\$medoids, append=F, row.names=F, col.names=F, sep="$separatorOut")
	write.table(file="$outputFile.sizes", fit\$size, append=F, row.names=F, col.names=F, sep="$separatorOut")
	write("$app_name : input data with cluster assigned to each row (last column) written in $outputFile (cluster centers in $outputFile.centers, cluster sizes in $outputFile.sizes and all-together info in $outputFile.info)", stderr())
	info=c()
	ncl=dim(fit\$medoids)[1]; ndi=dim(fit\$medoids)[2]
	for( i in 1:ncl){
		crap = ""
		for(j in 1:ndi){
			crap = paste(crap, " ", fit\$medoids[i,j], sep="")
		}
#		info[i] = paste("center[", i, "]=(", crap, ") withinss[", i, "]=(", fit\$withinss[i], ") size[", i, "]=(", fit\$size[i], ")", sep="")
		info[i] = paste("center[", i, "]=(", crap, ") size[", i, "]=(", fit\$clusinfo[i,1], ")", sep="")
	}
	write("$app_name : centers:", stderr())
	write(fit\$medoids, stderr())
	write("$app_name : other info:", stderr())
	write(fit\$clusinfo, stderr())
} else
if( ("$clusteringMethod" == "kmeans_hw") || ("$clusteringMethod" == "kmeans_lloyd") || ("$clusteringMethod" == "kmeans_forgy") || ("$clusteringMethod" == "kmeans_macqueen") ){
	if( "$findNumberOfClustersAndExit" == "TRUE" ){
		wss <- (nrow(dataIn)-1)*sum(apply(dataIn,2,var))
		for (i in 2:$numClusters){ wss[i] <- sum(kmeans(dataIn, centers=i, algorith=algorithm)\$withinss) }
		png(filename="$outputFile.png")
		plot(1:$numClusters, wss, type="b", xlab="Number of Clusters (method:$clusteringMethod, dist:$distanceMetric, iters:$maxIterations)", ylab="Within groups sum of squares")
		dev.off()
		write("$app_name : a plot of within groups sum of squares by number of clusters (max $numClusters) extracted is written in $outputFile.png", stderr())
		quit(status=0)
	}
	# K-Means Cluster Analysis
	algorith=""
	if("$clusteringMethod" == "kmeans_hw") algorithm="Hartigan-Wong"
	else if( "$clusteringMethod" == "kmeans_lloyd") algorithm="Lloyd"
	else if( "$clusteringMethod" == "kmeans_forgy") algorithm="Forgy"
	else if( "$clusteringMethod" == "kmeans_macqueen") algorithm="MacQueen"
	write(paste("$app_name : doing ",algorithm," (kmeans) clustering", sep=""), stderr())
	
	fit <- kmeans(dataIn, $numClusters, algorith=algorithm)
	# get cluster means
	aggregate(dataIn,by=list(fit\$cluster),FUN=mean)
	# append cluster assignment
	dataOut <- data.frame(dataIn, fit\$cluster)
	centers = fit\$centers
	sizes = fit\$size
	write.table(file="$outputFile.centers", fit\$centers, append=F, row.names=F, col.names=F, sep="$separatorOut")
	write.table(file="$outputFile.sizes", fit\$size, append=F, row.names=F, col.names=F, sep="$separatorOut")
	write("$app_name : input data with cluster assigned to each row (last column) written in $outputFile (cluster centers in $outputFile.centers, cluster sizes in $outputFile.sizes and all-together info in $outputFile.info)", stderr())
	info=c()
	ncl=dim(fit\$centers)[1]; ndi=dim(fit\$centers)[2]
	for( i in 1:ncl){
        	crap = ""
	        for(j in 1:ndi){
        	        crap = paste(crap, " ", fit\$centers[i,j], sep="")
	        }
        	info[i] = paste("center[", i, "]=(", crap, ") withinss[", i, "]=(", fit\$withinss[i], ") size[", i, "]=(", fit\$size[i], ")", sep="")
	}
	write("$app_name : centers:", stderr())
	write(fit\$centers, stderr())
	write("$app_name : sizes:", stderr())
	write(fit\$size, stderr())
} else
if( "$clusteringMethod" == "amap" ){
	# K-Means Cluster Analysis
	library(amap)

	if( "$findNumberOfClustersAndExit" == "TRUE" ){
		wss <- (nrow(dataIn)-1)*sum(apply(dataIn,2,var))
		for (i in 2:$numClusters){ wss[i] <- sum(Kmeans(dataIn, centers=$numClusters, iter.max=$maxIterations, method="$distanceMetric")\$withinss) }
		png(filename="$outputFile.png")
		plot(1:$numClusters, wss, type="b", xlab="Number of Clusters (method:$clusteringMethod, dist:$distanceMetric, iters:$maxIterations)", ylab="Within groups sum of squares")
		dev.off()
		write("$app_name : a plot of within groups sum of squares by number of clusters (max $numClusters) extracted is written in $outputFile.png", stderr())
		quit(status=0)
	}

	fit <- Kmeans(dataIn, centers=$numClusters, iter.max=$maxIterations, method="$distanceMetric")
	# get cluster means
	aggregate(dataIn,by=list(fit\$cluster),FUN=mean)
	# append cluster assignment
	dataOut <- data.frame(dataIn, fit\$cluster)
	centers = fit\$centers
	sizes = fit\$size
	info=c()
	ncl=dim(fit\$centers)[1]; ndi=dim(fit\$centers)[2]
	for( i in 1:ncl){
        	crap = ""
	        for(j in 1:ndi){
        	        crap = paste(crap, " ", fit\$centers[i,j], sep="")
	        }
        	info[i] = paste("center[", i, "]=(", crap, ") withinss[", i, "]=(", fit\$withinss[i], ") size[", i, "]=(", fit\$size[i], ")", sep="")
	}
	write("$app_name : centers:", stderr())
	write(fit\$centers, stderr())
	write("$app_name : sizes:", stderr())
	write(fit\$size, stderr())
} else {
	write("$app_name : error, unknown clustering method $clusteringMethod", stderr())
	quit(status=1)
}

if( "$header" == "TRUE" ){
	# add some headers to the output
	names_of_data = names(data)
	names_str = paste('"', names_of_data[1], '"', sep="")
	for(i in 2:length(names_of_data)){ quoted_name = paste('"', names_of_data[i], '"', sep=""); names_str = paste(names_str, quoted_name, sep="\t") }
	info = c(paste("names=(",names_str, ")", sep=""), info)
	# the cluster ID is assigned as the last column of dataout, now change the header of that column
	names(dataOut)[length(dataOut)] = "AssignedClusterID"
	centers = data.frame(1:ncl, centers); names(centers)[1] = "AssignedClusterID"
}

write.table(file="$outputFile.centers", centers, append=F, row.names=F, col.names=$header, quote=T, sep="$separatorOut")
write.table(file="$outputFile.sizes", sizes, append=F, row.names=F, col.names=F, quote=F, sep="$separatorOut")
write.table(file="$outputFile", dataOut, append=F, row.names=F, col.names=$header, quote=T, sep="$separatorOut")
write.table(file="$outputFile.info", info, append=F, row.names=F, col.names=F, quote=F, sep="$separatorOut")
write("$app_name : input data with cluster assigned to each row (last column) written in $outputFile (cluster centers in $outputFile.centers, cluster sizes in $outputFile.sizes and all-together info in $outputFile.info)", stderr())

quit(status=0)
EOC

echo "$app_name : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
        cat /tmp/rlog$rand.txt; echo "$appname : error executing R script ($Rscript_name)"
        echo "Command line:"
        echo "$0 ${*}"
        exit 1
fi      
echo "Command line:"
echo "$0 ${*}"
#if [ "$header" == "TRUE" ]; then
#	head -1 "$inputFile" > /tmp/$rand.1 && cat /tmp/$rand.1 "$outputFile.centers" > /tmp/$rand.2; mv /tmp/$rand.2 "$outputFile.centers"
#fi
echo "$0 : output files: $outputFile.info, $outputFile.centers, $outputFile.sizes and labelled inputs as $outputFile"
rm -f /tmp/$rand.* /tmp/rlog"$rand".txt "${Rscript_name}"
exit 0
