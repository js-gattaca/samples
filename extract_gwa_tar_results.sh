# Usage: ./extract_results.sh
# Does four things
# 		a) Check whether all the tar files are present (N = 270)
#		b) If all tar files are present then extract them otherwise make a file 'missing_tar_files' and exit
# 		c) After extraction, check the number of extracted files (N=25606) if some files are missing them make 'missing_files' and exit
#		d) combine results for different chromosomes and only keep with pvalue <= 0.001 in new directory gwa_results_pval_cutoff
# Needs a directory called 'tar_files' with all the results
# Makes a new directory 'gwa_results' to store the the extracted results
# ------------------------------------------------------------------------------ #
# remove the missing files if they already exist
if [ -e missing_tar_files ];then rm -f missing_tar_files; fi
if [ -e missing_files ];then rm -f missing_files; fi

echo 'Checking whether all tar files are present'
# count the number of tar files; it should be 270; if not then make a file 'missing_tar_files' with missing files and exit
n_tar_files=`find tar_files/ -type f -name "*.tar.gz" | wc -l` # count the number of tar files

if [ $n_tar_files -eq 270 ];then
	echo 'All tar files are present.'
else
	m_tar_files=$((270 - $n_tar_files))
	echo $m_tar_files 'tar files are missing'
	echo 'Writing names of missing tar files to missing_tar_files file....'
	
	# make a file with all the missing tar files
	for file in `cat /Users/js/Desktop/lng/common_files/all_tar_files`;do
		check='tar_files/'$file
		if [ ! -e $check ];then
			echo $file >> missing_tar_files
		fi
	done
	echo 'Exiting as some tar files are  missing'
	exit 0
fi

# if all tar files are present then extract them
echo 'Extracting tar files and saving the results in gwa_results directory'
mkdir gwa_results # make a directory 

for file in tar_files/*_results.tar.gz;do
	dir=`echo $(basename $file) | cut -d. -f 1`
	tar -zxf $file -C tar_files 
	mv tar_files/$dir/* gwa_results
	rm -r tar_files/$dir	
done

# count the number of extracted files it should be 25606 otherwise make a file with names of missing files
echo 'Checking whether all files are present after extraction.'
n_files=`find gwa_results/ -type f -name "*.results" | wc -l` #count the number of files
if [ $n_files -eq 25606 ];then
	echo 'All extracted files are present.'
else
	m_files=$((25606 - $n_files))
	echo $m_files 'files are missing'
	echo 'Writing names of missing files to missing_files file....'
	
	# which files are missing
	while read line;do
		file=$(basename $line)
		fileName=$(echo $file | cut -d. -f 1)
		checkFile=$fileName'_'*'.results'
		if [ ! -e gwa_results/${checkFile} ];then
			echo $line >> missing_files
		fi
	done < /Users/js/Desktop/lng/common_files/all_imputed_filtered_files
	echo 'Exiting as some files are missing..'
	exit 0
fi

# combining GWA results from different chromosomes and only keep with pvalue <= 0.001
echo 'Combining results by chromosome...'
python /Users/js/Desktop/lng/common_files/combining_imputed_GWA_all_results.py

echo 'Done'
