# combining GWA results from different split files

import os
import glob
import csv

# make a new directory
os.makedirs("gwa_results_chromo")

def combining_gwa(chromo, outfile):
	f_out = open(outfile, 'w')
	csv_file_writer = csv.writer(f_out)
	csv_file_writer.writerow(['snp', 'Pval', 'Beta', 'SE']) # print out the header

	chromo_files = glob.glob('gwa_results/lng_hyb_c' + str(chromo) + '_*_*.results') # gets all files for a chromsome as a list
	for read_one_file in chromo_files:
		f_in = open(read_one_file)
		one_file = csv.reader(f_in, delimiter = '\t')

		next(one_file) # skip the first row
		for row in one_file:
			if len(row) < 5: continue
			if row[4] != ''	:
				csv_file_writer.writerow([row[0], row[4], row[1], row[2]])
		
		f_in.close()
	f_out.close()

def main():
	for chromo in range(1,23):
		outfile = os.path.join('gwa_results_chromo', 'chr' + str(chromo) + '_gwa_out.csv')
		combining_gwa(chromo, outfile)

if __name__ == '__main__':
	main()