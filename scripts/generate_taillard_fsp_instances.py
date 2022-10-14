# generate_taillard_fsp_instances.py  ************************************************************
# Generates a problem instance file (txt) for each instance listed in Flow shop sequencing txt 
#  files from Taillard website.
# References: http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/ordonnancement.html
# http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/flowshop.dir/tai20_5.txt

import sys, getopt
import os
import os.path

def main(argv):

   instance_path = ''
   
   print('Taillard execution script\n')
   print('===============================\n')
   
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["instance_path="])
   except getopt.GetoptError:
      print('generate_taillard_fsp_instances.py --instance_path <instance_path>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('generate_taillard_fsp_instances.py --instance_path <instance_path>')
         sys.exit()
      elif opt in ("-p", "--instance_path"):
         instance_path = arg

   if(instance_path == ''):
      print('Please specify the fsp instances path! Aborting.')
      sys.exit()
   
   print('Instances dir is ' + instance_path + '\n')
      
   process_instance_dir(instance_path)
   
   
  
def process_instance_dir(instance_path):

	# List of taillard txt files to be processed (in the designated order)
	file_list = ['tai20_5.txt', 'tai20_10.txt', 'tai20_20.txt', 'tai50_5.txt', 'tai50_10.txt', 'tai50_20.txt', 'tai100_5.txt', 'tai100_10.txt', 'tai100_20.txt', 'tai200_10.txt', 'tai200_20.txt', 'tai500_20.txt']
	
	# For each input file, create 10 instance files, numbered in order
	txt_count = 1
	instance_count = 1
	prefix = 'tail0'
	folder = instance_path + '\\output'
	# create folder
	if not os.path.exists(folder):
		os.makedirs(folder)
	for input_path in file_list:
		print('Processing txt file ' + str(txt_count) + ' of ' + str(len(file_list)) + ': ' + str(input_path) + '...')
		input_file = open(instance_path + '\\' + input_path, "r")
		lines = input_file.readlines()
		print(str(lines))
		line_count = 0
		headers = []
		matrices = []
		while line_count < len(lines):
			if lines[line_count].find('number') >= 0:
				# ignore this line and capture the next line (number of jobs, number of machines)
				line_count += 1
				data = lines[line_count].split()
				n = data[0]
				m = data[1]
				seed = data[2]
				ub = data[3]
				lb = data[4]
				print('Captured header: ' + str(data))
				headers.append(n + '\t' + m + '\t' + seed + '\t' + ub + '\t' + lb + '\n')
				line_count += 1
			elif lines[line_count].find('processing') >= 0:
				# ignore this line and capture the next lines, until 'number' is found again
				line_count += 1
				matrix = ''
				while line_count < len(lines) and lines[line_count].find('number') < 0:
					matrix += lines[line_count].strip(' ')
					line_count += 1
				print('Captured matrix: ' + str(matrix))
				matrices.append(matrix)
			# end if
		# end while
		
		for i in range(0, len(matrices)):
			instance_file = open(folder + '\\' + prefix + str(instance_count) + '.in', "w")
			instance_file.write(str(headers[i]))
			instance_file.write(str(matrices[i]))
			instance_file.close()
			instance_count += 1
		
		input_file.close()
		txt_count += 1
		
	
	
def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
	
if __name__ == "__main__":
   main(sys.argv[1:])
