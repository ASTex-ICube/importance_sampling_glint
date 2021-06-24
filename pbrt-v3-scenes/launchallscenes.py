import subprocess

process = subprocess.run(['pbrt','fig1_left.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig1_left.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig1_middle.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig1_middle.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig1_right.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig1_right.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_1.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_1.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_2.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_2.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_3.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_3.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_4.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_4.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_5.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_5.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_6.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_6.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_7.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_7.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_8.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_8.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_9.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_9.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_10.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_10.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_11.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_11.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig3_12.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig3_12.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig4_left.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig4_left.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig4_middle.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig4_middle.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig4_right.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig4_right.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig5_bottom.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig5_bottom.log", "w")
output_file.write(process.stdout)
output_file.close()

process = subprocess.run(['pbrt','fig5_top.pbrt'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
output_file = open("fig5_top.log", "w")
output_file.write(process.stdout)
output_file.close()