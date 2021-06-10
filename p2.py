import matplotlib.pyplot as plt
import os
 
 
def run(array_size):
   run_command = "gcc -o simple -fopenmp b.c -lm"
   if os.system(run_command) == 0:
       r_c = "./simple "+str(array_size)
       os.system(r_c)
 
 
 
 
def read_file():
   serial = []
   parallel = []
   a_size = []
   f = open("time2.txt")
   count = 0
   for line in f:
       line = line[:-1]
       if count == 0:
           a_size.append(int(line))
       elif count == 1:
           parallel.append(float(line))
       elif count == 2:
           serial.append(float(line))
       count += 1
       if count  == 3:
           count = 0
   plt.plot(a_size,serial,label = "serial")
   plt.plot(a_size,parallel,label = "parallel")
   plt.legend()
   plt.xlabel('array size')
   plt.ylabel('time taken')
   filename_ = "openmp"
   plt.savefig(filename_+".png")
   plt.close()
 
 
def run_file(sizes):
   for size in sizes:
       run(size)
   read_file()
 
 
 
 
a_sizes = [i for i in range(1,1000,10)]
run_file(a_sizes)