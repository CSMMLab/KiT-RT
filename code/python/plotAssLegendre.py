# importing the required module
import matplotlib.pyplot as plt
import csv
import numpy as np

'''
# x axis values
x = [1, 2, 3]
# corresponding y axis values
y = [2, 4, 1]

# plotting the points
plt.plot(x, y)

# naming the x axis
plt.xlabel('x - axis')
# naming the y axis
plt.ylabel('y - axis')

# giving a title to my graph
plt.title('My first graph!')

# function to show the plot
plt.show()
'''

# Read csv data
my = []
moments = list(list())

with open('assLegendre.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        my.append( row[0])
        moments.append( row[1:])

#transform to array
momentsArr = np.array( moments)
myArr = np.array(my)

# plot stuff
print(len(momentsArr[0,:]))
#for idx in range (0, len(momentsArr[0,:])):
   # print(idx)
# print("\n")
print(momentsArr[:, 2])
print(myArr)
plt.plot(myArr, momentsArr[:,2], label=1 )
print("------------------\n")
axes = plt.gca()
    #axes.set_xlim([-1, 1])
    #axes.set_ylim([-2, 2])

#plt.legend()
plt.show()

