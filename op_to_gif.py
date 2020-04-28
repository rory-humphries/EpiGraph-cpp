import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio


files = []
for filename in os.listdir("./output/"):
    files += [filename]
    
data = pd.read_csv("./output/0.csv")
X = np.array(data["x"])
Y = np.array(data["y"])

base = np.zeros((max(Y)+1, max(X)+1))
bg = np.ones((max(Y)+1, max(X)+1))
pic = np.zeros((max(Y)+1, max(X)+1))
for i in data.values:
    base[int(i[1])][int(i[0])] = 1
base = base[::-1]

s=[]
i=[]
r=[]
x=[]

maxs = 0
maxi = 0
maxr = 0
maxx = 0
for filename in range(1, len(files)):
    df = pd.read_csv("./output/"+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I))
    maxr = max(maxr, max(df.R))
    maxx = max(maxx, max(df.X))
    s.append(sum(df.S))
    i.append(sum(df.I))
    r.append(sum(df.R))
    x.append(sum(df.X))
plt.plot(s, label = 'S')
plt.plot(i, label = 'I')
plt.plot(r, label = 'R')
plt.plot(x, label = 'X')
plt.title('Evolution of the SIRX model')
plt.ylabel('Number of individuals')
plt.xlabel('Time')
plt.savefig('evolution.png', bbox_inches='tight')

for filename in range(1, len(files)):
    df = pd.read_csv("./output/"+str(filename)+".csv")
    fig = plt.figure(figsize=(10,10))
    for i in df.values: 
        pic[int(i[1])][int(i[0])] = i[3]
    pic = np.ma.masked_where(pic < 1, pic)
    #print(max(df.I))
    base = np.ma.masked_where(base == 0, base)
    plt.imshow(bg[::-1], cmap = plt.get_cmap('inferno'),alpha = 0.5,origin='lower')
    plt.imshow(base, cmap = plt.get_cmap('inferno'),alpha = 0.7,interpolation="nearest",origin='lower')

    plt.imshow(pic[::-1], alpha = 1,  cmap = plt.get_cmap('inferno'), aspect="equal",interpolation="nearest", vmin = 1, vmax = 1000)
    #norm=matplotlib.colors.LogNorm() plt.scatter(df["East"],df["North"], marker = ".",s=1,c = df.I, , cmap = plt.get_cmap('inferno'))
    #norm=colors.LogNorm(),
    plt.title("Day " + str(filename))
    plt.colorbar().set_label("% of Infected Population per Square Kilometre")
    
    plt.grid(False)
    plt.savefig("./to_gif/"+ str(filename), bbox_inches='tight')
    #plt.show()
    plt.close()
    

    
    
files = []
for filename in os.listdir("./to_gif/"):
    if (".png" not in filename):
        continue
    files += [filename]
images = []

for filename in range(1,len(files)):
    images.append(imageio.imread("./to_gif/"+str(filename)+".png"))
imageio.mimsave('output.gif', images)
    
