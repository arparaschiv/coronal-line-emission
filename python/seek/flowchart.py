import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pylab import *



fig=figure(figsize=[9,8])
plt.rc('text',usetex='True')

#plt.figure(figsize = (12,9))

t1=r'Observed \\ \hphantom{ab} Stokes data\\ \hphantom{ab} $IQUV$ \\ \hphantom{ab} at each \\ \hphantom{ab}  elongation $y$'
t2=r'Preprocess $IQUV$\\  \hphantom{ab} integrate and \\ \hphantom{ab} fit components, \\ \hphantom{ab} $\div$  by max($I$)'
t3=r'Rotate $QU$ by  $\alpha$\\ \hphantom{ab}around LOS \\ \hphantom{ab} ($x$-axis)'
t4=r'Select database \\ \hphantom{ab} for given $y$'
t5=r'Find acceptable \textbf{b}\\ \hphantom{ab} vectors from\\ \hphantom{ab} matches in\\ \hphantom{ab} database'
t6=r'Rotate \textbf{b} by -$\alpha$ \\ \hphantom{ab} around LOS \\ \hphantom{ab} ($x$-axis)'
t7=r'Scale $B$ with $V$ \\ \hphantom{ab} return \textbf{B}(x,y,z)\\ \hphantom{ab} including all\\ \hphantom{ab} degenerate\\ \hphantom{ab}  solutions'
From = [t1,t2,t3,t4,t5,t6,t7]        
To =   [t2,t3,t4,t5,t6,t7,t7]

df = pd.DataFrame({ 'from':From,
                   'to':To})
# Define Node Positions
pos = {t1:(1,3),
       t2:(2,3),
       t3:(3,3),
       t4:(2,2),
       t5:(1,1),
       t6:(2,1),
       t7:(3,1)}

# Define Node Colors
NodeColors = {t1:[1,0,0],
              t2:[1,1,0],
              t3:[0,1,1],
              t4:[1,0,0],
              t5:[1,0,0],
              t6:[1,0,0],
              t7:[1,0,0]}

Labels = {}
i = 0
for a in From:
    Labels[a]=a
    i +=1
Labels[To[-1]]=To[-1]

# Build your graph. Note that we use the DiGraph function to create the graph! This adds arrows
G=nx.from_pandas_edgelist(df, 'from', 'to', create_using=nx.DiGraph() )

# Define the colormap and set nodes to circles, but the last one to a triangle
Circles = []
Triangle = []
Colors = []
for n in G.nodes:
    if n != 'Everyone$^{Dies}$':
        Circles.append(n)
    else:
        Triangle.append(n)
    Colors.append(NodeColors[n])

# By making a white node that is larger, I can make the arrow "start" beyond the node
nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Circles,
                       node_size=.95e4,
                       verts=((-0.3,0.2),(0.5,0.2),(0.4,-0.2),(-0.4,-0.2)),
                       node_shape='s',
                       node_color='gray',
                       alpha=1)

nodes = nx.draw_networkx_nodes(G, pos,  
                       nodelist = Circles,
                       node_size=.8e4,
                       verts=((-0.3,0.2),(0.5,0.2),(0.4,-0.2),(-0.4,-0.2)),
                       node_shape='s',
                       node_color='white')
#                       edgecolors='black')
#,
#                       alpha=0.5)


nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Triangle,
                       node_size=1.35e4,
                       node_shape='>')
#                       node_color='white',
#                       alpha=1)

nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Triangle,
                       node_size=1.1e4,
                       node_shape='>')
#                       node_color=Colors,
#                       edgecolors='black',
#                       alpha=0.5)


nx.draw_networkx_labels(G, pos, Labels, font_size=12)

# Again by making the node_size larer, I can have the arrows end before they actually hit the node
edges = nx.draw_networkx_edges(G, pos, node_size=1.8e4,
                               arrowstyle='->',width=2,arrowsizes=10)

plt.xlim(0,4.5)
plt.ylim(0,4)
plt.axis('off')
plt.tight_layout()

savefig("flowc.pdf")
import os

plt.close('all')
os.system("open flowc.pdf")

 
