import pypangraph as pp
import pandas as pd

def get_adj(neighbour, block):
    if (neighbour[1] and block[1]) or (not neighbour[1] and not block[1]):
        adj = (neighbour[0], 't')
    else:
        adj = (neighbour[0], 'h')
    return adj

graph = pp.Pangraph.from_json("/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/st131_cl416_community_0_subcommunity_503/pangraph.json")

path_dict = graph.to_path_dictionary()

adjacencies = {}

for plasmid in path_dict.keys():
    adjacencies[plasmid]=[]
    path = path_dict[plasmid]
    for i in range(len(path)-1):
        radj = get_adj(path[i+1], path[i])
        adjacencies[plasmid].append(((path[i][0],'h'),radj))
    radj = get_adj(path[0], path[len(path)-1])
    adjacencies[plasmid].append(((path[len(path)-1][0],'h'),radj))
    #print(adjacencies[plasmid])

all_adjacencies = list(set([el for plasmid in adjacencies.keys() for el in adjacencies[plasmid]]))
presence_absence = {}
for adjacency in all_adjacencies:
    presence_absence[adjacency]=[]
    for plasmid in adjacencies.keys():
        if adjacency in adjacencies[plasmid]:
            presence_absence[adjacency].append(True)
        else:
            presence_absence[adjacency].append(False)

df = pd.DataFrame(data=presence_absence)

for adjacency in all_adjacencies:
    if df[adjacency].all():
        print(adjacency)

bl_count = graph.to_blockcount_df()
block_PA = bl_count > 0

for plasmid in path_dict.keys():
    path = path_dict[plasmid]
    new_path = []
    new_blocks = {}
    new_block = 0
    for i in range(len(path)-1):
        radj = get_adj(path[i+1], path[i])
        if df[((path[i][0], 'h'),radj)].all():
            if path[i][0] not in new_blocks.keys():
                new_block = new_block + 1
                new_path.append((new_block, True))
            new_blocks[path[i+1][0]]=new_block
        else:
            new_path.append(path[i])
    print(new_path)

print(new_blocks)
    
                
