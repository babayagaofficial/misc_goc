import pickle
from typing import cast
from plasnet.communities import Communities
from plasnet.output_producer import OutputProducer
from plasnet.base_graph import BaseGraph
from plasnet.list_of_graphs import ListOfGraphs
from plasnet.Templates import Templates
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import json
from pathlib import Path


communities_pickle = "/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/communities.pkl"
communities = cast(Communities, Communities.load(communities_pickle))

communities = communities.get_graphs_sorted_by_size()

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
hubs = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/hub_plasmids.csv")

graph = communities[0]

meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

plasmids = meta[meta["ST"].isin([131,95,73,69])]["Plasmid_ID"]

'''
with open("/home/daria/Documents/projects/ABC/host_presence/exclude_in_503.txt") as f:
    to_exclude = f.read().split("\n")

plasmids = [plasmid for plasmid in plasmids if not plasmid in to_exclude]
'''

subgraph = graph.subgraph(plasmids)



#subcomm = "community_0_subcommunity_503"

#subcomm_plasmids = list(set(typing[typing["type"]==subcomm]["plasmid"]).intersection(set(plasmids)))

subcomm_plasmids = ["30134_6#119_2"]

subcomm_boundary = list(nx.edge_boundary(subgraph, subcomm_plasmids))

low_dist_neighbours = [edge[1] for edge in subcomm_boundary if subgraph.edges[edge]['td']<5  and subgraph.edges[edge]['sd']<=0.5]

neighbours = list(set(typing[typing["plasmid"].isin(low_dist_neighbours)]["type"]))

neighbours = [neighbour for neighbour in neighbours if len(typing[typing["type"]==neighbour]["plasmid"])>11]
print(neighbours)

subcomm_graph = nx.Graph(subgraph.subgraph(subcomm_plasmids+low_dist_neighbours))


for edge in subcomm_graph.edges:
    if subcomm_graph.edges[edge]['td']>=12:
        subcomm_graph.remove_edge(edge[0],edge[1])

        
print(subcomm_graph.nodes.data())


data = nx.cytoscape_data(subcomm_graph)
with open(f"network_c0_sc503.json", "w") as file:
    json.dump(data, file)


blub = subcomm_graph

visualisation_src = Templates.read_template("visualisation_template")

graph_as_cy_dict = nx.cytoscape_data(blub)
elements_as_cy_json = json.dumps(graph_as_cy_dict["elements"])

filters = (
            f'<label for="hide_hubs">'
            f"Hide hub plasmids (1 present)"
            f"</label>"
            f'<input type="checkbox" id="hide_hubs" name="hide_hubs"><br/>'
        )
custom_buttons = '<div><input type="submit" value="Redraw" onclick="redraw()"></div>'

final_html_lines = []
for line in visualisation_src:
    line = line.replace("<libs_relative_path>", "..")
    line = line.replace("<samples_selectors>", "")
    line = line.replace("<elements_tag>", elements_as_cy_json)
    line = line.replace("<movementThreshold>", str(len(blub.edges())))
    line = line.replace("<maxSimulationTime>", str(10000))
    line = line.replace("<filters_tag>", filters)
    line = line.replace("<custom_buttons_tag>", custom_buttons)
    final_html_lines.append(line)


html = "\n".join(final_html_lines)



#html, json_dict = blub.produce_visualisation()
html_path = Path("/home/daria/Documents/projects/ABC/st131/graphs/network_hub.html")
html_path.write_text(html)
