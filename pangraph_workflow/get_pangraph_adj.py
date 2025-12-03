import pypangraph as pp

graph = pp.Pangraph.from_json("/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/st131_cl416_community_0_subcommunity_503/pangraph.json")

with open("st131_cl416_community_0_subcommunity_503_adj.txt", "w") as f:
    f.write("#Species\tGene_1\tExt_1\tSpecies\tGene_2\tExt_2\tWeight\n")
    for plasmid in graph.paths.keys():
        path = graph.paths[plasmid].nodes
        for i in range(len(path)-1):
            if graph.nodes[path[i]]["strand"]:
                ext_1 = "t"
            else:
                ext_1 = "h"
            if graph.nodes[path[i+1]]["strand"]:
                ext_2 = "h"
            else:
                ext_2 = "t"
            gene_1 = graph.nodes[path[i]]["block_id"]
            gene_2 = graph.nodes[path[i+1]]["block_id"]
            f.write(f"{plasmid}\t{gene_1}\t{ext_1}\t{plasmid}\t{gene_2}\t{ext_2}\t1\n")
        if graph.nodes[path[len(path)-1]]["strand"]:
            ext_1 = "t"
        else:
            ext_1 = "h"
        if graph.nodes[path[0]]["strand"]:
            ext_2 = "h"
        else:
            ext_2 = "t"
        gene_1 = graph.nodes[path[len(path)-1]]["block_id"]
        gene_2 = graph.nodes[path[0]]["block_id"]
        f.write(f"{plasmid}\t{gene_1}\t{ext_1}\t{plasmid}\t{gene_2}\t{ext_2}\t1\n")

