import pyfastaq
from operator import itemgetter
import numpy as np

def random_colour_codes(num):
    colour_codes = []
    for i in range(256):
        for j in range(256):
            for k in range(256):
                colour_codes.append([i,j,k])
    rng = np.random.default_rng(seed=42)
    all_colours = rng.choice(colour_codes, num,replace=False)
    return all_colours


def read_in_unimog(path):
    seqs = {}
    with open(path, "r") as file:
        for line in file:
            if line[0] == ">":
                genome = line[1:].strip("\n")
            else:
                seq = line.split(" ")
                seq.remove(")\n")
                seqs[genome] = [int(el) for el in seq]
    return seqs

def read_in_accessory(path):
    with open(path) as file:
        accessory = file.read().split(" ")
    accessory.remove('')
    accessory = [int(el) for el in accessory]
    colours = {}
    abs_accessory = list(set([abs(el) for el in accessory]))
    colour_codes = random_colour_codes(len(abs_accessory))
    colours = {}
    for i in range(len(abs_accessory)):
        hexcode = [hex(colour_codes[i][j]).replace('0x','') for j in range(3)]
        colours[abs_accessory[i]] = '#'+hexcode[0]+hexcode[1]+hexcode[2]
        colours[-abs_accessory[i]] = '#'+hexcode[0]+hexcode[1]+hexcode[2]
    return accessory, colours

# coords = list of tuples [(x1, y1), (x2, y2) ...]
def svg_polygon(coords, fill_colour, border_colour, border_width = 1, opacity=-1):
    return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
        + 'fill="' + fill_colour + '" '

    if opacity != -1:
        return_string += 'fill-opacity="' + str(opacity) + '" '

    return_string += 'stroke="' + border_colour + '" ' \
                     + 'stroke-width="' + str(border_width) + '" ' \
                     + '/>'
    return return_string

def write_svg_contigs(sequences, colours, accessory, contig_height, y_space, match_height, filehandle):
    y_top = 0
    for genome in sequences.keys():
        y_bottom = y_top + contig_height
        print(contigs_svg(sequences[genome], y_top, y_bottom), file=filehandle)
        print(write_accessory(sequences[genome], accessory, colours, y_top), file=filehandle)
        y_top += contig_height + 2 * y_space + match_height

def contigs_svg(seq, y_top, y_bottom):
    gene_len = 3
    lines = []
    x_left = 0
    x_right = gene_len*len(seq)
    coords = [(x_left, y_top), (x_right, y_top), (x_right, y_bottom), (x_left, y_bottom)]
    lines.append(
        svg_polygon(
            coords,
            'black',
            'black',
            border_width=0
        )
    )
    return '\n'.join(lines)

def write_svg_matches_between_two_assemblies(seq_1, seq_2, y_top, match_height):
    # top assembly = the query in nucmer matches
    # bottom assembly = the reference in nucmer matches
    lines = []
    y_bottom = y_top + match_height
    gene_len = 3
    i=0
    while i < len(seq_1):
        strand = True
        length = 0
        j = i
        neighbours_shared=True
        if seq_1[i] in seq_2:
             loc = seq_2.index(seq_1[i])
             while neighbours_shared:  
                length =+ gene_len
                if seq_1[j+1]!=seq_2[j+1]:
                    neighbours_shared = False
                j=+1
        elif -seq_1[i] in seq_2:
            strand = False
            loc = seq_2.index(-seq_1[i])
            while neighbours_shared:
                length =+ gene_len
                if seq_1[j+1]!=-seq_2[j-1]:
                    neighbours_shared = False
                j=+1
        top_start = i*gene_len
        top_end = top_start + length
        bottom_start = loc*gene_len
        bottom_end = bottom_start + length
        coords = [(top_start, y_top), (top_end, y_top), (bottom_end, y_bottom), (bottom_start, y_bottom)]
        if strand:
            lines.append((0, abs(top_start - top_end), svg_polygon(coords, 'lightseagreen', 'lightseagreen', opacity=0.8, border_width=0)))
        else:
            lines.append((1, abs(top_start - top_end), svg_polygon(coords, 'peru', 'peru', opacity=0.8, border_width=0)))
    lines.sort(key=itemgetter(0, 1))
    return '\n'.join([x[-1] for x in lines])

def write_accessory(seq, accessory, colours, y_top):
    lines = []
    accessory_width = 3
    gene_len = 3
    for i in range(len(seq)):
        if seq[i] in accessory:
            y_bottom = y_top + accessory_width
            if seq[i]>0:
                x_left = i*gene_len
                x_right = x_left + 0.6*gene_len
            else:
                x_right = (i+1)*gene_len
                x_left = x_right - 0.6*gene_len
            coords = [(x_left, y_top), (x_right, y_top), (x_right, y_bottom), (x_left, y_bottom)]
            lines.append(
                svg_polygon(
                    coords,
                    colours[seq[i]],
                    'black',
                    border_width=0
                )
            )
    return '\n'.join(lines)


def write_svg_header(filehandle, width, height):
    print(r'''<?xml version="1.0" standalone="no"?>''', file=filehandle)
    print(r'''<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">''', file=filehandle)
    print(r'<svg width="' + str(width) + '" height="' + str(height) + '">', file=filehandle)

def write_all_svg_matches(sequences, contig_height, match_height, y_space, filehandle):
    y_top = contig_height + y_space
    genomes = list(sequences.keys())
    for i in range(len(genomes) - 1):
        top_assembly = sequences[genomes[i]]
        bottom_assembly = sequences[genomes[i+1]]
        print(write_svg_matches_between_two_assemblies(top_assembly, bottom_assembly, y_top, match_height), file=filehandle)
        y_top += contig_height + 2 * y_space + match_height

def run(outprefix, unimog, acc_file):

    sequences = read_in_unimog(unimog)
    accessory, colours = read_in_accessory(acc_file)
    print(colours)

    contig_height = 1.5
    match_height = 30
    y_space = 1
    svg_height = (len(sequences.keys()) - 1) * (contig_height + 2 * y_space + match_height) + y_space + contig_height
    svg_width = 400

    svg_file = outprefix + '.svg'
    svg_fh = pyfastaq.utils.open_file_write(svg_file)
    write_svg_header(svg_fh, svg_width, svg_height)
    write_svg_contigs(sequences, colours, accessory, contig_height, y_space, match_height, svg_fh)
    write_all_svg_matches(sequences, contig_height, match_height, y_space, svg_fh)
    print('</svg>', file=svg_fh)
    pyfastaq.utils.close(svg_fh)

genes = "/home/daria/Documents/projects/ABC/rate_estimation/genes/accessory_st73_cl455_community_0_subcommunity_494.txt"
unimogs = "/home/daria/Documents/projects/ABC/rate_estimation/unimogs/st73_cl455_community_0_subcommunity_494.unimog"
out = "gene_alignments/st73_cl455_community_0_subcommunity_494"

run(out, unimogs, genes)
