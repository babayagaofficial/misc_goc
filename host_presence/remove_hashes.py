

with open("presence_per_host.tsv", "r") as f:
	with open("presence_per_host_hashless.tsv", "w") as v:
		for line in f:
			v.write(line.replace("#","_"))