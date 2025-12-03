from scipy.spatial.distance import pdist, squareform
from scipy.cluster.vq import kmeans, whiten
from scipy import stats
from scipy.cluster.hierarchy import centroid
from dcj_downstream import *
from statistics import mean
import random
import itertools
import matplotlib.pyplot as plt


def summary_stats(subcomm, subcomm_df, subcomm_vals, dists_subcomm):
    sizes = []
    ratios = {}
    diffs = {}
    centres = {}
    for clade in subcomm:
        clade_name = clade.replace("_community_0_subcommunity_501","").replace("_community_0_subcommunity_503","")
        ints = list(subcomm_df[subcomm_df["clade"]==clade_name].index)
        sizes.append(len(ints))
        ratios[clade_name] = []
        diffs[clade_name] = []
        for i in ints:
            inner_dists = [dists_subcomm[i][j] for j in ints if j!=i]
            outer_dists = [dists_subcomm[i][j] for j in range(len(dists_subcomm[0])) if not j in ints]
            diffs[clade_name].append(abs(min(inner_dists)-min(outer_dists)))
            ratios[clade_name].append(min(inner_dists)/min(outer_dists))

        vals = [subcomm_vals[j] for j in ints]
        centre = [mean([vals[i][0] for i in range(len(vals))]),mean([vals[i][1] for i in range(len(vals))])]
        vals.append(centre)
        vals[0], vals[-1] = vals[-1], vals[0]
        pdists = pdist(vals)
        dists = [pdists[j - 2 // 2] for j in range(1,len(vals))]
        centres[clade_name] = dists
    return centres, diffs, ratios, sizes
        

clades_501 = ["st73_cl429_community_0_subcommunity_501", "st73_cl597_community_0_subcommunity_501", "st73_cl613_community_0_subcommunity_501"]
clades_503 = ["st73_cl568_community_0_subcommunity_503", "st95_cl378_community_0_subcommunity_503", "st95_cl400_community_0_subcommunity_503", "st95_cl461_community_0_subcommunity_503", "st131_cl416_community_0_subcommunity_503"]

subcomm_501 = pd.DataFrame(columns=["DCJ-Indel", "Fitch", "SNPs"])
subcomm_503 = pd.DataFrame(columns=["DCJ-Indel", "Fitch", "SNPs"])


for clade in clades_501:
    df = read_in(clade)
    subcomm_501 = pd.concat([subcomm_501,df])
    subcomm_501.index = [i for i in range(len(subcomm_501))]
for clade in clades_503:
    df = read_in(clade)
    subcomm_503 = pd.concat([subcomm_503,df])
    subcomm_503.index = [i for i in range(len(subcomm_503))]

subcomm_501_vals = subcomm_501[["DCJ-Indel","SNPs"]].to_numpy()
subcomm_503_vals = subcomm_503[["DCJ-Indel","SNPs"]].to_numpy()

dists_501 = squareform(pdist(subcomm_501_vals))
dists_503 = squareform(pdist(subcomm_503_vals))

centres_501, diffs_501, ratios_501, sizes_501 = summary_stats(clades_501,subcomm_501,subcomm_501_vals,dists_501)
centres_503, diffs_503, ratios_503, sizes_503 = summary_stats(clades_503,subcomm_503,subcomm_503_vals,dists_503)



'''     
actual_ratio = mean(list(itertools.chain(list(ratios_501.values())))[0])
actual_diff = mean(list(itertools.chain(list(diffs_501.values())))[0])

n = 100
nn= {i:{} for i in range(n)}
rng = random.seed(42)

for it in range(n):
    rand = [i for i in range(len(dists_501[0]))]
    random.shuffle(rand)
    ratios = np.ndarray((n,len(rand)))
    diffs = np.ndarray((n,len(rand)))
    centres = np.ndarray((n,len(rand)))
    for colour in range(len(sizes_501)):
        if colour == 0:
            ints = rand[:sizes_501[colour]]
        else:
            ints = rand[sizes_501[colour-1]:sizes_501[colour]]
        for i in ints:
            inner_dists = [dists_501[i][j] for j in ints if j!=i]
            outer_dists = [dists_501[i][j] for j in range(len(rand)) if not j in ints]
            diffs[it][i]=abs(min(inner_dists)-min(outer_dists))
            ratios[it][i] = min(inner_dists)/min(outer_dists)
        vals = [subcomm_501_vals[j] for j in ints]
        centre = [mean([vals[i][0] for i in range(len(vals))]),mean([vals[i][1] for i in range(len(vals))])]
        vals.append(centre)
        vals[0], vals[-1] = vals[-1], vals[0]
        pdists = pdist(vals)
        dists = [pdists[j - 2 // 2] for j in range(1,len(vals))]
        #centres[it] = dists


sampled_ratios = [mean(ratios[it]) for it in range(n)]

f = lambda x: 1 if x>=actual_diff else 0
p = sum([f(ratio) for ratio in sampled_ratios])/n #converges to likelihood that the ratio is greater than the actual ratio (smaller ratio-> more clustering)
print(p)





res = stats.fligner(ratios_501['st73_cl429'], ratios_501['st73_cl597'], ratios_501['st73_cl613'])

print(res.pvalue)

res = stats.fligner(centres_501['st73_cl429'], centres_501['st73_cl597'], centres_501['st73_cl613'])

print(res.pvalue)

res = stats.fligner(diffs_501['st73_cl429'], diffs_501['st73_cl597'], diffs_501['st73_cl613'])

print(res.pvalue)


def statistic(*samples):
    return stats.fligner(*samples).statistic
ref = stats.permutation_test(
    (ratios_501['st73_cl429'], ratios_501['st73_cl597'], ratios_501['st73_cl613']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)

ref = stats.permutation_test(
    (centres_501['st73_cl429'], centres_501['st73_cl597'], centres_501['st73_cl613']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)

ref = stats.permutation_test(
    (diffs_501['st73_cl429'], diffs_501['st73_cl597'], diffs_501['st73_cl613']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)

res = stats.fligner(ratios_503['st73_cl568'], ratios_503['st95_cl378'], ratios_503['st95_cl400'], ratios_503['st95_cl461'], ratios_503['st131_cl416'])

print(res.pvalue)

res = stats.fligner(centres_503['st73_cl568'], centres_503['st95_cl378'], centres_503['st95_cl400'], centres_503['st95_cl461'], centres_503['st131_cl416'])

print(res.pvalue)

res = stats.fligner(diffs_503['st73_cl568'], diffs_503['st95_cl378'], diffs_503['st95_cl400'], diffs_503['st95_cl461'], diffs_503['st131_cl416'])

print(res.pvalue)


def statistic(*samples):
    return stats.fligner(*samples).statistic
ref = stats.permutation_test(
    (ratios_503['st73_cl568'], ratios_503['st95_cl378'], ratios_503['st95_cl400'], ratios_503['st95_cl461'], ratios_503['st131_cl416']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)

ref = stats.permutation_test(
    (centres_503['st73_cl568'], centres_503['st95_cl378'], centres_503['st95_cl400'], centres_503['st95_cl461'], centres_503['st131_cl416']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)

ref = stats.permutation_test(
    (diffs_503['st73_cl568'], diffs_503['st95_cl378'], diffs_503['st95_cl400'], diffs_503['st95_cl461'], diffs_503['st131_cl416']), statistic,
    permutation_type='independent', alternative='greater'
)

print(ref.pvalue)
'''
