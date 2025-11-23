from io import StringIO
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chisquare, ks_2samp, pearsonr, spearmanr, uniform, gamma

from sklearn.metrics import mean_squared_error, mean_absolute_percentage_error, r2_score
from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_codons

alphabet = 'ACGT'
transitions = ['A>G', 'C>T', 'G>A', 'T>C']
transversions = ['A>C', 'A>T', 'C>A', 'C>G', 'G>C', 'G>T', 'T>A', 'T>G']
amino_acid_codes = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}

grantham_matrix_str = '''
FIRST	R	L	P	T	A	V	G	I	F	Y	C	H	Q	N	K	D	E	M	W
S	110	145	74	58	99	124	56	142	155	144	112	89	68	46	121	65	80	135	177
R	0	102	103	71	112	96	125	97	97	77	180	29	43	86	26	96	54	91	101
L	0	0	98	92	96	32	138	5	22	36	198	99	113	153	107	172	138	15	61
P	0	0	0	38	27	68	42	95	114	110	169	77	76	91	103	108	93	87	147
T	0	0	0	0	58	69	59	89	103	92	149	47	42	65	78	85	65	81	128
A	0	0	0	0	0	64	60	94	113	112	195	86	91	111	106	126	107	84	148
V	0	0	0	0	0	0	109	29	50	55	192	84	96	133	97	152	121	21	88
G	0	0	0	0	0	0	0	135	153	147	159	98	87	80	127	94	98	127	184
I	0	0	0	0	0	0	0	0	21	33	198	94	109	149	102	168	134	10	61
F	0	0	0	0	0	0	0	0	0	22	205	100	116	158	102	177	140	28	40
Y	0	0	0	0	0	0	0	0	0	0	194	83	99	143	85	160	122	36	37
C	0	0	0	0	0	0	0	0	0	0	0	174	154	139	202	154	170	196	215
H	0	0	0	0	0	0	0	0	0	0	0	0	24	68	32	81	40	87	115
Q	0	0	0	0	0	0	0	0	0	0	0	0	0	46	53	61	29	101	130
N	0	0	0	0	0	0	0	0	0	0	0	0	0	0	94	23	42	142	174
K	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	101	56	95	110
D	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	160	181
E	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	126	152
M	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	67
'''

grantham = pd.read_csv(StringIO(grantham_matrix_str), sep='\t', index_col=0).replace(0, np.nan)
grantham.index.name = 'aa1'
grantham_long = grantham.melt(ignore_index=False, var_name='aa2', value_name='grantham_distance').dropna().reset_index()
grantham_long['aa1'] = grantham_long['aa1'].map(amino_acid_codes)
grantham_long['aa2'] = grantham_long['aa2'].map(amino_acid_codes)
grantham_long = pd.concat([grantham_long, grantham_long.rename(columns={'aa1':'aa2', 'aa2':'aa1'})], ignore_index=True)
del grantham


def collect_possible_changes(gc, spectrum_dict=None):
    coda = CodonAnnotation(gc)
    nucls = alphabet
    i = 1
    data = []
    for cdn1 in possible_codons:
        aa1 = amino_acid_codes[coda.translate_codon(cdn1)]
        for pic in range(3):
            nuc1 = cdn1[pic]
            for nuc2 in nucls:
                if nuc1 == nuc2:
                    continue
                cdn2 = list(cdn1)
                cdn2[pic] = nuc2
                cdn2 = ''.join(cdn2)
                aa2 = amino_acid_codes[coda.translate_codon(cdn2)]
                is_syn = aa1 == aa2
                sbs = f'{nuc1}>{nuc2}'
                data.append((pic, cdn1, cdn2, aa1, aa2, is_syn, sbs))
                i += 1

    df_changes = pd.DataFrame(
        data, columns=['pic', 'cdn1', 'cdn2', 'aa1', 'aa2', 'is_syn', 'sbs'])
    
    if spectrum_dict is not None:
        df_changes['rate'] = df_changes['sbs'].map(spectrum_dict)    

    return df_changes


def nuc_spectrum_to_matrix(spec):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    M = np.zeros((4,4))
    for i1,n1 in enumerate(alphabet):
        for i2,n2 in enumerate(alphabet):
            if n1!=n2:
                M[i1,i2] = spec[f"{n1}>{n2}"]
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()    
    return M


def cdn_spectrum_to_matrix(df_changes: pd.DataFrame):
    cdn_sbs = df_changes.groupby(['cdn1', 'cdn2'])['rate'].sum().reset_index()
    Pcdn = cdn_sbs.pivot(index='cdn1', columns='cdn2', values='rate').fillna(0.)
    Pcdn = (Pcdn.T / Pcdn.sum(axis=1)).T
    return Pcdn


def get_equilibrium_probabilities(Pmatrix: np.ndarray):
    Q = Pmatrix.copy()
    Q = Q - np.diag(np.sum(Q, axis=1))

    values, vectors = np.linalg.eig(Q.T)
    pi = vectors[:, np.isclose(values, 0)].real.flatten()
    pi = pi / pi.sum()
    return pi


def simulate_markov_continual(
        transition_Pmatrix: np.ndarray, initial_vector: np.ndarray, 
        num_iterations: int, delta_t=0.01):
    pi = initial_vector.copy()
    data = [pi.copy()]
    Q = transition_Pmatrix.copy()
    Q = Q - np.diag(np.sum(Q, axis=1))
    for _ in range(num_iterations):
        pi_new = pi + delta_t * (pi @ Q)
        pi_new = pi_new / pi_new.sum()
        data.append(pi_new.copy())
        if np.linalg.norm(pi_new - pi) < 1e-8:
            break
        pi = pi_new
    return data


def get_equilibrium_freqs(spectrum: pd.DataFrame, rate_col='MutSpec', gc=1):
    coda = CodonAnnotation(gc)
    spectrum_dict = spectrum.set_index('Mut')[rate_col].to_dict()
    df_changes = collect_possible_changes(gc, spectrum_dict)
    
    M = cdn_spectrum_to_matrix(df_changes)
    eq_prob = get_equilibrium_probabilities(M.values)

    eq_freqs_cdn = pd.Series(dict(zip(possible_codons, eq_prob)))
    eq_freqs_cdn.name = 'eq_freq'
    eq_freqs_cdn.index.name = 'cdn'
    eq_freqs_cdn = eq_freqs_cdn.reset_index()
    eq_freqs_cdn['aa'] = eq_freqs_cdn['cdn']\
        .map(coda.translate_codon).map(amino_acid_codes)
    
    eq_freqs_aa = eq_freqs_cdn[eq_freqs_cdn.aa !='*'].groupby('aa')['eq_freq'].sum()
    eq_freqs_aa /= eq_freqs_aa.sum()
    eq_freqs_aa = eq_freqs_aa.sort_values(ascending=False).reset_index()
    
    return eq_freqs_cdn, eq_freqs_aa


def get_random_spectra(n=10, tstv_sampling='poisson'):
    """For normal spectrum ts/tv must be between 2 and 10
    - param tstv_sampling: 'poisson', 'random'
    """
    columns = [f'{n1}>{n2}' for n1 in alphabet for n2 in alphabet if n1 != n2]
    data = uniform.rvs(size=(n, 12))
    rnd_spectra = pd.DataFrame(data, columns=columns)
    sampled_tstv_ratios = (rnd_spectra[transitions].sum(axis=1) / \
                          rnd_spectra[transversions].sum(axis=1)).values
    
    if tstv_sampling == 'poisson':
        tstv_ratios = np.random.poisson(3.5, n) + np.random.random(n)
        multiplier = tstv_ratios / sampled_tstv_ratios
        rnd_spectra[transitions] *= multiplier[:, None]
    elif tstv_sampling == 'random':
        tstv_ratios = sampled_tstv_ratios
    else:
        raise ValueError(f"Unknown tstv_sampling: {tstv_sampling}")

    rnd_spectra = (rnd_spectra.T / rnd_spectra.sum(axis=1)).T
    rnd_spectra['tstv_ratio'] = tstv_ratios
    rnd_spectra.reset_index(names='replica', inplace=True)
    rnd_spectra_long = rnd_spectra.melt(
        ['replica', 'tstv_ratio'], var_name='Mut', value_name='MutSpec')
    return rnd_spectra_long


def get_random_spectrum(tstv_ratio=None):
    """For normal spectrum ts/tv must be between 2 and 20"""
    rnd_spectrum = pd.DataFrame({
        'Mut': [f'{n1}>{n2}' for n1 in alphabet for n2 in alphabet if n1 != n2],
        'MutSpec': uniform.rvs(size=12),
    })
    if tstv_ratio is not None:
        ts_frac = rnd_spectrum[rnd_spectrum['Mut'].isin(transitions)].MutSpec.sum()
        tv_frac = rnd_spectrum[~rnd_spectrum['Mut'].isin(transitions)].MutSpec.sum()
        multiplier = tstv_ratio / (ts_frac / tv_frac)
        rnd_spectrum.loc[rnd_spectrum['Mut'].isin(transitions), 'MutSpec'] *= multiplier

    rnd_spectrum['MutSpec'] = rnd_spectrum['MutSpec'] / rnd_spectrum['MutSpec'].sum()
    return rnd_spectrum


def get_equilibrium_freqs_rnd(gc=1, tstv_ratio=2.0):
    """Use random spectrum for aa substitutions matrix generation"""
    rnd_spectrum = get_random_spectrum(tstv_ratio=tstv_ratio)
    res = get_equilibrium_freqs(rnd_spectrum, rate_col='MutSpec', gc=gc)
    return rnd_spectrum, res


def prepare_exp_aa_subst(spectrum: pd.DataFrame, rate_col='rate', gc=1, save_path=None):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index('Mut')[rate_col].to_dict()

    df_changes['rate'] = df_changes['sbs'].map(spectrum_dict)

    ## Calculate expected AA substitutions matrix
    exp_aa_subst = df_changes[(df_changes.aa1 != '*')&(df_changes.aa2 != '*')]\
        .groupby(['aa1', 'aa2'])['rate'].sum().reset_index()
    
    if save_path:
        exp_aa_subst.to_csv(save_path, float_format='%g', index=False)
    exp_aa_subst_matrix = exp_aa_subst.pivot(index='aa1', columns='aa2', values='rate').fillna(0.)
    return exp_aa_subst, exp_aa_subst_matrix


def prepare_exp_aa_subst_codons(spectrum: pd.DataFrame, rate_col='rate', gc=1, save_path=None, codon_pos='all'): 
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index('Mut')[rate_col].to_dict()
    df_changes['rate'] = df_changes['sbs'].map(spectrum_dict)

    if codon_pos == 'all':
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*"')
    elif codon_pos == '12':
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*" & pic in [0, 1]')
    elif isinstance(codon_pos, int):
        pic=codon_pos - 1
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*" & pic == @pic')
    else:
        raise ValueError(f"Invalid codon_pos value: {codon_pos}. Allowed values are 'all', '12', 1, 2, or 3.")

    ## Calculate expected AA substitutions matrix
    exp_aa_subst = df_changes_flt.groupby(['aa1', 'aa2'])['rate'].sum().reset_index()
    
    if save_path:
        exp_aa_subst.to_csv(save_path, float_format='%g', index=False)
    exp_aa_subst_matrix = exp_aa_subst.pivot(index='aa1', columns='aa2', values='rate').fillna(0.)
    return exp_aa_subst, exp_aa_subst_matrix


def prepare_exp_cdn_subst(spectrum: pd.DataFrame, rate_col='rate', gc=1, save_path=None):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index('Mut')[rate_col].to_dict()

    df_changes['rate'] = df_changes['sbs'].map(spectrum_dict)

    ## Calculate expected AA substitutions matrix
    exp_cdn_subst = df_changes[(df_changes.cdn1 != '*')&(df_changes.cdn2 != '*')]\
        .groupby(['cdn1', 'cdn2'])['rate'].sum().reset_index()
    
    if save_path:
        exp_cdn_subst.to_csv(save_path, float_format='%g', index=False)
    exp_cdn_subst_matrix = exp_cdn_subst.pivot(index='cdn1', columns='cdn2', values='rate').fillna(0.)
    return exp_cdn_subst, exp_cdn_subst_matrix


def prepare_rnd_exp_aa_subst(gc=1, tstv_ratio=None, save_path=None):
    """Use random spectrum for aa substitutions matrix generation"""
    rnd_spectrum = get_random_spectrum(tstv_ratio=tstv_ratio)
    res = prepare_exp_aa_subst(
        rnd_spectrum, rate_col='MutSpec', gc=gc, save_path=save_path)
    return res


def prepare_aa_subst(obs_df: pd.DataFrame, exp_aa_subst: pd.DataFrame, ref_aa_freqs: dict):
    aa_subst = obs_df.groupby(['aa1', 'aa2'])['count'].sum().rename('nobs').reset_index()
    aa_subst = aa_subst[(aa_subst['aa1'] != aa_subst['aa2']) & 
                        (aa_subst['aa1'] != '*') & 
                        (aa_subst['aa2'] != '*')]
    aa_subst['aa1'] = aa_subst['aa1'].map(amino_acid_codes)
    aa_subst['aa2'] = aa_subst['aa2'].map(amino_acid_codes)
    ref_aa_total_cnt = sum([x for x in ref_aa_freqs.values() if x > 0])
    aa_subst['ref_aa1_freq'] = aa_subst['aa1'].map(ref_aa_freqs) / ref_aa_total_cnt
    aa_subst['nobs_scaled'] = (aa_subst['nobs'] / aa_subst['ref_aa1_freq']).replace(np.inf, np.nan)
    aa_subst['nobs_scaled'] = aa_subst['nobs_scaled'] / aa_subst['nobs_scaled'].sum() * aa_subst['nobs'].sum()
    aa_subst = aa_subst.merge(exp_aa_subst.rename(columns={'rate': 'rate_exp'}), 'right').fillna(0)
    aa_subst = aa_subst[aa_subst['aa1'] != aa_subst['aa2']]
    aa_subst['nexp'] = aa_subst['rate_exp'] / aa_subst['rate_exp'].sum() * aa_subst['nobs_scaled'].sum()
    aa_subst['diff'] = aa_subst['nobs_scaled'] - aa_subst['nexp']
    aa_subst['pe'] = aa_subst['diff'] / aa_subst['nobs_scaled'] * 100  # %
    aa_subst['nobs_freqs'] = aa_subst['nobs_scaled'] / aa_subst['nobs_scaled'].sum()
    aa_subst['nexp_freqs'] = aa_subst['nexp'] / aa_subst['nexp'].sum()
    aa_subst['obs_relative_freq'] = (aa_subst['nobs_freqs'] / aa_subst['nexp_freqs'])
    aa_subst = aa_subst.merge(grantham_long, 'left')
    return aa_subst


def calc_accuracy(y_true, y_pred):
    '''
    LEGACY
    
    Acc = Sum(TP1, TP2, …, TPn)/total_obs_cnt'''
    TP = np.min([y_true, y_pred], axis=0)
    total_obs_cnt = y_true.sum()
    acc = TP.sum() / total_obs_cnt
    return acc


def calc_f1(y_true, y_pred):
    '''
    LEGACY

    Acc = Sum(TP1, TP2, …, TPn)/total_obs_cnt
    '''
    f1_weighted = 0.
    f1_macro = 0.
    total = y_true.sum()
    n = len(y_true)

    for o, e in zip(y_true, y_pred):
        if o < e:
            tp = o
            fp = e - o
            fn = 0
        elif o > e:
            tp = e
            fp = 0
            fn = o - e
        else:
            tp = o
            fp = fn = 0

        f1 = (2 * tp) / (2 * tp + fp + fn)

        f1_macro += f1 / n
        f1_weighted += f1 * (o / total)

    return f1_macro, f1_weighted


def calc_slope(y_true, y_pred):
    """
    Calculate slope of the regression line
    """
    try:
        slope, intercept = np.polyfit(y_true, y_pred, 1)
    except:
        slope, intercept = np.nan, np.nan

    # slope2 = np.linalg.lstsq(
    #     np.vstack([y_true, np.ones_like(y_true)]).T,
    #     y_pred,
    # )[0]
    return slope, intercept


def weighted_average_percentage_error(y_true, y_pred):
    return mean_absolute_percentage_error(y_true, y_pred, sample_weight=y_true)


def calc_metrics(aa_subst: pd.DataFrame):
    aa_subst = aa_subst.dropna()
    y_true = aa_subst.nobs_scaled
    y_pred = aa_subst.nexp

    # 2. Kolmogorov-Smirnov test
    ks_stat, ks_p = ks_2samp(y_true, y_pred)

    # 3. Log-Likelihood
    log_likelihood = np.sum(aa_subst.nobs_freqs * np.log(aa_subst.nexp_freqs) + \
                                (1 - aa_subst.nobs_freqs) * np.log(1 - aa_subst.nexp_freqs))

    # 4. RMSE
    try:
        rmse = mean_squared_error(aa_subst.nobs_freqs, aa_subst.nexp_freqs) ** 0.5
    except:
        rmse = np.nan

    mask = y_true > 0
    mape = mean_absolute_percentage_error(y_true[mask], y_pred[mask])
    wape = weighted_average_percentage_error(y_true[mask], y_pred[mask])

    # Spearman's rank correlation coefficient
    spearman_corr, spearman_p = spearmanr(y_true, y_pred)
    pearson_corr, pearson_p = pearsonr(y_true, y_pred)

    # total number of ns mutations
    mut_count = np.sum(y_true)

    slope, intercept = calc_slope(aa_subst.nobs_freqs, aa_subst.nexp_freqs)

    r2 = r2_score(y_true, y_pred)

    corr_chem_vs_rel_freq = pearsonr(
        aa_subst.obs_relative_freq, 
        aa_subst.grantham_distance,
    )
    
    metrics = {
        'r2': r2,
        'mape': mape,
        'wape': wape,
        'slope': slope,
        'intercept': intercept,
        'spearman_corr': spearman_corr,
        'spearman_p': spearman_p,
        'pearson_corr': pearson_corr,
        'pearson_corr_squared': pearson_corr**2,
        'pearson_p': pearson_p,
        'ks_stat': ks_stat,
        'ks_p': ks_p,
        'rmse': rmse,
        'log_likelihood': log_likelihood,
        'mut_count': mut_count,
        'mut_type_count': aa_subst.nobs.ne(0).sum(),
        'corr_chem_vs_rel_freq': corr_chem_vs_rel_freq.correlation,
        'corr_chem_vs_rel_freq_pval': corr_chem_vs_rel_freq.pvalue,
    }
    return metrics


def plot_exp_heatmap(exp_aa_subst_matrix: pd.DataFrame, save_path: str, show=True, annot=False):
    """
    Plot expected amino acid substitution matrix.
    """
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()
    max_flow_value = max(freqs_to.max(), freqs_from.max())+0.1

    fig, axs = plt.subplots(2, 3, figsize=(11, 10), 
                            width_ratios=[0.1, 1, .05], height_ratios=[1, 0.1])
    sns.heatmap(exp_aa_subst_matrix, annot=annot, fmt=".2f", 
                ax=axs[0, 1], cbar_ax=axs[0, 2], cmap="light:r", 
                cbar_kws={'label': 'Substitution rate'}, 
                mask=exp_aa_subst_matrix==0,
    )
    axs[0, 1].set_ylabel('')
    axs[0, 1].set_xlabel('')
    axs[0, 1].set_xticks([])
    axs[0, 1].set_yticks([])
    # axs[0, 1].set_title('Expected substitution rates between amino acids')

    sns.barplot(freqs_from.reset_index(), y='aa1', x=0, ax=axs[0, 0],
                color='lightgray', edgecolor='black')
    axs[0, 0].set_ylabel('From', fontsize=14)
    axs[0, 0].set_xlabel('')
    axs[0, 0].spines['top'].set_visible(False)
    axs[0, 0].spines['right'].set_visible(False)
    axs[0, 0].spines['bottom'].set_visible(False)
    axs[0, 0].spines['left'].set_visible(False)
    axs[0, 0].invert_xaxis()
    axs[1, 1].set_ylim(0, max_flow_value)
    axs[0, 0].set_xticks([])
    axs[0, 0].set_yticklabels(axs[0, 0].get_yticklabels(), fontsize=12)

    sns.barplot(freqs_to.reset_index(), x='aa2', y=0, ax=axs[1, 1],
                color='lightgray', edgecolor='black')
    axs[1, 1].set_ylabel('')
    axs[1, 1].set_xlabel('To', fontsize=14)
    axs[1, 1].spines['top'].set_visible(False)
    axs[1, 1].spines['right'].set_visible(False)
    axs[1, 1].spines['bottom'].set_visible(False)
    axs[1, 1].spines['left'].set_visible(False)
    axs[1, 1].invert_yaxis()
    axs[1, 1].set_ylim(max_flow_value, 0)
    axs[1, 1].set_yticks([])
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), fontsize=12)

    axs[1, 0].remove()
    axs[1, 2].remove()

    plt.tight_layout()
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_aa_eq_freqs(exp_aa_subst_matrix, save_path: str, show=True, figsize=(12, 3)):
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()

    _freqs_to = (freqs_to / freqs_to.sum()).copy()
    _freqs_from = (freqs_from / freqs_from.sum()).copy()
    _freqs_to.index.name = 'aa'
    _freqs_from.index.name = 'aa'

    plt.figure(figsize=figsize)
    ax = sns.barplot((_freqs_to - _freqs_from).sort_values(ascending=False).reset_index(), 
                x='aa', y=0, color='lightgray', edgecolor=".65")

    for bar in ax.patches:
        if bar.get_height() > 0.05:
            bar.set_color('orangered')
            bar.set_edgecolor('0.02')
        elif bar.get_height() > 0.02:
            bar.set_color('lightcoral')
            bar.set_edgecolor('0.05')
        elif bar.get_height() < -0.05:
            bar.set_color('royalblue')
            bar.set_edgecolor('0.02')
        elif bar.get_height() < -0.02:
            bar.set_color('lightblue')
            bar.set_edgecolor('0.45')

    plt.legend([
        plt.Rectangle((0,0),1,1,fc="orangered", edgecolor = '0.45'), 
        plt.Rectangle((0,0),1,1,fc="lightcoral", edgecolor = '0.45'), 
        plt.Rectangle((0,0),1,1,fc='lightblue', edgecolor = '0.45'),
        plt.Rectangle((0,0),1,1,fc='royalblue', edgecolor = '0.45'),
        ],[
            'Top Gainers', 'Gainers', 'Loosers', 'Top Loosers'], 
        loc='upper right')

    # plt.title(f'Expected aa subst based assignment')
    plt.ylabel('Rates [To - From]', fontsize=13)
    plt.xlabel('')
    plt.xticks(fontsize=12)
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()


def plot_obs_vs_exp(
        aa_subst, save_path=None, ax=None, show=True, 
        text='', text_x=-2.2, text_y=-4.
    ):
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
    aa_subst = aa_subst.copy()
    aa_subst['nobs_freqs_log'] = np.log10(aa_subst['nobs_freqs'])
    aa_subst['nexp_freqs_log'] = np.log10(aa_subst['nexp_freqs'])

    rho, pval = spearmanr(aa_subst['nobs_freqs'], aa_subst['nexp_freqs'])
    k, b = calc_slope(aa_subst['nobs_freqs'], aa_subst['nexp_freqs'])

    if ax is None:
        plt.figure(figsize=(5, 5))
        ax = plt.gca()

    ax = sns.scatterplot(aa_subst, ax=ax,
                x='nobs_freqs_log', y='nexp_freqs_log', 
                color='blue', alpha=0.5)
    ax = sns.regplot(aa_subst[aa_subst['nobs_freqs']>0], 
                color='blue', scatter=False,
                x='nobs_freqs_log', y='nexp_freqs_log', 
                ax=ax, label=f'y={k:.2f}x+{b:.1g}')

    # ticks = np.log10(np.array([1e-4, 1e-3, 1e-2, 1e-1]))
    ticks = np.linspace(-4, -1, 4)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    formatter = lambda x, pos: '10$^{' + str(int(x)) + '}$'
    ax.get_xaxis().set_major_formatter(formatter)
    ax.get_yaxis().set_major_formatter(formatter)
    # ax.get_xaxis().set_minor_locator(LogLocator(-10, 'auto'))

    ax.text(text_x, text_y, 
            f"{text}\nSpearman: {rho:.2f}\np-value: {pval:.2g}", 
            fontsize=10)
    ax.plot([-4, -1], [-4, -1], color='black', linestyle='--', label='y=x')
    ax.set_xlabel('Observed AA frequencies')
    ax.set_ylabel('Expected AA frequencies')
    ax.legend(loc='upper left', fontsize=10)
    if save_path:
        plt.savefig(save_path)
    if show:
        plt.show()


def plot_subst_freqs(aa_subst, title=''):
    aa_subst = aa_subst.copy()
    aa_subst['nobs_freqs_log'] = np.log10(aa_subst['nobs_freqs'])
    aa_subst['nexp_freqs_log'] = np.log10(aa_subst['nexp_freqs'])

    y_true, y_pred = aa_subst['nobs_freqs'], aa_subst['nexp_freqs']
    cor_res = pearsonr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    print(r2)
    print(f"Pearson correlation: {cor_res.correlation:.3f} (p-value: {cor_res.pvalue:.3g})")

    plt.figure(figsize=(8, 8))
    ax = sns.regplot(aa_subst[aa_subst['nobs_freqs']>0], 
                color='blue', scatter_kws={'alpha':0.5, 's':50},
                y='nobs_freqs_log', x='nexp_freqs_log')

    ticks = np.linspace(-4, -1, 4)
    ticks_minor = np.log10(np.concat([
        np.linspace(10**-4, 10**-3, 10),
        np.linspace(10**-3, 10**-2, 10)[1:],
        np.linspace(10**-2, 10**-1, 10)[1:],
    ]))
    ax.set_xticks(ticks, ticks, size=14)
    ax.set_yticks(ticks, ticks, size=14)
    ax.set_xticks(ticks_minor, minor=True)
    ax.set_yticks(ticks_minor, minor=True)
    formatter = lambda x, pos: '10$^{' + str(int(x)) + '}$'
    ax.get_xaxis().set_major_formatter(formatter)
    ax.get_yaxis().set_major_formatter(formatter)
    # ax.get_xaxis().set_minor_locator(LogLocator(-10, 'auto'))
    # ax.get_xaxis().set_minor_locator(LogLocator())

    # plt.text(-2, -4., 
    #          f"r={cor_res.correlation:.1f} (p={cor_res.pvalue:.1g})", 
    #          fontsize=10)
    plt.plot([-4, -1], [-4, -1], color='black', linestyle='--',)
    plt.ylabel('Observed AA substitution frequencies', fontsize=16)
    plt.xlabel('Predicted AA substitution frequencies', fontsize=16)
    # plt.grid()
    # plt.ylabel('Наблюдаемые частоты замещений аминокислот', fontsize=14)
    # plt.xlabel('Ожидаемые частоты замещений аминокислот', fontsize=14)
    plt.title(title, fontsize=16)
    # plt.legend(title=f"spearmanr={cor_res.correlation:.2f} (p={cor_res.pvalue:.1g})", title_fontsize=14)
    plt.legend(title=f"$R^2$ = {r2:.2f}", title_fontsize=16)
    # plt.savefig('./figures/obs_exp_aa_freqs_20A.pdf')
    return ax


def calculate_truncated_mean(alpha, scale, lower, upper):
    """
    Calculates the conditional expectation E[X | lower < X < upper] 
    for a Gamma distribution analytically.
    """
    # The mean of a truncated Gamma is related to the CDF of Gamma(alpha + 1)
    # E[X] = alpha * scale * (F(u | a+1) - F(l | a+1)) / (F(u | a) - F(l | a))
    
    # Denominator: Probability mass in this bin
    prob_mass = gamma.cdf(upper, a=alpha, scale=scale) - gamma.cdf(lower, a=alpha, scale=scale)
    
    if prob_mass == 0:
        return 0.0
        
    # Numerator components using alpha + 1
    num_upper = gamma.cdf(upper, a=alpha+1, scale=scale)
    num_lower = gamma.cdf(lower, a=alpha+1, scale=scale)
    
    expected_value = (alpha * scale) * (num_upper - num_lower) / prob_mass
    return expected_value


def categorize_site_rates_robust_plus_invariant(
        mut_counts, n_categories=6, hotspot_percentile=99.0, plot=False):
    """
    Categorizes sites into:
    - Category 0: Invariable/Zero observed mutations
    - Categories 1 to N-1: Gamma distributed rates (Normal sites)
    - Category N: Hotspots (Extreme outliers)
    """    
    # 1. Normalize relative to the GLOBAL mean (including zeros)
    rates_observed = np.array(mut_counts, dtype=float)
    global_mean = np.mean(rates_observed)
    rates_observed = rates_observed / global_mean
    
    # 2. Identify Zeros (Invariant + Stochastic Zeros)
    mask_zero = rates_observed == 0
    rates_zero = rates_observed[mask_zero]
    
    # Data remaining for Gamma/Hotspot analysis
    rates_nonzero = rates_observed[~mask_zero]
    
    # 3. Identify Hotspots within the NON-ZERO data
    # We calculate percentile based on the non-zero distribution to capture the true tail
    cutoff_value = np.percentile(rates_nonzero, hotspot_percentile)
    
    mask_normal = rates_nonzero <= cutoff_value
    mask_hotspot = rates_nonzero > cutoff_value
    
    rates_normal = rates_nonzero[mask_normal]
    rates_hotspot = rates_nonzero[mask_hotspot]
    
    print(f"Summary: {len(rates_zero)} zero sites, {len(rates_normal)} gamma sites, {len(rates_hotspot)} hotspots.")

    # 4. Fit Gamma to NORMAL (non-zero, non-hotspot) sites
    # IMPORTANT: floc=0 constraints the fit to start at 0, but since we removed exact zeros,
    # the fit will naturally approximate the shape of the low-rate tail.
    alpha, loc, scale = gamma.fit(rates_normal, floc=0)
    
    # 5. Discretize Normal Rates
    # We reserve Cat 0 for zeros, and Cat N for hotspots.
    # Available bins for Gamma = n_categories - 2
    n_gamma_cats = n_categories - 2
    
    quantiles = np.linspace(0, 1, n_gamma_cats + 1)
    gamma_thresholds = gamma.ppf(quantiles, a=alpha, scale=scale)
    gamma_thresholds[-1] = cutoff_value # Cap at hotspot cutoff
    
    # 6. Calculate Rates
    category_rates = []
    
    # --- Cat 0: Zeros ---
    category_rates.append(0.0)
    
    # --- Cat 1 to N-1: Gamma ---
    for i in range(n_gamma_cats):
        l, u = gamma_thresholds[i], gamma_thresholds[i+1]
        # Use the analytical mean function defined previously
        mean_rate = calculate_truncated_mean(alpha, scale, l, u)
        category_rates.append(mean_rate)
        
    # --- Cat N: Hotspot ---
    if len(rates_hotspot) > 0:
        hotspot_rate = np.mean(rates_hotspot)
    else:
        hotspot_rate = cutoff_value * 1.5
    category_rates.append(hotspot_rate)
    
    category_rates = np.array(category_rates)
    
    # 7. Assign Categories
    # Construct full thresholds: [0, ...gamma_thresholds..., inf]
    # We need to handle the exact zeros separately or they might get confused with the first bin
    
    site_categories = np.zeros(len(rates_observed), dtype=int)
    
    # Indices of non-zero sites
    nonzero_indices = np.where(~mask_zero)[0]
    
    # Thresholds for digitization (exclude 0, start from first gamma threshold)
    # We use the gamma thresholds + infinity
    digitization_bins = np.concatenate([gamma_thresholds, [np.inf]])
    
    # Digitize returns 1 for the first bin, so we add 0 to shift (Cat 0 is already 0)
    # But wait: Cat 0 is Zeros. Cat 1 is first Gamma.
    # np.digitize will return 1 for values < gamma_thresholds[0] (which is 0).
    # We need to map the non-zero values carefully.
    
    # Digitize the non-zero rates
    # bins: [0.1, 0.5, 1.2, 5.0, inf]
    # value 0.3 -> index 1. value 6.0 -> index 4.
    gamma_hotspot_cats = np.digitize(rates_nonzero, digitization_bins[1:]) + 1
    
    site_categories[nonzero_indices] = gamma_hotspot_cats

    # 7. Plotting
    if plot:
        plt.figure(figsize=(12, 6))
        
        # Histogram of all data (zoom in on relevant part)
        max_plot_x = np.percentile(rates_observed, 99.5) * 1.2
        bins = np.linspace(0, max_plot_x, 100)
        
        plt.hist(rates_observed, bins=bins, density=True, alpha=0.5, color='gray', label="Observed Rates")
        
        # Plot the Fitted Gamma (scaled to the proportion of normal data)
        x = np.linspace(0, max_plot_x, 500)
        # We multiply PDF by proportion of normal sites because the histogram includes hotspots
        prop_normal = len(rates_normal) / len(rates_observed)
        plt.plot(x, gamma.pdf(x, a=alpha, scale=scale) * prop_normal, 
                 'r-', lw=2, label=f"Gamma Fit (Normal Data Only)")
        
        # Plot thresholds
        for th in gamma_thresholds[1:-1]:
            plt.axvline(th, color='blue', linestyle='--', alpha=0.6)
            
        plt.axvline(cutoff_value, color='red', linestyle='-', lw=2, label="Hotspot Cutoff")
        
        # Visualize the discrete rates
        for rate in category_rates:
            plt.scatter(rate, 0, color='black', zorder=5, s=50, marker='x')

        plt.xlabel("Relative Site Rate")
        plt.ylabel("Density")
        plt.title(f"Site Categorization: {n_gamma_cats} Gamma + 1 Hotspot")
        plt.legend()
        plt.show()


    return site_categories, category_rates, alpha