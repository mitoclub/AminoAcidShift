import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chisquare, ks_2samp, pearsonr, spearmanr, uniform
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
    
    metrics = {
        'r2': r2,
        'mape': mape,
        'wape': wape,
        'slope': slope,
        'intercept': intercept,
        'spearman_corr': spearman_corr,
        'spearman_p': spearman_p,
        'pearson_corr': pearson_corr,
        'pearson_p': pearson_p,
        'ks_stat': ks_stat,
        'ks_p': ks_p,
        'rmse': rmse,
        'log_likelihood': log_likelihood,
        'mut_count': mut_count,
        'mut_type_count': aa_subst.nobs.ne(0).sum(),
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
