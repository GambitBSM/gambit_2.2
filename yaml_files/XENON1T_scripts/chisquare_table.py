"""
Chi-squared table for solar and DM ALPs
"""

import best_fit_loglike_arxiv_v2 as lnl

def chisq(l):
    return -2. * l

def row(model, tritium):
    """
    @returns Single row in tex format
    """
    if tritium:
        bkg = "xe1t_3h"
    else:
        bkg = "xe1t"

    if model == "dm":
        if tritium:
            signal = "dm_xe1t_3h_alp"
            name = r"DM ALP + \ce{^3H}"
        else:
            signal = "dm_xe1t_alp"
            name = "DM ALP"
    elif model == "alp":
        if tritium:
            signal = "xe1t_3h_alp"
            name = r"Solar ALP + \ce{^3H}"
        else:
            signal = "xe1t_alp"
            name = "Solar ALP"
    elif model == "bkg":
        if tritium:
            signal = "xe1t_3h"
            name = r"Background + \ce{^3H}"
        else:
            signal = "xe1t"
            name = "Background"

    bkg_xe1t = chisq(getattr(lnl, bkg))
    signal_xe1t = chisq(getattr(lnl, signal))

    bkg_xe1t_r = chisq(getattr(lnl, "{}_r".format(bkg)))
    signal_xe1t_r = chisq(getattr(lnl, "{}_r".format(signal)))

    bkg_xe1t_r_wd = chisq(getattr(lnl, "{}_r_wd".format(bkg)))
    signal_xe1t_r_wd = chisq(getattr(lnl, "{}_r_wd".format(signal)))

    data = [signal_xe1t, bkg_xe1t - signal_xe1t,
            signal_xe1t_r, bkg_xe1t_r - signal_xe1t_r,
            signal_xe1t_r_wd, bkg_xe1t_r_wd -  signal_xe1t_r_wd]
    data_str = [r"\num{{{}}}".format(n) for n in data]
    return " & ".join([name] + data_str)


if __name__ == "__main__":

    rows = []
    for tritium in [False, True]:
        for model in ["bkg", "alp", "dm"]:
            rows.append(row(model, tritium))

    combined = "\\\\\n".join(rows)
    print(combined)
