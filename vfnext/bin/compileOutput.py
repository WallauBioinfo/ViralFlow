#!/usr/bin/env python
"""Compile per-sample ViralFlow results into batch-level CSV files.

Run once at the end of an ILLUMINA run, by the `compileOutputs` process. It is
not used in NANOPORE mode.

Input layout
------------
`--dataDir` is the run's `outDir`. The script walks it looking for directories
named `<sample>_results`, and treats the part before `_results` as the sample
code (`cod`). A sample is only processed if its consensus FASTA
`<cod>.depth<depth>.fa` is present; `--depth` must therefore match the `depth`
the run used, or nothing is found.

Two virus tags, and how they differ
-----------------------------------
`-virus_tag` selects between two behaviours. `sars-cov2` expects the lineage
tools to have run, `custom` does not:

| | sars-cov2 | custom |
|---|---|---|
| Per-sample files expected | + Pangolin CSV, + NextClade CSV | intrahost, mutations, Picard metrics, wgs only |
| Batch CSVs written | + `pango.csv`, + `nextclade.csv` | the rest only |
| `major_summary.csv` / `minor_summary.csv` | written, split on whether the sample has minor variants | not written |
| `lineage_summary.csv` | written from the Pangolin calls | not written |
| `short_summary.csv` | n/a | written: depth and coverage per sample |

Pangolin output is looked up under two names, `<cod>.fa.pango.out.csv` first and
`<cod>.all.fa.pango.out.csv` as a fallback, because runPangolin picks the name
based on whether the sample had minor variants to append.

Degradation, not failure
------------------------
Missing per-sample files are recorded in `errors_detected.csv` and the sample is
skipped; the run still produces whatever could be compiled. Each batch CSV is
also written independently, so an absent input silently drops that one file with
a warning rather than aborting. The single hard failure is finding no
`*_results` directory with a consensus at all, which exits 1.

Outputs
-------
Written to `--outputDir`: `seqbatch.fa` (all consensus sequences), plus
`mutations.csv`, `reads_count.csv`, `wgs.csv`, `errors_detected.csv`, and the
tag-specific files listed above.

Logging
-------
Messages go to stdout through `logging`, which is where the `print` calls this
replaced already sent them, so `.command.out` reads as before. `--log-level`
controls verbosity; `DEBUG` adds the per-sample progress counter and dumps the
lineage table.
"""

import argparse
import logging
import os
import sys

import pandas as pd
from Bio import SeqIO

__author__ = "Antonio Marinho"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Antonio Marinho"
__email__ = "amarinhosn@pm.me"
__date__ = "2022/06/07"
__username__ = "amarinhosn"

LOGGER = logging.getLogger("compileOutput")


def configure_logging(level):
    """Send logs to stdout, where the print calls this replaced already went."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(levelname)s: %(message)s",
        stream=sys.stdout,
    )


# --- FUNCTIONS ---------------------------------------------------------------
# --- compile output
def get_consensus_seq(fasta_path, cod):
    """
    get sequence line from a single sequence fasta file
    """
    with open(fasta_path, "r") as fasta_fl:
        for line in fasta_fl:
            if line.startswith(">"):
                # sanity check
                assert cod in line
            else:
                return line


def parse_fabc(fbc_path):
    """
    parse depth data on 'fa.bc' file
    """
    # load file
    fl = open(fbc_path, "r")
    # get only pos and depth per line and return data as dataframe
    pos_lst = []
    depth_lst = []
    for line in fl:
        l_dt = line.split("\t")
        pos_lst.append(int(l_dt[1]))
        depth_lst.append(int(l_dt[3]))

    df = pd.DataFrame.from_dict({"pos": [pos_lst], "depth": [depth_lst]})
    return df


def __check_if_outfls_exist(cod, out_fls_lst, mut_fl):
    """
    Check if expected files exist
    """
    # for a list of expected files, check if they exist at specified path
    # if not, then print a warning, store missing error data on a dictionary
    # and return it
    err_dct_lst = []
    for fl in out_fls_lst:
        try:
            assert os.path.isfile(fl)
        except AssertionError:
            srvc = fl.split(".")[-2]
            if fl == mut_fl:
                srvc = "Mutations"
            LOGGER.warning("missing %s output for %s [expects %s]", srvc, cod, fl)
            type = "ERROR"
            kind = f"No {srvc} output [expected={fl}]"
            err_dct = {"cod": cod, "type": type, "kind": kind}
            err_dct_lst.append(err_dct)
            continue
    return err_dct_lst


def __parse_xsv(cod, xsv_path, sep=","):
    """
    get dataframe from tabular data (ex: csv and tsv)
    """
    df = pd.read_csv(xsv_path, sep=sep)
    df["cod"] = [cod for i in range(0, len(df))]
    return df


def __parse_metrics(cod, metrics_path):
    """
    get read counts from picard output
    """
    with open(metrics_path, "r") as metrics_fl:
        for line in metrics_fl:
            if line.startswith("#"):
                continue
            # get reads summary
            if line.startswith("PAIR") or line.startswith("UNPAIRED"):
                d_line = line.split("\t")
                dct = {
                    "cod": cod,
                    "total_reads": int(d_line[1]),
                    "mepf_reads_aligned": int(d_line[5]),
                    "pct_pf_reads_aligned": float(d_line[6]),
                }
                break

    df = pd.DataFrame(dct, index=[0])
    return df


def get_mut(row):
    """
    get mutation data on a single string
    """
    return row["REF"] + str(row["POS"]) + row["ALT"]


def compile_output_fls(data_dir, out_dir, depth, virus_tag):
    """
    Compile Viralflow output data files on subdirs
    Parameters
    ==========
    data_dir:<path>
        path for dir containing results output subdirs
    out_dir:<path>
        output dir to write compiled data
    depth:<depth>
        depth prefix used
    """
    # create outdir if does not exist
    try:
        assert os.path.exists(out_dir)
    except AssertionError:
        LOGGER.info("creating output dir")
        os.system("mkdir " + out_dir)
    if out_dir.endswith("/") is False:
        out_dir += "/"
    LOGGER.info("compiling output files")
    # create multifasta file
    multifas_fl = open(out_dir + "seqbatch.fa", "w")

    # create list to store data
    pango_df_lst = []
    nxtcd_df_lst = []
    mtrcs_df_lst = []
    ithst_df_lst = []
    wgs_df_lst = []
    mut_df_lst = []
    err_dct_lst = []
    skip_lst = []
    dpth_lst = []

    # iterate over subdirs on data dir
    c = 0
    for path, subdirs, files in os.walk(data_dir):
        for name in files:
            # check only files at *.results subdirs
            sub_dir = path.split("/")[-1]

            if sub_dir.endswith("_results") is False:
                continue

            # get consensus
            cod = sub_dir.replace("_results", "")  # sub_dir.split("_")[0]
            prfx = path + f"/{cod}."
            if cod in skip_lst:
                continue

            # look for consensus sequence fasta, if yes, then check for
            # 1) expected files
            # 2) parse and compiled results in single files
            if name.endswith(f"depth{depth}.fa"):
                c += 1
                out_fls_lst = []
                LOGGER.debug("%s samples processed", c)
                fasta_path = os.path.join(path, name)
                # sanity
                assert cod in name
                # write sequence to new
                seq = get_consensus_seq(fasta_path, cod)
                if len(seq) == 1:
                    type = "ERROR"
                    kind = "No consensus sequence obtained"
                    LOGGER.warning("No consensus sequence obtained for %s", cod)
                    err_dct_lst.append({"cod": cod, "type": type, "kind": kind})
                    skip_lst.append(cod)
                    continue
                # write to multifasta
                multifas_fl.write(">" + cod + "\n")
                multifas_fl.write(seq + "\n")
                # continue

                # -- SET EXPECTED OUTPUT FILES --------------------------------
                # depth data
                depth_fl = f"{prfx}depth{depth}.fa.bc"
                if virus_tag == "sars-cov2":
                    # pangolin
                    try:
                        pango_fl = f"{prfx}fa.pango.out.csv"
                        assert os.path.isfile(pango_fl)
                    except AssertionError:
                        pango_fl = f"{prfx}all.fa.pango.out.csv"
                    out_fls_lst.append(pango_fl)

                    # nextclade
                    try:
                        nxtcl_fl = f"{prfx}depth{depth}.fa.nextclade.csv"
                        assert os.path.isfile(nxtcl_fl)
                    except AssertionError:
                        nxtcl_fl = f"{prfx}depth{depth}.all.fa.nextclade.csv"

                    out_fls_lst.append(nxtcl_fl)
                    # chromosomes
                    # chrms_fl = f"{path}/chromosomes.report"

                # intrahost
                ithst_fl = f"{prfx}depth{depth}.fa.bc.intrahost.tsv"
                out_fls_lst.append(ithst_fl)
                # mutation
                mut_fl = f"{prfx}tsv"
                out_fls_lst.append(mut_fl)
                # picard
                pcrd_mtrcs_fl = f"{path}/metrics.alignment_summary_metrics"
                out_fls_lst.append(pcrd_mtrcs_fl)

                # wgs
                wgs_fl = f"{path}/wgs"
                out_fls_lst.append(wgs_fl)

                # --- CHECK IF FILES EXIST ------------------------------------
                # out_fls_lst = [pango_fl, ithst_fl, mut_fl, pcrd_mtrcs_fl]#,depth_fl, chrms_fl]
                err_dct_ = __check_if_outfls_exist(cod, out_fls_lst, mut_fl)
                # if any file is missing, skip current sample
                if len(err_dct_) > 0:
                    for err in err_dct_:
                        err_dct_lst.append(err)
                    break

                # --- PARSE DATA ----------------------------------------------
                if virus_tag == "sars-cov2":
                    # pangolin
                    pango_df_lst.append(__parse_xsv(cod, pango_fl, sep=","))
                    # chromossomes
                    # chrms_df_lst.append(__parse_xsv(cod, chrms_fl, sep="\t"))
                    # nextclade
                    nxtcd_df_lst.append(__parse_xsv(cod, nxtcl_fl, sep=";"))

                # intrahost
                ithst_df_lst.append(__parse_xsv(cod, ithst_fl, sep="\t"))
                # metrics (picard)
                mtrcs_df_lst.append(__parse_metrics(cod, pcrd_mtrcs_fl))
                # Mutations
                mut_df = pd.read_csv(mut_fl, sep="\t")
                val_mut_df = mut_df.loc[mut_df["PASS"].eq(True)]
                # wgs
                wgs_df_lst.append(__parse_wgs(wgs_fl, cod))

                # check if no mutation data
                if len(val_mut_df) == 0:
                    LOGGER.warning("no mutation data for %s", cod)
                    mut_lst = []
                if len(val_mut_df) > 0:
                    mut_lst = val_mut_df.apply(get_mut, axis=1).values

                dct = {"cod": cod, "mut": mut_lst}
                mut_df_lst.append(pd.DataFrame(dct))
                # depth data
                dpth_df = parse_fabc(depth_fl)
                dpth_df["cod"] = [cod]
                dpth_lst.append(dpth_df)

    LOGGER.info("Total %s samples processed", c)
    # mount compiled dataframes
    if c == 0:
        LOGGER.error("No ViralFlow output files found")
        sys.exit(1)
    if c > 0:
        LOGGER.info("writing compiled data")
        if virus_tag == "sars-cov2":
            try:
                assert len(pango_df_lst) > 0
                all_pango_df = pd.concat(pango_df_lst, ignore_index=True)
                all_pango_df.to_csv(out_dir + "/pango.csv", index=False)
                LOGGER.info("wrote pango.csv")
            except AssertionError:
                LOGGER.warning("No data from pangolin")

            try:
                assert len(nxtcd_df_lst) > 0
                all_nxtcd_df = pd.concat(nxtcd_df_lst, ignore_index=True)
                all_nxtcd_df.to_csv(out_dir + "/nextclade.csv", index=False)
                LOGGER.info("wrote nextclade.csv")
            except AssertionError:
                LOGGER.warning("No data from Nextclade")

        try:
            assert len(mut_df_lst) > 0
            all_mut_df = pd.concat(mut_df_lst, ignore_index=True)
            all_mut_df.to_csv(out_dir + "/mutations.csv", index=False)
            LOGGER.info("wrote mutations.csv")
        except AssertionError:
            LOGGER.warning("No mutation data")

        try:
            assert len(mtrcs_df_lst) > 0
            all_metrics_df = pd.concat(mtrcs_df_lst, ignore_index=True)
            all_metrics_df.to_csv(out_dir + "/reads_count.csv", index=False)
            LOGGER.info("wrote reads_count.csv")
        except AssertionError:
            LOGGER.warning("No reads_count data")
        try:
            assert len(wgs_df_lst) > 0
            all_wgs_df = pd.concat(wgs_df_lst, ignore_index=True)
            all_wgs_df.to_csv(out_dir + "/wgs.csv")
            LOGGER.info("wrote wgs.csv")
        except AssertionError:
            LOGGER.warning("No wgs data")

        errors_df = pd.DataFrame(err_dct_lst)
        errors_df.to_csv(out_dir + "/errors_detected.csv", index=False)
        LOGGER.info("wrote errors_detected.csv")

    LOGGER.info("DONE")


# get lineage summary


def load_short_summary_df(wgs_df, pang_df):
    """ """
    # get slices
    wgs_slice = wgs_df[["cod", "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE"]]
    pang_slice = pang_df[["cod", "taxon", "lineage", "scorpio_call"]]
    # merge dataframes
    short_summary = pd.merge(wgs_slice, pang_slice, right_on="cod", left_on="cod")
    short_summary = short_summary.rename(
        columns={
            "MEAN_COVERAGE": "mean_depth_coverage",
            "SD_COVERAGE": "sd_depth_coverage",
            "MEDIAN_COVERAGE": "median_depth_coverage",
        }
    )
    # return df
    return short_summary


# --- coverage ------------------------------------------------------------------------
def computeCoverage(seq):
    """
    (N - total bases) / total_bases
    """
    total_N = sum([1 for i in seq if i == "N"])
    total_bases = len(seq)
    return (total_bases - total_N) / total_bases


def loadCoverageDF(multifasta_path):
    """
    Compute coverage for each sequence and return sample cod, sequence and
    coverage.
    Parameters
    ----------
    multifasta_path:<str>
        path for fasta file containing samples consensus sequence
    Returns
    -------
    pd.DataFrame
    """
    dct_lst = []
    for record in SeqIO.parse(multifasta_path, "fasta"):
        cod = record.id
        seq = record.seq
        cov = "{:.4f}".format(computeCoverage(seq))
        dct_lst.append({"cod": cod, "coverage_breadth": cov})
    return pd.DataFrame(dct_lst)


# -------------------------------------------------------------------------------------
def __parse_wgs(wgs_flpth, cod):

    with open(wgs_flpth, "r") as wgs_fl:
        for line in wgs_fl:
            # skip lines
            if line.startswith("#") or line.startswith("\n"):
                continue
            # get column names
            if line.startswith("GENOME_TERRITORY"):
                columns_name = line.replace("\n", "").split("\t")
                continue
            # get values (wgs has only one line of data before histogram data)
            else:
                values = line.replace("\n", "").split("\t")
                break
    assert len(columns_name) == len(values)
    dct = {}
    dct["cod"] = cod

    for i, col_i in enumerate(columns_name):
        dct[col_i] = [values[i]]

    df = pd.DataFrame.from_dict(dct)
    return df


def get_lineages_summary(wgs_csv, outdir, multifasta, virus_tag, pango_csv=None):
    """
    count lineages on a given pangolin dataframe
    """

    def isMinor(row):
        """
        locate minor row on short summary
        """
        if row["taxon"].endswith("_minor"):
            return True
        else:
            return False

    # Lineage summary is currently only supported for sars-cov2, for other viruses only
    # mean depth and covarage will be provided

    if virus_tag == "sars-cov2":
        # load pango df
        LOGGER.info("compute lineage summary")
        if doFileExists(pango_csv):
            pango_df = pd.read_csv(pango_csv, index_col=False)

            LOGGER.info("%s total samples", len(pango_df))
            lineage_df = pango_df["lineage"].value_counts()
            lineage_df = lineage_df.rename_axis("lineage")
            lineage_df = lineage_df.rename("count")
            # lineage_df =  lineage_df.Series.rename(index='lineage')
            lineage_df.to_csv(outdir + "/lineage_summary.csv", index=True)
            LOGGER.info("wrote %slineage_summary.csv", outdir)
            LOGGER.debug("lineage table:\n%s", lineage_df)
        else:
            LOGGER.warning(
                "%s was not found, no lineage summary will be written", pango_csv
            )

    # short summary
    LOGGER.info("generating short summary [sample, depth, coverage, lineage]")

    # if wgs csv does not exist, no short summary can be writen
    if doFileExists(wgs_csv):
        wgs_df = pd.read_csv(wgs_csv)
        wgs_slice = wgs_df[["cod", "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE"]]

        if virus_tag == "sars-cov2":
            short_summary_df = load_short_summary_df(wgs_df, pango_df)
            # split minors from majors data
            minor_btable = short_summary_df.apply(isMinor, axis=1)
            minors_df = short_summary_df.loc[minor_btable]
            major_df = short_summary_df.loc[~minor_btable]
            # write csvs
            major_df.to_csv(outdir + "major_summary.csv", index=False)
            LOGGER.info("wrote %smajor_summary.csv", outdir)
            minors_df.to_csv(outdir + "minor_summary.csv", index=False)
            LOGGER.info("wrote %sminor_summary.csv", outdir)
            # check for empty dataframes
            if len(minors_df) == 0:
                LOGGER.info("No minor sequences available")
            if len(major_df) == 0:
                LOGGER.warning("No major sequence available")

            if doFileExists(multifasta):
                cov_df = loadCoverageDF(multifasta)
                short_summary_df = short_summary_df.merge(cov_df, on="cod")
                short_summary_df.to_csv(f"{outdir}short_summary.csv")
            else:
                LOGGER.warning(
                    "%s was not found, no short_summary will be written", multifasta
                )

        if virus_tag == "custom":
            short_summary_df = wgs_slice
            if doFileExists(multifasta):
                cov_df = loadCoverageDF(multifasta)
                short_summary_df = short_summary_df.merge(cov_df, on="cod")
                short_summary_df = short_summary_df.rename(
                    columns={
                        "MEAN_COVERAGE": "mean_depth_coverage",
                        "SD_COVERAGE": "sd_depth_coverage",
                        "MEDIAN_COVERAGE": "median_depth_coverage",
                    }
                )

                if "lineage" not in short_summary_df.columns:
                    short_summary_df["lineage"] = None

                short_summary_df.to_csv(f"{outdir}short_summary.csv")

            else:
                LOGGER.warning(
                    "%s was not found, no short_summary will be written", multifasta
                )
    else:
        LOGGER.warning("%s was not found, no lineage summary will be written", wgs_csv)


def doFileExists(file_path):
    return os.path.exists(file_path)


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This scripts compile outputs of samples processed on viralow on single files"
    )
    parser.add_argument(
        "-dD",
        "--dataDir",
        required=True,
        help="directory containing viralflow output files",
    )
    parser.add_argument(
        "-oD", "--outputDir", required=True, help="directory to output compiled files"
    )
    parser.add_argument(
        "-dp",
        "--depth",
        type=int,
        default=5,
        help="Minimum depth value to consider set as input to ViralFlow as minor variant, default = 5",
    )
    parser.add_argument("-virus_tag", type=str, help="viral tag provided")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="verbosity; DEBUG adds the per-sample counter and the lineage table",
    )

    args = parser.parse_args()
    configure_logging(args.log_level)
    # add check for virus tag
    valid_virus = ["sars-cov2", "custom"]
    assert args.virus_tag in valid_virus, f"{args.virus_tag} not a valid virus tag"

    compile_output_fls(args.dataDir, args.outputDir, args.depth, args.virus_tag)
    if args.virus_tag == "sars-cov2":
        get_lineages_summary(
            "./wgs.csv",
            args.outputDir,
            "./seqbatch.fa",
            args.virus_tag,
            pango_csv="./pango.csv",
        )
    if args.virus_tag == "custom":
        get_lineages_summary(
            "./wgs.csv", args.outputDir, "./seqbatch.fa", args.virus_tag
        )
