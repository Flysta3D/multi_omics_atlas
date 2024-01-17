from typing import  Optional, Tuple, Union
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

#####################################################################
# Perform gene list enrichment analysis using GSEApy                #
# Reference: https://github.com/zqfang/GSEApy                       #
#            https://maayanlab.cloud/modEnrichr/                    #
#            https://maayanlab.cloud/FlyEnrichr/help#background&q=0 #
#####################################################################


def basic_sets(
    organism: Literal["human", "mouse", "yeast", "fly", "fish", "worm"] = "fly",
) -> list:
    """
    View the enrichr libraries provided by Enrichr.
    See Also: https://maayanlab.cloud/modEnrichr/ and https://github.com/zqfang/GSEApy/blob/master/gseapy/parser.py.
    Args:
        organism: Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
    Returns:
        A list of enrichr libraries from selected organism database.
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError(
            "You need to install the package `gseapy`."
            "\ninstall gseapy via `pip install gseapy`"
        )

    return gp.get_library_name(organism=organism)


def basic_enrichr(
    genes: Union[str, list],
    gene_sets: Union[str, list, tuple] = "GO_Biological_Process_2018",
    organism: Literal["human", "mouse", "yeast", "fly", "fish", "worm"] = "fly",
    outdir: str = "./enrichr",
    **kwargs,
):
    """
    Perform gene list enrichment analysis with default enrichr libraries using GSEApy.
    See Also: https://maayanlab.cloud/modEnrichr/ and https://gseapy.readthedocs.io/en/latest/gseapy_example.html#2.-Enrichr-Example.
    Args:
        genes: Flat file with list of genes, one gene id per row, or a python list object.
               The input `identifier` should be the same type to `gene_sets`.
        gene_sets: str, list, tuple of Enrichr Library name(s).
                str: 'KEGG_2016'
                list: ['KEGG_2016','KEGG_2013']
                Use comma to separate each other, e.g. "KEGG_2016,huMAP,GO_Biological_Process_2018"
                Please view the enrichr libraries through `basic_sets`.
        organism: Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
        outdir: Output file directory
        **kwargs: Other parameters used in gseapy.enrichr.
    Returns:
        An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError(
            "You need to install the package `gseapy`."
            "\ninstall gseapy via `pip install gseapy`"
        )

    basic_enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism=organism,
        outdir=outdir,
        no_plot=True,
        verbose=True,
        **kwargs,
    )

    return basic_enr


def custom_enrichr(
    genes: Union[str, list],
    gene_sets: Union[str, dict],
    background: Union[int, str, list],
    outdir: str = "./enrichr",
    **kwargs,
):
    """
    Perform gene list enrichment analysis with custom gene sets using GSEApy.
    See Also: https://maayanlab.cloud/modEnrichr/ and https://gseapy.readthedocs.io/en/latest/gseapy_example.html#2.-Enrichr-Example.
    Args:
        genes: Flat file with list of genes, one gene id per row, or a python list object.
               The input `identifier` should be the same type to `gene_sets`.
        gene_sets: Custom defined gene_sets (dict, or gmt file).
                dict, gene_sets={‘A’:[‘gene1’, ‘gene2’,…],
                                 ’B’:[‘gene2’, ‘gene4’,…], …}
                gmt: “genes.gmt”
        background: The background is a lookup table of expected ranks and variances for each term in the library.
                    See Also: https://maayanlab.cloud/FlyEnrichr/help#background&q=0.
                    There are 3 ways to set this argument:
                *. 1.(Recommended) Input a list of background genes. The background gene list is defined by your
                     experiment. e.g. the expressed genes in your RNA-seq. The gene identifier in gmt/dict should be the
                     same type to the background genes.
                *. 2.Specify a number, e.g. the number of total expressed genes. This works, but not recommend. It
                     assumes that all your genes could be found in background. If genes exist in gmt but not included
                     in background, they will affect the significance of the statistical test.
                *. 3.(Default) Set a Biomart dataset name. The background will be all annotated genes from the BioMart
                     datasets you’ve chosen. The program will try to retrieve the background information automatically.
        outdir: Output file directory
        **kwargs: Other parameters used in gseapy.enrichr.
    Returns:
        An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError(
            "You need to install the package `gseapy`."
            "\ninstall gseapy via `pip install gseapy`"
        )

    custom_enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        background=background,
        outdir=outdir,
        no_plot=True,
        verbose=True,
        **kwargs,
    )

    return custom_enr
