from pathlib import Path
import pandas as pd
import json
import altair as alt
# Import plotting functions from plotting
from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit, bowtie2_alignment_plot
from utils import parse_bowtielog, parse_fastp_report

# function to read in svg code
def return_svg(svg: str):
    with open(svg, "r") as f:
        return f.read()
    
# function to create the report
def html_template_report(
    sample_name: str,
    out_path: str,
    total_reads: int,
    number_aligned: int,
    number_unaligned: int,
    bowtie_plot: str,
    fastp_df: str,
    megahit_histogram: str,
    kaiju_raw: str,
    kraken_raw: str,
    kaiju_and_cat: str,
    cat_kaiju_df: str,
    svg: str,
) -> None:
    """
    Creates html report
    """
  
    html = f"""
    <!DOCTYPE html>
    <html>
        <head>
            <title>Report of {sample_name} </title>
            <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
            <script src="https://cdn.jsdelivr.net/npm/vega-lite@4"></script>
            <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
            
            <style>
                body {{
                font-family: Arial, sans-serif;
                font-size: 16px;
                line-height: 1.5;
                }}
                header {{
                background: linear-gradient(to right, #000000, #333333);
                color: white;
                padding: 40px;
                text-align: center;
                box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19);
                }}
                h1 {{
                font-size: 3em;
                margin: 0;
                text-align: center;
                }}
                h2 {{
                font-size: 2em;
                margin-bottom:30px;
                text-align: center;
                }}
                section {{
                margin: 80px 0;
                }}
                hr {{
                border: 0;
                height: 1px;
                background: #333;
                background-image: linear-gradient(to right, #ccc, #333, #ccc);
                }}
                table {{
                border-collapse: collapse;
                }}
                th, td {{
                border: 1px solid #ccc;
                padding: 10px;
                text-align: left;
                }}
                th {{
                background-color: #eee;
                }}
            </style>

        </head>
        
        <body>
        
            <header>
                <h1>
                Pandemic Prepardeness Report 
                </h1>
                <h2>
                Report of {sample_name}
                </h2>
            </header>
            
            <h2>
                Information about Reads
            </h2>
            <table style="margin-bottom:30px;">
                <tr>
                    <th>Total reads pairs</th>
                    <th>Reads aligned to human genome</th>
                    <th>Reads NOT aligned to human genome</th>
                </tr>
                <tr>
                    <td>{total_reads:,}</td>
                    <td>{number_aligned:,} ({number_aligned / total_reads * 100:.2f}%)</td>
                    <td>{number_unaligned:,} ({number_unaligned / total_reads * 100:.2f}%)</td>
                </tr>
            </table>

            <div id="aligned" style="display:flex;width:100%;height:100%;"></div>
            <script type="text/javascript">
                vegaEmbed("#aligned", {bowtie_plot});
            </script>
            
            <h2>
            Information from fastp
            </h2>
            <div style="margin: auto;">
            {fastp_df}
            </div>
            
            <hr />
            
            <!-- KRAKEN RAW -->
            <h2> 
            Raw reads classfied with Kraken and Kaiju
            </h2>
            
            <div style="border: 1px solid black; padding: 10px; margin: 10px; width: 1000px; border-radius: 10px; box-shadow: 5px 5px 5px grey; background-color: #cccccc; display:flex;justify-content:center;">
            <ul>
            <li>KRAKEN2 and KAIJU are both popular bioinformatics tools for accurately and efficiently identifying and categorizing large numbers of short DNA sequences.</li>
            <li>KRAKEN2 is based on a probabilistic model and uses a pre-constructed database of genomic sequences, while KAIJU uses a database of reference sequences and an efficient indexing method.</li>
            <li>KRAKEN2 and KAIJU have been shown to produce highly accurate results in a short amount of time, but KRAKEN2 may require more computational resources.</li>
            <li>Both KRAKEN2 and KAIJU are open-source and freely available for use, making them valuable resources for researchers in the field of genomics.</li>
            </ul>
            </div>

            <h3>
            Kraken classification
            </h3>
            <div id="kraken_raw" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#kraken_raw", {kraken_raw});
            </script>
            
            <!-- KAIJU RAW -->
            <h3>
            Kaiju classification
            </h3>
            <div id="kaiju_raw" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#kaiju_raw", {kaiju_raw});
            </script>
            
            <hr />
            
            <h2>
            Contig information
            </h2>
            <h3>
            Histogram of megahit contigs
            </h3>
            <div style="border: 1px solid black; padding: 10px; margin: 10px; width: 1000px; border-radius: 10px; box-shadow: 5px 5px 5px grey; background-color: #cccccc; display:flex;justify-content:center;">
            <ul>
            <li>MEGAHIT is a popular and highly-efficient bioinformatics tool for assembling large and complex genomes.</li>
            <li>It is based on a De Bruijn graph approach and is capable of accurately assembling millions of reads into a complete genome in a short amount of time.</li>
            <li>MEGAHIT has been widely used in a variety of genomic studies and has been shown to produce high-quality assemblies that are comparable to those generated by other leading assemblers.</li>
            </ul>
            </div>


            <!-- PLOT MEGAHIT -->
            
            <div id="megahit_histo" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#megahit_histo", {megahit_histogram});
            </script>
            
            <!-- PLOT Kaiju and Cat -->
            <h3>
            Contigs classified with Kaiju and CAT
            </h3>
            
            <div style="border: 1px solid black; padding: 10px; margin: 10px; width: 1000px; border-radius: 10px; box-shadow: 5px 5px 5px grey; background-color: #cccccc; display:flex;justify-content:center;">
            <ul>
            <li>CAT (Classification and Annotation Tool) is a popular bioinformatics tool that uses a combination of machine learning algorithms and a pre-constructed database of genomic sequences to classify and annotate DNA sequences.</li>
            <li>The BLASTP aspect of CAT involves using the Basic Local Alignment Search Tool (BLASTP) to compare a query sequence to the database of reference sequences and identify similar sequences.</li>
            <li>BLASTP is a widely-used algorithm for sequence alignment and is known for its speed and accuracy in identifying similar sequences.</li>
          </ul>
          </div>
            
            <div id="kaiju_and_cat" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#kaiju_and_cat", {kaiju_and_cat});
            </script>
            
            <h3>
            Table containing information about contigs
            </h3>
            {cat_kaiju_df}
            <hr />
            
            
        </body>
    </html>
    """
    
    output = Path(out_path) / f"{sample_name}-report.html"
    with open(output, "w") as f:
        print(html, file=f)
        
        
# function that writes the report using the right samples
def create_report(
    sample: str,
    out_path: str,
) -> None:
    """
    Generate the report
    """
    # Sample and sample name
    sample = Path(sample)
    sample_name = sample.parts[-1]
    
    # Number of bars to include in the figures:
    number = 10
    
    # Number of reads and number of reads aligned to reference genome (from bowtie2logfile)
    bowtie2log = list(sample.rglob("*bowtie_raw.log"))[0]
    
    total_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie2log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned
    
    # Bowtie2 alignment plot:
    bowtie_plot = bowtie2_alignment_plot.plot_alignment(bowtie2log).to_json()
    
    # fastp dataframe
    fastp_report = list(sample.rglob("*fastp/*.html"))[0]
    fastp_df = parse_fastp_report.parse_fastp(fastp_report).to_html(classes=["center-table"])
    
    # Raw bracken and kaiju report
    cleaned_bracken_report = list(sample.rglob("*bracken_raw.csv"))[0]
    cleaned_kaiju_report = list(sample.rglob("*kaiju_raw.csv"))[0]
    
    # Raw bracken and kaiju plots
    bracken_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, number=number,virus_only=True
    )

    bracken_domain_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, level="domain", virus_only=False
    )

    kaiju_raw_plot = kaiju_raw.bar_chart_kaiju_raw(file=cleaned_kaiju_report).to_json()

    species_and_domain_bracken = (
        alt.hconcat(bracken_bar_plot, bracken_domain_bar_plot)
        .resolve_scale(color="independent")
    ).to_json()
    
    
    # Contigs (Megahit)
    megahit_csv = list(sample.rglob("megahit/*.csv"))[0]
    megahit_histogram = contig_quality.megahit_contig_histogram(file=megahit_csv).to_json()
    
    # Contigs (CAT and Kaiju)
    # files
    kaiju_megahit_report = list(sample.rglob("*megahit.out"))[0]
    cat_megahit_out = list(sample.rglob("*contigs_names.txt"))[0]
    cat_kaiju_csv = list(sample.rglob("*cat_kaiju_merged.csv"))[0]

    # plots
    kaiju_bar_plot = kaiju_megahit.bar_chart_kaiju_megahit(file=kaiju_megahit_report)
    cat_bar_plot = cat_megahit.bar_chart_cat_megahit(file=cat_megahit_out)
    kaiju_and_cat = (
        alt.hconcat(kaiju_bar_plot, cat_bar_plot)
        .resolve_scale(color="independent")
        .to_json()
    )

    # cat and kaiju dataframe
    cat_kaiju_df = pd.read_csv(cat_kaiju_csv)[["name", "taxon_id", "length", "last_level_kaiju", "last_level_cat"]].head(10).to_html()


    # test svg
    svg = return_svg("visualization (1).svg")
    
    # generate the html report
    html_template_report(
        sample_name=sample_name, 
        out_path=out_path, 
        total_reads=total_reads,
        number_aligned=number_aligned,
        number_unaligned=number_unaligned,
        kraken_raw=species_and_domain_bracken, 
        kaiju_raw=kaiju_raw_plot, 
        svg=svg,
        bowtie_plot=bowtie_plot,
        fastp_df=fastp_df,
        megahit_histogram=megahit_histogram,
        kaiju_and_cat=kaiju_and_cat,
        cat_kaiju_df=cat_kaiju_df,
    )

# read in the data (testing sample11 for now)
sample_folder = Path("../virusclassification_nextflow/results/")
samples = [x for x in sample_folder.iterdir() if x.is_dir() and not x.stem.startswith(".")]

for sample in samples:
    # create report
    create_report(sample=sample, out_path=".")